// External sort-merge of annotator *.out tables into per-sample CSVs.
// Streams seq/taxID/features only so Kraken k-mer LCA maps never sit in RAM.

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

#if __cplusplus >= 201703L
#include <filesystem>
namespace fs = std::filesystem;
#else
#error "C++17 is required"
#endif

namespace {

constexpr size_t kDefaultChunk = 500000;
constexpr size_t kIoBuf = 1 << 20;

enum class ToolKind { Kaiju, Kraken, Kraken2, Custom };

std::string trim_cr(std::string s) {
    if (!s.empty() && s.back() == '\r') {
        s.pop_back();
    }
    return s;
}

std::string to_lower(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (true) {
        size_t pos = line.find('\t', start);
        if (pos == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

struct Schema {
    std::vector<std::string> tax_ids;
    std::vector<std::string> feat_ids;
    int seq_col = 0;
    std::vector<int> tax_cols;
    std::vector<int> feat_cols;
};

struct Record {
    std::string seq;
    std::vector<std::string> tax;
    std::vector<std::string> feat;
};

void resize_record(Record* rec, const Schema& sch) {
    rec->tax.assign(sch.tax_ids.size(), "0");
    rec->feat.assign(sch.feat_ids.size(), "");
}

std::string drop_prefix_ci(const std::string& name, const std::string& prefix) {
    std::string low = to_lower(name);
    std::string pre = to_lower(prefix);
    if (low == pre) {
        return "";
    }
    if (low.size() > pre.size() && low.compare(0, pre.size(), pre) == 0 && name[pre.size()] == '_') {
        return name.substr(pre.size() + 1);
    }
    return name;
}

bool is_skip_header(const std::string& h) {
    std::string l = to_lower(h);
    return l == "classified" || l == "taxa" || l == "k-mer" || l == "kmer" || l == "k_mer" ||
           l == "sample" || l == "true" || l == "read_type";
}

bool is_acgt_kmer_id(const std::string& name) {
    if (name.empty() || name.size() > 8) {
        return false;
    }
    for (unsigned char ch : name) {
        char u = static_cast<char>(std::toupper(ch));
        if (u != 'A' && u != 'C' && u != 'G' && u != 'T') {
            return false;
        }
    }
    return true;
}

bool is_latent_id(const std::string& name) {
    std::string l = to_lower(name);
    if (l.size() < 2 || l[0] != 'z') {
        return false;
    }
    for (size_t i = 1; i < l.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(l[i]))) {
            return false;
        }
    }
    return true;
}

bool is_allowed_feature_header(const std::string& h) {
    std::string l = to_lower(h);
    if (l == "feat" || l.rfind("feat_", 0) == 0 || l.rfind("feat.", 0) == 0) {
        return true;
    }
    if (l == "length" || l == "gc" || l == "gc_content" || l == "mag_id" || l == "map_mag_id") {
        return true;
    }
    return is_acgt_kmer_id(h) || is_latent_id(h);
}

bool is_gtdb_taxid(const std::string& text) {
    return text.find("__") != std::string::npos;
}

bool is_nonnegative_taxid(const std::string& text) {
    if (text.empty()) {
        return false;
    }
    std::string body = text;
    if (body.size() > 2 && body.compare(body.size() - 2, 2, ".0") == 0) {
        body.resize(body.size() - 2);
    }
    for (unsigned char ch : body) {
        if (!std::isdigit(ch)) {
            return false;
        }
    }
    return true;
}

bool is_unclassified_taxid(const std::string& text) {
    std::string l = to_lower(text);
    return l.empty() || l == "0" || l == "unclassified" || l == "unclassified_root" || l == "nan" ||
           l == "none" || l == "na" || l == "<na>" || l == "null" || l == "other";
}

void require_taxid_token(const std::string& raw) {
    if (is_unclassified_taxid(raw) || is_nonnegative_taxid(raw) || is_gtdb_taxid(raw)) {
        return;
    }
    throw std::runtime_error(
        "Invalid taxid '" + raw +
        "'. Taxonomy values must be non-negative integers, unclassified, or GTDB tokens");
}

bool is_tax_header(const std::string& h) {
    std::string l = to_lower(h);
    if (l.rfind("feat", 0) == 0) {
        return false;
    }
    if (l == "taxa" || l == "taxon" || l == "taxonomy") {
        return false;
    }
    if (l == "taxid" || l == "tax_id" || l == "tax") {
        return true;
    }
    return l.rfind("taxid_", 0) == 0 || l.rfind("tax_", 0) == 0;
}

bool looks_like_header(const std::vector<std::string>& f) {
    if (f.empty()) {
        return false;
    }
    std::string a = to_lower(f[0]);
    while (!a.empty() && (a[0] == '#' || a[0] == '@')) {
        a.erase(a.begin());
    }
    return a == "seq" || a == "read_id" || a == "readid" || a == "sequenceid" ||
           a == "anonymous_read_id";
}

std::string tax_header_id(const std::string& h) {
    std::string id = drop_prefix_ci(h, "taxID");
    if (id != h) {
        return id;
    }
    id = drop_prefix_ci(h, "taxid");
    if (id != h) {
        return id;
    }
    id = drop_prefix_ci(h, "tax_id");
    if (id != h) {
        return id;
    }
    return drop_prefix_ci(h, "tax");
}

std::string feat_header_id(const std::string& h) {
    std::string id = drop_prefix_ci(h, "feat");
    return id == h ? h : id;
}

Schema schema_from_header(const std::vector<std::string>& header) {
    Schema s;
    s.seq_col = 0;
    for (size_t i = 0; i < header.size(); ++i) {
        std::string l = to_lower(header[i]);
        if (l == "seq" || l == "read_id" || l == "readid" || l == "sequenceid" ||
            l == "anonymous_read_id") {
            s.seq_col = static_cast<int>(i);
            continue;
        }
        if (is_skip_header(header[i])) {
            continue;
        }
        if (is_tax_header(header[i])) {
            s.tax_ids.push_back(tax_header_id(header[i]));
            s.tax_cols.push_back(static_cast<int>(i));
        } else if (is_allowed_feature_header(header[i])) {
            s.feat_ids.push_back(feat_header_id(header[i]));
            s.feat_cols.push_back(static_cast<int>(i));
        } else {
            throw std::runtime_error(
                "Column '" + header[i] +
                "' is not a taxID column or an allowed feature column");
        }
    }
    return s;
}

Schema native_schema(ToolKind kind) {
    Schema s;
    s.seq_col = 1;
    s.tax_ids.push_back("");
    s.tax_cols.push_back(2);
    if (kind == ToolKind::Kraken || kind == ToolKind::Kraken2) {
        s.feat_ids.push_back("length");
        s.feat_cols.push_back(3);
    }
    return s;
}

Schema unheadered_custom_schema(size_t nfields) {
    Schema s;
    s.seq_col = 0;
    if (nfields >= 2) {
        s.tax_ids.push_back("");
        s.tax_cols.push_back(1);
    }
    if (nfields == 3) {
        s.feat_ids.push_back("length");
        s.feat_cols.push_back(2);
    } else if (nfields > 3) {
        throw std::runtime_error(
            "Unheadered custom table has extra columns; name feature columns "
            "(feat_*, length, GC, MAG_ID, ACGT k-mer ids, or z<n>)");
    }
    return s;
}

void write_schema(std::ostream& out, const Schema& s) {
    out << "#schema\tTAX\t" << s.tax_ids.size();
    for (const auto& id : s.tax_ids) {
        out << '\t' << id;
    }
    out << "\tFEAT\t" << s.feat_ids.size();
    for (const auto& id : s.feat_ids) {
        out << '\t' << id;
    }
    out << '\n';
}

bool parse_schema_line(const std::string& line, Schema* s) {
    if (line.size() < 7 || line.compare(0, 7, "#schema") != 0) {
        return false;
    }
    auto f = split_tab(line);
    s->tax_ids.clear();
    s->feat_ids.clear();
    size_t i = 1;
    if (i < f.size() && to_lower(f[i]) == "tax") {
        ++i;
    }
    if (i >= f.size()) {
        return true;
    }
    int n_tax = 0;
    try {
        n_tax = std::stoi(f[i++]);
    } catch (...) {
        return false;
    }
    for (int k = 0; k < n_tax && i < f.size(); ++k) {
        s->tax_ids.push_back(f[i++]);
    }
    if (i < f.size() && to_lower(f[i]) == "feat") {
        ++i;
    }
    if (i >= f.size()) {
        return true;
    }
    int n_feat = 0;
    try {
        n_feat = std::stoi(f[i++]);
    } catch (...) {
        return false;
    }
    for (int k = 0; k < n_feat && i < f.size(); ++k) {
        s->feat_ids.push_back(f[i++]);
    }
    return true;
}

std::string combined_col(const std::string& prefix, const std::string& tool, size_t i,
                         const std::string& id) {
    std::string base = prefix + "_" + tool + "_" + std::to_string(i);
    if (id.empty()) {
        return base;
    }
    return base + "_" + id;
}


struct Options {
    std::string input_dir;
    std::string output_dir;
    int split_n = 1;
    size_t chunk_rows = kDefaultChunk;
    std::string tmp_dir;
    bool keep_tmp = false;
    std::string truth_table;
    std::string add_file;
    std::string add_name;
    std::string base_csv;
    std::string out_csv;
    std::string retrue_csv;
};

std::vector<std::string> split_underscore(const std::string& s) {
    std::vector<std::string> parts;
    size_t start = 0;
    while (true) {
        size_t pos = s.find('_', start);
        if (pos == std::string::npos) {
            parts.push_back(s.substr(start));
            break;
        }
        parts.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return parts;
}

std::string csv_escape(const std::string& s) {
    bool quote = s.find_first_of(",\"\n\r") != std::string::npos;
    if (!quote) {
        return s;
    }
    std::string out;
    out.push_back('"');
    for (char c : s) {
        if (c == '"') {
            out.push_back('"');
        }
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}

std::string normalize_taxid(std::string t) {
    require_taxid_token(t);
    if (is_unclassified_taxid(t)) {
        return "0";
    }
    if (t.size() > 2 && t.compare(t.size() - 2, 2, ".0") == 0 && is_nonnegative_taxid(t)) {
        t.resize(t.size() - 2);
    }
    return t;
}

std::string first_length_token(const std::string& raw) {
    if (raw.empty()) {
        return "";
    }
    size_t bar = raw.find('|');
    return bar == std::string::npos ? raw : raw.substr(0, bar);
}

std::string extract_kraken2_taxid(const std::string& taxa) {
    const std::string key = "taxid ";
    for (size_t i = 0; i + key.size() <= taxa.size(); ++i) {
        bool match = true;
        for (size_t j = 0; j < key.size(); ++j) {
            char a = static_cast<char>(std::tolower(static_cast<unsigned char>(taxa[i + j])));
            if (a != key[j]) {
                match = false;
                break;
            }
        }
        if (!match) {
            continue;
        }
        size_t k = i + key.size();
        std::string digits;
        while (k < taxa.size() && std::isdigit(static_cast<unsigned char>(taxa[k]))) {
            digits.push_back(taxa[k]);
            ++k;
        }
        if (!digits.empty()) {
            return digits;
        }
    }
    return "0";
}

std::string extract_read_type(const std::string& seq) {
    const std::string key = "read_type:";
    for (size_t i = 0; i < seq.size(); ++i) {
        bool match = true;
        for (size_t j = 0; j < key.size(); ++j) {
            if (i + j >= seq.size()) {
                match = false;
                break;
            }
            char a = static_cast<char>(std::tolower(static_cast<unsigned char>(seq[i + j])));
            if (a != key[j]) {
                match = false;
                break;
            }
        }
        if (!match) {
            continue;
        }
        size_t k = i + key.size();
        std::string token;
        while (k < seq.size()) {
            unsigned char c = static_cast<unsigned char>(seq[k]);
            if (!(std::isalnum(c) || seq[k] == '_' || seq[k] == '+' || seq[k] == '-')) {
                break;
            }
            token.push_back(static_cast<char>(std::tolower(c)));
            ++k;
        }
        if (!token.empty()) {
            return token;
        }
    }
    return "";
}

std::string extract_true_taxid(const std::string& seq) {
    const std::string key = "taxid:";
    for (size_t i = 0; i < seq.size(); ++i) {
        bool match = true;
        for (size_t j = 0; j < key.size(); ++j) {
            if (i + j >= seq.size()) {
                match = false;
                break;
            }
            char a = static_cast<char>(std::tolower(static_cast<unsigned char>(seq[i + j])));
            if (a != key[j]) {
                match = false;
                break;
            }
        }
        if (!match) {
            continue;
        }
        size_t k = i + key.size();
        std::string digits;
        while (k < seq.size() && std::isdigit(static_cast<unsigned char>(seq[k]))) {
            digits.push_back(seq[k]);
            ++k;
        }
        if (!digits.empty()) {
            return digits;
        }
    }
    std::string lower = to_lower(seq);
    static const std::pair<const char*, const char*> prefixes[] = {
        {"scer.fna", "4932"},
        {"ecoli.fna", "562"},
        {"hsap.fna", "9606"},
        {"phix.fna", "2886930"},
    };
    for (const auto& p : prefixes) {
        std::string pref = p.first;
        if (lower.rfind(pref, 0) == 0 || lower.find(std::string("|") + pref) != std::string::npos) {
            return p.second;
        }
    }
    return "";
}

struct TruthIndex {
    std::unordered_map<std::string, std::string> exact;
    std::vector<std::string> prefixes;
    std::unordered_map<std::string, std::string> prefix_tax;

    void add(const std::string& key, std::string tax) {
        tax = normalize_taxid(std::move(tax));
        if (key.empty() || tax.empty() || tax == "0") {
            return;
        }
        exact[key] = tax;
        if (key.size() > 2 && key[key.size() - 2] == '/') {
            exact[key.substr(0, key.size() - 2)] = tax;
        }
        if (key.find('/') == std::string::npos && key.size() >= 4) {
            prefixes.push_back(key);
            prefix_tax[key] = tax;
        }
    }

    void freeze() {
        std::sort(prefixes.begin(), prefixes.end());
        prefixes.erase(std::unique(prefixes.begin(), prefixes.end()), prefixes.end());
    }

    std::string lookup(const std::string& seq) const {
        auto hit = exact.find(seq);
        if (hit != exact.end()) {
            return hit->second;
        }
        std::string token = seq;
        size_t sp = token.find_first_of(" \t");
        if (sp != std::string::npos) {
            token.resize(sp);
            hit = exact.find(token);
            if (hit != exact.end()) {
                return hit->second;
            }
        }
        if (token.size() > 2 && token[token.size() - 2] == '/') {
            hit = exact.find(token.substr(0, token.size() - 2));
            if (hit != exact.end()) {
                return hit->second;
            }
        }
        std::string best;
        for (const auto& key : prefixes) {
            if (key.size() > token.size() || key.size() <= best.size()) {
                continue;
            }
            if (token.compare(0, key.size(), key) == 0) {
                best = key;
            }
        }
        if (!best.empty()) {
            auto pt = prefix_tax.find(best);
            if (pt != prefix_tax.end()) {
                return pt->second;
            }
        }
        return "";
    }
};

std::string strip_meta_prefix(const std::string& s) {
    if (!s.empty() && s[0] == '#') {
        return s.substr(1);
    }
    if (s.size() >= 2 && s[0] == '@' && s[1] == '@') {
        return s.substr(2);
    }
    return s;
}

std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    for (char c : line) {
        if (c == ' ' || c == '\t') {
            if (!cur.empty()) {
                fields.push_back(cur);
                cur.clear();
            }
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) {
        fields.push_back(cur);
    }
    return fields;
}

int col_index(const std::vector<std::string>& header, std::initializer_list<const char*> names) {
    for (size_t i = 0; i < header.size(); ++i) {
        std::string h = to_lower(header[i]);
        while (!h.empty() && (h.front() == '@' || h.front() == '#')) {
            h.erase(h.begin());
        }
        for (const char* n : names) {
            if (h == n) {
                return static_cast<int>(i);
            }
        }
    }
    return -1;
}

bool load_truth_table(const fs::path& path, TruthIndex* idx) {
    std::ifstream in(path, std::ios::in | std::ios::binary);
    if (!in) {
        std::cerr << "Truth table not found: " << path << "\n";
        return false;
    }
    std::string line;
    int seq_cols[6];
    int n_seq = 0;
    int tax_col = -1;
    bool header_done = false;
    auto remember_seq = [&](int c) {
        if (c < 0) {
            return;
        }
        for (int i = 0; i < n_seq; ++i) {
            if (seq_cols[i] == c) {
                return;
            }
        }
        if (n_seq < 4) {
            seq_cols[n_seq++] = c;
        }
    };
    while (std::getline(in, line)) {
        line = trim_cr(line);
        if (line.empty() || (line[0] == '@' && (line.size() < 2 || line[1] != '@'))) {
            continue;
        }
        if (!header_done) {
            std::string low = to_lower(line);
            bool looks_header = line[0] == '#' || (line.size() >= 2 && line[0] == '@' && line[1] == '@') ||
                                low.find("tax_id") != std::string::npos || low.find("taxid") != std::string::npos ||
                                low.find("anonymous_read") != std::string::npos || low.find("sequenceid") != std::string::npos;
            if (looks_header) {
                auto header = split_tab(strip_meta_prefix(line));
                if (header.size() < 2) {
                    header = split_ws(strip_meta_prefix(line));
                }
                remember_seq(col_index(header, {"anonymous_read_id", "sequenceid", "seq"}));
                remember_seq(col_index(header, {"read_id"}));
                remember_seq(col_index(header, {"contig_id"}));
                tax_col = col_index(header, {"tax_id", "taxid", "true"});
                if (tax_col < 0) {
                    for (size_t i = 0; i < header.size(); ++i) {
                        if (to_lower(header[i]).find("tax") != std::string::npos) {
                            tax_col = static_cast<int>(i);
                            break;
                        }
                    }
                }
                header_done = true;
                if (tax_col >= 0 && n_seq > 0) {
                    continue;
                }
            }
            if (!header_done) {
                header_done = true;
                tax_col = 1;
                seq_cols[0] = 0;
                n_seq = 1;
            }
        }
        auto fields = split_tab(line);
        if (tax_col >= 0 && static_cast<int>(fields.size()) <= tax_col) {
            fields = split_ws(line);
        }
        if (tax_col < 0 || tax_col >= static_cast<int>(fields.size())) {
            continue;
        }
        std::string tax = fields[static_cast<size_t>(tax_col)];
        if (n_seq == 0) {
            idx->add(fields[0], tax);
            continue;
        }
        for (int i = 0; i < n_seq; ++i) {
            int c = seq_cols[i];
            if (c >= 0 && c < static_cast<int>(fields.size())) {
                idx->add(fields[static_cast<size_t>(c)], tax);
            }
        }
    }
    idx->freeze();
    return true;
}

std::string true_taxid_for(const std::string& seq, const TruthIndex* truth) {
    if (truth != nullptr) {
        std::string hit = truth->lookup(seq);
        if (!hit.empty()) {
            return hit;
        }
    }
    return extract_true_taxid(seq);
}

std::string match_tool(const std::string& filename) {
    if (filename.size() < 4 || filename.compare(filename.size() - 4, 4, ".out") != 0) {
        return "";
    }
    size_t last_dot = filename.rfind('.');
    size_t prev_dot = filename.rfind('.', last_dot == std::string::npos ? 0 : last_dot - 1);
    if (prev_dot == std::string::npos || last_dot == std::string::npos) {
        return "";
    }
    std::string arg = filename.substr(prev_dot + 1, last_dot - prev_dot - 1);
    if (arg.rfind("custom_", 0) == 0) {
        return arg.substr(7);
    }
    if (arg == "kraken1") {
        return "kraken1";
    }
    if (arg == "krakenuniq" || arg == "krakenu" || arg == "krakenunique") {
        return "krakenuniq";
    }
    if (arg == "metaphlan" || arg == "metaphlan4" || arg == "mpa" || arg == "mp4") {
        return "metaphlan";
    }
    if (arg == "dummy9606" || arg == "constant9606" || arg == "constant" || arg == "dummy" ||
        arg == "random") {
        return arg;
    }
    return arg;
}

std::string sample_name_from(const std::string& filename, int split_n) {
    // ``{sample}_{run}.{tool}.out`` or ``{sample}_{run}.custom_{tool}.out``.
    // Sample ids like ``sample_0`` or ``1_full`` contain underscores, so taking
    // the first ``-s`` token is wrong. Drop the ``.{tool}.out`` suffix, then
    // the last ``_`` component (the run name).
    size_t last_dot = filename.rfind('.');
    size_t prev_dot = (last_dot == std::string::npos || last_dot == 0)
        ? std::string::npos
        : filename.rfind('.', last_dot - 1);
    if (prev_dot != std::string::npos) {
        std::string stem = filename.substr(0, prev_dot);
        size_t us = stem.rfind('_');
        if (us != std::string::npos && us > 0) {
            return stem.substr(0, us);
        }
        if (!stem.empty()) {
            return stem;
        }
    }
    auto parts = split_underscore(filename);
    if (split_n < 1) {
        split_n = 1;
    }
    if (static_cast<int>(parts.size()) < split_n) {
        return filename;
    }
    std::string out = parts[0];
    for (int i = 1; i < split_n; ++i) {
        out.push_back('_');
        out += parts[static_cast<size_t>(i)];
    }
    return out;
}

ToolKind tool_kind(const std::string& tool) {
    if (tool == "kraken2") {
        return ToolKind::Kraken2;
    }
    if (tool == "kraken" || tool == "kraken1" || tool == "krakenuniq" || tool == "krakenu" ||
        tool == "krakenunique") {
        return ToolKind::Kraken;
    }
    if (tool == "kaiju") {
        return ToolKind::Kaiju;
    }
    return ToolKind::Custom;
}

bool fill_from_fields(const std::vector<std::string>& f, const Schema& sch, Record* rec) {
    rec->seq.clear();
    resize_record(rec, sch);
    if (sch.seq_col < 0 || sch.seq_col >= static_cast<int>(f.size())) {
        if (f.empty()) {
            return false;
        }
        rec->seq = f[0];
    } else {
        rec->seq = f[static_cast<size_t>(sch.seq_col)];
    }
    for (size_t i = 0; i < sch.tax_cols.size() && i < rec->tax.size(); ++i) {
        int c = sch.tax_cols[i];
        if (c >= 0 && c < static_cast<int>(f.size())) {
            rec->tax[i] = normalize_taxid(f[static_cast<size_t>(c)]);
        }
    }
    for (size_t i = 0; i < sch.feat_cols.size() && i < rec->feat.size(); ++i) {
        int c = sch.feat_cols[i];
        if (c >= 0 && c < static_cast<int>(f.size())) {
            std::string val = f[static_cast<size_t>(c)];
            if (i < sch.feat_ids.size() && sch.feat_ids[i] == "length") {
                val = first_length_token(val);
            }
            rec->feat[i] = val;
        }
    }
    for (char& ch : rec->seq) {
        if (ch == '\t') {
            ch = ' ';
        }
    }
    // Custom files often copy the FASTQ token, including a mate mark. Native
    // tools and fastq_record_id drop only a trailing /1 or /2; r/10 stays.
    if (rec->seq.size() >= 2 && rec->seq[rec->seq.size() - 2] == '/' &&
        (rec->seq.back() == '1' || rec->seq.back() == '2')) {
        rec->seq.resize(rec->seq.size() - 2);
    }
    return !rec->seq.empty();
}

bool parse_line(const std::string& line, ToolKind kind, const Schema& sch, Record* rec) {
    if (line.empty() || line[0] == '#') {
        return false;
    }
    auto f = split_tab(line);
    if (kind == ToolKind::Custom) {
        if (f.size() < 1) {
            return false;
        }
        return fill_from_fields(f, sch, rec);
    }
    if (f.size() < 3) {
        return false;
    }
    resize_record(rec, sch);
    rec->seq = f[1];
    if (kind == ToolKind::Kraken2) {
        rec->tax[0] = normalize_taxid(extract_kraken2_taxid(f[2]));
        if (!rec->feat.empty() && f.size() >= 4) {
            rec->feat[0] = first_length_token(f[3]);
        }
    } else if (kind == ToolKind::Kraken) {
        rec->tax[0] = normalize_taxid(f[2]);
        if (!rec->feat.empty() && f.size() >= 4) {
            rec->feat[0] = first_length_token(f[3]);
        }
    } else {
        rec->tax[0] = normalize_taxid(f[2]);
    }
    for (char& ch : rec->seq) {
        if (ch == '\t') {
            ch = ' ';
        }
    }
    return !rec->seq.empty();
}

void write_tsv_record(std::ostream& out, const Record& r) {
    out << r.seq;
    for (const auto& t : r.tax) {
        out << '\t' << t;
    }
    for (const auto& feat : r.feat) {
        out << '\t' << feat;
    }
    out << '\n';
}

bool read_tsv_record(std::istream& in, const Schema& sch, Record* rec) {
    std::string line;
    if (!std::getline(in, line)) {
        return false;
    }
    line = trim_cr(line);
    if (line.empty() || (line.size() >= 7 && line.compare(0, 7, "#schema") == 0)) {
        return read_tsv_record(in, sch, rec);
    }
    auto f = split_tab(line);
    rec->seq = f.empty() ? "" : f[0];
    resize_record(rec, sch);
    size_t i = 1;
    for (size_t k = 0; k < rec->tax.size(); ++k) {
        rec->tax[k] = i < f.size() ? normalize_taxid(f[i++]) : "0";
    }
    for (size_t k = 0; k < rec->feat.size(); ++k) {
        rec->feat[k] = i < f.size() ? f[i++] : "";
    }
    return true;
}

class BufferedIn {
public:
    explicit BufferedIn(const fs::path& path) : file_(path, std::ios::in | std::ios::binary) {
        std::string line;
        if (std::getline(file_, line)) {
            line = trim_cr(line);
            if (!parse_schema_line(line, &schema_)) {
                file_.clear();
                file_.seekg(0);
                schema_.tax_ids = {""};
                schema_.feat_ids = {"length"};
            }
        }
    }
    bool ok() const { return static_cast<bool>(file_); }
    const Schema& schema() const { return schema_; }
    bool next(Record* rec) { return read_tsv_record(file_, schema_, rec); }

private:
    std::ifstream file_;
    Schema schema_;
};

fs::path write_run(const fs::path& dir, size_t idx, std::vector<Record>& recs, const Schema& sch) {
    std::stable_sort(recs.begin(), recs.end(), [](const Record& a, const Record& b) {
        return a.seq < b.seq;
    });
    fs::path path = dir / ("run_" + std::to_string(idx) + ".tsv");
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    std::vector<char> buf(kIoBuf);
    out.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    write_schema(out, sch);
    for (const Record& r : recs) {
        write_tsv_record(out, r);
    }
    return path;
}

fs::path merge_runs(const std::vector<fs::path>& runs, const fs::path& out_path, const Schema& sch) {
    struct Item {
        Record rec;
        size_t src;
        bool operator>(const Item& o) const {
            if (rec.seq != o.rec.seq) {
                return rec.seq > o.rec.seq;
            }
            return src > o.src;
        }
    };
    std::vector<BufferedIn> streams;
    streams.reserve(runs.size());
    for (const auto& p : runs) {
        streams.emplace_back(p);
    }
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> heap;
    for (size_t i = 0; i < streams.size(); ++i) {
        Record r;
        if (streams[i].next(&r)) {
            heap.push(Item{std::move(r), i});
        }
    }
    std::ofstream out(out_path, std::ios::out | std::ios::binary | std::ios::trunc);
    std::vector<char> buf(kIoBuf);
    out.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    write_schema(out, sch);
    while (!heap.empty()) {
        Item cur = heap.top();
        heap.pop();
        write_tsv_record(out, cur.rec);
        Record nxt;
        if (streams[cur.src].next(&nxt)) {
            heap.push(Item{std::move(nxt), cur.src});
        }
    }
    return out_path;
}

fs::path sort_annotator(const fs::path& input, ToolKind kind, const fs::path& tmp,
                        size_t chunk_rows, size_t file_index) {
    fs::path work = tmp / ("ann_" + std::to_string(file_index));
    fs::create_directories(work);
    std::ifstream in(input, std::ios::in | std::ios::binary);
    std::vector<char> buf(kIoBuf);
    in.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));

    Schema sch = native_schema(kind);
    bool header_checked = kind != ToolKind::Custom;
    std::vector<Record> chunk;
    chunk.reserve(std::min(chunk_rows, static_cast<size_t>(65536)));
    std::vector<fs::path> runs;
    size_t run_i = 0;
    std::string line;
    while (std::getline(in, line)) {
        line = trim_cr(line);
        if (line.empty()) {
            continue;
        }
        if (!header_checked) {
            header_checked = true;
            auto fields = split_tab(line);
            if (looks_like_header(fields)) {
                sch = schema_from_header(fields);
                continue;
            }
            sch = unheadered_custom_schema(fields.size());
        }
        Record rec;
        if (!parse_line(line, kind, sch, &rec)) {
            continue;
        }
        chunk.push_back(std::move(rec));
        if (chunk.size() >= chunk_rows) {
            runs.push_back(write_run(work, run_i++, chunk, sch));
            chunk.clear();
        }
    }
    if (!chunk.empty()) {
        runs.push_back(write_run(work, run_i++, chunk, sch));
        chunk.clear();
        chunk.shrink_to_fit();
    }
    fs::path sorted = tmp / ("sorted_" + std::to_string(file_index) + ".tsv");
    if (runs.empty()) {
        std::ofstream out(sorted, std::ios::trunc);
        write_schema(out, sch);
        return sorted;
    }
    if (runs.size() == 1) {
        fs::rename(runs[0], sorted);
        return sorted;
    }
    merge_runs(runs, sorted, sch);
    return sorted;
}

struct JoinStream {
    BufferedIn in;
    Record rec;
    bool alive = false;
    explicit JoinStream(const fs::path& path) : in(path) { alive = in.next(&rec); }
};

bool share_seq(const fs::path& left, const fs::path& right) {
    BufferedIn a(left);
    BufferedIn b(right);
    Record ra;
    Record rb;
    bool have_a = a.next(&ra);
    bool have_b = b.next(&rb);
    while (have_a && have_b) {
        if (ra.seq == rb.seq) {
            return true;
        }
        if (ra.seq < rb.seq) {
            have_a = a.next(&ra);
        } else {
            have_b = b.next(&rb);
        }
    }
    return false;
}

// Natives always open rows. A custom file whose read ids share nothing with
// those natives must not outer-join its own ids in (that doubles the table
// with zero-filled rows). With no native file, the first tool is the anchor.
std::vector<char> row_sources(const std::vector<fs::path>& sorted, const std::vector<std::string>& tools) {
    std::vector<char> introduce(sorted.size(), 1);
    std::vector<size_t> natives;
    for (size_t i = 0; i < tools.size(); ++i) {
        if (tool_kind(tools[i]) != ToolKind::Custom) {
            natives.push_back(i);
        }
    }
    if (!natives.empty()) {
        for (size_t i = 0; i < tools.size(); ++i) {
            if (tool_kind(tools[i]) != ToolKind::Custom) {
                continue;
            }
            bool hit = false;
            for (size_t native : natives) {
                if (share_seq(sorted[i], sorted[native])) {
                    hit = true;
                    break;
                }
            }
            if (!hit) {
                introduce[i] = 0;
                std::cerr << "WARNING: " << tools[i]
                          << " read ids do not overlap the other annotators; "
                          << "its rows were not added to the table\n";
            }
        }
        return introduce;
    }
    for (size_t i = 1; i < sorted.size(); ++i) {
        if (!share_seq(sorted[0], sorted[i])) {
            introduce[i] = 0;
            std::cerr << "WARNING: " << tools[i]
                      << " read ids do not overlap the other annotators; "
                      << "its rows were not added to the table\n";
        }
    }
    return introduce;
}

void merge_sample(const std::vector<fs::path>& sorted, const std::vector<std::string>& tools,
                  const fs::path& csv_path, const TruthIndex* truth) {
    const std::vector<char> introduce = row_sources(sorted, tools);
    std::vector<JoinStream> streams;
    streams.reserve(sorted.size());
    for (const auto& p : sorted) {
        streams.emplace_back(p);
    }
    std::ofstream out(csv_path, std::ios::out | std::ios::binary | std::ios::trunc);
    std::vector<char> buf(kIoBuf);
    out.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));

    out << "seq";
    for (size_t i = 0; i < streams.size(); ++i) {
        const Schema& sch = streams[i].in.schema();
        for (const auto& id : sch.tax_ids) {
            out << ',' << csv_escape(combined_col("taxID", tools[i], i, id));
        }
        for (const auto& id : sch.feat_ids) {
            out << ',' << csv_escape(combined_col("feat", tools[i], i, id));
        }
    }
    out << ",length,true,read_type\n";

    uint64_t n = 0;
    while (true) {
        std::string min_seq;
        bool any = false;
        for (size_t i = 0; i < streams.size(); ++i) {
            if (!introduce[i] || !streams[i].alive) {
                continue;
            }
            if (!any || streams[i].rec.seq < min_seq) {
                min_seq = streams[i].rec.seq;
                any = true;
            }
        }
        if (!any) {
            break;
        }
        for (size_t i = 0; i < streams.size(); ++i) {
            if (introduce[i]) {
                continue;
            }
            auto& s = streams[i];
            while (s.alive && s.rec.seq < min_seq) {
                s.alive = s.in.next(&s.rec);
            }
        }
        std::vector<std::vector<Record>> copies(streams.size());
        for (size_t i = 0; i < streams.size(); ++i) {
            auto& s = streams[i];
            while (s.alive && s.rec.seq == min_seq) {
                copies[i].push_back(s.rec);
                s.alive = s.in.next(&s.rec);
            }
        }
        size_t rows = 1;
        for (size_t i = 0; i < copies.size(); ++i) {
            if (introduce[i] && copies[i].size() > rows) {
                rows = copies[i].size();
            }
        }
        for (size_t r = 0; r < rows; ++r) {
            std::string length;
            out << csv_escape(min_seq);
            for (size_t i = 0; i < streams.size(); ++i) {
                const Schema& sch = streams[i].in.schema();
                const Record* rec = nullptr;
                if (!copies[i].empty()) {
                    rec = &copies[i][r < copies[i].size() ? r : 0];
                }
                for (size_t k = 0; k < sch.tax_ids.size(); ++k) {
                    out << ',' << csv_escape(rec && k < rec->tax.size() ? rec->tax[k] : "0");
                }
                for (size_t k = 0; k < sch.feat_ids.size(); ++k) {
                    const std::string feat = rec && k < rec->feat.size() ? rec->feat[k] : "";
                    if (length.empty() && sch.feat_ids[k] == "length" && !feat.empty()) {
                        length = feat;
                    }
                    out << ',' << csv_escape(feat);
                }
            }
            out << ',' << csv_escape(length) << ',' << csv_escape(true_taxid_for(min_seq, truth))
                << ',' << csv_escape(extract_read_type(min_seq)) << '\n';
            ++n;
        }
    }
    std::cerr << "Exported " << csv_path.filename().string() << " (" << n << " rows)\n";
}

struct SavedTool {
    std::string filename;
    std::string tool;
    std::vector<std::string> tax_ids;
    std::vector<std::string> feat_ids;
};

fs::path tools_sidecar(const fs::path& csv) { return fs::path(csv.string() + ".tools"); }

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    bool quote = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (quote) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    cur.push_back('"');
                    ++i;
                } else {
                    quote = false;
                }
            } else {
                cur.push_back(c);
            }
        } else if (c == '"') {
            quote = true;
        } else if (c == ',') {
            fields.push_back(std::move(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    fields.push_back(std::move(cur));
    return fields;
}

std::vector<std::string> split_pipe(const std::string& text) {
    std::vector<std::string> parts;
    if (text.empty()) {
        return parts;
    }
    size_t start = 0;
    while (true) {
        size_t pos = text.find('|', start);
        if (pos == std::string::npos) {
            parts.push_back(text.substr(start));
            break;
        }
        parts.push_back(text.substr(start, pos - start));
        start = pos + 1;
    }
    return parts;
}

std::string join_pipe(const std::vector<std::string>& parts) {
    std::string out;
    for (size_t i = 0; i < parts.size(); ++i) {
        if (i) {
            out.push_back('|');
        }
        out += parts[i];
    }
    return out;
}

std::vector<SavedTool> read_sidecar(const fs::path& csv) {
    std::ifstream in(tools_sidecar(csv));
    std::vector<SavedTool> tools;
    std::string line;
    while (std::getline(in, line)) {
        line = trim_cr(line);
        if (line.empty()) {
            continue;
        }
        auto f = split_tab(line);
        if (f.size() < 4) {
            throw std::runtime_error("Bad tools sidecar for " + csv.string());
        }
        SavedTool tool;
        tool.filename = f[0];
        tool.tool = f[1];
        int n_tax = std::stoi(f[2]);
        int n_feat = std::stoi(f[3]);
        auto take_ids = [](int n, const std::string& text) {
            if (n <= 0) {
                return std::vector<std::string>{};
            }
            if (text.empty()) {
                if (n != 1) {
                    throw std::runtime_error("Tools sidecar id list is short");
                }
                return std::vector<std::string>{""};
            }
            auto parts = split_pipe(text);
            if (static_cast<int>(parts.size()) != n) {
                throw std::runtime_error("Tools sidecar id list length mismatch");
            }
            return parts;
        };
        auto tax = take_ids(n_tax, f.size() > 4 ? f[4] : "");
        auto feat = take_ids(n_feat, f.size() > 5 ? f[5] : "");
        tool.tax_ids = std::move(tax);
        tool.feat_ids = std::move(feat);
        tools.push_back(std::move(tool));
    }
    return tools;
}

void write_sidecar(const fs::path& csv, const std::vector<SavedTool>& tools) {
    fs::path dest = tools_sidecar(csv);
    fs::path tmp = fs::path(dest.string() + ".tmp");
    std::ofstream out(tmp, std::ios::trunc);
    for (const auto& tool : tools) {
        out << tool.filename << '\t' << tool.tool << '\t' << tool.tax_ids.size() << '\t'
            << tool.feat_ids.size() << '\t' << join_pipe(tool.tax_ids) << '\t' << join_pipe(tool.feat_ids)
            << '\n';
    }
    out.close();
    fs::rename(tmp, dest);
}

SavedTool tool_from_sorted(const fs::path& sorted, const std::string& filename, const std::string& tool) {
    BufferedIn in(sorted);
    if (!in.ok()) {
        throw std::runtime_error("Cannot read sorted annotator " + sorted.string());
    }
    SavedTool saved;
    saved.filename = filename;
    saved.tool = tool;
    saved.tax_ids = in.schema().tax_ids;
    saved.feat_ids = in.schema().feat_ids;
    return saved;
}

fs::path explode_tool(const fs::path& csv, const SavedTool& tool, size_t col0, const fs::path& dest) {
    std::ifstream in(csv);
    std::string header;
    if (!std::getline(in, header)) {
        throw std::runtime_error("Empty annotation CSV " + csv.string());
    }
    std::ofstream out(dest, std::ios::trunc);
    Schema sch;
    sch.tax_ids = tool.tax_ids;
    sch.feat_ids = tool.feat_ids;
    write_schema(out, sch);
    std::string line;
    while (std::getline(in, line)) {
        line = trim_cr(line);
        if (line.empty()) {
            continue;
        }
        auto fields = split_csv_line(line);
        if (fields.empty()) {
            continue;
        }
        Record rec;
        rec.seq = fields[0];
        rec.tax.assign(tool.tax_ids.size(), "0");
        rec.feat.assign(tool.feat_ids.size(), "");
        for (size_t k = 0; k < tool.tax_ids.size(); ++k) {
            size_t col = col0 + k;
            rec.tax[k] = col < fields.size() ? fields[col] : "0";
        }
        for (size_t k = 0; k < tool.feat_ids.size(); ++k) {
            size_t col = col0 + tool.tax_ids.size() + k;
            rec.feat[k] = col < fields.size() ? fields[col] : "";
        }
        write_tsv_record(out, rec);
    }
    return dest;
}

int run_append(const Options& opt) {
    fs::path out_csv(opt.out_csv);
    fs::path tmp = opt.tmp_dir.empty() ? (out_csv.parent_path() / ".combine_tmp") : fs::path(opt.tmp_dir);
    fs::create_directories(out_csv.parent_path());
    fs::create_directories(tmp);
    std::string tool_name = match_tool(opt.add_name);
    if (tool_name.empty()) {
        std::cerr << "Not an annotator output: " << opt.add_name << "\n";
        return 1;
    }
    int rc = 0;
    try {
        fs::path sorted_new = sort_annotator(opt.add_file, tool_kind(tool_name), tmp, opt.chunk_rows, 0);
        SavedTool added = tool_from_sorted(sorted_new, opt.add_name, tool_name);
        std::vector<SavedTool> tools;
        std::vector<fs::path> sorted;
        if (!opt.base_csv.empty() && fs::exists(opt.base_csv)) {
            auto original = read_sidecar(opt.base_csv);
            size_t col = 1;
            size_t index = 1;
            std::vector<std::pair<SavedTool, fs::path>> exploded;
            for (const auto& tool : original) {
                size_t width = tool.tax_ids.size() + tool.feat_ids.size();
                if (tool.filename == opt.add_name) {
                    col += width;
                    continue;
                }
                fs::path part = tmp / ("base_" + std::to_string(index++) + ".tsv");
                explode_tool(opt.base_csv, tool, col, part);
                col += width;
                exploded.push_back({tool, part});
            }
            exploded.push_back({added, sorted_new});
            std::sort(exploded.begin(), exploded.end(),
                      [](const auto& a, const auto& b) { return a.first.filename < b.first.filename; });
            tools.clear();
            for (const auto& item : exploded) {
                tools.push_back(item.first);
                sorted.push_back(item.second);
            }
        } else {
            tools.push_back(added);
            sorted.push_back(sorted_new);
        }
        std::vector<std::string> names;
        names.reserve(tools.size());
        for (const auto& tool : tools) {
            names.push_back(tool.tool);
        }
        TruthIndex truth;
        const TruthIndex* truth_ptr = nullptr;
        if (!opt.truth_table.empty()) {
            if (!load_truth_table(opt.truth_table, &truth)) {
                return 1;
            }
            truth_ptr = &truth;
        }
        fs::path tmp_csv = fs::path(out_csv.string() + ".tmp");
        merge_sample(sorted, names, tmp_csv, truth_ptr);
        fs::rename(tmp_csv, out_csv);
        write_sidecar(out_csv, tools);
    } catch (const std::exception& ex) {
        std::cerr << "combine_annotations append failed: " << ex.what() << "\n";
        rc = 1;
    }
    if (!opt.keep_tmp) {
        std::error_code ec;
        fs::remove_all(tmp, ec);
    }
    return rc;
}

int run_retrue(const Options& opt) {
    fs::path csv(opt.retrue_csv);
    TruthIndex truth;
    if (!load_truth_table(opt.truth_table, &truth)) {
        return 1;
    }
    std::ifstream in(csv);
    std::string header;
    if (!std::getline(in, header)) {
        std::cerr << "Empty annotation CSV " << csv << "\n";
        return 1;
    }
    header = trim_cr(header);
    auto cols = split_csv_line(header);
    int true_col = -1;
    for (size_t i = 0; i < cols.size(); ++i) {
        if (to_lower(cols[i]) == "true") {
            true_col = static_cast<int>(i);
            break;
        }
    }
    if (true_col < 0 || cols.empty()) {
        std::cerr << "Annotation CSV has no true column: " << csv << "\n";
        return 1;
    }
    fs::path tmp = fs::path(csv.string() + ".tmp");
    std::ofstream out(tmp, std::ios::trunc);
    out << header << '\n';
    std::string line;
    while (std::getline(in, line)) {
        line = trim_cr(line);
        if (line.empty()) {
            continue;
        }
        auto fields = split_csv_line(line);
        if (static_cast<int>(fields.size()) <= true_col) {
            fields.resize(static_cast<size_t>(true_col) + 1);
        }
        fields[static_cast<size_t>(true_col)] = true_taxid_for(fields[0], &truth);
        for (size_t i = 0; i < fields.size(); ++i) {
            if (i) {
                out << ',';
            }
            out << csv_escape(fields[i]);
        }
        out << '\n';
    }
    in.close();
    out.close();
    fs::rename(tmp, csv);
    return 0;
}

void usage() {
    std::cerr
        << "Usage: samovar_combine_annotations -i DIR -o DIR [-s N] [--chunk-rows N] [--tmp DIR] [--truth-table FILE]\n"
        << "       samovar_combine_annotations --add FILE --name NAME --out CSV [--base CSV] [--truth-table FILE]\n"
        << "       samovar_combine_annotations --retrue CSV --truth-table FILE\n"
        << "  --truth-table  CAMI-style or seq\\ttaxid map used for the true column (else parse genome headers)\n"
        << "  --add          Append one annotator output to a per-sample annotation CSV\n";
}

bool parse_args(int argc, char** argv, Options* opt) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "-i" || a == "--input_dir") {
            opt->input_dir = need(a.c_str());
        } else if (a == "-o" || a == "--output_dir") {
            opt->output_dir = need(a.c_str());
        } else if (a == "-s" || a == "--split_sample_name") {
            opt->split_n = std::stoi(need(a.c_str()));
        } else if (a == "--chunk-rows") {
            opt->chunk_rows = static_cast<size_t>(std::stoull(need(a.c_str())));
            if (opt->chunk_rows < 1000) {
                opt->chunk_rows = 1000;
            }
        } else if (a == "--tmp") {
            opt->tmp_dir = need(a.c_str());
        } else if (a == "--truth-table" || a == "--truth_table" || a == "-t") {
            opt->truth_table = need(a.c_str());
        } else if (a == "--add") {
            opt->add_file = need(a.c_str());
        } else if (a == "--name") {
            opt->add_name = need(a.c_str());
        } else if (a == "--base") {
            opt->base_csv = need(a.c_str());
        } else if (a == "--out") {
            opt->out_csv = need(a.c_str());
        } else if (a == "--retrue") {
            opt->retrue_csv = need(a.c_str());
        } else if (a == "--keep-tmp") {
            opt->keep_tmp = true;
        } else if (a == "-h" || a == "--help") {
            usage();
            std::exit(0);
        } else {
            std::cerr << "Unknown argument: " << a << "\n";
            usage();
            return false;
        }
    }
    if (!opt->retrue_csv.empty()) {
        return !opt->truth_table.empty();
    }
    if (!opt->add_file.empty() || !opt->out_csv.empty() || !opt->add_name.empty()) {
        return !opt->add_file.empty() && !opt->out_csv.empty() && !opt->add_name.empty();
    }
    return !opt->input_dir.empty() && !opt->output_dir.empty();
}

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse_args(argc, argv, &opt)) {
        usage();
        return 2;
    }
    if (!opt.retrue_csv.empty()) {
        return run_retrue(opt);
    }
    if (!opt.add_file.empty()) {
        return run_append(opt);
    }
    fs::path in_dir(opt.input_dir);
    fs::path out_dir(opt.output_dir);
    if (!fs::exists(in_dir) || !fs::is_directory(in_dir)) {
        std::cerr << "Input directory not found: " << in_dir << "\n";
        return 1;
    }
    fs::create_directories(out_dir);

    fs::path tmp = opt.tmp_dir.empty() ? (out_dir / ".combine_tmp") : fs::path(opt.tmp_dir);
    fs::create_directories(tmp);

    struct FileInfo {
        fs::path path;
        std::string sample;
        std::string tool;
        std::string name;
    };
    std::vector<FileInfo> files;
    for (const auto& ent : fs::directory_iterator(in_dir)) {
        if (!ent.is_regular_file()) {
            continue;
        }
        std::string name = ent.path().filename().string();
        if (name.size() < 4 || name.compare(name.size() - 4, 4, ".out") != 0) {
            continue;
        }
        std::string tool = match_tool(name);
        if (tool.empty()) {
            continue;
        }
        files.push_back(FileInfo{ent.path(), sample_name_from(name, opt.split_n), tool, name});
    }
    std::sort(files.begin(), files.end(),
              [](const FileInfo& a, const FileInfo& b) { return a.name < b.name; });

    std::unordered_map<std::string, std::vector<FileInfo>> by_sample;
    for (const auto& f : files) {
        by_sample[f.sample].push_back(f);
    }
    if (by_sample.empty()) {
        std::cerr << "No parseable *.out files in " << in_dir << "\n";
        return 1;
    }

    std::vector<std::string> samples;
    samples.reserve(by_sample.size());
    for (const auto& kv : by_sample) {
        samples.push_back(kv.first);
    }
    std::sort(samples.begin(), samples.end());

    TruthIndex truth;
    const TruthIndex* truth_ptr = nullptr;
    if (!opt.truth_table.empty()) {
        if (!load_truth_table(opt.truth_table, &truth)) {
            return 1;
        }
        truth_ptr = &truth;
        std::cerr << "Loaded ground-truth table " << opt.truth_table << " (" << truth.exact.size()
                  << " keys)\n";
    }

    int rc = 0;
    try {
        for (const auto& sample : samples) {
            auto& group = by_sample[sample];
            std::cerr << "Merging sample " << sample << " (" << group.size() << " annotators)\n";
            fs::path sample_tmp = tmp / sample;
            fs::create_directories(sample_tmp);
            std::vector<fs::path> sorted;
            std::vector<std::string> tools;
            for (size_t i = 0; i < group.size(); ++i) {
                sorted.push_back(sort_annotator(group[i].path, tool_kind(group[i].tool), sample_tmp,
                                                opt.chunk_rows, i));
                tools.push_back(group[i].tool);
            }
            merge_sample(sorted, tools, out_dir / (sample + ".annotation.csv"), truth_ptr);
        }
        std::cerr << "True annotations extracted\n";
    } catch (const std::exception& ex) {
        std::cerr << "combine_annotations failed: " << ex.what() << "\n";
        rc = 1;
    }

    if (!opt.keep_tmp) {
        std::error_code ec;
        fs::remove_all(tmp, ec);
    }
    return rc;
}
