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
        } else {
            s.feat_ids.push_back(feat_header_id(header[i]));
            s.feat_cols.push_back(static_cast<int>(i));
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
        for (size_t i = 2; i < nfields; ++i) {
            s.feat_ids.push_back(std::to_string(i - 2));
            s.feat_cols.push_back(static_cast<int>(i));
        }
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
    if (t.empty() || t == "NA" || t == "NaN" || t == "nan" || t == "None") {
        return "0";
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
    std::sort(recs.begin(), recs.end(), [](const Record& a, const Record& b) {
        return a.seq < b.seq;
    });
    fs::path path = dir / ("run_" + std::to_string(idx) + ".tsv");
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    std::vector<char> buf(kIoBuf);
    out.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    write_schema(out, sch);
    std::string prev;
    for (const Record& r : recs) {
        if (r.seq == prev) {
            continue;
        }
        write_tsv_record(out, r);
        prev = r.seq;
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
    std::string prev;
    while (!heap.empty()) {
        Item cur = heap.top();
        heap.pop();
        if (cur.rec.seq != prev) {
            write_tsv_record(out, cur.rec);
            prev = cur.rec.seq;
        }
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

void merge_sample(const std::vector<fs::path>& sorted, const std::vector<std::string>& tools,
                  const fs::path& csv_path, const TruthIndex* truth) {
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
        for (auto& s : streams) {
            if (!s.alive) {
                continue;
            }
            if (!any || s.rec.seq < min_seq) {
                min_seq = s.rec.seq;
                any = true;
            }
        }
        if (!any) {
            break;
        }
        std::string length;
        std::vector<std::vector<std::string>> tax_row(streams.size());
        std::vector<std::vector<std::string>> feat_row(streams.size());
        for (size_t i = 0; i < streams.size(); ++i) {
            const Schema& sch = streams[i].in.schema();
            tax_row[i].assign(sch.tax_ids.size(), "0");
            feat_row[i].assign(sch.feat_ids.size(), "");
        }
        for (size_t i = 0; i < streams.size(); ++i) {
            auto& s = streams[i];
            if (!s.alive || s.rec.seq != min_seq) {
                continue;
            }
            tax_row[i] = s.rec.tax;
            feat_row[i] = s.rec.feat;
            const Schema& sch = s.in.schema();
            for (size_t k = 0; k < sch.feat_ids.size() && k < s.rec.feat.size(); ++k) {
                if (length.empty() && sch.feat_ids[k] == "length" && !s.rec.feat[k].empty()) {
                    length = s.rec.feat[k];
                }
            }
            s.alive = s.in.next(&s.rec);
            while (s.alive && s.rec.seq == min_seq) {
                s.alive = s.in.next(&s.rec);
            }
        }
        out << csv_escape(min_seq);
        for (size_t i = 0; i < streams.size(); ++i) {
            for (const auto& t : tax_row[i]) {
                out << ',' << csv_escape(t);
            }
            for (const auto& feat : feat_row[i]) {
                out << ',' << csv_escape(feat);
            }
        }
        out << ',' << csv_escape(length) << ',' << csv_escape(true_taxid_for(min_seq, truth))
            << ',' << csv_escape(extract_read_type(min_seq)) << '\n';
        ++n;
    }
    std::cerr << "Exported " << csv_path.filename().string() << " (" << n << " rows)\n";
}

void usage() {
    std::cerr
        << "Usage: samovar_combine_annotations -i DIR -o DIR [-s N] [--chunk-rows N] [--tmp DIR] [--truth-table FILE]\n"
        << "  --truth-table  CAMI-style or seq\\ttaxid map used for the true column (else parse genome headers)\n";
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
    return !opt->input_dir.empty() && !opt->output_dir.empty();
}

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse_args(argc, argv, &opt)) {
        usage();
        return 2;
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
