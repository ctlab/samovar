// External sort-merge of annotator *.out tables into per-sample CSVs.
// Streams seq/taxID/length only so Kraken k-mer columns never sit in RAM.

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

struct Record {
    std::string seq;
    std::string taxid;
    std::string length;
};

struct Options {
    std::string input_dir;
    std::string output_dir;
    int split_n = 1;
    size_t chunk_rows = kDefaultChunk;
    std::string tmp_dir;
    bool keep_tmp = false;
};

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

enum class ToolKind { Kaiju, Kraken, Kraken2, Custom };

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

bool parse_line(const std::string& line, ToolKind kind, Record* rec) {
    if (line.empty()) {
        return false;
    }
    auto f = split_tab(line);
    rec->seq.clear();
    rec->taxid = "0";
    rec->length.clear();
    if (kind == ToolKind::Custom) {
        if (f.size() < 2) {
            return false;
        }
        rec->seq = f[0];
        rec->taxid = normalize_taxid(f[1]);
        if (f.size() >= 3) {
            rec->length = first_length_token(f[2]);
        }
        return !rec->seq.empty();
    }
    if (f.size() < 3) {
        return false;
    }
    rec->seq = f[1];
    if (kind == ToolKind::Kraken2) {
        rec->taxid = normalize_taxid(extract_kraken2_taxid(f[2]));
        if (f.size() >= 4) {
            rec->length = first_length_token(f[3]);
        }
    } else if (kind == ToolKind::Kraken) {
        rec->taxid = normalize_taxid(f[2]);
        if (f.size() >= 4) {
            rec->length = first_length_token(f[3]);
        }
    } else {
        rec->taxid = normalize_taxid(f[2]);
    }
    for (char& c : rec->seq) {
        if (c == '\t') {
            c = ' ';
        }
    }
    return !rec->seq.empty();
}

void write_tsv_record(std::ostream& out, const Record& r) {
    out << r.seq << '\t' << r.taxid << '\t' << r.length << '\n';
}

bool read_tsv_record(std::istream& in, Record* rec) {
    std::string line;
    if (!std::getline(in, line)) {
        return false;
    }
    line = trim_cr(line);
    if (line.empty()) {
        return read_tsv_record(in, rec);
    }
    auto f = split_tab(line);
    rec->seq = f.empty() ? "" : f[0];
    rec->taxid = f.size() > 1 ? normalize_taxid(f[1]) : "0";
    rec->length = f.size() > 2 ? f[2] : "";
    return true;
}

class BufferedIn {
public:
    explicit BufferedIn(const fs::path& path) : file_(path, std::ios::in | std::ios::binary) {}
    bool ok() const { return static_cast<bool>(file_); }
    bool next(Record* rec) { return read_tsv_record(file_, rec); }

private:
    std::ifstream file_;
};

fs::path write_run(const fs::path& dir, size_t idx, std::vector<Record>& recs) {
    std::sort(recs.begin(), recs.end(), [](const Record& a, const Record& b) {
        return a.seq < b.seq;
    });
    fs::path path = dir / ("run_" + std::to_string(idx) + ".tsv");
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    std::vector<char> buf(kIoBuf);
    out.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
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

fs::path merge_runs(const std::vector<fs::path>& runs, const fs::path& out_path) {
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

    std::vector<Record> chunk;
    chunk.reserve(std::min(chunk_rows, static_cast<size_t>(65536)));
    std::vector<fs::path> runs;
    size_t run_i = 0;
    std::string line;
    while (std::getline(in, line)) {
        line = trim_cr(line);
        Record rec;
        if (!parse_line(line, kind, &rec)) {
            continue;
        }
        chunk.push_back(std::move(rec));
        if (chunk.size() >= chunk_rows) {
            runs.push_back(write_run(work, run_i++, chunk));
            chunk.clear();
        }
    }
    if (!chunk.empty()) {
        runs.push_back(write_run(work, run_i++, chunk));
        chunk.clear();
        chunk.shrink_to_fit();
    }
    fs::path sorted = tmp / ("sorted_" + std::to_string(file_index) + ".tsv");
    if (runs.empty()) {
        std::ofstream(sorted, std::ios::trunc).close();
        return sorted;
    }
    if (runs.size() == 1) {
        fs::rename(runs[0], sorted);
        return sorted;
    }
    merge_runs(runs, sorted);
    return sorted;
}

struct JoinStream {
    BufferedIn in;
    Record rec;
    bool alive = false;
    explicit JoinStream(const fs::path& path) : in(path) {
        alive = in.next(&rec);
    }
};

void merge_sample(const std::vector<fs::path>& sorted, const std::vector<std::string>& tools,
                  const fs::path& csv_path) {
    std::vector<JoinStream> streams;
    streams.reserve(sorted.size());
    for (const auto& p : sorted) {
        streams.emplace_back(p);
    }
    std::ofstream out(csv_path, std::ios::out | std::ios::binary | std::ios::trunc);
    std::vector<char> buf(kIoBuf);
    out.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));

    out << "seq";
    for (size_t i = 0; i < tools.size(); ++i) {
        out << ",taxID_" << tools[i] << '_' << i;
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
        std::vector<std::string> tax(tools.size(), "0");
        std::string length;
        for (size_t i = 0; i < streams.size(); ++i) {
            auto& s = streams[i];
            if (!s.alive || s.rec.seq != min_seq) {
                continue;
            }
            tax[i] = normalize_taxid(s.rec.taxid);
            if (length.empty() && !s.rec.length.empty()) {
                length = s.rec.length;
            }
            s.alive = s.in.next(&s.rec);
            while (s.alive && s.rec.seq == min_seq) {
                s.alive = s.in.next(&s.rec);
            }
        }
        out << csv_escape(min_seq);
        for (const auto& t : tax) {
            out << ',' << csv_escape(t);
        }
        out << ',' << csv_escape(length) << ',' << csv_escape(extract_true_taxid(min_seq))
            << ',' << csv_escape(extract_read_type(min_seq)) << '\n';
        ++n;
    }
    std::cerr << "Exported " << csv_path.filename().string() << " (" << n << " rows)\n";
}

void usage() {
    std::cerr
        << "Usage: samovar_combine_annotations -i DIR -o DIR [-s N] [--chunk-rows N] [--tmp DIR]\n";
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
            merge_sample(sorted, tools, out_dir / (sample + ".annotation.csv"));
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
