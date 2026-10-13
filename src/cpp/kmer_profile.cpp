// ACGT dinucleotide (k=2) profile for the Feature contract.
// One output row per read id. R1 and R2 counts are summed.
// --kmer-max N keeps the N k-mers with the largest across-read variance
// (ties break toward the canonical AA..TT order). N is capped at 16.
//
// -t threads splits reads across an OpenMP team. Counting is independent
// per read, so the team stays busy on short reads as long as there are
// at least as many reads as threads. Gzip input is inflated with pigz
// when it is on PATH, otherwise zlib, and the count still runs in parallel.

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr int KMER_N = 16;
constexpr uint32_t COUNT_MAX = 0xffffffffu;

const char* KMER_IDS[KMER_N] = {
    "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
    "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT",
};

struct Table {
    std::vector<std::string> ids;
    std::vector<uint32_t> counts;
};

int worker_count(int requested, size_t nitems) {
    int n = requested < 1 ? 1 : requested;
    if (nitems > 0 && static_cast<size_t>(n) > nitems) {
        n = static_cast<int>(nitems);
    }
    return n < 1 ? 1 : n;
}

int base4(unsigned char c) {
    switch (c) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return -1;
    }
}

void count_seq(const char* seq, size_t n, uint32_t* counts) {
    int prev = -1;
    for (size_t i = 0; i < n; ++i) {
        int b = base4(static_cast<unsigned char>(seq[i]));
        if (b < 0) {
            prev = -1;
            continue;
        }
        if (prev >= 0) {
            uint32_t& cell = counts[prev * 4 + b];
            if (cell < COUNT_MAX) {
                ++cell;
            }
        }
        prev = b;
    }
}

void replace_all(std::string& text, const char* from, const char* to) {
    const size_t nfrom = std::strlen(from);
    const size_t nto = std::strlen(to);
    size_t pos = 0;
    while ((pos = text.find(from, pos)) != std::string::npos) {
        text.replace(pos, nfrom, to);
        pos += nto;
    }
}

std::string read_id(const char* header, size_t n) {
    while (n && (*header == ' ' || *header == '\t')) {
        ++header;
        --n;
    }
    if (n && *header == '@') {
        ++header;
        --n;
    }
    size_t len = 0;
    while (len < n && header[len] != ' ' && header[len] != '\t') {
        ++len;
    }
    std::string id(header, len);
    replace_all(id, "/1", "");
    replace_all(id, "/2", "");
    return id;
}

void add_counts(uint32_t* dst, const uint32_t* src) {
    for (int k = 0; k < KMER_N; ++k) {
        uint64_t sum = static_cast<uint64_t>(dst[k]) + src[k];
        dst[k] = sum > COUNT_MAX ? COUNT_MAX : static_cast<uint32_t>(sum);
    }
}

bool load_plain(const std::string& path, std::vector<char>& out) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return false;
    }
    in.seekg(0, std::ios::end);
    std::streamoff len = in.tellg();
    if (len < 0) {
        return false;
    }
    in.seekg(0, std::ios::beg);
    out.resize(static_cast<size_t>(len));
    if (len > 0) {
        in.read(out.data(), len);
        if (!in) {
            return false;
        }
    }
    return true;
}

bool load_zlib(const std::string& path, std::vector<char>& out) {
    gzFile handle = gzopen(path.c_str(), "rb");
    if (!handle) {
        return false;
    }
    char buf[1 << 20];
    int n = 0;
    while ((n = gzread(handle, buf, sizeof(buf))) > 0) {
        out.insert(out.end(), buf, buf + n);
    }
    int err = 0;
    gzerror(handle, &err);
    gzclose(handle);
    return n == 0 && err == Z_OK;
}

bool load_pigz(const std::string& path, int threads, std::vector<char>& out) {
    int fds[2];
    if (pipe(fds) != 0) {
        return false;
    }
    pid_t pid = fork();
    if (pid < 0) {
        close(fds[0]);
        close(fds[1]);
        return false;
    }
    if (pid == 0) {
        close(fds[0]);
        if (dup2(fds[1], STDOUT_FILENO) < 0) {
            _exit(127);
        }
        close(fds[1]);
        std::string t = std::to_string(threads < 1 ? 1 : threads);
        execlp("pigz", "pigz", "-dc", "-p", t.c_str(), "--", path.c_str(), static_cast<char*>(nullptr));
        _exit(127);
    }
    close(fds[1]);
    char buf[1 << 20];
    while (true) {
        ssize_t n = read(fds[0], buf, sizeof(buf));
        if (n == 0) {
            break;
        }
        if (n < 0) {
            if (errno == EINTR) {
                continue;
            }
            break;
        }
        out.insert(out.end(), buf, buf + n);
    }
    close(fds[0]);
    int status = 0;
    while (waitpid(pid, &status, 0) < 0) {
        if (errno != EINTR) {
            out.clear();
            return false;
        }
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        out.clear();
        return false;
    }
    return true;
}

bool is_gzip(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return false;
    }
    unsigned char magic[2] = {0, 0};
    in.read(reinterpret_cast<char*>(magic), 2);
    return in.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
}

bool load_bytes(const std::string& path, int threads, std::vector<char>& out) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (path.empty() || !fs::is_regular_file(path, ec) || fs::file_size(path, ec) == 0) {
        return true;
    }
    if (!is_gzip(path)) {
        return load_plain(path, out);
    }
    if (threads > 1 && load_pigz(path, threads, out)) {
        return true;
    }
    out.clear();
    return load_zlib(path, out);
}

std::vector<size_t> line_starts(const std::vector<char>& text, int threads) {
    const size_t n = text.size();
    std::vector<size_t> starts;
    if (n == 0) {
        return starts;
    }
    const int workers = worker_count(threads, n);
    std::vector<std::vector<size_t>> local(static_cast<size_t>(workers));
#ifdef _OPENMP
#pragma omp parallel num_threads(workers)
#endif
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        const size_t begin = (n * static_cast<size_t>(tid)) / static_cast<size_t>(workers);
        const size_t end = (n * static_cast<size_t>(tid + 1)) / static_cast<size_t>(workers);
        std::vector<size_t> mine;
        mine.reserve((end - begin) / 40 + 1);
        for (size_t i = begin; i < end; ++i) {
            if (i == 0 || text[i - 1] == '\n') {
                mine.push_back(i);
            }
        }
        local[static_cast<size_t>(tid)] = std::move(mine);
    }
    size_t total = 0;
    for (const auto& part : local) {
        total += part.size();
    }
    starts.reserve(total);
    for (const auto& part : local) {
        starts.insert(starts.end(), part.begin(), part.end());
    }
    return starts;
}

void line_span(const std::vector<char>& text, const std::vector<size_t>& starts, size_t line,
               const char*& ptr, size_t& len) {
    size_t begin = starts[line];
    size_t end = (line + 1 < starts.size()) ? starts[line + 1] : text.size();
    if (end > begin && text[end - 1] == '\n') {
        --end;
    }
    if (end > begin && text[end - 1] == '\r') {
        --end;
    }
    ptr = text.data() + begin;
    len = end - begin;
}

Table collapse(std::vector<std::string> ids, std::vector<uint32_t> counts) {
    Table out;
    const size_t n = ids.size();
    std::unordered_map<std::string, size_t> index;
    index.reserve(n * 2 + 1);
    out.ids.reserve(n);
    out.counts.reserve(n * KMER_N);
    for (size_t i = 0; i < n; ++i) {
        if (ids[i].empty()) {
            continue;
        }
        auto it = index.find(ids[i]);
        if (it == index.end()) {
            index.emplace(ids[i], out.ids.size());
            out.ids.push_back(std::move(ids[i]));
            out.counts.insert(out.counts.end(), counts.begin() + static_cast<std::ptrdiff_t>(i * KMER_N),
                              counts.begin() + static_cast<std::ptrdiff_t>((i + 1) * KMER_N));
        } else {
            add_counts(&out.counts[it->second * KMER_N], &counts[i * KMER_N]);
        }
    }
    return out;
}

Table count_fastq(const std::string& path, int threads) {
    std::vector<char> text;
    if (!load_bytes(path, threads, text)) {
        throw std::runtime_error("cannot read FASTQ: " + path);
    }
    std::vector<size_t> starts = line_starts(text, threads);
    const size_t nrec = starts.size() / 4;
    std::vector<std::string> ids(nrec);
    std::vector<uint32_t> counts(nrec * KMER_N, 0);
    const int workers = worker_count(threads, nrec);
    const long long nrec_ll = static_cast<long long>(nrec);
#ifdef _OPENMP
#pragma omp parallel for num_threads(workers) schedule(static)
#endif
    for (long long r = 0; r < nrec_ll; ++r) {
        const char* header = nullptr;
        const char* seq = nullptr;
        size_t hlen = 0;
        size_t slen = 0;
        const size_t rec = static_cast<size_t>(r);
        line_span(text, starts, rec * 4, header, hlen);
        line_span(text, starts, rec * 4 + 1, seq, slen);
        ids[rec] = read_id(header, hlen);
        if (!ids[rec].empty()) {
            count_seq(seq, slen, &counts[rec * KMER_N]);
        }
    }
    return collapse(std::move(ids), std::move(counts));
}

void merge_into(Table& dst, const Table& src) {
    std::unordered_map<std::string, size_t> index;
    index.reserve(dst.ids.size() * 2 + 1);
    for (size_t i = 0; i < dst.ids.size(); ++i) {
        index.emplace(dst.ids[i], i);
    }
    for (size_t j = 0; j < src.ids.size(); ++j) {
        auto it = index.find(src.ids[j]);
        if (it == index.end()) {
            index.emplace(src.ids[j], dst.ids.size());
            dst.ids.push_back(src.ids[j]);
            dst.counts.insert(dst.counts.end(), src.counts.begin() + static_cast<std::ptrdiff_t>(j * KMER_N),
                              src.counts.begin() + static_cast<std::ptrdiff_t>((j + 1) * KMER_N));
        } else {
            add_counts(&dst.counts[it->second * KMER_N], &src.counts[j * KMER_N]);
        }
    }
}

struct Rank {
    unsigned __int128 disp;
    int idx;
};

std::vector<int> select_kmers(const Table& table, int kmer_max) {
    const size_t n = table.ids.size();
    std::vector<Rank> ranks(KMER_N);
    for (int k = 0; k < KMER_N; ++k) {
        unsigned __int128 sum = 0;
        unsigned __int128 sumsq = 0;
        for (size_t i = 0; i < n; ++i) {
            unsigned __int128 c = table.counts[i * KMER_N + static_cast<size_t>(k)];
            sum += c;
            sumsq += c * c;
        }
        unsigned __int128 disp = 0;
        if (n > 1) {
            disp = static_cast<unsigned __int128>(n) * sumsq - sum * sum;
        }
        ranks[static_cast<size_t>(k)] = Rank{disp, k};
    }
    int keep = KMER_N;
    if (kmer_max >= 0) {
        keep = std::min(kmer_max, KMER_N);
    }
    if (kmer_max >= 0) {
        std::sort(ranks.begin(), ranks.end(), [](const Rank& a, const Rank& b) {
            if (a.disp != b.disp) {
                return a.disp > b.disp;
            }
            return a.idx < b.idx;
        });
        ranks.resize(static_cast<size_t>(keep));
        std::sort(ranks.begin(), ranks.end(), [](const Rank& a, const Rank& b) { return a.idx < b.idx; });
    }
    std::vector<int> cols;
    cols.reserve(ranks.size());
    for (const Rank& rank : ranks) {
        cols.push_back(rank.idx);
    }
    return cols;
}

void write_table(const std::string& path, const Table& table, const std::vector<int>& cols) {
    namespace fs = std::filesystem;
    fs::path dest(path);
    std::error_code ec;
    if (!dest.parent_path().empty()) {
        fs::create_directories(dest.parent_path(), ec);
    }
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("cannot write " + path);
    }
    out << "seq";
    for (int col : cols) {
        out << '\t' << KMER_IDS[col];
    }
    out << '\n';
    for (size_t i = 0; i < table.ids.size(); ++i) {
        out << table.ids[i];
        const uint32_t* row = &table.counts[i * KMER_N];
        for (int col : cols) {
            out << '\t' << row[col];
        }
        out << '\n';
    }
    if (!out) {
        throw std::runtime_error("cannot write " + path);
    }
}

void usage() {
    std::cerr
        << "Usage: samovar_kmer -i R1 [-I R2] -o OUT [-t THREADS] [--kmer-max N]\n"
        << "  Count ACGT dinucleotides per read. --kmer-max keeps that many\n"
        << "  k-mers with the highest across-read variance (at most 16).\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::string r1;
    std::string r2;
    std::string dest;
    int threads = 1;
    int kmer_max = -1;
    bool kmer_max_set = false;
    try {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            auto need = [&](const char* flag) -> std::string {
                if (i + 1 >= argc) {
                    throw std::runtime_error(std::string("missing value for ") + flag);
                }
                return argv[++i];
            };
            if (arg == "-h" || arg == "--help") {
                usage();
                return 0;
            }
            if (arg == "-i" || arg == "-1") {
                r1 = need("-i");
            } else if (arg == "-I" || arg == "-2") {
                r2 = need("-I");
            } else if (arg == "-o" || arg == "--output") {
                dest = need("-o");
            } else if (arg == "-t" || arg == "--threads") {
                threads = std::stoi(need("-t"));
            } else if (arg == "-d") {
                (void)need("-d");
            } else if (arg == "--kmer-max") {
                kmer_max = std::stoi(need("--kmer-max"));
                kmer_max_set = true;
            } else if (arg.rfind("--kmer-max=", 0) == 0) {
                kmer_max = std::stoi(arg.substr(std::strlen("--kmer-max=")));
                kmer_max_set = true;
            } else {
                throw std::runtime_error("unknown argument: " + arg);
            }
        }
        if (r1.empty() || dest.empty()) {
            throw std::runtime_error("-i and -o are required");
        }
        if (threads < 1) {
            throw std::runtime_error("-t must be >= 1");
        }
        if (kmer_max_set && kmer_max < 0) {
            throw std::runtime_error("--kmer-max must be >= 0");
        }
        if (!kmer_max_set) {
            kmer_max = -1;
        }
    } catch (const std::exception& exc) {
        std::cerr << "samovar_kmer: " << exc.what() << "\n";
        usage();
        return 2;
    }

#ifdef _OPENMP
    omp_set_dynamic(0);
#endif

    try {
        Table table = count_fastq(r1, threads);
        if (!r2.empty()) {
            merge_into(table, count_fastq(r2, threads));
        }
        write_table(dest, table, select_kmers(table, kmer_max));
    } catch (const std::exception& exc) {
        std::cerr << "samovar_kmer: " << exc.what() << "\n";
        return 1;
    }
    return 0;
}
