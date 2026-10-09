// Pairwise Spearman rank correlation for sampled per-read columns.
// One matrix per sample. Rows with a missing value are dropped per pair.
// Ties use midranks. A constant column yields 0 off the diagonal.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

double parse_double(const std::string& token) {
    if (token.empty()) {
        return kNaN;
    }
    char* end = nullptr;
    double value = std::strtod(token.c_str(), &end);
    if (end == token.c_str()) {
        return kNaN;
    }
    return value;
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    for (char ch : line) {
        if (ch == '\r') {
            continue;
        }
        if (ch == '\t') {
            fields.push_back(cur);
            cur.clear();
            continue;
        }
        cur.push_back(ch);
    }
    fields.push_back(cur);
    return fields;
}

void average_ranks(const std::vector<double>& values, std::vector<double>& ranks) {
    const int n = static_cast<int>(values.size());
    ranks.resize(static_cast<size_t>(n));
    std::vector<int> order(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        order[static_cast<size_t>(i)] = i;
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return values[static_cast<size_t>(a)] < values[static_cast<size_t>(b)];
    });
    int i = 0;
    while (i < n) {
        int j = i;
        while (j + 1 < n &&
               values[static_cast<size_t>(order[static_cast<size_t>(j + 1)])] ==
                   values[static_cast<size_t>(order[static_cast<size_t>(i)])]) {
            ++j;
        }
        const double avg = 0.5 * (static_cast<double>(i + 1) + static_cast<double>(j + 1));
        for (int k = i; k <= j; ++k) {
            ranks[static_cast<size_t>(order[static_cast<size_t>(k)])] = avg;
        }
        i = j + 1;
    }
}

double pearson(const std::vector<double>& left, const std::vector<double>& right) {
    const int n = static_cast<int>(left.size());
    if (n < 2) {
        return 0.0;
    }
    double mean_l = 0.0;
    double mean_r = 0.0;
    for (int i = 0; i < n; ++i) {
        mean_l += left[static_cast<size_t>(i)];
        mean_r += right[static_cast<size_t>(i)];
    }
    mean_l /= static_cast<double>(n);
    mean_r /= static_cast<double>(n);
    double num = 0.0;
    double den_l = 0.0;
    double den_r = 0.0;
    for (int i = 0; i < n; ++i) {
        const double dl = left[static_cast<size_t>(i)] - mean_l;
        const double dr = right[static_cast<size_t>(i)] - mean_r;
        num += dl * dr;
        den_l += dl * dl;
        den_r += dr * dr;
    }
    if (den_l <= 0.0 || den_r <= 0.0) {
        return 0.0;
    }
    return num / std::sqrt(den_l * den_r);
}

double spearman_pair(const std::vector<double>& left, const std::vector<double>& right, int rows) {
    std::vector<double> xs;
    std::vector<double> ys;
    xs.reserve(static_cast<size_t>(rows));
    ys.reserve(static_cast<size_t>(rows));
    for (int r = 0; r < rows; ++r) {
        const double a = left[static_cast<size_t>(r)];
        const double b = right[static_cast<size_t>(r)];
        if (std::isnan(a) || std::isnan(b)) {
            continue;
        }
        xs.push_back(a);
        ys.push_back(b);
    }
    if (xs.size() < 2) {
        return 0.0;
    }
    bool vary_x = false;
    bool vary_y = false;
    for (size_t i = 1; i < xs.size(); ++i) {
        vary_x = vary_x || xs[i] != xs[0];
        vary_y = vary_y || ys[i] != ys[0];
    }
    if (!vary_x || !vary_y) {
        return 0.0;
    }
    std::vector<double> rank_x;
    std::vector<double> rank_y;
    average_ranks(xs, rank_x);
    average_ranks(ys, rank_y);
    const double rho = pearson(rank_x, rank_y);
    if (std::isnan(rho)) {
        return 0.0;
    }
    return rho;
}

struct Sample {
    std::string name;
    int rows = 0;
    std::vector<std::string> columns;
    std::vector<std::vector<double>> data;
};

std::vector<double> pairwise(const Sample& sample, bool against_last) {
    const int n = static_cast<int>(sample.columns.size());
    std::vector<double> matrix(static_cast<size_t>(n * n), 0.0);
    for (int i = 0; i < n; ++i) {
        matrix[static_cast<size_t>(i * n + i)] = 1.0;
    }
    if (n < 2 || sample.rows < 2) {
        return matrix;
    }
    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(static_cast<size_t>(n * (n - 1) / 2));
    if (against_last) {
        const int last = n - 1;
        for (int i = 0; i < last; ++i) {
            pairs.emplace_back(i, last);
        }
    } else {
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                pairs.emplace_back(i, j);
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int p = 0; p < static_cast<int>(pairs.size()); ++p) {
        const int i = pairs[static_cast<size_t>(p)].first;
        const int j = pairs[static_cast<size_t>(p)].second;
        const double rho = spearman_pair(sample.data[static_cast<size_t>(i)],
                                          sample.data[static_cast<size_t>(j)],
                                          sample.rows);
        matrix[static_cast<size_t>(i * n + j)] = rho;
        matrix[static_cast<size_t>(j * n + i)] = rho;
    }
    return matrix;
}

void write_matrix(const Sample& sample, const std::vector<double>& matrix) {
    const int n = static_cast<int>(sample.columns.size());
    std::cout << "SAMPLE\t" << sample.name << "\n";
    for (int j = 0; j < n; ++j) {
        if (j) {
            std::cout << '\t';
        }
        std::cout << sample.columns[static_cast<size_t>(j)];
    }
    std::cout << '\n';
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j) {
                std::cout << '\t';
            }
            std::cout << matrix[static_cast<size_t>(i * n + j)];
        }
        std::cout << '\n';
    }
}

bool read_sample(Sample* sample) {
    std::string line;
    if (!std::getline(std::cin, line)) {
        return false;
    }
    if (line.empty()) {
        return false;
    }
    {
        const std::vector<std::string> fields = split_tab(line);
        if (fields.size() < 3 || fields[0] != "SAMPLE") {
            throw std::runtime_error("expected SAMPLE header");
        }
        sample->name = fields[1];
        sample->rows = std::atoi(fields[2].c_str());
    }
    if (!std::getline(std::cin, line)) {
        throw std::runtime_error("missing column header");
    }
    sample->columns = split_tab(line);
    const int cols = static_cast<int>(sample->columns.size());
    sample->data.assign(static_cast<size_t>(cols), std::vector<double>(static_cast<size_t>(sample->rows), kNaN));
    for (int r = 0; r < sample->rows; ++r) {
        if (!std::getline(std::cin, line)) {
            throw std::runtime_error("truncated sample rows");
        }
        std::vector<std::string> fields = split_tab(line);
        if (static_cast<int>(fields.size()) < cols) {
            fields.resize(static_cast<size_t>(cols));
        }
        for (int c = 0; c < cols; ++c) {
            sample->data[static_cast<size_t>(c)][static_cast<size_t>(r)] = parse_double(fields[static_cast<size_t>(c)]);
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    int threads = 1;
    bool against_last = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--threads" && i + 1 < argc) {
            threads = std::max(1, std::atoi(argv[++i]));
        } else if (arg == "--against-last") {
            against_last = true;
        }
    }
#ifdef _OPENMP
    omp_set_num_threads(threads);
#else
    (void)threads;
#endif
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string line;
    if (!std::getline(std::cin, line)) {
        return 0;
    }
    int n_samples = 0;
    {
        std::istringstream header(line);
        std::string tag;
        header >> tag >> n_samples;
        if (tag != "SAMPLES") {
            std::cerr << "expected SAMPLES header\n";
            return 2;
        }
    }
    std::cout.setf(std::ios::fmtflags(0), std::ios::floatfield);
    std::cout.precision(17);
    try {
        for (int s = 0; s < n_samples; ++s) {
            Sample sample;
            if (!read_sample(&sample)) {
                std::cerr << "missing sample " << s << "\n";
                return 2;
            }
            write_matrix(sample, pairwise(sample, against_last));
        }
    } catch (const std::exception& exc) {
        std::cerr << exc.what() << "\n";
        return 2;
    }
    return 0;
}
