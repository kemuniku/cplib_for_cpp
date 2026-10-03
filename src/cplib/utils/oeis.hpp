#pragma once
#include <cplib/math/fractions.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <curl/curl.h>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <string_view>

namespace cplib {
// This optional network utility requires Boost headers and linking with -lcurl.
enum OEISResultKind {
    oeisSequence,
    oeisMatrixRowMajor,
    oeisMatrixAntidiagonal,
    oeisNumerator,
    oeisDenominator
};

struct OEISResult {
    std::string id, name, data, query, url;
    OEISResultKind kind = oeisSequence;
};

namespace detail::oeis {
template <class R> std::string makeQuery(const R &v) {
    std::ostringstream s;
    bool first = true;
    for (auto x : v) {
        if (!first)
            s << ',';
        first = false;
        if constexpr (std::is_integral_v<decltype(x)> && sizeof(x) == 1)
            s << +x;
        else
            s << x;
    }
    return s.str();
}

inline std::vector<OEISResult> parse(const std::string &body, const std::string &query,
                                     OEISResultKind kind) {
    boost::property_tree::ptree root;
    std::istringstream input(body);
    boost::property_tree::read_json(input, root);
    auto results = root.get_child_optional("results");
    const auto &entries = results ? *results : root;
    std::vector<OEISResult> out;
    for (auto &[key, entry] : entries) {
        if (!key.empty())
            continue;
        auto number = entry.get_optional<Int>("number");
        if (!number)
            continue;
        std::ostringstream id;
        id << 'A' << std::setw(6) << std::setfill('0') << *number;
        out.push_back({id.str(), entry.get<std::string>("name", ""),
                       entry.get<std::string>("data", ""), query, "https://oeis.org/" + id.str(),
                       kind});
    }
    return out;
}

class Client {
    struct CurlInit {
        CurlInit() {
            if (curl_global_init(CURL_GLOBAL_DEFAULT) != CURLE_OK)
                throw std::runtime_error("curl initialization failed");
        }

        ~CurlInit() {
            curl_global_cleanup();
        }
    };

    static CURL *makeHandle() {
        static CurlInit initialization;
        (void)initialization;
        return curl_easy_init();
    }

    std::unique_ptr<CURL, decltype(&curl_easy_cleanup)> handle{makeHandle(), curl_easy_cleanup};

    static std::size_t write(char *p, std::size_t size, std::size_t count, void *data) noexcept {
        try {
            static_cast<std::string *>(data)->append(p, size * count);
            return size * count;
        } catch (...) {
            return 0;
        }
    }

public:
    Client() {
        if (!handle)
            throw std::runtime_error("curl allocation failed");
    }

    std::string operator()(const std::string &query) {
        auto h = handle.get();
        std::unique_ptr<char, decltype(&curl_free)> escaped(
            curl_easy_escape(h, query.c_str(), query.size()), curl_free);
        if (!escaped)
            throw std::runtime_error("URL encoding failed");
        std::string url = "https://oeis.org/search?q=" + std::string(escaped.get()) + "&fmt=json",
                    body;
        curl_easy_setopt(h, CURLOPT_URL, url.c_str());
        curl_easy_setopt(h, CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(h, CURLOPT_TIMEOUT, 30L);
        curl_easy_setopt(h, CURLOPT_USERAGENT, "cppcplib OEIS lookup");
        curl_easy_setopt(h, CURLOPT_WRITEFUNCTION, &write);
        curl_easy_setopt(h, CURLOPT_WRITEDATA, &body);
        auto code = curl_easy_perform(h);
        if (code != CURLE_OK)
            throw std::runtime_error(curl_easy_strerror(code));
        long status = 0;
        curl_easy_getinfo(h, CURLINFO_RESPONSE_CODE, &status);
        if (status != 200)
            throw std::runtime_error("OEIS request failed: " + std::to_string(status));
        return body;
    }
};

template <class Fetch> auto search(Fetch &fetch, const std::string &query, OEISResultKind kind) {
    return query.empty() ? std::vector<OEISResult>{} : parse(fetch(query), query, kind);
}

inline void append(std::vector<OEISResult> &a, std::vector<OEISResult> b) {
    a.insert(a.end(), std::make_move_iterator(b.begin()), std::make_move_iterator(b.end()));
}
}

// Fetch(query) returns a JSON response body; injectable for offline use/testing.
// 整数列をOEISで検索する。通信・JSON解析を除く処理時間はO(n)。
template <std::integral T, class Fetch>
auto searchOEIS(const std::vector<T> &values, Fetch &&fetch) {
    return detail::oeis::search(fetch, detail::oeis::makeQuery(values), oeisSequence);
}

// 二次元整数配列を行方向と反対角線方向で検索する。
template <std::integral T, class Fetch>
auto searchOEIS(const std::vector<std::vector<T>> &matrix, Fetch &&fetch) {
    std::vector<T> rowMajor, diagonal;
    Int width = 0;
    for (auto &row : matrix) {
        width = std::max<Int>(width, row.size());
        rowMajor.insert(rowMajor.end(), row.begin(), row.end());
    }
    std::vector<OEISResult> out;
    if (rowMajor.empty())
        return out;
    for (Int d = 0; d < Int(matrix.size()) + width - 1; ++d)
        for (Int r = 0; r < Int(matrix.size()); ++r)
            if (d >= r && d - r < Int(matrix[r].size()))
                diagonal.push_back(matrix[r][d - r]);
    detail::oeis::append(
        out, detail::oeis::search(fetch, detail::oeis::makeQuery(rowMajor), oeisMatrixRowMajor));
    detail::oeis::append(out, detail::oeis::search(fetch, detail::oeis::makeQuery(diagonal),
                                                   oeisMatrixAntidiagonal));
    return out;
}

// 有理数列を既約分数に直し、分子列と分母列を別々に検索する。
template <std::integral T, class Fetch>
auto searchOEIS(const std::vector<Fraction<T>> &values, Fetch &&fetch) {
    std::vector<T> num, den;
    for (auto x : values) {
        if (!x.den)
            throw std::invalid_argument("OEIS rational search requires nonzero denominators");
        x.reduce();
        num.push_back(x.num);
        den.push_back(x.den);
    }
    std::vector<OEISResult> out;
    detail::oeis::append(out,
                         detail::oeis::search(fetch, detail::oeis::makeQuery(num), oeisNumerator));
    detail::oeis::append(
        out, detail::oeis::search(fetch, detail::oeis::makeQuery(den), oeisDenominator));
    return out;
}

template <class R> auto searchOEIS(const R &values) {
    detail::oeis::Client client;
    return searchOEIS(values, client);
}

// OEIS検索結果を項目情報が分かる文字列にする。
inline std::string to_string(const OEISResult &item) {
    static constexpr const char *names[] = {"整数列", "2次元配列（行方向の連結）",
                                            "2次元配列（反対角線方向）", "有理数列（分子）",
                                            "有理数列（分母）"};
    std::string out = item.id;
    if (!item.name.empty())
        out += ": " + item.name;
    out += "\n  検索形式: " + std::string(names[item.kind]) + "\n  URL: " + item.url;
    if (!item.data.empty()) {
        out += "\n  OEISの項: ";
        std::istringstream in(item.data);
        std::string term;
        Int count = 0;
        while (std::getline(in, term, ',')) {
            if (count == 12) {
                out += ", …";
                break;
            }
            if (count++)
                out += ", ";
            auto a = term.find_first_not_of(" \t\r\n"), b = term.find_last_not_of(" \t\r\n");
            if (a != std::string::npos)
                out += term.substr(a, b - a + 1);
        }
    }
    return out;
}

// OEIS検索結果の列を見やすい複数行の文字列にする。
inline std::string to_string(const std::vector<OEISResult> &items) {
    if (items.empty())
        return "OEIS検索結果はありませんでした。";
    std::string out = "OEIS検索結果（" + std::to_string(items.size()) + "件）";
    for (std::size_t i = 0; i < items.size(); ++i)
        out += "\n\n" + std::to_string(i + 1) + ". " + to_string(items[i]);
    return out;
}

inline std::ostream &operator<<(std::ostream &out, const OEISResult &item) {
    return out << to_string(item);
}
}
