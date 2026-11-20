#pragma once

#include <bits/stdc++.h>
#include <cplib/tmpl/fastio.hpp>

namespace cplib {
using namespace std;

// ★★★★★ いわゆるQCFium、おまじない的につけとくと速い
#ifndef LOCAL_TEST
#pragma GCC target("avx")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
#endif // LOCAL_TEST

// ★★★★★ 型名を短くする、ちょっと嬉しい
using ll = long long;
using i128 = __int128_t;
using str = string;
using pii = pair<int, int>;
using pll = pair<ll, ll>;
using vi = vector<int>;
using vvi = vector<vi>;
using vvvi = vector<vvi>;
using vl = vector<ll>;
using vvl = vector<vl>;
using vvvl = vector<vvl>;
using vb = vector<bool>;
using vvb = vector<vb>;
using vvvb = vector<vvb>;
using vc = vector<char>;
using vvc = vector<vc>;
using vvvc = vector<vvc>;
using vd = vector<double>;
using vvd = vector<vd>;
using vvvd = vector<vvd>;
using vs = vector<string>;
using vvs = vector<vector<string>>;
using vvvs = vector<vector<vector<string>>>;

// ★★★ 多次元vector初期化用の関数、ちょっと癖があるけど短く書ける
// テンプレなし：vector dp(n+1, vector(m+1, vector<ll>(k+1, 0)));
// テンプレあり：auto dp = vvv<ll>(n+1, m+1, k+1, 0);
// 指定値でh×wの二次元vectorを作る。初期化する一時配列を含め O((h+1)(w+1))。
template <typename T> vector<vector<T>> vv(int h, int w, T val = T()) {
    return vector(h, vector<T>(w, val));
}

// 指定値でh1×h2×h3の三次元vectorを作る。一時配列を含め O((h1+1)(h2+1)(h3+1))。
template <typename T> vector<vector<vector<T>>> vvv(int h1, int h2, int h3, T val = T()) {
    return vector(h1, vector(h2, vector<T>(h3, val)));
}

// 指定値で四次元vectorを作る。一時配列を含め O((h1+1)(h2+1)(h3+1)(h4+1))。
template <typename T>
vector<vector<vector<vector<T>>>> vvvv(int h1, int h2, int h3, int h4, T val = T()) {
    return vector(h1, vector(h2, vector(h3, vector<T>(h4, val))));
}

// ★★ いわゆるheapq、C++のデフォルトpriority_queueは大きい順（Python、Nimとは逆）に注意
template <class T> using priority_queue_min = priority_queue<T, vector<T>, greater<T>>;

// ★ 定数系、お好みで
constexpr double PI = 3.14159265358979323;

inline float EPS = 1e-8;
inline double EPSL = 1e-10;

// ★★★★ 入出力高速化、おまじない
struct Nyan {
    // 標準入出力の設定を変更する。O(1) とストリーム設定の費用。
    Nyan() {
        cin.tie(nullptr);
        ios::sync_with_stdio(false);
        cout << fixed << setprecision(18);
    }
};

inline Nyan nyan;

// マクロ類
// わりと全部使う
// all: sort(all(a)) とか lower_bount(all(a), x) とか
// rep: オーバーロードしているので、引数の個数で挙動が変わる、基本Pythonのrangeに近い感覚で使えるようにしてるはず
//      rep(5) -> 5回ループ（新たにループ変数は作らない）
//      rep(i, 5) -> i=0,1,...,4
//      rep(i, 5) -> i=0,1,...,4
//      rep(i, 1, 6) -> i=1,...,4
//      rep(i, 1, 6, 2) -> i=1,3,5
//      rep(i, 10, -1, -1) -> i=10,9,.., 0
// smod, sdiv: python-like mod, div
// uniq: 重複削除、ソートされるのに注意
// vl a = {1, 3, 2, 5, 2, 3}; uniq(a); // a = {1, 2, 3, 5}
#define all(a) (a).begin(), (a).end()
#define len(x) ((ll)(x).size())
#define sz(x) ((ll)(x).size())
#define rep1(n) for (ll dummy_iter = 0LL; dummy_iter < n; ++dummy_iter) // 0,1,...,n-1
#define rep2(i, n)                                                                                 \
    for (ll i = 0LL, i##_counter = 0LL; i##_counter < ll(n);                                       \
         ++(i##_counter), (i) = i##_counter) // i=0,1,...,n-1
#define rep3(i, s, t)                                                                              \
    for (ll i = ll(s), i##_counter = ll(s); i##_counter < ll(t);                                   \
         ++(i##_counter), (i) = (i##_counter)) // i=s,s+1,...,t-1
#define rep4(i, s, t, step)                                                                        \
    for (ll i##_counter = step > 0 ? ll(s) : -ll(s), i##_end = step > 0 ? ll(t) : -ll(t),          \
            i##_step = abs(step), i = ll(s);                                                       \
         i##_counter < i##_end;                                                                    \
         i##_counter += i##_step, i = step > 0 ? i##_counter : -i##_counter) // i=s,s+step,...,<t
#define overload4(a, b, c, d, e, ...) e
#define rep(...) overload4(__VA_ARGS__, rep4, rep3, rep2, rep1)(__VA_ARGS__)
#define repe(a, v) for (auto &&a : (v))      // iterate over all elements in v
#define repc(a, v) for (const auto &a : (v)) // iterate over all elements in v
#define smod(n, m) ((((n) % (m)) + (m)) % (m))
#define sdiv(n, m) (((n) - smod(n, m)) / (m))
#define uniq(a)                                                                                    \
    {                                                                                              \
        sort(all(a));                                                                              \
        (a).erase(unique(all(a)), (a).end());                                                      \
    }

// 整数範囲に対する写像・絞り込みを公開する。時間は範囲長と各コールバックの費用に比例する。
using detail::template_support::mapIt;
using detail::template_support::filterIt;

// 真ならYes、偽ならNoを出力する。計算量 O(1) と出力費用。
inline void yes(bool value = true) {
    print(value ? "Yes" : "No");
}

// ★★ Yes, No なくても困らない
// 真ならYes、偽ならNoを出力し0を返す。O(1) と出力費用。
inline int Yes(bool b = true) {
    cout << (b ? "Yes\n" : "No\n");
    return 0;
};

// 真ならYES、偽ならNOを出力し0を返す。O(1) と出力費用。
inline int YES(bool b = true) {
    cout << (b ? "YES\n" : "NO\n");
    return 0;
};

// 真ならNo、偽ならYesを出力し0を返す。O(1) と出力費用。
inline int No(bool b = true) {
    return Yes(!b);
};

// 真ならNO、偽ならYESを出力し0を返す。O(1) と出力費用。
inline int NO(bool b = true) {
    return YES(!b);
};

// ★★★★ max, min, sum のvector向けオーバーロード、デフォルトがちょっと使いにくいので
// 空でない配列の最大値を返す。O(N)。
template <typename T, size_t N> T max(array<T, N> &a) {
    return *max_element(all(a));
};

// 空でない配列の最小値を返す。O(N)。
template <typename T, size_t N> T min(array<T, N> &a) {
    return *min_element(all(a));
};

// 空でないvectorの最大値を返す。要素数Nとして O(N)。
template <typename T> T max(vector<T> &a) {
    return *max_element(all(a));
};

// 空でないvectorの最小値を返す。要素数Nとして O(N)。
template <typename T> T min(vector<T> &a) {
    return *min_element(all(a));
};

// 半開区間[l,r)をコピーする。O(r-l+1)。
template <typename T> vector<T> vec_slice(const vector<T> &a, int l, int r) {
    vector<T> rev;
    rep(i, l, r) rev.push_back(a[i]);
    return rev;
};

// 半開区間[l,r)をコピーする。O(r-l+1)。
template <typename T> vector<T> vec_slice(const vector<T> &&a, int l, int r) {
    vector<T> rev;
    rep(i, l, r) rev.push_back(a[i]);
    return rev;
};

// 要素の和にzeroを加えて返す。要素数Nとして O(N)。
template <typename T> T sum(vector<T> &a, T zero = T(0)) {
    T rev = zero;
    rep(i, sz(a)) rev += a[i];
    return rev;
};

// 要素の和にzeroを加えて返す。要素数Nとして O(N)。
template <typename T> T sum(vector<T> &&a, T zero = T(0)) {
    T rev = zero;
    rep(i, sz(a)) rev += a[i];
    return rev;
};

// ★★★ vector の各要素を1増やす／減らす、グラフ系の入力受けでちょっと嬉しい
// vector<ll> a = {1, 2, 4}; ++a; // a = {2, 3, 5}
// 全要素を1減らす。要素数Nとして O(N)。
template <class T> inline vector<T> &operator--(vector<T> &v) {
    repe(x, v)-- x;
    return v;
}

// 全要素を1増やす。要素数Nとして O(N)。
template <class T> inline vector<T> &operator++(vector<T> &v) {
    repe(x, v)++ x;
    return v;
}

// ★★★★ 整数pow/sqrt
// 整数aの非負整数n乗を返す。O(log(n+1))。
inline ll powll(ll a, ll n) {
    ll res = 1;
    while (n > 0) {
        if (n & 1)
            res = (res * a);
        if (n > 1)
            a = (a * a);
        n >>= 1;
    }
    return res;
}

// aの非負整数n乗を法modで求める。O(log(n+1))。
inline ll powm32(ll a, ll n, int mod) {
    ll res = 1;
    while (n > 0) {
        if (n & 1)
            res = (res * a) % mod;
        if (n > 1)
            a = (a * a) % mod;
        n >>= 1;
    }
    return res % mod;
}

// aの非負整数n乗を法modで求める。O(log(n+1))。
inline ll powm64(i128 a, i128 n, ll mod) {
    i128 res = 1;
    while (n > 0) {
        if (n & 1)
            res = (res * a) % mod;
        if (n > 1)
            a = (a * a) % mod;
        n >>= 1;
    }
    return res % mod;
}

// 非負整数xの平方根を切り捨てる。平方根の近似後の補正回数をKとして O(K+1)。
inline ll sqrtll(ll x) {
    assert(x >= 0);
    ll rev = sqrt(x);
    while (rev * rev > x)
        --rev;
    while ((rev + 1) * (rev + 1) <= x)
        ++rev;
    return rev;
}

// ★★★★ chmax, chmin
using detail::template_support::chmax;
using detail::template_support::chmin;

// ★★★★★ vector を直接cinできるようにする
// map, set, multiset とかを直接 cout できるようにする
template <class T, class U> inline istream &operator>>(istream &is, pair<T, U> &p);
template <class T> inline istream &operator>>(istream &is, vector<T> &v);
template <class T, class U> inline ostream &operator<<(ostream &os, const pair<T, U> &p);
template <class T> inline ostream &operator<<(ostream &os, const vector<T> &v);
template <typename T, typename S> ostream &operator<<(ostream &os, const map<T, S> &mp);
template <typename T> ostream &operator<<(ostream &os, const set<T> &st);
template <typename T> ostream &operator<<(ostream &os, const multiset<T> &st);
template <typename T> ostream &operator<<(ostream &os, const unordered_set<T> &st);
template <typename T> ostream &operator<<(ostream &os, deque<T> q);
template <class T, class Container, class Compare>
ostream &operator<<(ostream &os, priority_queue<T, Container, Compare> pq);

// 入出力演算子は要素を順に読み書きする。pairは2要素、それ以外はN要素の入出力費用と O(N)。
// dequeはコピーして走査し、priority_queueはコピーと取り出しに O(N log(N+1)) を要する。
// overload operators
// pairの2要素を順に読み取る。計算量は各要素の入力費用の合計。
template <class T, class U> inline istream &operator>>(istream &is, pair<T, U> &p) {
    is >> p.first >> p.second;
    return is;
}

// vectorの全要素を順に読み取る。要素数 N に対し O(N) と各要素の入力費用。
template <class T> inline istream &operator>>(istream &is, vector<T> &v) {
    repe(x, v) is >> x;
    return is;
}

// pairの2要素を空白区切りで出力する。計算量は各要素の出力費用の合計。
template <class T, class U> inline ostream &operator<<(ostream &os, const pair<T, U> &p) {
    os << p.first << " " << p.second;
    return os;
}

// vectorの全要素を空白区切りで出力する。要素数 N に対し O(N) と各要素の出力費用。
template <class T> inline ostream &operator<<(ostream &os, const vector<T> &v) {
    rep(i, sz(v)) {
        os << v.at(i);
        if (i != sz(v) - 1)
            os << " ";
    }
    return os;
}

// mapのキーと値を順に出力する。要素数 N に対し O(N) と各要素の出力費用。
template <typename T, typename S> ostream &operator<<(ostream &os, const map<T, S> &mp) {
    for (auto &[key, val] : mp) {
        os << key << ":" << val << " ";
    }
    return os;
}

// setの要素を順に出力する。要素数 N に対し O(N) と各要素の出力費用。
template <typename T> ostream &operator<<(ostream &os, const set<T> &st) {
    auto itr = st.begin();
    for (int i = 0; i < (int)st.size(); i++) {
        os << *itr << (i + 1 != (int)st.size() ? " " : "");
        itr++;
    }
    return os;
}

// multisetの要素を順に出力する。要素数 N に対し O(N) と各要素の出力費用。
template <typename T> ostream &operator<<(ostream &os, const multiset<T> &st) {
    auto itr = st.begin();
    for (int i = 0; i < (int)st.size(); i++) {
        os << *itr << (i + 1 != (int)st.size() ? " " : "");
        itr++;
    }
    return os;
}

// unordered_setの要素を列挙して出力する。要素数 N に対し O(N) と各要素の出力費用。
template <typename T> ostream &operator<<(ostream &os, const unordered_set<T> &st) {
    ll cnt = 0;
    for (auto &e : st) {
        os << e << (++cnt != (int)st.size() ? " " : "");
    }
    return os;
}

// dequeをコピーし先頭から出力する。要素数 N に対し O(N) とコピー・出力費用。
template <typename T> ostream &operator<<(ostream &os, deque<T> q) {
    while (q.size()) {
        os << q.front();
        q.pop_front();
        if (q.size())
            os << " ";
    }
    return os;
}

// priority_queueをコピーして優先順に出力する。要素数 N に対し O(N log(N+1)) とコピー・出力費用。
template <class T, class Container, class Compare>
ostream &operator<<(ostream &os, priority_queue<T, Container, Compare> pq) {
    while (pq.size()) {
        os << pq.top() << " ";
        pq.pop();
    }
    return os;
}

#define dout(x) cout << fixed << setprecision(10) << x << endl

#define read1(a) cin >> a;
#define read2(a, b) cin >> a >> b;
#define read3(a, b, c) cin >> a >> b >> c;
#define read4(a, b, c, d) cin >> a >> b >> c >> d;
#define read5(a, b, c, d, e) cin >> a >> b >> c >> d >> e;
#define read6(a, b, c, d, e, f) cin >> a >> b >> c >> d >> e >> f;
#define read7(a, b, c, d, e, f, g) cin >> a >> b >> c >> d >> e >> f >> g;
#define read8(a, b, c, d, e, f, g, h) cin >> a >> b >> c >> d >> e >> f >> g >> h;

#define overload_read(a, b, c, d, e, f, g, h, i, ...) i
#define read(...)                                                                                  \
    overload_read(__VA_ARGS__, read8, read7, read6, read5, read4, read3, read2, read1)(__VA_ARGS__)

#ifdef LOCAL_TEST
#define inner_output1(a)                                                                           \
    cout << a << endl;                                                                             \
    cerr << "[OUTPUT #" << __LINE__ << "] " << a << endl;
#define inner_output2(a, b)                                                                        \
    cout << a << " " << b << endl;                                                                 \
    cerr << "[OUTPUT #" << __LINE__ << "] " << a << " " << b << endl;
#define inner_output3(a, b, c)                                                                     \
    cout << a << " " << b << " " << c << endl;                                                     \
    cerr << "[OUTPUT #" << __LINE__ << "] " << a << " " << b << " " << c << endl;
#define inner_output4(a, b, c, d)                                                                  \
    cout << a << " " << b << " " << c << " " << d << endl;                                         \
    cerr << "[OUTPUT #" << __LINE__ << "] " << a << " " << b << " " << c << " " << d << endl;
#define inner_output5(a, b, c, d, e)                                                               \
    cout << a << " " << b << " " << c << " " << d << " " << e << endl;                             \
    cerr << "[OUTPUT #" << __LINE__ << "] " << a << " " << b << " " << c << " " << d << " " << e   \
         << endl;
#define inner_output6(a, b, c, d, e, f)                                                            \
    cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << endl;                 \
    cerr << "[OUTPUT #" << __LINE__ << "] " << a << " " << b << " " << c << " " << d << " " << e   \
         << " " << f << endl;
#define inner_output7(a, b, c, d, e, f, g)                                                         \
    cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << endl;     \
    cerr << "[OUTPUT #" << __LINE__ << "] " << a << " " << b << " " << c << " " << d << " " << e   \
         << " " << f << " " << g << endl;
#define inner_output8(a, b, c, d, e, f, g, h)                                                      \
    cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h  \
         << endl;                                                                                  \
    cerr << "[OUTPUT #" << __LINE__ << "] " << a << " " << b << " " << c << " " << d << " " << e   \
         << " " << f << " " << g << " " << h << endl;

#else

#define inner_output1(a) cout << a << endl;
#define inner_output2(a, b) cout << a << " " << b << endl;
#define inner_output3(a, b, c) cout << a << " " << b << " " << c << endl;
#define inner_output4(a, b, c, d) cout << a << " " << b << " " << c << " " << d << endl;
#define inner_output5(a, b, c, d, e)                                                               \
    cout << a << " " << b << " " << c << " " << d << " " << e << endl;
#define inner_output6(a, b, c, d, e, f)                                                            \
    cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << endl;
#define inner_output7(a, b, c, d, e, f, g)                                                         \
    cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << endl;
#define inner_output8(a, b, c, d, e, f, g, h)                                                      \
    cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h  \
         << endl;

#endif
#define overload_inner_output(a, b, c, d, e, f, g, h, i, ...) i
#define out(...)                                                                                   \
    overload_inner_output(__VA_ARGS__, inner_output8, inner_output7, inner_output6, inner_output5, \
                          inner_output4, inner_output3, inner_output2, inner_output1)(__VA_ARGS__)

#define ii(...)                                                                                    \
    ll __VA_ARGS__;                                                                                \
    read(__VA_ARGS__)
#define si(...)                                                                                    \
    string __VA_ARGS__;                                                                            \
    read(__VA_ARGS__)
#define ci(...)                                                                                    \
    char __VA_ARGS__;                                                                              \
    read(__VA_ARGS__)
#define di(...)                                                                                    \
    double __VA_ARGS__;                                                                            \
    read(__VA_ARGS__)
#define li(name, size)                                                                             \
    ;                                                                                              \
    vector<ll> name(size);                                                                         \
    read(name)
#define lli(name, H, W)                                                                            \
    ;                                                                                              \
    vector name(H, vector<ll>(W));                                                                 \
    rep(i, H) cin >> name[i];

#ifdef LOCAL_TEST
#define inner_debug1(a) cerr << "[DEBUG#" << __LINE__ << "] " << #a << " = " << a << endl;
#define inner_debug2(a, b)                                                                         \
    cerr << "[DEBUG#" << __LINE__ << "] " << #a << " = " << a << ", " << #b << " = " << b << endl;
#define inner_debug3(a, b, c)                                                                      \
    cerr << "[DEBUG#" << __LINE__ << "] " << #a << " = " << a << ", " << #b << " = " << b << ", "  \
         << #c << " = " << c << endl;
#define inner_debug4(a, b, c, d)                                                                   \
    cerr << "[DEBUG#" << __LINE__ << "] " << #a << " = " << a << ", " << #b << " = " << b << ", "  \
         << #c << " = " << c << ", " << #d << " = " << d << endl;
#define inner_debug5(a, b, c, d, e)                                                                \
    cerr << "[DEBUG#" << __LINE__ << "] " << #a << " = " << a << ", " << #b << " = " << b << ", "  \
         << #c << " = " << c << ", " << #d << " = " << d << ", " << #e << " = " << e << endl;
#define inner_debug6(a, b, c, d, e, f)                                                             \
    cerr << "[DEBUG#" << __LINE__ << "] " << #a << " = " << a << ", " << #b << " = " << b << ", "  \
         << #c << " = " << c << ", " << #d << " = " << d << ", " << #e << " = " << e << ", " << #f \
         << " = " << f << endl;
#define inner_debug7(a, b, c, d, e, f, g)                                                          \
    cerr << "[DEBUG#" << __LINE__ << "] " << #a << " = " << a << ", " << #b << " = " << b << ", "  \
         << #c << " = " << c << ", " << #d << " = " << d << ", " << #e << " = " << e << ", " << #f \
         << " = " << f << ", " << #g << " = " << g << endl;
#define inner_debug8(a, b, c, d, e, f, g, h)                                                       \
    cerr << "[DEBUG#" << __LINE__ << "] " << #a << " = " << a << ", " << #b << " = " << b << ", "  \
         << #c << " = " << c << ", " << #d << " = " << d << ", " << #e << " = " << e << ", " << #f \
         << " = " << f << ", " << #g << " = " << g << ", " << #h << " = " << h << endl;

#define overload_inner_debug(a, b, c, d, e, f, g, h, i, ...) i
#define debug(...)                                                                                 \
    overload_inner_debug(__VA_ARGS__, inner_debug8, inner_debug7, inner_debug6, inner_debug5,      \
                         inner_debug4, inner_debug3, inner_debug2, inner_debug1)(__VA_ARGS__)

#else
#define debug(...) ;
#endif // LOCAL_TEST

// 非零整数の下位から連続する0ビット数を返す。O(1)。
inline ll ctz(ll x) {
    return __builtin_ctzll(x);
}

// 非零整数の上位から連続する0ビット数を返す。O(1)。
inline ll clz(ll x) {
    return __builtin_clzll(x);
}

// 立っているビット数を返す。O(1)。
inline ll popcount(ll x) {
    return __builtin_popcountll(x);
}

// xが半開区間[a,b)に入るか判定する。O(1)。
inline bool inrange(ll x, ll a, ll b) {
    return a <= x && x < b;
}

// vectorで最初に一致する添字を返す。なければ-1。要素数Nとして O(N)。
template <typename T> inline ll findll(vector<T> &v, T x) {
    auto tmp = find(all(v), x);
    if (tmp == v.end()) {
        return -1;
    } else {
        return distance(v.begin(), tmp);
    }
}

// 文字列で最初に一致する添字を返す。なければ-1。長さNとして O(N)。
inline ll findll(string &s, char x) {
    auto tmp = find(all(s), x);
    if (tmp == s.end()) {
        return -1;
    } else {
        return distance(s.begin(), tmp);
    }
}

// xが非負、yが正のときx/yの切り上げを返す。O(1)。
inline ll ceildiv(ll x, ll y) {
    return (x + y - 1) / y;
}

#define allit(a, pred)                                                                             \
    [&] {                                                                                          \
        repc(it, a) {                                                                              \
            if (!(pred))                                                                           \
                return false;                                                                      \
        }                                                                                          \
        return true;                                                                               \
    }()
#define anyit(a, pred)                                                                             \
    [&] {                                                                                          \
        repc(it, a) {                                                                              \
            if ((pred))                                                                            \
                return true;                                                                       \
        }                                                                                          \
        return false;                                                                              \
    }()

#define mapit(a, pred)                                                                             \
    ([&]() {                                                                                       \
        decltype(a)::value_type it;                                                                \
        vector<decltype(pred)> result_mapit(a.size());                                             \
        rep(idx, a.size()) {                                                                       \
            decltype(a)::value_type &it = a[idx];                                                  \
            result_mapit[idx] = pred;                                                              \
        }                                                                                          \
        return result_mapit;                                                                       \
    })()

#define filterit(a, pred)                                                                          \
    ([&]() {                                                                                       \
        decltype(a) result_filterit;                                                               \
        rep(idx, a.size()) {                                                                       \
            decltype(a)::value_type &it = a[idx];                                                  \
            if (pred) {                                                                            \
                result_filterit.push_back(it);                                                     \
            }                                                                                      \
        }                                                                                          \
        return result_filterit;                                                                    \
    })()

#define applyit(a, pred)                                                                           \
    {rep(idx, a.size()){decltype(a)::value_type &it = a[idx];                                      \
    a[idx] = pred;                                                                                 \
    }                                                                                              \
    }

#define countit(a, pred)                                                                           \
    ([&]() {                                                                                       \
        ll result_countit = 0;                                                                     \
        rep(idx, a.size()) {                                                                       \
            decltype(a)::value_type &it = a[idx];                                                  \
            if (pred) {                                                                            \
                result_countit++;                                                                  \
            }                                                                                      \
        }                                                                                          \
        return result_countit;                                                                     \
    })()

#define sortByIt(a, pred)                                                                          \
    {                                                                                              \
        sort(all(a), [&](const decltype(a)::value_type &left_value,                                \
                         const decltype(a)::value_type &right_value) {                             \
            auto it = left_value;                                                                  \
            auto x_value = pred;                                                                   \
            it = right_value;                                                                      \
            auto y_value = pred;                                                                   \
            return x_value < y_value;                                                              \
        });                                                                                        \
    }

// 最小要素の最初の添字を返す。空なら0。要素数Nとして O(N+1)。
inline ll minIndex(vector<ll> &a) {
    ll minIndex = 0;
    rep(i, 1, sz(a)) {
        if (a[i] < a[minIndex]) {
            minIndex = i;
        }
    }
    return minIndex;
}

// 最大要素の最初の添字を返す。空なら0。要素数Nとして O(N+1)。
inline ll maxIndex(vector<ll> &a) {
    ll maxIndex = 0;
    rep(i, 1, sz(a)) {
        if (a[i] > a[maxIndex]) {
            maxIndex = i;
        }
    }
    return maxIndex;
}

// N要素のコピーを昇順にソートする。O(N log(N+1))。
template <typename T> vector<T> sorted(vector<T> X) {
    sort(all(X));
    return X;
}

// 文字cで文字列を分割する。長さNとして O(N+1)。
inline vector<string> split(const string &s, char c) {
    vector<string> res;
    res.push_back("");
    repc(a, s) {
        if (a == c) {
            res.push_back("");
        } else {
            res.back() += a;
        }
    }
    return res;
}

} // namespace cplib
