#include<bits/stdc++.h>
using namespace std;
using ll = long long;using vl = vector<ll>;
#define rep(i,N) for(ll i = 0;i < N;i++)
#define rep2(i,l,r) for(ll i = l;i < r;i++)
#define all(x) x.begin(),x.end()
#define len(x) ((ll)(x).size())
template <class T> inline bool chmax(T& M, const T& x) { if (M < x) { M = x; return true; } return false; }
template <class T> inline bool chmin(T& m, const T& x) { if (m > x) { m = x; return true; } return false; }
const ll INF = 1LL<<61; // 2.3e18
template <class T> using priority_queue_min = priority_queue<T, vector<T>, greater<T>>;
template <class T> void vin(vector<T>& v){rep(i,len(v)) cin >> v[i];}
template <class T> void vout(vector<T> v){rep(i,len(v)) cout << v[i] << ((i != len(v)-1) ? ' ' : '\n');}
struct Nyan { Nyan() { cin.tie(nullptr); ios::sync_with_stdio(false); cout << fixed << setprecision(18); } } nyan;