#include <cplib/geometry/base.hpp>
#include <cplib/geometry/angle.hpp>
#include <cplib/geometry/argsort.hpp>
#include <cplib/geometry/distance.hpp>
#include <cplib/geometry/polygon.hpp>
#include <cplib/geometry/projection.hpp>
#include <random>
#include <set>
#include <unordered_set>
using namespace cplib;
int main(){
    using P=Point<Int>;P o{0,0},x{2,0};auto l=initLine(o,x);assert(ccw(l,P{0,1})==COUNTER_CLOCKWISE&&ccw(l,P{0,-1})==CLOCKWISE&&ccw(l,P{-1,0})==ONLINE_BACK&&ccw(l,P{3,0})==ONLINE_FRONT&&ccw(l,o)==ON_SEGMENT&&ccw(l,o,true)==ONLINE_BACK&&ccw(l,x,true)==ONLINE_FRONT);assert(dot(P{2,3},P{4,5})==23&&cross(P{2,3},P{4,5})==-2);assert(manhattan(P{-2,3})==5);assert(to_string(P{2,3})=="(2, 3)");std::unordered_set<P> set;set.insert({2,3});assert(set.contains({2,3}));
    std::vector<P> oct={{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};for(Int i=0;i<8;++i)assert(angle(P{1,0},oct[i])==(i<=4?i:i-8));assert(is_parallel(P{1,1},P{-2,-2})&&is_orthogonal(P{1,1},P{1,-1}));
    std::vector<std::pair<Int,Int>> angles;for(auto p:oct)angles.emplace_back(p.x,p.y);auto shuffled=angles;std::reverse(shuffled.begin(),shuffled.end());assert(argsorted(shuffled)==angles);assert(argcmp({std::numeric_limits<Int>::max(),1},{1,std::numeric_limits<Int>::max()})<0);
    using F=Fraction<Int>;auto fl=initLine(F(1),F(1),F(-1));auto vertical=initLine(F(1),F(0),F(0));assert((cross_point(fl,vertical)==Point<F>{F(0),F(1)}));auto fp=projection(fl,Point<F>{F(0),F(0)});assert((fp==Point<F>{F(1,2),F(1,2)}));assert((reflection(fl,Point<F>{F(0),F(0)})==Point<F>{F(1),F(1)}));assert(norm(Point<F>{F(0),F(0)},fl)==F(1,2));
    auto square=initPolygon<Int>({{0,0},{4,0},{4,4},{0,4}});assert(area(square)==32&&is_convex(square)&&is_convex_ccw(square));for(Int i=-1;i<=5;++i)for(Int j=-1;j<=5;++j){assert(contains(square,P{i,j})==(0<=i&&i<=4&&0<=j&&j<=4));assert(contains(square,P{i,j},true)==(0<i&&i<4&&0<j&&j<4));}auto rational=initPolygon<F>({{F(0),F(0)},{F(1),F(0)},{F(0),F(1)}});assert(area(rational)==F(1,2));
    std::mt19937 rng(39);for(int trial=0;trial<1000;++trial){auto point=[&](){return P{Int(rng()%21)-10,Int(rng()%21)-10};};P a=point(),b=point(),c=point(),d=point();if(a==b||c==d)continue;auto s=initSegment(a,b),t=initSegment(c,d);auto orient=[](P a,P b,P c){return cross(b-a,c-a);};auto sign=[](Int x){return (x>0)-(x<0);};bool boxes=std::max(std::min(a.x,b.x),std::min(c.x,d.x))<=std::min(std::max(a.x,b.x),std::max(c.x,d.x))&&std::max(std::min(a.y,b.y),std::min(c.y,d.y))<=std::min(std::max(a.y,b.y),std::max(c.y,d.y));bool want=boxes&&sign(orient(a,b,c))*sign(orient(a,b,d))<=0&&sign(orient(c,d,a))*sign(orient(c,d,b))<=0;assert(intersect(s,t)==want);using D=Point<double>;auto dl=initLine(D{double(a.x),double(a.y)},D{double(b.x),double(b.y)});D dp{double(c.x),double(c.y)};auto proj=projection(dl,dp);assert(std::abs(cross(vector(dl),proj-dl.s))<1e-7);assert(std::abs(dot(vector(dl),dp-proj))<1e-7);assert(std::abs(distance(dp,dl)*distance(dp,dl)-norm(dp,dl))<1e-7);assert(reflection(dl,reflection(dl,dp))==dp);}
    // strictの同一線分判定は元コードの挙動を保持する。
    auto same=initSegment(P{0,0},P{4,0});assert(intersect(same,same)&&!intersect(same,same,true));assert(intersect(same,initSegment(P{1,0},P{3,0}),true));
    for(int trial=0;trial<300;++trial){std::set<std::pair<Int,Int>> unique;for(int i=0;i<30;++i)unique.emplace(Int(rng()%101)-50,Int(rng()%101)-50);std::vector<P> points;for(auto [x,y]:unique)points.push_back({x,y});auto hull=convex_hull(points);assert(is_convex_ccw(hull)&&area(hull)>0);for(auto p:points)assert(contains(hull,p));std::reverse(points.begin(),points.end());assert(convex_hull(points).v==hull.v);}
    auto flat=convex_hull(std::vector<P>{{0,0},{1,0},{2,0}});assert((flat.v==std::vector<P>{{2,0},{0,0}}));assert(convex_hull(std::vector<P>{{0,0},{1,0},{2,0}},false).v.size()==4);
    assert(geometry_eq(1.0,1.0+1e-12)&&!geometry_eq(-1e9,-1e9+0.01));
}
