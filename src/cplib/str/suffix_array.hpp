#pragma once
#include <cplib/common.hpp>
namespace cplib {
namespace detail {
template<class Range> Int sa_symbol(const Range& s,Int i){if constexpr(std::is_same_v<typename Range::value_type,char>)return static_cast<unsigned char>(s[i]);else return s[i];}
// SA-IS。LMS位置表は半分、通常の作業配列は32ビットで保持する。O(n+upper)。
template<class I,class Range> std::vector<I> sa_is_impl(const Range& s,Int upper){
    Int n=s.size();if(!n)return {};if(n==1)return {0};if(n==2)return sa_symbol(s,0)<sa_symbol(s,1)?std::vector<I>{0,1}:std::vector<I>{1,0};
    std::vector<std::uint8_t> isS(n);std::vector<I> bucket(upper+2);++bucket[sa_symbol(s,n-1)+1];
    for(Int i=n-2;i>=0;--i){if(s[i]==s[i+1])isS[i]=isS[i+1];else if(sa_symbol(s,i)<sa_symbol(s,i+1))isS[i]=1;++bucket[sa_symbol(s,i)+1];}
    for(Int c=0;c<=upper;++c)bucket[c+1]+=bucket[c];
    std::vector<I> sa(n),buf(upper+1);
    auto induce=[&](const std::vector<I>& lms){
        std::fill(sa.begin(),sa.end(),I(-1));for(Int j=0;j<=upper;++j)buf[j]=bucket[j+1];
        for(Int j=Int(lms.size())-1;j>=0;--j){I d=lms[j];sa[--buf[sa_symbol(s,d)]]=d;}
        for(Int j=0;j<=upper;++j)buf[j]=bucket[j];
        sa[buf[sa_symbol(s,n-1)]++]=I(n-1);
        for(Int j=0;j<n;++j){I v=sa[j];if(v>=1 && sa_symbol(s,v-1)>=sa_symbol(s,v))sa[buf[sa_symbol(s,v-1)]++]=v-1;}
        for(Int j=0;j<=upper;++j)buf[j]=bucket[j+1];
        for(Int j=n-1;j>=0;--j){I v=sa[j];if(v>=1 && isS[v-1])sa[--buf[sa_symbol(s,v-1)]]=v-1;}
    };
    std::vector<I> lmsMap((n+1)/2),lms;lms.reserve(n/2);
    for(Int i=1;i<n;++i)if(!isS[i-1] && isS[i]){lmsMap[i>>1]=lms.size();lms.push_back(i);isS[i]=2;}
    Int m=lms.size();induce(lms);
    if(m>1){
        std::vector<I> sortedLms;sortedLms.reserve(m);for(I v:sa)if(v>0 && isS[v]==2)sortedLms.push_back(v);
        std::vector<I> recS(m);Int recUpper=0;recS[lmsMap[sortedLms[0]>>1]]=0;
        for(Int i=1;i<m;++i){
            I left=sortedLms[i-1],right=sortedLms[i];I endLeft=lmsMap[left>>1]+1<m?lms[lmsMap[left>>1]+1]:I(n),endRight=lmsMap[right>>1]+1<m?lms[lmsMap[right>>1]+1]:I(n);
            bool same=endLeft<n && endRight<n && endLeft-left==endRight-right;
            if(same){Int length=endLeft-left;for(Int d=0;d<length;++d)if(s[left+d]!=s[right+d]){same=false;break;}if(same && s[endLeft]!=s[endRight])same=false;}
            if(!same)++recUpper;
            recS[lmsMap[right>>1]]=recUpper;
        }
        if(recUpper+1<m){auto recSa=sa_is_impl<I>(recS,recUpper);for(Int i=0;i<m;++i)sortedLms[i]=lms[recSa[i]];}
        induce(sortedLms);
    }
    return sa;
}
template<class Range> std::vector<Int> sa_is(const Range& s,Int upper){if(s.size()<=std::size_t(std::numeric_limits<std::int32_t>::max())){auto sa=sa_is_impl<std::int32_t>(s,upper);return {sa.begin(),sa.end()};}return sa_is_impl<Int>(s,upper);}
template<class I,class Range> std::vector<Int> lcp_impl(const Range& s,std::span<const Int> sa){
    Int n=s.size();std::vector<I> phi(n);phi[sa[0]]=-1;for(Int i=1;i<n;++i)phi[sa[i]]=sa[i-1];Int h=0;
    for(Int i=0;i<n;++i){Int j=phi[i];if(j<0){h=0;continue;}Int limit=n-std::max(i,j);while(h<limit && s[i+h]==s[j+h])++h;phi[i]=h;if(h>0)--h;}
    std::vector<Int> result(n-1);for(Int i=0;i<n-1;++i)result[i]=phi[sa[i+1]];return result;
}
}
// 値域0..upperの整数列をSA-ISで処理する。O(n+upper)。
template<class Range> std::vector<Int> suffix_array(const Range& s,Int upper){assert(upper>=0);for(auto v:s){assert(v>=0&&v<=upper);(void)v;}return detail::sa_is(s,upper);}
// バイト列はO(n+256)、他の型は座標圧縮を含めO(n log n)。
template<class Range> std::vector<Int> suffix_array(const Range& s){
    if constexpr(std::is_same_v<typename Range::value_type,char>||std::is_same_v<typename Range::value_type,unsigned char>)return detail::sa_is(s,255);
    else {Int n=s.size();if(!n)return {};std::vector<Int> idx(n),compressed(n);std::iota(idx.begin(),idx.end(),0);std::sort(idx.begin(),idx.end(),[&](Int l,Int r){return s[l]<s[r];});Int upper=0;for(Int i=0;i<n;++i){if(i>0&&s[idx[i-1]]!=s[idx[i]])++upper;compressed[idx[i]]=upper;}return detail::sa_is(compressed,upper);}
}
// 隣接接尾辞のLCPを求める。O(n)。
template<class Range> std::vector<Int> lcp_array(const Range& s,std::span<const Int> sa){assert(sa.size()==s.size());if(s.size()<=1)return {};for(Int v:sa){assert(v>=0&&v<Int(s.size()));(void)v;}if(s.size()<=std::size_t(std::numeric_limits<std::int32_t>::max()))return detail::lcp_impl<std::int32_t>(s,sa);return detail::lcp_impl<Int>(s,sa);}
}
