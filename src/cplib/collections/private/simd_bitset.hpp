#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <cplib/collections/private/bitset_avx2_impl.hpp>
#include <cplib/collections/private/bitset_avx512_impl.hpp>
#include <cplib/collections/private/bitset_search_impl.hpp>
#include <cplib/collections/private/bitset_avx512_shift_assign.hpp>
#include <cplib/collections/private/bitset_avx512_fuse_shift.hpp>
#include <stdexcept>
#include <compare>
#include <string_view>
namespace cplib::detail {
struct SimdBitSetAccess;
struct BitSetKernel2 {
template<class... A> static auto avxAnd(A... a){return bitset_avx2::cplib_bs_and(a...);}
template<class... A> static auto avxAndNot(A... a){return bitset_avx2::cplib_bs_andnot(a...);}
template<class... A> static auto avxOr(A... a){return bitset_avx2::cplib_bs_or(a...);}
template<class... A> static auto avxXor(A... a){return bitset_avx2::cplib_bs_xor(a...);}
template<class... A> static auto avxNot(A... a){return bitset_avx2::cplib_bs_not(a...);}
template<class... A> static auto avxShl(A... a){return bitset_avx2::cplib_bs_shl(a...);}
template<class... A> static auto avxShr(A... a){return bitset_avx2::cplib_bs_shr(a...);}
template<class... A> static auto avxPopcount(A... a){return bitset_avx2::cplib_bs_popcount(a...);}
template<class... A> static auto avxAndPopcount(A... a){return bitset_avx2::cplib_bs_andpopcount(a...);}
template<class... A> static auto avxOrPopcount(A... a){return bitset_avx2::cplib_bs_orpopcount(a...);}
template<class... A> static auto avxXorPopcount(A... a){return bitset_avx2::cplib_bs_xorpopcount(a...);}
template<class... A> static auto avxFromBools(A... a){return bitset_avx2::cplib_bs_from_bools(a...);}
template<class... A> static auto avxSelectAssign(A... a){return bitset_avx2::cplib_bs_select(a...);}
template<class... A> static auto avxOrAndAssign(A... a){return bitset_avx2::cplib_bs_orand(a...);}
template<class... A> static auto avxAndOrAssign(A... a){return bitset_avx2::cplib_bs_andor(a...);}
template<class... A> static auto avxXorAndAssign(A... a){return bitset_avx2::cplib_bs_xorand(a...);}
template<class... A> static auto avxMajority(A... a){return bitset_avx2::cplib_bs_majority(a...);}
template<class... A> static auto avxXnorAssign(A... a){return bitset_avx2::cplib_bs_xnor(a...);}
template<class... A> static auto avxPopcountRange(A... a){return bitset_avx2::cplib_bs_popcount_range(a...);}
template<class... A> static auto avxAndpopcountRange(A... a){return bitset_avx2::cplib_bs_andpopcount_range(a...);}
template<class... A> static auto avxOrpopcountRange(A... a){return bitset_avx2::cplib_bs_orpopcount_range(a...);}
template<class... A> static auto avxXorpopcountRange(A... a){return bitset_avx2::cplib_bs_xorpopcount_range(a...);}
template<class... A> static auto avxFromStringChar(A... a){return bitset_avx2::cplib_bs_from_string_char(a...);}
template<class... A> static auto avxFromStringEqual(A... a){return bitset_avx2::cplib_bs_from_string_equal(a...);}
template<class... A> static auto avxIntersects(A... a){return bitset_avx2::cplib_bs_intersects(a...);}
template<class... A> static auto avxSubset(A... a){return bitset_avx2::cplib_bs_subset(a...);}
template<class... A> static auto avxSetRange(A... a){return bitset_avx2::cplib_bs_set_range(a...);}
template<class... A> static auto avxClearRange(A... a){return bitset_avx2::cplib_bs_clear_range(a...);}
template<class... A> static auto avxFlipRange(A... a){return bitset_avx2::cplib_bs_flip_range(a...);}
template<class... A> static auto avxCmp(A... a){return bitset_avx2::cplib_bs_cmp(a...);}
template<class... A> static auto avxAll(A... a){return bitset_avx2::cplib_bs_all(a...);}
template<class... A> static auto avxAny(A... a){return bitset_avx2::cplib_bs_any(a...);}
static std::size_t searchNext(const UInt* x,std::size_t lo,std::size_t hi){return bitset_search::cplib_bs_search_next_4(x,lo,hi);}
static std::size_t searchPrev(const UInt* x,std::size_t lo,std::size_t hi){return bitset_search::cplib_bs_search_prev_4(x,lo,hi);}
};
struct BitSetKernel512 {
template<class... A> static auto avxAnd(A... a){return bitset_avx512::cplib_bs512_and(a...);}
template<class... A> static auto avxAdd(A... a){return bitset_avx512::cplib_bs512_add(a...);}
template<class... A> static auto avxAndNot(A... a){return bitset_avx512::cplib_bs512_andnot(a...);}
template<class... A> static auto avxOr(A... a){return bitset_avx512::cplib_bs512_or(a...);}
template<class... A> static auto avxXor(A... a){return bitset_avx512::cplib_bs512_xor(a...);}
template<class... A> static auto avxNot(A... a){return bitset_avx512::cplib_bs512_not(a...);}
template<class... A> static auto avxShl(A... a){return bitset_avx512::cplib_bs512_shl(a...);}
template<class... A> static auto avxShr(A... a){return bitset_avx512::cplib_bs512_shr(a...);}
template<class... A> static auto avxPopcount(A... a){return bitset_avx512::cplib_bs512_popcount(a...);}
template<class... A> static auto avxAndPopcount(A... a){return bitset_avx512::cplib_bs512_andpopcount(a...);}
template<class... A> static auto avxOrPopcount(A... a){return bitset_avx512::cplib_bs512_orpopcount(a...);}
template<class... A> static auto avxXorPopcount(A... a){return bitset_avx512::cplib_bs512_xorpopcount(a...);}
template<class... A> static auto avxFromBools(A... a){return bitset_avx512::cplib_bs512_from_bools(a...);}
template<class... A> static auto avxSelectAssign(A... a){return bitset_avx512::cplib_bs512_select(a...);}
template<class... A> static auto avxOrAndAssign(A... a){return bitset_avx512::cplib_bs512_orand(a...);}
template<class... A> static auto avxAndOrAssign(A... a){return bitset_avx512::cplib_bs512_andor(a...);}
template<class... A> static auto avxXorAndAssign(A... a){return bitset_avx512::cplib_bs512_xorand(a...);}
template<class... A> static auto avxMajority(A... a){return bitset_avx512::cplib_bs512_majority(a...);}
template<class... A> static auto avxXnorAssign(A... a){return bitset_avx512::cplib_bs512_xnor(a...);}
template<class... A> static auto avxAndAssignPopcount(A... a){return bitset_avx512::cplib_bs512_and_update_count(a...);}
template<class... A> static auto avxOrAssignPopcount(A... a){return bitset_avx512::cplib_bs512_or_update_count(a...);}
template<class... A> static auto avxXorAssignPopcount(A... a){return bitset_avx512::cplib_bs512_xor_update_count(a...);}
template<class... A> static auto avxAndNotAssignPopcount(A... a){return bitset_avx512::cplib_bs512_andnot_update_count(a...);}
template<class... A> static auto avxSelectAssignPopcount(A... a){return bitset_avx512::cplib_bs512_select_update_count(a...);}
template<class... A> static auto avxOrAndAssignPopcount(A... a){return bitset_avx512::cplib_bs512_orand_update_count(a...);}
template<class... A> static auto avxAndOrAssignPopcount(A... a){return bitset_avx512::cplib_bs512_andor_update_count(a...);}
template<class... A> static auto avxXorAndAssignPopcount(A... a){return bitset_avx512::cplib_bs512_xorand_update_count(a...);}
template<class... A> static auto avxXnorAssignPopcount(A... a){return bitset_avx512::cplib_bs512_xnor_update_count(a...);}
template<class... A> static auto avxPopcountRange(A... a){return bitset_avx512::cplib_bs512_popcount_range(a...);}
template<class... A> static auto avxAndpopcountRange(A... a){return bitset_avx512::cplib_bs512_andpopcount_range(a...);}
template<class... A> static auto avxOrpopcountRange(A... a){return bitset_avx512::cplib_bs512_orpopcount_range(a...);}
template<class... A> static auto avxXorpopcountRange(A... a){return bitset_avx512::cplib_bs512_xorpopcount_range(a...);}
template<class... A> static auto avxFromStringChar(A... a){return bitset_avx512::cplib_bs512_from_string_char(a...);}
template<class... A> static auto avxFromStringEqual(A... a){return bitset_avx512::cplib_bs512_from_string_equal(a...);}
template<class... A> static auto avxIntersects(A... a){return bitset_avx512::cplib_bs512_intersects(a...);}
template<class... A> static auto avxSubset(A... a){return bitset_avx512::cplib_bs512_subset(a...);}
template<class... A> static auto avxSetRange(A... a){return bitset_avx512::cplib_bs512_set_range(a...);}
template<class... A> static auto avxClearRange(A... a){return bitset_avx512::cplib_bs512_clear_range(a...);}
template<class... A> static auto avxFlipRange(A... a){return bitset_avx512::cplib_bs512_flip_range(a...);}
template<class... A> static auto avxCmp(A... a){return bitset_avx512::cplib_bs512_cmp(a...);}
template<class... A> static auto avxAll(A... a){return bitset_avx512::cplib_bs512_all(a...);}
template<class... A> static auto avxAny(A... a){return bitset_avx512::cplib_bs512_any(a...);}
static std::size_t searchNext(const UInt* x,std::size_t lo,std::size_t hi){return bitset_search::cplib_bs_search_next(x,lo,hi);}
static std::size_t searchPrev(const UInt* x,std::size_t lo,std::size_t hi){return bitset_search::cplib_bs_search_prev(x,lo,hi);}
};
}
namespace cplib {
// 動的長と固定長は別ストレージ。AVX512版は元カーネルのCPU判定でAVX2へ切り替える。
template<int Width,Int Size=-1> class SimdBitSet {
    friend struct detail::SimdBitSetAccess;
    static_assert((Width==2||Width==512)&&Size>=-1);
    using Kernel=std::conditional_t<Width==2,detail::BitSetKernel2,detail::BitSetKernel512>;
    using Storage=std::conditional_t<Size<0,std::vector<UInt>,std::array<UInt,Size<0?0:(Size+63)/64>>;
    Storage bits_{};Int size_=Size<0?0:Size;
    void same(const SimdBitSet& y)const{if(size_!=y.size_)throw std::invalid_argument("BitSet sizes must match");}
    void check(Int i)const{if(i<0||i>=size_)throw std::out_of_range("BitSet index out of bounds");}
    void range(Int l,Int r)const{if(l<0||l>r||r>size_)throw std::out_of_range("BitSet range out of bounds");}
    void trim(){if(size_&63)bits_.back()&=(UInt(1)<<(size_&63))-1;}
    SimdBitSet empty()const{if constexpr(Size<0)return SimdBitSet(size_);else return {};}
    void initialize_bools(const std::vector<bool>& v){if(Int(v.size())>size_)throw std::invalid_argument("initial value is longer than BitSet size");if constexpr(Size<0){for(std::size_t i=0;i<v.size();++i)if(v[i])bits_[i>>6]|=UInt(1)<<(i&63);}else{std::vector<unsigned char> bytes(v.size());for(std::size_t i=0;i<v.size();++i)bytes[i]=v[i];if(size_)Kernel::avxFromBools(bits_.data(),bytes.data(),bytes.size(),bits_.size());}}
public:
    SimdBitSet()=default;explicit SimdBitSet(Int n) requires(Size<0):size_(n){if(n<0)throw std::invalid_argument("BitSet size must be non-negative");bits_.resize((n>>6)+bool(n&63));}
    explicit SimdBitSet(const std::vector<bool>& v){if constexpr(Size<0){size_=v.size();bits_.resize((size_+63)/64);}initialize_bools(v);}
    SimdBitSet(const std::vector<bool>& v,Int n) requires(Size<0):SimdBitSet(n){initialize_bools(v);}
    explicit SimdBitSet(std::span<const bool> v) requires(Size>=0){static_assert(sizeof(bool)==1);if(Int(v.size())>size_)throw std::invalid_argument("initial value is longer than BitSet size");if(size_)Kernel::avxFromBools(bits_.data(),v.data(),v.size(),bits_.size());}
    void fromString(std::string_view s,char match){if(Int(s.size())>size_)throw std::invalid_argument("source string is longer than BitSet size");if(size_)Kernel::avxFromStringChar(bits_.data(),s.data(),nullptr,static_cast<unsigned char>(match),s.size(),bits_.size());}
    void fromString(std::string_view s,std::string_view reference){if(Int(s.size())>size_||s.size()!=reference.size())throw std::invalid_argument("source/reference size mismatch");if(size_)Kernel::avxFromStringEqual(bits_.data(),s.data(),reference.data(),0,s.size(),bits_.size());}
    Int len()const{return size_;}std::span<const UInt> words()const{return bits_;}
    bool operator[](Int i)const{check(i);return (bits_[i>>6]>>(i&63))&1;}void set(Int i,bool v){check(i);if(v)bits_[i>>6]|=UInt(1)<<(i&63);else bits_[i>>6]&=~(UInt(1)<<(i&63));}void flip(Int i){check(i);bits_[i>>6]^=UInt(1)<<(i&63);}
    struct Reference {SimdBitSet* owner;Int i;operator bool()const{return std::as_const(*owner)[i];}Reference& operator=(bool v){owner->set(i,v);return *this;}Reference& operator=(Int v){if(v==0||v==1)owner->set(i,v==1);return *this;}Reference& operator=(int v){return *this=Int(v);}Reference& operator=(const Reference& v){return *this=bool(v);}};
    Reference operator[](Int i){return {this,i};}CPLIB_BACKWARDS_INDEX_OVERLOADS
    void andInto(const SimdBitSet& x,const SimdBitSet& y){same(x);same(y);if(size_)Kernel::avxAnd(bits_.data(),x.bits_.data(),y.bits_.data(),bits_.size());}
    SimdBitSet operator&(const SimdBitSet& y)const{same(y);auto out=empty();out.andInto(*this,y);return out;}
    SimdBitSet& operator&=(const SimdBitSet& y){andInto(*this,y);return *this;}
    void orInto(const SimdBitSet& x,const SimdBitSet& y){same(x);same(y);if(size_)Kernel::avxOr(bits_.data(),x.bits_.data(),y.bits_.data(),bits_.size());}
    SimdBitSet operator|(const SimdBitSet& y)const{same(y);auto out=empty();out.orInto(*this,y);return out;}
    SimdBitSet& operator|=(const SimdBitSet& y){orInto(*this,y);return *this;}
    void xorInto(const SimdBitSet& x,const SimdBitSet& y){same(x);same(y);if(size_)Kernel::avxXor(bits_.data(),x.bits_.data(),y.bits_.data(),bits_.size());}
    SimdBitSet operator^(const SimdBitSet& y)const{same(y);auto out=empty();out.xorInto(*this,y);return out;}
    SimdBitSet& operator^=(const SimdBitSet& y){xorInto(*this,y);return *this;}
    void andNotInto(const SimdBitSet& x,const SimdBitSet& y){same(x);same(y);if(size_)Kernel::avxAndNot(bits_.data(),x.bits_.data(),y.bits_.data(),bits_.size());}
    void andNotAssign(const SimdBitSet& y){andNotInto(*this,y);}
    void selectAssign(const SimdBitSet& y,const SimdBitSet& z){same(y);same(z);if(size_)Kernel::avxSelectAssign(bits_.data(),bits_.data(),y.bits_.data(),z.bits_.data(),bits_.size());}
    void orAndAssign(const SimdBitSet& y,const SimdBitSet& z){same(y);same(z);if(size_)Kernel::avxOrAndAssign(bits_.data(),bits_.data(),y.bits_.data(),z.bits_.data(),bits_.size());}
    void andOrAssign(const SimdBitSet& y,const SimdBitSet& z){same(y);same(z);if(size_)Kernel::avxAndOrAssign(bits_.data(),bits_.data(),y.bits_.data(),z.bits_.data(),bits_.size());}
    void xorAndAssign(const SimdBitSet& y,const SimdBitSet& z){same(y);same(z);if(size_)Kernel::avxXorAndAssign(bits_.data(),bits_.data(),y.bits_.data(),z.bits_.data(),bits_.size());}
    SimdBitSet majority(const SimdBitSet& y,const SimdBitSet& z)const{same(y);same(z);auto out=empty();if(size_)Kernel::avxMajority(out.bits_.data(),bits_.data(),y.bits_.data(),z.bits_.data(),bits_.size());return out;}
    void xnorInto(const SimdBitSet& x,const SimdBitSet& y){same(x);same(y);if(size_){Kernel::avxXnorAssign(bits_.data(),x.bits_.data(),y.bits_.data(),x.bits_.data(),bits_.size());trim();}}
    SimdBitSet xnor(const SimdBitSet& y)const{auto out=empty();out.xnorInto(*this,y);return out;}void xnorAssign(const SimdBitSet& y){xnorInto(*this,y);}
    SimdBitSet operator~()const{auto out=empty();if(size_){Kernel::avxNot(out.bits_.data(),bits_.data(),bits_.size());out.trim();}return out;}
    SimdBitSet operator<<(Int shift)const{if(shift<0)throw std::invalid_argument("negative shift");auto out=empty();if(shift<size_){Kernel::avxShl(out.bits_.data(),bits_.data(),bits_.size(),shift);out.trim();}return out;}
    SimdBitSet operator>>(Int shift)const{if(shift<0)throw std::invalid_argument("negative shift");auto out=empty();if(shift<size_)Kernel::avxShr(out.bits_.data(),bits_.data(),bits_.size(),shift);return out;}
    Int popcount()const{return size_?Int(Kernel::avxPopcount(bits_.data(),bits_.data(),bits_.size())):0;}
    Int popcountrange(Int l,Int r)const{range(l,r);return l<r?Int(Kernel::avxPopcountRange(bits_.data(),bits_.data(),l,r)):0;}
    Int andpopcount(const SimdBitSet& y)const{same(y);return size_?Int(Kernel::avxAndPopcount(bits_.data(),y.bits_.data(),bits_.size())):0;}
    Int andpopcountrange(const SimdBitSet& y,Int l,Int r)const{same(y);range(l,r);return l<r?Int(Kernel::avxAndpopcountRange(bits_.data(),y.bits_.data(),l,r)):0;}
    Int orpopcount(const SimdBitSet& y)const{same(y);return size_?Int(Kernel::avxOrPopcount(bits_.data(),y.bits_.data(),bits_.size())):0;}
    Int orpopcountrange(const SimdBitSet& y,Int l,Int r)const{same(y);range(l,r);return l<r?Int(Kernel::avxOrpopcountRange(bits_.data(),y.bits_.data(),l,r)):0;}
    Int xorpopcount(const SimdBitSet& y)const{same(y);return size_?Int(Kernel::avxXorPopcount(bits_.data(),y.bits_.data(),bits_.size())):0;}
    Int xorpopcountrange(const SimdBitSet& y,Int l,Int r)const{same(y);range(l,r);return l<r?Int(Kernel::avxXorpopcountRange(bits_.data(),y.bits_.data(),l,r)):0;}
    Int xnorpopcount(const SimdBitSet& y)const{return size_-xorpopcount(y);}
    bool intersects(const SimdBitSet& y)const{same(y);return size_&&Kernel::avxIntersects(bits_.data(),y.bits_.data(),bits_.size());}
    bool isSubsetOf(const SimdBitSet& y)const{same(y);return !size_||Kernel::avxSubset(bits_.data(),y.bits_.data(),bits_.size());}
    Int nextSetBit(Int start)const{if(start<0||start>size_)throw std::out_of_range("BitSet index");if(start==size_)return -1;Int i=start>>6;UInt word=bits_[i]&(~UInt(0)<<(start&63));if(word)return i*64+std::countr_zero(word);std::size_t j=Kernel::searchNext(bits_.data(),i+1,bits_.size());return j<bits_.size()?Int(j*64+std::countr_zero(bits_[j])):-1;}
    Int prevSetBit(Int start)const{if(start< -1||start>=size_)throw std::out_of_range("BitSet index");if(start==-1)return -1;Int i=start>>6;UInt word=bits_[i]&(~UInt(0)>>(63-(start&63)));if(word)return i*64+63-std::countl_zero(word);std::size_t j=Kernel::searchPrev(bits_.data(),0,i);return j<std::size_t(i)?Int(j*64+63-std::countl_zero(bits_[j])):-1;}
    Int lowestBit()const{return nextSetBit(0);}
    void setRange(Int l,Int r){range(l,r);if(l<r)Kernel::avxSetRange(bits_.data(),l,r);}
    void clearRange(Int l,Int r){range(l,r);if(l<r)Kernel::avxClearRange(bits_.data(),l,r);}
    void flipRange(Int l,Int r){range(l,r);if(l<r)Kernel::avxFlipRange(bits_.data(),l,r);}
    void clear(){clearRange(0,size_);}void fill(){setRange(0,size_);}void flipAll(){flipRange(0,size_);}
    int cmp(const SimdBitSet& y)const{same(y);return size_?Kernel::avxCmp(bits_.data(),y.bits_.data(),bits_.size()):0;}bool lexLess(const SimdBitSet& y)const{return cmp(y)<0;}auto operator<=>(const SimdBitSet& y)const{return cmp(y)<=>0;}bool operator==(const SimdBitSet& y)const{return size_==y.size_&&bits_==y.bits_;}
    bool all()const{return !size_||Kernel::avxAll(bits_.data(),size_);}bool any()const{return size_&&Kernel::avxAny(bits_.data(),size_);}
    std::string to_string()const{std::string out(size_,'0');for(Int i=0;i<size_;++i)if((*this)[i])out[size_-1-i]='1';return out;}
    struct Iterator {const SimdBitSet* owner;std::size_t word;UInt remaining;void skip(){if(!remaining&&word<owner->bits_.size()){word=Kernel::searchNext(owner->bits_.data(),word+1,owner->bits_.size());if(word<owner->bits_.size())remaining=owner->bits_[word];}}Int operator*()const{return word*64+std::countr_zero(remaining);}Iterator& operator++(){remaining&=remaining-1;skip();return *this;}bool operator==(const Iterator& other)const{return owner==other.owner&&word==other.word&&remaining==other.remaining;}};
    Iterator begin()const{Iterator it{this,0,bits_.empty()?0:bits_[0]};it.skip();return it;}Iterator end()const{return {this,bits_.size(),0};}const SimdBitSet& items()const{return *this;}
    Int andAssignPopcount(const SimdBitSet& y) requires(Width==512){same(y);return size_?Int(Kernel::avxAndAssignPopcount(bits_.data(),y.bits_.data(),bits_.data(),size_)):0;}
    Int orAssignPopcount(const SimdBitSet& y) requires(Width==512){same(y);return size_?Int(Kernel::avxOrAssignPopcount(bits_.data(),y.bits_.data(),bits_.data(),size_)):0;}
    Int xorAssignPopcount(const SimdBitSet& y) requires(Width==512){same(y);return size_?Int(Kernel::avxXorAssignPopcount(bits_.data(),y.bits_.data(),bits_.data(),size_)):0;}
    Int andNotAssignPopcount(const SimdBitSet& y) requires(Width==512){same(y);return size_?Int(Kernel::avxAndNotAssignPopcount(bits_.data(),y.bits_.data(),bits_.data(),size_)):0;}
    Int selectAssignPopcount(const SimdBitSet& y,const SimdBitSet& z) requires(Width==512){same(y);same(z);return size_?Int(Kernel::avxSelectAssignPopcount(bits_.data(),y.bits_.data(),z.bits_.data(),size_)):0;}
    Int orAndAssignPopcount(const SimdBitSet& y,const SimdBitSet& z) requires(Width==512){same(y);same(z);return size_?Int(Kernel::avxOrAndAssignPopcount(bits_.data(),y.bits_.data(),z.bits_.data(),size_)):0;}
    Int andOrAssignPopcount(const SimdBitSet& y,const SimdBitSet& z) requires(Width==512){same(y);same(z);return size_?Int(Kernel::avxAndOrAssignPopcount(bits_.data(),y.bits_.data(),z.bits_.data(),size_)):0;}
    Int xorAndAssignPopcount(const SimdBitSet& y,const SimdBitSet& z) requires(Width==512){same(y);same(z);return size_?Int(Kernel::avxXorAndAssignPopcount(bits_.data(),y.bits_.data(),z.bits_.data(),size_)):0;}
    Int xnorAssignPopcount(const SimdBitSet& y) requires(Width==512){same(y);return size_?Int(Kernel::avxXnorAssignPopcount(bits_.data(),y.bits_.data(),bits_.data(),size_)):0;}
    SimdBitSet operator+(const SimdBitSet& y)const requires(Width==512){same(y);auto out=empty();if(size_){Kernel::avxAdd(out.bits_.data(),bits_.data(),y.bits_.data(),bits_.size());out.trim();}return out;}
    SimdBitSet& operator+=(const SimdBitSet& y) requires(Width==512){same(y);if(size_){Kernel::avxAdd(bits_.data(),bits_.data(),y.bits_.data(),bits_.size());trim();}return *this;}
    SimdBitSet& operator<<=(Int k) requires(Width==512){if(k<0)throw std::invalid_argument("negative shift");if(k>=size_)clear();else if(k){detail::bitset_shift_assign::cplib_assign_left(bits_.data(),bits_.size(),k);trim();}return *this;}
    SimdBitSet& operator>>=(Int k) requires(Width==512){if(k<0)throw std::invalid_argument("negative shift");if(k>=size_)clear();else if(k)detail::bitset_shift_assign::cplib_assign_right(bits_.data(),bits_.size(),k);return *this;}
    void orShiftLeftAssign(Int k) requires(Width==512){if(k<0)throw std::invalid_argument("negative shift");if(k>0&&k<size_){detail::bitset_fuse_shift::cplib_fuse_left(bits_.data(),bits_.size(),k);trim();}}
    void orShiftRightAssign(Int k) requires(Width==512){if(k<0)throw std::invalid_argument("negative shift");if(k>0&&k<size_)detail::bitset_fuse_shift::cplib_fuse_right(bits_.data(),bits_.size(),k);}
    bool lastBit()const requires(Width==512){return size_&&(*this)[size_-1];}
};
template<int W,Int N> auto len(const SimdBitSet<W,N>& x){return x.len();}
template<int W,Int N> auto popcount(const SimdBitSet<W,N>& x){return x.popcount();}
template<int W,Int N> auto lowestBit(const SimdBitSet<W,N>& x){return x.lowestBit();}
template<int W,Int N> auto all(const SimdBitSet<W,N>& x){return x.all();}
template<int W,Int N> auto any(const SimdBitSet<W,N>& x){return x.any();}
template<int W,Int N> auto lastBit(const SimdBitSet<W,N>& x){return x.lastBit();}
template<int W,Int N> auto to_string(const SimdBitSet<W,N>& x){return x.to_string();}
template<int W,Int N> auto cmp(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.cmp(y);}
template<int W,Int N> auto lexLess(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.lexLess(y);}
template<int W,Int N> auto andpopcount(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.andpopcount(y);}
template<int W,Int N> auto orpopcount(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.orpopcount(y);}
template<int W,Int N> auto xorpopcount(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.xorpopcount(y);}
template<int W,Int N> auto xnorpopcount(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.xnorpopcount(y);}
template<int W,Int N> auto intersects(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.intersects(y);}
template<int W,Int N> auto isSubsetOf(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.isSubsetOf(y);}
template<int W,Int N> auto xnor(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.xnor(y);}
template<int W,Int N> void andInto(SimdBitSet<W,N>& dst,const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){dst.andInto(x,y);}
template<int W,Int N> void orInto(SimdBitSet<W,N>& dst,const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){dst.orInto(x,y);}
template<int W,Int N> void xorInto(SimdBitSet<W,N>& dst,const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){dst.xorInto(x,y);}
template<int W,Int N> void andNotInto(SimdBitSet<W,N>& dst,const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){dst.andNotInto(x,y);}
template<int W,Int N> void xnorInto(SimdBitSet<W,N>& dst,const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){dst.xnorInto(x,y);}
template<int W,Int N> auto selectAssign(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,const SimdBitSet<W,N>& z){return x.selectAssign(y,z);}
template<int W,Int N> auto orAndAssign(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,const SimdBitSet<W,N>& z){return x.orAndAssign(y,z);}
template<int W,Int N> auto andOrAssign(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,const SimdBitSet<W,N>& z){return x.andOrAssign(y,z);}
template<int W,Int N> auto xorAndAssign(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,const SimdBitSet<W,N>& z){return x.xorAndAssign(y,z);}
template<int W,Int N> auto selectAssignPopcount(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,const SimdBitSet<W,N>& z){return x.selectAssignPopcount(y,z);}
template<int W,Int N> auto orAndAssignPopcount(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,const SimdBitSet<W,N>& z){return x.orAndAssignPopcount(y,z);}
template<int W,Int N> auto andOrAssignPopcount(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,const SimdBitSet<W,N>& z){return x.andOrAssignPopcount(y,z);}
template<int W,Int N> auto xorAndAssignPopcount(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,const SimdBitSet<W,N>& z){return x.xorAndAssignPopcount(y,z);}
template<int W,Int N> auto andNotAssign(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.andNotAssign(y);}
template<int W,Int N> auto xnorAssign(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.xnorAssign(y);}
template<int W,Int N> auto andAssignPopcount(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.andAssignPopcount(y);}
template<int W,Int N> auto orAssignPopcount(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.orAssignPopcount(y);}
template<int W,Int N> auto xorAssignPopcount(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.xorAssignPopcount(y);}
template<int W,Int N> auto andNotAssignPopcount(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.andNotAssignPopcount(y);}
template<int W,Int N> auto xnorAssignPopcount(SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y){return x.xnorAssignPopcount(y);}
template<int W,Int N> void orShiftLeftAssign(SimdBitSet<W,N>& x,Int k){x.orShiftLeftAssign(k);}
template<int W,Int N> void orShiftRightAssign(SimdBitSet<W,N>& x,Int k){x.orShiftRightAssign(k);}
template<int W,Int N> void flip(SimdBitSet<W,N>& x,Int k){x.flip(k);}
template<int W,Int N> Int nextSetBit(const SimdBitSet<W,N>& x,Int k){return x.nextSetBit(k);}
template<int W,Int N> Int prevSetBit(const SimdBitSet<W,N>& x,Int k){return x.prevSetBit(k);}
template<int W,Int N> void clear(SimdBitSet<W,N>& x){x.clear();}
template<int W,Int N> void fill(SimdBitSet<W,N>& x){x.fill();}
template<int W,Int N> void flipAll(SimdBitSet<W,N>& x){x.flipAll();}
template<int W,Int N> void setRange(SimdBitSet<W,N>& x,Int l,Int r){x.setRange(l,r);}
template<int W,Int N> void clearRange(SimdBitSet<W,N>& x,Int l,Int r){x.clearRange(l,r);}
template<int W,Int N> void flipRange(SimdBitSet<W,N>& x,Int l,Int r){x.flipRange(l,r);}
template<int W,Int N> Int popcountrange(const SimdBitSet<W,N>& x,Int l,Int r){return x.popcountrange(l,r);}
template<int W,Int N> Int andpopcountrange(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,Int l,Int r){return x.andpopcountrange(y,l,r);}
template<int W,Int N> Int orpopcountrange(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,Int l,Int r){return x.orpopcountrange(y,l,r);}
template<int W,Int N> Int xorpopcountrange(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,Int l,Int r){return x.xorpopcountrange(y,l,r);}
template<int W,Int N> auto majority(const SimdBitSet<W,N>& x,const SimdBitSet<W,N>& y,const SimdBitSet<W,N>& z){return x.majority(y,z);}
}

namespace cplib::detail {
struct SimdBitSetAccess {
 template<int W,Int N> static UInt* data(SimdBitSet<W,N>& x){return x.bits_.data();}
 template<int W,Int N> static void trim(SimdBitSet<W,N>& x){x.trim();}
};
}
