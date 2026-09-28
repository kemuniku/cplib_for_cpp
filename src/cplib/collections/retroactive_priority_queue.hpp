#pragma once
#include <cplib/collections/private/retroactive_common.hpp>
namespace cplib {
template<class T> class RetroactivePriorityQueue {
 struct Node{Int balance=0,minPrefix=0,remaining=-1,removed=-1;};
 Int capacity,size=1,count_=0,dummyRemoved=0;SortOrder order;std::vector<QueueDebugKind> operations;std::vector<T> values;std::vector<bool> remaining;std::vector<Node> tree;T total{},pushTotal{};
 bool before(Int a,Int b)const{if(a==0)return false;if(b==0)return true;if(values[a]<values[b])return order==Ascending;if(values[b]<values[a])return order==Descending;return a<b;}
 Int choose(Int a,Int b,bool remains)const{if(a<0)return b;if(b<0)return a;return before(a,b)==remains?a:b;}
 void pull(Int i){auto a=tree[i*2],b=tree[i*2+1];tree[i]={a.balance+b.balance,std::min(a.minPrefix,a.balance+b.minPrefix),choose(a.remaining,b.remaining,true),choose(a.removed,b.removed,false)};}
 void refresh(Int t){Node node;if(!t){node.balance=dummyRemoved;if(dummyRemoved<capacity)node.remaining=0;if(dummyRemoved>0)node.removed=0;}else if(operations[t]==qdkPush){if(remaining[t])node.remaining=t;else{node.removed=t;node.balance=1;}}else if(operations[t]==qdkPop)node.balance=-1;node.minPrefix=node.balance;Int i=size+t;tree[i]=node;while(i>1){i/=2;pull(i);}}
 Int findBridge(Int node,Int l,Int r,Int t,Int prefix,bool first)const{if(first?r<t:l>=t)return -1;if(prefix+tree[node].minPrefix>0)return -1;if(r-l==1)return r;Int m=(l+r)/2,rp=prefix+tree[node*2].balance;Int out=first?findBridge(node*2,l,m,t,prefix,first):findBridge(node*2+1,m,r,t,rp,first);if(out<0)out=first?findBridge(node*2+1,m,r,t,rp,first):findBridge(node*2,l,m,t,prefix,first);return out;}
 Int previousBridge(Int t)const{return std::max<Int>(0,findBridge(1,0,size,t,0,false));}Int nextBridge(Int t)const{return findBridge(1,0,size,t,0,true);}
 Int candidate(Int l,Int r,bool remains)const{Int out=-1;for(l+=size,r+=size;l<r;l/=2,r/=2){if(l&1){out=choose(out,remains?tree[l].remaining:tree[l].removed,remains);++l;}if(r&1){--r;out=choose(out,remains?tree[r].remaining:tree[r].removed,remains);}}return out;}
 void changeRemaining(Int t,bool remains,QueueDelta<Int,T>& delta){assert(t>=0);if(!t)dummyRemoved+=remains?-1:1;else{remaining[t]=remains;QueueValue<Int,T> entry{t-1,values[t]};if(remains){++count_;if constexpr(std::is_arithmetic_v<T>)total+=values[t];delta.added.push_back(entry);}else{--count_;if constexpr(std::is_arithmetic_v<T>)total-=values[t];auto it=std::find_if(delta.added.begin(),delta.added.end(),[&](auto& x){return x.time==entry.time;});if(it!=delta.added.end())delta.added.erase(it);else delta.removed.push_back(entry);}}refresh(t);}
 void eraseOperation(Int t,QueueDelta<Int,T>& d){if(operations[t]==qdkNone)return;if(operations[t]==qdkPush){if constexpr(std::is_arithmetic_v<T>)pushTotal=detail::queue_sub(pushTotal,values[t]);if(remaining[t])changeRemaining(t,false,d);else changeRemaining(candidate(0,nextBridge(t+1),true),false,d);}else changeRemaining(candidate(previousBridge(t),capacity+1,false),true,d);operations[t]=qdkNone;values[t]=T{};refresh(t);}
public:
 explicit RetroactivePriorityQueue(Int n,SortOrder ord=Ascending):capacity(n),order(ord){assert(n>=0);while(size<n+1)size*=2;operations.resize(n+1);values.resize(n+1);remaining.resize(n+1);tree.resize(size*2);refresh(0);}
 QueueDelta<Int,T> erase(Int t){assert(0<=t&&t<capacity);QueueDelta<Int,T> d;eraseOperation(t+1,d);return d;}
 QueueDelta<Int,T> setPush(Int t,T value){assert(0<=t&&t<capacity);Int i=t+1;QueueDelta<Int,T> d;eraseOperation(i,d);Int x=candidate(previousBridge(i),capacity+1,false);operations[i]=qdkPush;values[i]=value;if constexpr(std::is_arithmetic_v<T>)pushTotal=detail::queue_add(pushTotal,value);if(x<0||before(x,i))changeRemaining(i,true,d);else{remaining[i]=false;refresh(i);changeRemaining(x,true,d);}return d;}
 QueueDelta<Int,T> setPop(Int t){assert(0<=t&&t<capacity);Int i=t+1;QueueDelta<Int,T> d;if(operations[i]==qdkPop)return d;eraseOperation(i,d);changeRemaining(candidate(0,nextBridge(i),true),false,d);operations[i]=qdkPop;refresh(i);return d;}
 Int len()const{return count_;}T sum()const requires std::is_arithmetic_v<T>{return total;}T poppedSum()const requires std::is_arithmetic_v<T>{return detail::queue_sub(pushTotal,total);}
 std::optional<T> peek()const{Int i=tree[1].remaining;return i<=0?std::optional<T>{}:values[i];}
 bool isRemaining(Int t)const{assert(0<=t&&t<capacity);return operations[t+1]==qdkPush&&remaining[t+1];}
 std::vector<QueueDebugEntry<Int,T>> debugOperations()const{std::vector<QueueDebugEntry<Int,T>> out;for(Int t=0;t<capacity;++t){QueueDebugEntry<Int,T> e;e.time=t;e.kind=operations[t+1];if(e.kind==qdkPush)e.value=values[t+1];out.push_back(e);}return out;}
 auto debugTimeline()const{auto out=debugOperations();replayQueueDebug(out,order);return out;}std::string debugDump()const{return formatQueueDebug(debugTimeline());}
};
template<class T> auto initRetroactivePriorityQueue(Int n,SortOrder order=Ascending){return RetroactivePriorityQueue<T>(n,order);}
}
