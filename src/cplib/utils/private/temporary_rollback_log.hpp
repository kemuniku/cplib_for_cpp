#pragma once
#include <cplib/utils/private/auto_rollback.hpp>
#include <unordered_map>

namespace cplib {
struct TemporaryCheckpoint {
    Int position, parentStart;
};

struct TemporaryIndexCache {
    std::vector<Int> stamps;
    const void *owner = nullptr;
    const void *base = nullptr;
    std::size_t size = 0;
    Int length = 0;

    // 再利用する添字判定領域の長さをO(1)で返す。
    Int indexCapacity() const {
        return stamps.size();
    }
};

class TemporaryRollbackLog {
    struct Location {
        void *address;
        std::size_t size;
        bool operator==(const Location &) const = default;
    };

    struct Hash {
        std::size_t operator()(Location x) const {
            return std::hash<void *>{}(x.address) ^ (x.size * 0x9e3779b9);
        }
    };

    struct Entry {
        Location location;
        Int previous;
        std::size_t offset;
        Int *stamp;
    };

    std::vector<Entry> entries;
    std::unordered_map<Location, Int, Hash> latest;
    Int scopeStart = 0;
    std::vector<unsigned char> saved;
    std::vector<TemporaryIndexCache *> caches;

    template <class T> void save(T *p, Int previous, Int *stamp = nullptr) {
        static_assert(std::is_trivially_copyable_v<T>);
        auto offset = saved.size();
        saved.resize(offset + sizeof(T));
        std::memcpy(saved.data() + offset, p, sizeof(T));
        entries.push_back({{p, sizeof(T)}, previous, offset, stamp});
    }

public:
    // 空の一時変更履歴を作る。計算量 O(1)。
    TemporaryRollbackLog() = default;
    TemporaryRollbackLog(const TemporaryRollbackLog &) = delete;
    TemporaryRollbackLog &operator=(const TemporaryRollbackLog &) = delete;

    // 残る変更をすべて復元する。計算量は復元バイト数と履歴・キャッシュ数に比例。
    ~TemporaryRollbackLog() {
        restore(0);
    }

    // 現在保存している復元処理の数をO(1)で返す。
    Int size() const {
        return entries.size();
    }

    // 同じスコープでは同じアドレス・サイズの最初の値だけを保存する。期待償却 O(sizeof(T))。
    template <class T> void remember(T *p) {
        Location key{p, sizeof(T)};
        auto it = latest.find(key);
        Int prev = it == latest.end() ? -1 : it->second;
        if (prev >= scopeStart)
            return;
        Int i = size();
        save(p, prev);
        latest[key] = i;
    }

    // 添字用キャッシュで重複を判定し変更前の値を記録する。初回はキャッシュ長Lの確保 O(L)、以後は値の保存に償却 O(sizeof(T))。不一致時はrememberに委ねる。
    template <class T> void rememberIndexed(T *p, T *base, Int length, TemporaryIndexCache &cache) {
        if (!cache.owner) {
            if (length > Int(cache.stamps.size()))
                cache.stamps.resize(length);
            cache.owner = this;
            cache.base = base;
            cache.size = sizeof(T);
            cache.length = length;
            caches.push_back(&cache);
        }
        auto offset = reinterpret_cast<std::uintptr_t>(p) - reinterpret_cast<std::uintptr_t>(base);
        if (cache.owner != this || cache.base != base || cache.size != sizeof(T) ||
            cache.length != length || offset / sizeof(T) >= UInt(length) || offset % sizeof(T)) {
            remember(p);
            return;
        }
        auto &stamp = cache.stamps[offset / sizeof(T)];
        Int prev = stamp - 1;
        if (prev >= scopeStart)
            return;
        Int i = size();
        save(p, prev, &stamp);
        stamp = i + 1;
    }

    // 変更前の値を記録し参照を返す。計算量はrememberと同じ。
    template <class T> T &write(T &x) {
        remember(&x);
        return x;
    }

    // 変更前の値を記録してから代入する。計算量はrememberと代入の費用。
    template <class T, class V> void set(T &x, V &&v) {
        write(x) = std::forward<V>(v);
    }

    // 保存順の逆順で復元する。範囲が重なる値も最初の状態へ戻る。
    // 保存用vectorの容量を保ち、次の記録で再利用する。
    // 復元する履歴数K・合計サイズB、解放するキャッシュ数Cとして期待 O(K+B+C)。
    void restore(Int position) {
        assert(0 <= position && position <= size());
        while (size() > position) {
            auto e = entries.back();
            entries.pop_back();
            std::memcpy(e.location.address, saved.data() + e.offset, e.location.size);
            saved.resize(e.offset);
            if (e.stamp)
                *e.stamp = e.previous + 1;
            else if (e.previous < 0)
                latest.erase(e.location);
            else
                latest[e.location] = e.previous;
        }
        if (entries.empty()) {
            for (auto p : caches)
                p->owner = nullptr;
            caches.clear();
        }
    }

    // 入れ子の開始位置を保存し、新しいスコープで重複を判定する。O(1)。
    TemporaryCheckpoint beginTemporary() {
        TemporaryCheckpoint c{size(), scopeStart};
        scopeStart = size();
        return c;
    }

    // 入れ子の変更を復元し、親スコープの重複判定へ戻す。計算量はrestoreと同じ。
    void endTemporary(TemporaryCheckpoint c) {
        restore(c.position);
        scopeStart = c.parentStart;
    }
};
}
