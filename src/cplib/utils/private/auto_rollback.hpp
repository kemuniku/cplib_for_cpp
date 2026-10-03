#pragma once
#include <cplib/common.hpp>
#include <cstring>
#include <type_traits>
#include <functional>

namespace cplib {
// C++ cannot rewrite arbitrary assignments: use history.write(x) or set(x, value).
// Recorded objects must remain at the same address until rollback.
class RollbackLog {
    struct Entry {
        void *address;
        std::size_t offset, size;
    };

    std::vector<Entry> entries;
    std::vector<unsigned char> saved;

public:
    RollbackLog() = default;
    RollbackLog(const RollbackLog &) = delete;
    RollbackLog &operator=(const RollbackLog &) = delete;

    ~RollbackLog() {
        rollback(0);
    }

    Int snapshot() const {
        return entries.size();
    }

    // 安定したアドレスにあるxの変更前の値を保存する。時間・領域はsizeof(T)に比例する。
    template <class T> T &write(T &x) {
        static_assert(std::is_trivially_copyable_v<T>);
        auto offset = saved.size();
        saved.resize(offset + sizeof(T));
        std::memcpy(saved.data() + offset, &x, sizeof(T));
        entries.push_back({&x, offset, sizeof(T)});
        return x;
    }

    template <class T, class V> void set(T &x, V &&v) {
        write(x) = std::forward<V>(v);
    }

    // 指定位置以降の変更を逆順に復元する。復元する値の合計サイズに比例する時間。
    void rollback(Int position) {
        assert(0 <= position && position <= snapshot());
        while (snapshot() > position) {
            auto e = entries.back();
            std::memcpy(e.address, saved.data() + e.offset, e.size);
            saved.resize(e.offset);
            entries.pop_back();
        }
    }
};

namespace detail {
template <class Solver, class Apply, class Answer>
void run_auto_rollback(const Solver &solver, Apply apply, Answer answer) {
    RollbackLog log;
    std::vector<Int> checkpoints;
    solver.run(
        [&](auto &&...args) {
            checkpoints.push_back(log.snapshot());
            std::invoke(apply, log, std::forward<decltype(args)>(args)...);
        },
        [&] {
            log.rollback(checkpoints.back());
            checkpoints.pop_back();
        },
        answer);
}
}
}
