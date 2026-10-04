#pragma once
#include <cplib/tmpl/fastio.hpp>
#include <unordered_map>
#include <memory>
#include <cstring>
#include <functional>

namespace cplib {
namespace detail {
template <class T> inline unsigned char replayInputTypeToken;
}

// 型ごとの連続バッファと型の連続区間を記録する。1要素あたり期待償却O(1)。
template <bool OnlyInts = false> class ReplayInputState {
    struct Run {
        const void *kind;
        Int start;
    };

    struct BufferBase {
        virtual ~BufferBase() = default;
    };

    template <class T> struct Buffer : BufferBase {
        std::vector<T> values;
    };

    std::vector<UInt> entries_;
    std::vector<Run> runs_;
    std::unordered_map<const void *, std::unique_ptr<BufferBase>> buffers_;
    Int position_ = 0, runPosition_ = 0, depth_ = 0;

    template <class T> Buffer<T> &buffer(const void *kind) {
        auto it = buffers_.find(kind);
        if (it == buffers_.end())
            it = buffers_.emplace(kind, std::make_unique<Buffer<T>>()).first;
        return *static_cast<Buffer<T> *>(it->second.get());
    }

public:
    // 値を記録・再生する。小さなtrivially copyable型は直接格納し、ほかは型別バッファに保存する。
    // 再生する値の型は記録時と一致させること。
    template <class T, class Source> T read(Source &&source) {
        static_assert(!OnlyInts || std::is_same_v<T, Int>);
        constexpr bool inlineValue = std::is_trivially_copyable_v<T> && sizeof(T) <= sizeof(UInt);
        const void *kind = &detail::replayInputTypeToken<T>;
        T value{};
        if (position_ < Int(entries_.size())) {
            UInt entry = entries_[position_];
            if constexpr (!OnlyInts) {
                if (runPosition_ + 1 < Int(runs_.size()) &&
                    runs_[runPosition_ + 1].start == position_)
                    ++runPosition_;
                if (runs_[runPosition_].kind != kind)
                    throw std::invalid_argument("replayed input type mismatch");
            }
            if constexpr (inlineValue)
                std::memcpy(&value, &entry, sizeof(T));
            else
                value = buffer<T>(kind).values[entry];
            ++position_;
        } else {
            if (!depth_ && !entries_.empty()) {
                entries_.clear();
                if constexpr (!OnlyInts) {
                    runs_.clear();
                    runPosition_ = 0;
                    buffers_.clear();
                }
                position_ = 0;
            }
            value = std::invoke(std::forward<Source>(source));
            if (depth_) {
                if constexpr (!OnlyInts) {
                    if (runs_.empty() || runs_.back().kind != kind)
                        runs_.push_back({kind, Int(entries_.size())});
                    runPosition_ = runs_.size() - 1;
                }
                UInt entry = 0;
                if constexpr (inlineValue)
                    std::memcpy(&entry, &value, sizeof(T));
                else {
                    auto &b = buffer<T>(kind);
                    entry = b.values.size();
                    b.values.push_back(value);
                }
                entries_.push_back(entry);
                ++position_;
            }
        }
        return value;
    }

    // 例外で抜ける場合も読み取り位置を戻す。先読みはネストできる。
    template <class Body> decltype(auto) peekInput(Body &&body) {
        struct Guard {
            ReplayInputState *owner;
            Int position, run;

            ~Guard() {
                owner->position_ = position;
                owner->runPosition_ = run;
                --owner->depth_;
            }
        } guard{this, position_, runPosition_};

        ++depth_;
        return std::invoke(std::forward<Body>(body));
    }
};

template <class Source = cplib::InputSource, bool OnlyInts = false> class ReplayableInput {
    Source source_;
    ReplayInputState<OnlyInts> state_;

public:
    explicit ReplayableInput(Source source = {}) : source_(std::move(source)) {
    }

    template <class T = Int> T input() {
        return state_.template read<T>([&] { return source_.template input<T>(); });
    }

    template <class T = Int> std::vector<T> input(Int n) {
        assert(n >= 0);
        std::vector<T> out(n);
        for (auto &v : out)
            v = input<T>();
        return out;
    }

    Int ii() {
        return input<Int>();
    }

    auto lii(Int n) {
        return input<Int>(n);
    }

    std::string si()
        requires(!OnlyInts)
    {
        return input<std::string>();
    }

    template <class Body> decltype(auto) peekInput(Body &&body) {
        return state_.peekInput([&]() -> decltype(auto) {
            if constexpr (std::invocable<Body, ReplayableInput &>)
                return std::invoke(std::forward<Body>(body), *this);
            else
                return std::invoke(std::forward<Body>(body));
        });
    }
};

template <bool OnlyInts = false, class Source, class Body>
decltype(auto) replayableInput(Source source, Body &&body) {
    ReplayableInput<Source, OnlyInts> input(std::move(source));
    return std::invoke(std::forward<Body>(body), input);
}

template <bool OnlyInts = false, class Body> decltype(auto) replayableInput(Body &&body) {
    return replayableInput<OnlyInts>(cplib::InputSource{}, std::forward<Body>(body));
}

template <class Source, bool OnlyInts, class Body>
decltype(auto) peekInput(ReplayableInput<Source, OnlyInts> &input, Body &&body) {
    return input.peekInput(std::forward<Body>(body));
}
}
