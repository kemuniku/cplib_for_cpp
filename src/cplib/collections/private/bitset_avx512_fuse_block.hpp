#pragma once
#include <cplib/collections/private/simd_bitset.hpp>
#include <cplib/collections/private/bitset_avx512_fuse_arithmetic.hpp>
#include <functional>
#include <optional>

namespace cplib {
class FuseBlock;

// 式はSSA節点の参照だけを保持し、ビット列の一時配列を作らない。
struct FuseExpr {
    FuseBlock *owner{};
    int id = -1;
};

class FuseVariable {
    FuseBlock *owner_;
    int variable_;

public:
    FuseVariable(FuseBlock *owner, int variable) : owner_(owner), variable_(variable) {
    }

    operator FuseExpr() const;
    FuseVariable &operator=(FuseExpr value);

    FuseVariable &operator=(const FuseVariable &value) {
        return *this = static_cast<FuseExpr>(value);
    }

    FuseVariable &operator&=(FuseExpr value);
    FuseVariable &operator|=(FuseExpr value);
    FuseVariable &operator^=(FuseExpr value);
    FuseVariable &operator+=(FuseExpr value);
    FuseVariable &operator<<=(Int shift);
    FuseVariable &operator>>=(Int shift);
    void low(bool value);
};

class FuseBlock {
    friend class FuseVariable;

    enum Kind {
        input,
        bit_and,
        bit_or,
        bit_xor,
        bit_not,
        add,
        shift_left,
        shift_right,
        low_bit,
        ternary
    };

    struct Node {
        Kind kind;
        int a = -1, b = -1, c = -1;
        Int shift = 0;
        bool seed = false;
        const SimdBitSet<512> *source = nullptr;
    };

    struct Variable {
        SimdBitSet<512> *target;
        int id;
        bool written = false;
    };

    struct Capture {
        bool *target;
        int id;
    };

    struct Instruction {
        int id;
        Node node;
    };

    std::vector<Node> nodes_;
    std::vector<Variable> variables_;
    std::vector<Capture> captures_;
    std::vector<int> uses_;
    std::vector<Instruction> scalar_plan_, vector_plan_;

    int push(Node node) {
        nodes_.push_back(node);
        return int(nodes_.size()) - 1;
    }

    void check(FuseExpr expr) const {
        if (expr.owner != this || expr.id < 0)
            throw std::invalid_argument("fusion expressions belong to different blocks");
    }

    static bool logical(Kind kind) {
        return kind >= bit_and && kind <= bit_not;
    }

    int variable(SimdBitSet<512> &value) {
        for (int j = 0; j < int(variables_.size()); ++j)
            if (variables_[j].target == &value)
                return j;
        int id = push(Node{input, -1, -1, -1, 0, false, &value});
        variables_.push_back({&value, id, false});
        return int(variables_.size()) - 1;
    }

    void compile() {
        uses_.assign(nodes_.size(), 0);
        std::function<void(int)> visit = [&](int id) {
            if (uses_[id]++)
                return;
            auto node = nodes_[id];
            if (node.a >= 0)
                visit(node.a);
            if (node.b >= 0)
                visit(node.b);
        };
        for (auto v : variables_)
            if (v.written)
                visit(v.id);
        for (auto c : captures_)
            visit(c.id);
        scalar_plan_.clear();
        vector_plan_.clear();
        for (int id = 0; id < int(nodes_.size()); ++id)
            if (uses_[id])
                scalar_plan_.push_back({id, nodes_[id]});
        std::vector<bool> emitted(nodes_.size());
        std::function<void(int)> emit = [&](int id) {
            if (emitted[id])
                return;
            emitted[id] = true;
            auto node = nodes_[id];
            if (logical(node.kind)) {
                std::vector<int> atoms;
                std::function<void(int)> collect = [&](int j) {
                    if (j != id && (uses_[j] > 1 || !logical(nodes_[j].kind))) {
                        if (std::find(atoms.begin(), atoms.end(), j) == atoms.end())
                            atoms.push_back(j);
                    } else {
                        collect(nodes_[j].a);
                        if (nodes_[j].b >= 0)
                            collect(nodes_[j].b);
                    }
                };
                collect(id);
                if (atoms.size() <= 3) {
                    std::function<bool(int, int)> truth = [&](int j, int row) {
                        auto it = std::find(atoms.begin(), atoms.end(), j);
                        if (it != atoms.end())
                            return bool(row & (1 << (2 - (it - atoms.begin()))));
                        bool a = truth(nodes_[j].a, row);
                        if (nodes_[j].kind == bit_not)
                            return !a;
                        bool b = truth(nodes_[j].b, row);
                        if (nodes_[j].kind == bit_and)
                            return a && b;
                        if (nodes_[j].kind == bit_or)
                            return a || b;
                        return a != b;
                    };
                    unsigned mask = 0;
                    for (int row = 0; row < 8; ++row)
                        if (truth(id, row))
                            mask |= 1u << row;
                    for (int atom : atoms)
                        emit(atom);
                    while (atoms.size() < 3)
                        atoms.push_back(atoms[0]);
                    vector_plan_.push_back(
                        {id, Node{ternary, atoms[0], atoms[1], atoms[2], Int(mask)}});
                    return;
                }
            }
            if (node.a >= 0)
                emit(node.a);
            if (node.b >= 0)
                emit(node.b);
            vector_plan_.push_back({id, node});
        };
        for (auto v : variables_)
            if (v.written)
                emit(v.id);
        for (auto c : captures_)
            emit(c.id);
    }

    // 未対応の大きな内部シフト・異なる長さは通常演算の意味を保つ。
    void eager() {
        std::vector<std::optional<SimdBitSet<512>>> values(nodes_.size());
        for (auto ins : scalar_plan_) {
            auto node = ins.node;
            auto &out = values[ins.id];
            auto get = [&](int id) -> const SimdBitSet<512> & { return *values[id]; };
            switch (node.kind) {
            case input:
                out = *node.source;
                break;
            case bit_and:
                out = get(node.a) & get(node.b);
                break;
            case bit_or:
                out = get(node.a) | get(node.b);
                break;
            case bit_xor:
                out = get(node.a) ^ get(node.b);
                break;
            case bit_not:
                out = ~get(node.a);
                break;
            case add:
                out = get(node.a) + get(node.b);
                break;
            case shift_left:
                out = get(node.a) << node.shift;
                break;
            case shift_right:
                out = get(node.a) >> node.shift;
                break;
            case low_bit:
                out = get(node.a);
                out->set(0, node.seed);
                break;
            default:
                break;
            }
        }
        for (auto c : captures_)
            *c.target = values[c.id]->lastBit();
        for (auto v : variables_)
            if (v.written)
                *v.target = *values[v.id];
    }

    void scalar(std::size_t n, Int bits, const std::vector<const UInt *> &sources) {
        using namespace detail::bitset_fuse_arithmetic;
        std::vector<UInt> values(nodes_.size()), previous(nodes_.size());
        std::vector<unsigned> carries(nodes_.size());
        for (std::size_t i = 0; i < n; ++i) {
            for (auto ins : scalar_plan_) {
                auto node = ins.node;
                UInt a = node.a < 0 ? 0 : values[node.a], b = node.b < 0 ? 0 : values[node.b],
                     value = 0;
                switch (node.kind) {
                case input:
                    value = sources[ins.id][i];
                    break;
                case bit_and:
                    value = a & b;
                    break;
                case bit_or:
                    value = a | b;
                    break;
                case bit_xor:
                    value = a ^ b;
                    break;
                case bit_not:
                    value = ~a;
                    break;
                case add:
                    value = cplib_fuse_add64(a, b, &carries[ins.id]);
                    break;
                case shift_left:
                case shift_right:
                    if (node.kind == shift_right || node.shift >= 64)
                        value = cplib_fuse_shift0(sources[node.a], n, i, node.shift,
                                                  node.kind == shift_left);
                    else {
                        value = (a << node.shift) | (previous[ins.id] >> (64 - node.shift));
                        previous[ins.id] = a;
                    }
                    break;
                case low_bit:
                    value = i ? a : (node.seed ? a | 1 : a & ~UInt(1));
                    break;
                default:
                    break;
                }
                values[ins.id] = value;
            }
            for (auto v : variables_)
                if (v.written)
                    detail::SimdBitSetAccess::data(*v.target)[i] = values[v.id];
        }
        for (auto c : captures_)
            *c.target = n && ((values[c.id] >> ((bits - 1) & 63)) & 1);
    }

    __attribute__((target("avx512f"))) void vector512(std::size_t n, Int bits,
                                                      const std::vector<const UInt *> &sources);
    __attribute__((target("avx2"))) void vector256(std::size_t n, Int bits,
                                                   const std::vector<const UInt *> &sources);

public:
    FuseVariable var(SimdBitSet<512> &value) {
        return {this, variable(value)};
    }

    FuseExpr read(const SimdBitSet<512> &value) {
        for (auto v : variables_)
            if (v.target == &value)
                return {this, v.id};
        for (int id = 0; id < int(nodes_.size()); ++id)
            if (nodes_[id].kind == input && nodes_[id].source == &value)
                return {this, id};
        return {this, push(Node{input, -1, -1, -1, 0, false, &value})};
    }

    FuseExpr binary(char op, FuseExpr a, FuseExpr b) {
        check(a);
        check(b);
        Kind kind = op == '&' ? bit_and : op == '|' ? bit_or : op == '^' ? bit_xor : add;
        return {this, push(Node{kind, a.id, b.id})};
    }

    FuseExpr negate(FuseExpr a) {
        check(a);
        return {this, push(Node{bit_not, a.id})};
    }

    FuseExpr shift(FuseExpr a, Int k, bool left) {
        check(a);
        if (k < 0)
            throw std::invalid_argument("negative shift");
        if (!k)
            return a;
        return {this, push(Node{left ? shift_left : shift_right, a.id, -1, -1, k})};
    }

    FuseExpr low(FuseExpr a, bool value) {
        check(a);
        return {this, push(Node{low_bit, a.id, -1, -1, 0, value})};
    }

    void assign(SimdBitSet<512> &dst, FuseExpr value) {
        check(value);
        int v = variable(dst);
        variables_[v].id = value.id;
        variables_[v].written = true;
    }

    void lastBit(bool &dst, FuseExpr value) {
        check(value);
        captures_.push_back({&dst, value.id});
    }

    // m節点、nビットに対してO(m ceil(n/64))、通常経路の追加領域O(m)。
    // シフトした入力が出力と同一の場合だけ、その入力を一度複製する。
    void run(bool single = false) {
        compile();
        if (scalar_plan_.empty())
            return;
        // 自己シフトと x |= x << k / x >> k は元と同じ追加領域O(1)の専用カーネル。
        if (single && captures_.empty()) {
            Variable *output = nullptr;
            int count = 0;
            for (auto &v : variables_)
                if (v.written) {
                    output = &v;
                    ++count;
                }
            if (count == 1) {
                auto root = nodes_[output->id];
                auto self = [&](int id) {
                    return id >= 0 && nodes_[id].kind == input &&
                           nodes_[id].source == output->target;
                };
                if ((root.kind == shift_left || root.kind == shift_right) && self(root.a)) {
                    if (root.kind == shift_left)
                        *output->target <<= root.shift;
                    else
                        *output->target >>= root.shift;
                    return;
                }
                if (root.kind == bit_or) {
                    int shifted = self(root.a) ? root.b : self(root.b) ? root.a : -1;
                    if (shifted >= 0) {
                        auto node = nodes_[shifted];
                        if ((node.kind == shift_left || node.kind == shift_right) && self(node.a)) {
                            if (node.kind == shift_left)
                                output->target->orShiftLeftAssign(node.shift);
                            else
                                output->target->orShiftRightAssign(node.shift);
                            return;
                        }
                    }
                }
            }
        }
        Int bits = -1;
        bool supported = true, arithmetic = false;
        for (auto ins : scalar_plan_) {
            auto node = ins.node;
            if (node.kind == input) {
                if (bits < 0)
                    bits = node.source->len();
                else if (bits != node.source->len())
                    supported = false;
            }
            if (node.kind == add)
                arithmetic = true;
            if (node.kind == shift_left || node.kind == shift_right) {
                if (nodes_[node.a].kind != input && (node.kind == shift_right || node.shift >= 64))
                    supported = false;
            }
            if (node.kind == low_bit && bits == 0)
                supported = false;
        }
        for (auto v : variables_)
            if (v.written && v.target->len() != bits)
                supported = false;
        if (!supported) {
            eager();
            return;
        }
        std::size_t n = (bits + 63) / 64;
        std::vector<const UInt *> sources(nodes_.size());
        std::vector<std::vector<UInt>> snapshots(nodes_.size());
        for (auto ins : scalar_plan_)
            if (ins.node.kind == input)
                sources[ins.id] = ins.node.source->words().data();
        for (auto ins : scalar_plan_)
            if ((ins.node.kind == shift_right ||
                 (ins.node.kind == shift_left && ins.node.shift >= 64)) &&
                nodes_[ins.node.a].kind == input) {
                int id = ins.node.a;
                for (auto v : variables_)
                    if (v.written && v.target == nodes_[id].source && n && snapshots[id].empty()) {
                        snapshots[id].assign(sources[id], sources[id] + n);
                        sources[id] = snapshots[id].data();
                    }
            }
        if (__builtin_cpu_supports("avx512f"))
            vector512(n, bits, sources);
        else if (single && !arithmetic && __builtin_cpu_supports("avx2"))
            vector256(n, bits, sources);
        else
            scalar(n, bits, sources);
        for (auto v : variables_)
            if (v.written)
                detail::SimdBitSetAccess::trim(*v.target);
    }
};

inline FuseVariable::operator FuseExpr() const {
    return {owner_, owner_->variables_[variable_].id};
}

inline FuseVariable &FuseVariable::operator=(FuseExpr value) {
    owner_->check(value);
    auto &v = owner_->variables_[variable_];
    v.id = value.id;
    v.written = true;
    return *this;
}

inline FuseExpr operator&(FuseExpr a, FuseExpr b) {
    return a.owner->binary('&', a, b);
}

inline FuseExpr operator|(FuseExpr a, FuseExpr b) {
    return a.owner->binary('|', a, b);
}

inline FuseExpr operator^(FuseExpr a, FuseExpr b) {
    return a.owner->binary('^', a, b);
}

inline FuseExpr operator+(FuseExpr a, FuseExpr b) {
    return a.owner->binary('+', a, b);
}

inline FuseExpr operator~(FuseExpr a) {
    return a.owner->negate(a);
}

inline FuseExpr operator<<(FuseExpr a, Int k) {
    return a.owner->shift(a, k, true);
}

inline FuseExpr operator>>(FuseExpr a, Int k) {
    return a.owner->shift(a, k, false);
}

inline FuseVariable &FuseVariable::operator&=(FuseExpr value) {
    return *this = static_cast<FuseExpr>(*this) & value;
}

inline FuseVariable &FuseVariable::operator|=(FuseExpr value) {
    return *this = static_cast<FuseExpr>(*this) | value;
}

inline FuseVariable &FuseVariable::operator^=(FuseExpr value) {
    return *this = static_cast<FuseExpr>(*this) ^ value;
}

inline FuseVariable &FuseVariable::operator+=(FuseExpr value) {
    return *this = static_cast<FuseExpr>(*this) + value;
}

inline FuseVariable &FuseVariable::operator<<=(Int k) {
    return *this = static_cast<FuseExpr>(*this) << k;
}

inline FuseVariable &FuseVariable::operator>>=(Int k) {
    return *this = static_cast<FuseExpr>(*this) >> k;
}

inline void FuseVariable::low(bool value) {
    *this = owner_->low(static_cast<FuseExpr>(*this), value);
}
}

namespace cplib {
namespace detail::fuse_detail {
__attribute__((target("avx512f"))) inline __m512i ternary_logic(__m512i a, __m512i b, __m512i c,
                                                                unsigned mask) {
    switch (mask) {
    case 0:
        return _mm512_ternarylogic_epi64(a, b, c, 0);
    case 1:
        return _mm512_ternarylogic_epi64(a, b, c, 1);
    case 2:
        return _mm512_ternarylogic_epi64(a, b, c, 2);
    case 3:
        return _mm512_ternarylogic_epi64(a, b, c, 3);
    case 4:
        return _mm512_ternarylogic_epi64(a, b, c, 4);
    case 5:
        return _mm512_ternarylogic_epi64(a, b, c, 5);
    case 6:
        return _mm512_ternarylogic_epi64(a, b, c, 6);
    case 7:
        return _mm512_ternarylogic_epi64(a, b, c, 7);
    case 8:
        return _mm512_ternarylogic_epi64(a, b, c, 8);
    case 9:
        return _mm512_ternarylogic_epi64(a, b, c, 9);
    case 10:
        return _mm512_ternarylogic_epi64(a, b, c, 10);
    case 11:
        return _mm512_ternarylogic_epi64(a, b, c, 11);
    case 12:
        return _mm512_ternarylogic_epi64(a, b, c, 12);
    case 13:
        return _mm512_ternarylogic_epi64(a, b, c, 13);
    case 14:
        return _mm512_ternarylogic_epi64(a, b, c, 14);
    case 15:
        return _mm512_ternarylogic_epi64(a, b, c, 15);
    case 16:
        return _mm512_ternarylogic_epi64(a, b, c, 16);
    case 17:
        return _mm512_ternarylogic_epi64(a, b, c, 17);
    case 18:
        return _mm512_ternarylogic_epi64(a, b, c, 18);
    case 19:
        return _mm512_ternarylogic_epi64(a, b, c, 19);
    case 20:
        return _mm512_ternarylogic_epi64(a, b, c, 20);
    case 21:
        return _mm512_ternarylogic_epi64(a, b, c, 21);
    case 22:
        return _mm512_ternarylogic_epi64(a, b, c, 22);
    case 23:
        return _mm512_ternarylogic_epi64(a, b, c, 23);
    case 24:
        return _mm512_ternarylogic_epi64(a, b, c, 24);
    case 25:
        return _mm512_ternarylogic_epi64(a, b, c, 25);
    case 26:
        return _mm512_ternarylogic_epi64(a, b, c, 26);
    case 27:
        return _mm512_ternarylogic_epi64(a, b, c, 27);
    case 28:
        return _mm512_ternarylogic_epi64(a, b, c, 28);
    case 29:
        return _mm512_ternarylogic_epi64(a, b, c, 29);
    case 30:
        return _mm512_ternarylogic_epi64(a, b, c, 30);
    case 31:
        return _mm512_ternarylogic_epi64(a, b, c, 31);
    case 32:
        return _mm512_ternarylogic_epi64(a, b, c, 32);
    case 33:
        return _mm512_ternarylogic_epi64(a, b, c, 33);
    case 34:
        return _mm512_ternarylogic_epi64(a, b, c, 34);
    case 35:
        return _mm512_ternarylogic_epi64(a, b, c, 35);
    case 36:
        return _mm512_ternarylogic_epi64(a, b, c, 36);
    case 37:
        return _mm512_ternarylogic_epi64(a, b, c, 37);
    case 38:
        return _mm512_ternarylogic_epi64(a, b, c, 38);
    case 39:
        return _mm512_ternarylogic_epi64(a, b, c, 39);
    case 40:
        return _mm512_ternarylogic_epi64(a, b, c, 40);
    case 41:
        return _mm512_ternarylogic_epi64(a, b, c, 41);
    case 42:
        return _mm512_ternarylogic_epi64(a, b, c, 42);
    case 43:
        return _mm512_ternarylogic_epi64(a, b, c, 43);
    case 44:
        return _mm512_ternarylogic_epi64(a, b, c, 44);
    case 45:
        return _mm512_ternarylogic_epi64(a, b, c, 45);
    case 46:
        return _mm512_ternarylogic_epi64(a, b, c, 46);
    case 47:
        return _mm512_ternarylogic_epi64(a, b, c, 47);
    case 48:
        return _mm512_ternarylogic_epi64(a, b, c, 48);
    case 49:
        return _mm512_ternarylogic_epi64(a, b, c, 49);
    case 50:
        return _mm512_ternarylogic_epi64(a, b, c, 50);
    case 51:
        return _mm512_ternarylogic_epi64(a, b, c, 51);
    case 52:
        return _mm512_ternarylogic_epi64(a, b, c, 52);
    case 53:
        return _mm512_ternarylogic_epi64(a, b, c, 53);
    case 54:
        return _mm512_ternarylogic_epi64(a, b, c, 54);
    case 55:
        return _mm512_ternarylogic_epi64(a, b, c, 55);
    case 56:
        return _mm512_ternarylogic_epi64(a, b, c, 56);
    case 57:
        return _mm512_ternarylogic_epi64(a, b, c, 57);
    case 58:
        return _mm512_ternarylogic_epi64(a, b, c, 58);
    case 59:
        return _mm512_ternarylogic_epi64(a, b, c, 59);
    case 60:
        return _mm512_ternarylogic_epi64(a, b, c, 60);
    case 61:
        return _mm512_ternarylogic_epi64(a, b, c, 61);
    case 62:
        return _mm512_ternarylogic_epi64(a, b, c, 62);
    case 63:
        return _mm512_ternarylogic_epi64(a, b, c, 63);
    case 64:
        return _mm512_ternarylogic_epi64(a, b, c, 64);
    case 65:
        return _mm512_ternarylogic_epi64(a, b, c, 65);
    case 66:
        return _mm512_ternarylogic_epi64(a, b, c, 66);
    case 67:
        return _mm512_ternarylogic_epi64(a, b, c, 67);
    case 68:
        return _mm512_ternarylogic_epi64(a, b, c, 68);
    case 69:
        return _mm512_ternarylogic_epi64(a, b, c, 69);
    case 70:
        return _mm512_ternarylogic_epi64(a, b, c, 70);
    case 71:
        return _mm512_ternarylogic_epi64(a, b, c, 71);
    case 72:
        return _mm512_ternarylogic_epi64(a, b, c, 72);
    case 73:
        return _mm512_ternarylogic_epi64(a, b, c, 73);
    case 74:
        return _mm512_ternarylogic_epi64(a, b, c, 74);
    case 75:
        return _mm512_ternarylogic_epi64(a, b, c, 75);
    case 76:
        return _mm512_ternarylogic_epi64(a, b, c, 76);
    case 77:
        return _mm512_ternarylogic_epi64(a, b, c, 77);
    case 78:
        return _mm512_ternarylogic_epi64(a, b, c, 78);
    case 79:
        return _mm512_ternarylogic_epi64(a, b, c, 79);
    case 80:
        return _mm512_ternarylogic_epi64(a, b, c, 80);
    case 81:
        return _mm512_ternarylogic_epi64(a, b, c, 81);
    case 82:
        return _mm512_ternarylogic_epi64(a, b, c, 82);
    case 83:
        return _mm512_ternarylogic_epi64(a, b, c, 83);
    case 84:
        return _mm512_ternarylogic_epi64(a, b, c, 84);
    case 85:
        return _mm512_ternarylogic_epi64(a, b, c, 85);
    case 86:
        return _mm512_ternarylogic_epi64(a, b, c, 86);
    case 87:
        return _mm512_ternarylogic_epi64(a, b, c, 87);
    case 88:
        return _mm512_ternarylogic_epi64(a, b, c, 88);
    case 89:
        return _mm512_ternarylogic_epi64(a, b, c, 89);
    case 90:
        return _mm512_ternarylogic_epi64(a, b, c, 90);
    case 91:
        return _mm512_ternarylogic_epi64(a, b, c, 91);
    case 92:
        return _mm512_ternarylogic_epi64(a, b, c, 92);
    case 93:
        return _mm512_ternarylogic_epi64(a, b, c, 93);
    case 94:
        return _mm512_ternarylogic_epi64(a, b, c, 94);
    case 95:
        return _mm512_ternarylogic_epi64(a, b, c, 95);
    case 96:
        return _mm512_ternarylogic_epi64(a, b, c, 96);
    case 97:
        return _mm512_ternarylogic_epi64(a, b, c, 97);
    case 98:
        return _mm512_ternarylogic_epi64(a, b, c, 98);
    case 99:
        return _mm512_ternarylogic_epi64(a, b, c, 99);
    case 100:
        return _mm512_ternarylogic_epi64(a, b, c, 100);
    case 101:
        return _mm512_ternarylogic_epi64(a, b, c, 101);
    case 102:
        return _mm512_ternarylogic_epi64(a, b, c, 102);
    case 103:
        return _mm512_ternarylogic_epi64(a, b, c, 103);
    case 104:
        return _mm512_ternarylogic_epi64(a, b, c, 104);
    case 105:
        return _mm512_ternarylogic_epi64(a, b, c, 105);
    case 106:
        return _mm512_ternarylogic_epi64(a, b, c, 106);
    case 107:
        return _mm512_ternarylogic_epi64(a, b, c, 107);
    case 108:
        return _mm512_ternarylogic_epi64(a, b, c, 108);
    case 109:
        return _mm512_ternarylogic_epi64(a, b, c, 109);
    case 110:
        return _mm512_ternarylogic_epi64(a, b, c, 110);
    case 111:
        return _mm512_ternarylogic_epi64(a, b, c, 111);
    case 112:
        return _mm512_ternarylogic_epi64(a, b, c, 112);
    case 113:
        return _mm512_ternarylogic_epi64(a, b, c, 113);
    case 114:
        return _mm512_ternarylogic_epi64(a, b, c, 114);
    case 115:
        return _mm512_ternarylogic_epi64(a, b, c, 115);
    case 116:
        return _mm512_ternarylogic_epi64(a, b, c, 116);
    case 117:
        return _mm512_ternarylogic_epi64(a, b, c, 117);
    case 118:
        return _mm512_ternarylogic_epi64(a, b, c, 118);
    case 119:
        return _mm512_ternarylogic_epi64(a, b, c, 119);
    case 120:
        return _mm512_ternarylogic_epi64(a, b, c, 120);
    case 121:
        return _mm512_ternarylogic_epi64(a, b, c, 121);
    case 122:
        return _mm512_ternarylogic_epi64(a, b, c, 122);
    case 123:
        return _mm512_ternarylogic_epi64(a, b, c, 123);
    case 124:
        return _mm512_ternarylogic_epi64(a, b, c, 124);
    case 125:
        return _mm512_ternarylogic_epi64(a, b, c, 125);
    case 126:
        return _mm512_ternarylogic_epi64(a, b, c, 126);
    case 127:
        return _mm512_ternarylogic_epi64(a, b, c, 127);
    case 128:
        return _mm512_ternarylogic_epi64(a, b, c, 128);
    case 129:
        return _mm512_ternarylogic_epi64(a, b, c, 129);
    case 130:
        return _mm512_ternarylogic_epi64(a, b, c, 130);
    case 131:
        return _mm512_ternarylogic_epi64(a, b, c, 131);
    case 132:
        return _mm512_ternarylogic_epi64(a, b, c, 132);
    case 133:
        return _mm512_ternarylogic_epi64(a, b, c, 133);
    case 134:
        return _mm512_ternarylogic_epi64(a, b, c, 134);
    case 135:
        return _mm512_ternarylogic_epi64(a, b, c, 135);
    case 136:
        return _mm512_ternarylogic_epi64(a, b, c, 136);
    case 137:
        return _mm512_ternarylogic_epi64(a, b, c, 137);
    case 138:
        return _mm512_ternarylogic_epi64(a, b, c, 138);
    case 139:
        return _mm512_ternarylogic_epi64(a, b, c, 139);
    case 140:
        return _mm512_ternarylogic_epi64(a, b, c, 140);
    case 141:
        return _mm512_ternarylogic_epi64(a, b, c, 141);
    case 142:
        return _mm512_ternarylogic_epi64(a, b, c, 142);
    case 143:
        return _mm512_ternarylogic_epi64(a, b, c, 143);
    case 144:
        return _mm512_ternarylogic_epi64(a, b, c, 144);
    case 145:
        return _mm512_ternarylogic_epi64(a, b, c, 145);
    case 146:
        return _mm512_ternarylogic_epi64(a, b, c, 146);
    case 147:
        return _mm512_ternarylogic_epi64(a, b, c, 147);
    case 148:
        return _mm512_ternarylogic_epi64(a, b, c, 148);
    case 149:
        return _mm512_ternarylogic_epi64(a, b, c, 149);
    case 150:
        return _mm512_ternarylogic_epi64(a, b, c, 150);
    case 151:
        return _mm512_ternarylogic_epi64(a, b, c, 151);
    case 152:
        return _mm512_ternarylogic_epi64(a, b, c, 152);
    case 153:
        return _mm512_ternarylogic_epi64(a, b, c, 153);
    case 154:
        return _mm512_ternarylogic_epi64(a, b, c, 154);
    case 155:
        return _mm512_ternarylogic_epi64(a, b, c, 155);
    case 156:
        return _mm512_ternarylogic_epi64(a, b, c, 156);
    case 157:
        return _mm512_ternarylogic_epi64(a, b, c, 157);
    case 158:
        return _mm512_ternarylogic_epi64(a, b, c, 158);
    case 159:
        return _mm512_ternarylogic_epi64(a, b, c, 159);
    case 160:
        return _mm512_ternarylogic_epi64(a, b, c, 160);
    case 161:
        return _mm512_ternarylogic_epi64(a, b, c, 161);
    case 162:
        return _mm512_ternarylogic_epi64(a, b, c, 162);
    case 163:
        return _mm512_ternarylogic_epi64(a, b, c, 163);
    case 164:
        return _mm512_ternarylogic_epi64(a, b, c, 164);
    case 165:
        return _mm512_ternarylogic_epi64(a, b, c, 165);
    case 166:
        return _mm512_ternarylogic_epi64(a, b, c, 166);
    case 167:
        return _mm512_ternarylogic_epi64(a, b, c, 167);
    case 168:
        return _mm512_ternarylogic_epi64(a, b, c, 168);
    case 169:
        return _mm512_ternarylogic_epi64(a, b, c, 169);
    case 170:
        return _mm512_ternarylogic_epi64(a, b, c, 170);
    case 171:
        return _mm512_ternarylogic_epi64(a, b, c, 171);
    case 172:
        return _mm512_ternarylogic_epi64(a, b, c, 172);
    case 173:
        return _mm512_ternarylogic_epi64(a, b, c, 173);
    case 174:
        return _mm512_ternarylogic_epi64(a, b, c, 174);
    case 175:
        return _mm512_ternarylogic_epi64(a, b, c, 175);
    case 176:
        return _mm512_ternarylogic_epi64(a, b, c, 176);
    case 177:
        return _mm512_ternarylogic_epi64(a, b, c, 177);
    case 178:
        return _mm512_ternarylogic_epi64(a, b, c, 178);
    case 179:
        return _mm512_ternarylogic_epi64(a, b, c, 179);
    case 180:
        return _mm512_ternarylogic_epi64(a, b, c, 180);
    case 181:
        return _mm512_ternarylogic_epi64(a, b, c, 181);
    case 182:
        return _mm512_ternarylogic_epi64(a, b, c, 182);
    case 183:
        return _mm512_ternarylogic_epi64(a, b, c, 183);
    case 184:
        return _mm512_ternarylogic_epi64(a, b, c, 184);
    case 185:
        return _mm512_ternarylogic_epi64(a, b, c, 185);
    case 186:
        return _mm512_ternarylogic_epi64(a, b, c, 186);
    case 187:
        return _mm512_ternarylogic_epi64(a, b, c, 187);
    case 188:
        return _mm512_ternarylogic_epi64(a, b, c, 188);
    case 189:
        return _mm512_ternarylogic_epi64(a, b, c, 189);
    case 190:
        return _mm512_ternarylogic_epi64(a, b, c, 190);
    case 191:
        return _mm512_ternarylogic_epi64(a, b, c, 191);
    case 192:
        return _mm512_ternarylogic_epi64(a, b, c, 192);
    case 193:
        return _mm512_ternarylogic_epi64(a, b, c, 193);
    case 194:
        return _mm512_ternarylogic_epi64(a, b, c, 194);
    case 195:
        return _mm512_ternarylogic_epi64(a, b, c, 195);
    case 196:
        return _mm512_ternarylogic_epi64(a, b, c, 196);
    case 197:
        return _mm512_ternarylogic_epi64(a, b, c, 197);
    case 198:
        return _mm512_ternarylogic_epi64(a, b, c, 198);
    case 199:
        return _mm512_ternarylogic_epi64(a, b, c, 199);
    case 200:
        return _mm512_ternarylogic_epi64(a, b, c, 200);
    case 201:
        return _mm512_ternarylogic_epi64(a, b, c, 201);
    case 202:
        return _mm512_ternarylogic_epi64(a, b, c, 202);
    case 203:
        return _mm512_ternarylogic_epi64(a, b, c, 203);
    case 204:
        return _mm512_ternarylogic_epi64(a, b, c, 204);
    case 205:
        return _mm512_ternarylogic_epi64(a, b, c, 205);
    case 206:
        return _mm512_ternarylogic_epi64(a, b, c, 206);
    case 207:
        return _mm512_ternarylogic_epi64(a, b, c, 207);
    case 208:
        return _mm512_ternarylogic_epi64(a, b, c, 208);
    case 209:
        return _mm512_ternarylogic_epi64(a, b, c, 209);
    case 210:
        return _mm512_ternarylogic_epi64(a, b, c, 210);
    case 211:
        return _mm512_ternarylogic_epi64(a, b, c, 211);
    case 212:
        return _mm512_ternarylogic_epi64(a, b, c, 212);
    case 213:
        return _mm512_ternarylogic_epi64(a, b, c, 213);
    case 214:
        return _mm512_ternarylogic_epi64(a, b, c, 214);
    case 215:
        return _mm512_ternarylogic_epi64(a, b, c, 215);
    case 216:
        return _mm512_ternarylogic_epi64(a, b, c, 216);
    case 217:
        return _mm512_ternarylogic_epi64(a, b, c, 217);
    case 218:
        return _mm512_ternarylogic_epi64(a, b, c, 218);
    case 219:
        return _mm512_ternarylogic_epi64(a, b, c, 219);
    case 220:
        return _mm512_ternarylogic_epi64(a, b, c, 220);
    case 221:
        return _mm512_ternarylogic_epi64(a, b, c, 221);
    case 222:
        return _mm512_ternarylogic_epi64(a, b, c, 222);
    case 223:
        return _mm512_ternarylogic_epi64(a, b, c, 223);
    case 224:
        return _mm512_ternarylogic_epi64(a, b, c, 224);
    case 225:
        return _mm512_ternarylogic_epi64(a, b, c, 225);
    case 226:
        return _mm512_ternarylogic_epi64(a, b, c, 226);
    case 227:
        return _mm512_ternarylogic_epi64(a, b, c, 227);
    case 228:
        return _mm512_ternarylogic_epi64(a, b, c, 228);
    case 229:
        return _mm512_ternarylogic_epi64(a, b, c, 229);
    case 230:
        return _mm512_ternarylogic_epi64(a, b, c, 230);
    case 231:
        return _mm512_ternarylogic_epi64(a, b, c, 231);
    case 232:
        return _mm512_ternarylogic_epi64(a, b, c, 232);
    case 233:
        return _mm512_ternarylogic_epi64(a, b, c, 233);
    case 234:
        return _mm512_ternarylogic_epi64(a, b, c, 234);
    case 235:
        return _mm512_ternarylogic_epi64(a, b, c, 235);
    case 236:
        return _mm512_ternarylogic_epi64(a, b, c, 236);
    case 237:
        return _mm512_ternarylogic_epi64(a, b, c, 237);
    case 238:
        return _mm512_ternarylogic_epi64(a, b, c, 238);
    case 239:
        return _mm512_ternarylogic_epi64(a, b, c, 239);
    case 240:
        return _mm512_ternarylogic_epi64(a, b, c, 240);
    case 241:
        return _mm512_ternarylogic_epi64(a, b, c, 241);
    case 242:
        return _mm512_ternarylogic_epi64(a, b, c, 242);
    case 243:
        return _mm512_ternarylogic_epi64(a, b, c, 243);
    case 244:
        return _mm512_ternarylogic_epi64(a, b, c, 244);
    case 245:
        return _mm512_ternarylogic_epi64(a, b, c, 245);
    case 246:
        return _mm512_ternarylogic_epi64(a, b, c, 246);
    case 247:
        return _mm512_ternarylogic_epi64(a, b, c, 247);
    case 248:
        return _mm512_ternarylogic_epi64(a, b, c, 248);
    case 249:
        return _mm512_ternarylogic_epi64(a, b, c, 249);
    case 250:
        return _mm512_ternarylogic_epi64(a, b, c, 250);
    case 251:
        return _mm512_ternarylogic_epi64(a, b, c, 251);
    case 252:
        return _mm512_ternarylogic_epi64(a, b, c, 252);
    case 253:
        return _mm512_ternarylogic_epi64(a, b, c, 253);
    case 254:
        return _mm512_ternarylogic_epi64(a, b, c, 254);
    case 255:
        return _mm512_ternarylogic_epi64(a, b, c, 255);
    default:
        return _mm512_setzero_si512();
    }
}
}

__attribute__((target("avx2"))) inline void
FuseBlock::vector256(std::size_t n, Int bits, const std::vector<const UInt *> &sources) {
    using namespace detail::bitset_fuse_arithmetic;

    struct alignas(32) Value {
        __m256i value;
    };

    std::vector<Value> values(nodes_.size()), previous(nodes_.size());
    std::vector<unsigned> carries(nodes_.size());
    for (std::size_t i = 0; i < n; i += 4) {
        for (auto ins : scalar_plan_) {
            auto node = ins.node;
            __m256i a = node.a < 0 ? _mm256_setzero_si256() : values[node.a].value,
                    b = node.b < 0 ? _mm256_setzero_si256() : values[node.b].value,
                    value = _mm256_setzero_si256();
            switch (node.kind) {
            case input:
                value = cplib_fuse_load256(sources[ins.id], n, i);
                break;
            case bit_and:
                value = _mm256_and_si256(a, b);
                break;
            case bit_or:
                value = _mm256_or_si256(a, b);
                break;
            case bit_xor:
                value = _mm256_xor_si256(a, b);
                break;
            case bit_not:
                value = _mm256_xor_si256(a, _mm256_set1_epi64x(-1));
                break;
            case shift_left:
            case shift_right:
                if (node.kind == shift_right || node.shift >= 64)
                    value = cplib_fuse_shift256(sources[node.a], n, i, node.shift,
                                                node.kind == shift_left);
                else {
                    __m256i bridge = _mm256_permute2x128_si256(previous[ins.id].value, a, 0x21);
                    __m256i adjacent = _mm256_alignr_epi8(a, bridge, 8);
                    value = _mm256_or_si256(
                        _mm256_sll_epi64(a, _mm_cvtsi32_si128(int(node.shift))),
                        _mm256_srl_epi64(adjacent, _mm_cvtsi32_si128(int(64 - node.shift))));
                    previous[ins.id].value = a;
                }
                break;
            case low_bit: {
                alignas(32) UInt first[4] = {};
                first[0] = i == 0;
                __m256i mask = _mm256_loadu_si256((const __m256i *)first);
                value = node.seed ? _mm256_or_si256(a, mask) : _mm256_andnot_si256(mask, a);
                break;
            }
            default:
                break;
            }
            values[ins.id].value = value;
        }
        for (auto v : variables_)
            if (v.written) {
                auto *dst = detail::SimdBitSetAccess::data(*v.target) + i;
                if (n - i >= 4)
                    _mm256_storeu_si256((__m256i *)dst, values[v.id].value);
                else {
                    alignas(32) UInt tmp[4];
                    _mm256_storeu_si256((__m256i *)tmp, values[v.id].value);
                    std::copy_n(tmp, n - i, dst);
                }
            }
    }
    for (auto c : captures_) {
        alignas(32) UInt tmp[4];
        _mm256_storeu_si256((__m256i *)tmp, values[c.id].value);
        *c.target = n && ((tmp[(n - 1) % 4] >> ((bits - 1) & 63)) & 1);
    }
}

__attribute__((target("avx512f"))) inline void
FuseBlock::vector512(std::size_t n, Int bits, const std::vector<const UInt *> &sources) {
    using namespace detail::bitset_fuse_arithmetic;

    struct alignas(64) Value {
        __m512i value;
    };

    std::vector<Value> values(nodes_.size()), previous(nodes_.size());
    std::vector<unsigned> carries(nodes_.size());
    for (std::size_t i = 0; i < n; i += 8) {
        for (auto ins : vector_plan_) {
            auto node = ins.node;
            __m512i a = node.a < 0 ? _mm512_setzero_si512() : values[node.a].value,
                    b = node.b < 0 ? _mm512_setzero_si512() : values[node.b].value,
                    value = _mm512_setzero_si512();
            switch (node.kind) {
            case input:
                value = cplib_fuse_load512(sources[ins.id], n, i);
                break;
            case bit_and:
                value = _mm512_and_si512(a, b);
                break;
            case bit_or:
                value = _mm512_or_si512(a, b);
                break;
            case bit_xor:
                value = _mm512_xor_si512(a, b);
                break;
            case bit_not:
                value = _mm512_xor_si512(a, _mm512_set1_epi64(-1));
                break;
            case add:
                value = cplib_fuse_add512(a, b, &carries[ins.id]);
                break;
            case ternary:
                value = detail::fuse_detail::ternary_logic(a, b, values[node.c].value,
                                                           unsigned(node.shift));
                break;
            case shift_left:
            case shift_right:
                if (node.kind == shift_right || node.shift >= 64)
                    value = cplib_fuse_shift512(sources[node.a], n, i, node.shift,
                                                node.kind == shift_left);
                else {
                    __m512i adjacent = _mm512_alignr_epi64(a, previous[ins.id].value, 7);
                    value = _mm512_or_si512(
                        _mm512_sll_epi64(a, _mm_cvtsi32_si128(int(node.shift))),
                        _mm512_srl_epi64(adjacent, _mm_cvtsi32_si128(int(64 - node.shift))));
                    previous[ins.id].value = a;
                }
                break;
            case low_bit: {
                alignas(64) UInt first[8] = {};
                first[0] = i == 0;
                __m512i mask = _mm512_loadu_si512((const void *)first);
                value = node.seed ? _mm512_or_si512(a, mask) : _mm512_andnot_si512(mask, a);
                break;
            }
            default:
                break;
            }
            values[ins.id].value = value;
        }
        __mmask8 valid = __mmask8((1u << std::min<std::size_t>(8, n - i)) - 1);
        for (auto v : variables_)
            if (v.written)
                _mm512_mask_storeu_epi64(detail::SimdBitSetAccess::data(*v.target) + i, valid,
                                         values[v.id].value);
    }
    for (auto c : captures_) {
        alignas(64) UInt tmp[8];
        _mm512_storeu_si512((void *)tmp, values[c.id].value);
        *c.target = n && ((tmp[(n - 1) % 8] >> ((bits - 1) & 63)) & 1);
    }
}
}
