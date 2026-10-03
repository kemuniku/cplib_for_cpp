#pragma once
#include <cplib/utils/private/temporary_rollback_log.hpp>

namespace cplib {
// Pass the log to the update callback; the body receives the logged update and log.
// updateへ渡すlogに記録した変更を、bodyの終了・例外時に復元する。
// 記録対象の変更先は復元まで生存し、アドレスを保つこと。
template <class Update, class Body> decltype(auto) withAutoRollback(Update update, Body body) {
    RollbackLog log;
    auto tracked = [&](auto &&...args) -> decltype(auto) {
        return std::invoke(update, log, std::forward<decltype(args)>(args)...);
    };
    return std::invoke(body, tracked, log);
}

// bodyへ渡すlogに記録した変更を、終了・例外時に開始前へ復元する。
// 同じlogを使う入れ子では、各スコープ内の最初の値だけを保存する。
template <class Body> auto Temporary(TemporaryRollbackLog &log, Body body) {
    auto c = log.beginTemporary();

    struct Guard {
        TemporaryRollbackLog &log;
        TemporaryCheckpoint c;

        ~Guard() {
            log.endTemporary(c);
        }
    } guard{log, c};

    return std::invoke(body, log);
}

template <class Body> auto Temporary(Body body) {
    TemporaryRollbackLog log;
    return Temporary(log, std::move(body));
}
}
