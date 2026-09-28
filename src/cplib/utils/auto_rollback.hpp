#pragma once
#include <cplib/utils/private/temporary_rollback_log.hpp>
namespace cplib {
// Pass the log to the update callback; the body receives the logged update and log.
template<class Update,class Body> decltype(auto) withAutoRollback(Update update,Body body){RollbackLog log;auto tracked=[&](auto&&... args)->decltype(auto){return std::invoke(update,log,std::forward<decltype(args)>(args)...);};return std::invoke(body,tracked,log);}
template<class Body> auto Temporary(TemporaryRollbackLog& log,Body body){auto c=log.beginTemporary();struct Guard{TemporaryRollbackLog& log;TemporaryCheckpoint c;~Guard(){log.endTemporary(c);}} guard{log,c};return std::invoke(body,log);}
template<class Body> auto Temporary(Body body){TemporaryRollbackLog log;return Temporary(log,std::move(body));}
}
