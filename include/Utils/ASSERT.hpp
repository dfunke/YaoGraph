#pragma once

#pragma once

#include <exception>
#include <string>

#ifndef NDEBUG

#include <csignal>

#endif

class AssertionException : public std::exception {

public:
    AssertionException() { };

    AssertionException(const std::string &expr, const std::string &file,
                       const uint line)
        : _what(file + ":" + std::to_string(line) + ": failed assertion `" +
                expr + "'") { }

    AssertionException(const std::string &expr, const std::string &file,
                       const uint line, const std::string & msg)
        : _what(file + ":" + std::to_string(line) + ": failed assertion `" +
                expr + "' - " + msg) { }

    const char *what() const noexcept { return _what.c_str(); }

private:
    std::string _what;
};

#ifdef NDEBUG

#define ASSERT(e) ((void)(0))
#define MASSERT(e, msg) ((void)(0))
#define RASSERT(e) ((void)(0))

#else // NDEBUG

#define ASSERT(e)                                                              \
  ((void)((e) ? 0 : throw AssertionException(#e, __FILE__, __LINE__)))

#define MASSERT(e, msg)                                                              \
  ((void)((e) ? 0 : throw AssertionException(#e, __FILE__, __LINE__, msg)))

#define RASSERT(e)                                                              \
  ((void)((e) ? 0 : raise(SIGINT)))

#endif
