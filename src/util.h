//===-- util.h -------*- C++ -*-===//
///
/// \file
/// Contains a macro that wraps printf to conditionally print
///
//===----------------------------------------------------------------------===//

#define DEBUG 1
#define debug_print(fmt, ...)                                                  \
  do {                                                                         \
    if (DEBUG)                                                                 \
      fprintf(stderr, fmt, __VA_ARGS__);                                       \
  } while (0)
