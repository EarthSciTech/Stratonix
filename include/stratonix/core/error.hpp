/// @file error.hpp
/// @brief Error type used inside the library. The C API converts these to error codes.
#pragma once
#include <stdexcept>
#include <string>

namespace stratonix {

/// Error categories. Values are stable: they are mirrored in the C API (stx_status).
enum class ErrorCode : int {
  Ok = 0,
  InvalidArgument = 1,
  Io = 2,
  Parse = 3,
  Validation = 4,
  Build = 5,
  CheckFailed = 6,
  NotImplemented = 7,
  Unsupported = 8,
  Internal = 9,
};

/// Library exception: a code, an optional description field path and a message.
class Error : public std::runtime_error {
 public:
  Error(ErrorCode code, std::string message, std::string field = {})
      : std::runtime_error(field.empty() ? message : field + ": " + message),
        code_(code), field_(std::move(field)) {}
  ErrorCode code() const noexcept { return code_; }
  const std::string& field() const noexcept { return field_; }

 private:
  ErrorCode code_;
  std::string field_;
};

}  // namespace stratonix
