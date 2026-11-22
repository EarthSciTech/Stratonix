// Command implementations and a small argument parser (CLI11 is not available in every build
// environment; this keeps the same command structure).
#pragma once
#include <map>
#include <set>
#include <string>
#include <vector>

#include "stratonix/core/error.hpp"

namespace stratonix::cli {

/// Arguments: positionals, `--key value`, `--key=value`, `-o value`, and boolean `--flag`.
class Args {
 public:
  explicit Args(const std::vector<std::string>& a) {
    static const std::set<std::string> flags = {"help", "h", "zip", "reference-kernels", "no-figures", "verbose", "json"};
    for (std::size_t i = 0; i < a.size(); ++i) {
      std::string s = a[i];
      if (s.size() > 1 && s[0] == '-' && !(s.size() > 1 && (std::isdigit(static_cast<unsigned char>(s[1])) || s[1] == '.'))) {
        std::string k = s.substr(s[1] == '-' ? 2 : 1);
        auto eq = k.find('=');
        if (eq != std::string::npos) { kv_[k.substr(0, eq)] = k.substr(eq + 1); continue; }
        if (flags.count(k)) { flags_.insert(k); continue; }
        if (i + 1 >= a.size()) throw Error(ErrorCode::InvalidArgument, "option --" + k + " needs a value");
        kv_[k] = a[++i];
      } else {
        pos_.push_back(s);
      }
    }
  }
  bool flag(const std::string& k) const { return flags_.count(k) > 0; }
  bool has(const std::string& k) const { return kv_.count(k) > 0; }
  std::string get(const std::string& k, const std::string& def = {}) const {
    auto it = kv_.find(k);
    return it == kv_.end() ? def : it->second;
  }
  const std::vector<std::string>& positional() const { return pos_; }
  std::string pos(std::size_t i, const std::string& what) const {
    if (i >= pos_.size()) throw Error(ErrorCode::InvalidArgument, "missing argument: " + what);
    return pos_[i];
  }
  std::string out_dir() const {
    std::string o = get("o", get("output"));
    if (o.empty()) throw Error(ErrorCode::InvalidArgument, "give an output folder with -o <folder>");
    return o;
  }

 private:
  std::vector<std::string> pos_;
  std::map<std::string, std::string> kv_;
  std::set<std::string> flags_;
};

int cmd_validate(Args&);
int cmd_build(Args&);
int cmd_section(Args&);
int cmd_state(Args&);
int cmd_export(Args&);
int cmd_package(Args&);
int cmd_ensemble(Args&);
int cmd_fold_frame(Args&);
int cmd_plugins(Args&);
int cmd_bench(Args&);
int cmd_suite(Args&);
int cmd_verify(Args&);
int cmd_schema(Args&);
int cmd_version(Args&);

}  // namespace stratonix::cli
