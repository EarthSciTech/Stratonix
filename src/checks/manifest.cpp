#include "stratonix/checks/manifest.hpp"

#include <algorithm>
#include <chrono>
#include <ctime>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/registry/registry.hpp"

namespace stratonix {

std::string utc_timestamp() {
  auto now = std::chrono::system_clock::now();
  std::time_t t = std::chrono::system_clock::to_time_t(now);
  std::tm tm{};
#if defined(_WIN32)
  gmtime_s(&tm, &t);
#else
  gmtime_r(&t, &tm);
#endif
  char buf[32];
  std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", &tm);
  return buf;
}

Manifest Manifest::from_description(const desc::Description& d, const std::string& command) {
  Manifest m;
  m.command = command;
  m.model_name = d.model.name;
  m.model_version = d.model.version;
  m.evidence_basis = d.model.status;
  m.description_path = d.file.filename().generic_string();
  m.description_sha256 = d.file_sha256;
  m.resolved_description = d.resolved;
  for (auto& in : d.inputs) {
    ManifestFile f;
    f.path = in.path;
    f.role = "input:" + in.role;
    f.sharing = in.sharing;
    std::error_code ec;
    if (std::filesystem::exists(in.resolved, ec)) {
      f.sha256 = sha256_file(in.resolved);
      f.bytes = std::filesystem::file_size(in.resolved, ec);
    } else {
      f.sha256 = "missing";
    }
    m.inputs.push_back(f);
  }
  for (auto* e : registry::Registry::global().all())
    m.plugins.push_back(std::string(registry::kind_name(e->kind)) + "/" + e->name + "@" + e->version);
  m.timestamp_utc = utc_timestamp();
  return m;
}

void Manifest::add_outputs(const std::filesystem::path& root, const std::string& role_default) {
  std::vector<std::filesystem::path> files;
  for (auto& e : std::filesystem::recursive_directory_iterator(root))
    if (e.is_regular_file() && e.path().filename() != "manifest.json") files.push_back(e.path());
  std::sort(files.begin(), files.end(), [&](auto& a, auto& b) {
    return portable_path(a, root) < portable_path(b, root);
  });
  for (auto& p : files) {
    ManifestFile f;
    f.path = portable_path(p, root);
    f.sha256 = sha256_file(p);
    f.bytes = std::filesystem::file_size(p);
    f.role = role_default;
    outputs.push_back(f);
  }
}

std::string Manifest::to_json() const {
  using desc::Node;
  Node r = Node::make_map();
  r.set("tool", Node::make_string(tool));
  r.set("command", Node::make_string(command));
  auto& bi = build_info();
  Node b = Node::make_map();
  b.set("version", Node::make_string(bi.version));
  b.set("git_commit", Node::make_string(bi.git_commit));
  b.set("compiler", Node::make_string(bi.compiler));
  b.set("platform", Node::make_string(bi.platform));
  b.set("build_type", Node::make_string(bi.build_type));
  r.set("stratonix", b);
  r.set("timestamp_utc", Node::make_string(timestamp_utc));
  Node mm = Node::make_map();
  mm.set("name", Node::make_string(model_name));
  mm.set("version", Node::make_string(model_version));
  mm.set("evidence_basis", Node::make_string(evidence_basis));
  r.set("model", mm);
  Node dd = Node::make_map();
  dd.set("path", Node::make_string(description_path));
  dd.set("sha256", Node::make_string(description_sha256));
  dd.set("resolved", resolved_description);
  r.set("description", dd);
  auto files = [](const std::vector<ManifestFile>& v) {
    Node s = Node::make_seq();
    for (auto& f : v) {
      Node x = Node::make_map();
      x.set("path", Node::make_string(f.path));
      x.set("sha256", Node::make_string(f.sha256));
      x.set("bytes", Node::make_number(static_cast<double>(f.bytes)));
      x.set("role", Node::make_string(f.role));
      if (!f.sharing.empty()) x.set("sharing", Node::make_string(f.sharing));
      s.seq.push_back(x);
    }
    return s;
  };
  r.set("inputs", files(inputs));
  Node pl = Node::make_seq();
  for (auto& p : plugins) pl.seq.push_back(Node::make_string(p));
  r.set("plugins", pl);
  Node ex = Node::make_map();
  for (auto& [k, v] : extra) ex.set(k, Node::make_string(v));
  r.set("notes", ex);
  r.set("outputs", files(outputs));
  return desc::to_json(r);
}

void Manifest::write(const std::filesystem::path& file) const { write_file(file, to_json()); }

std::vector<std::string> verify_manifest(const std::filesystem::path& manifest_file) {
  std::vector<std::string> problems;
  auto doc = desc::parse_document(read_file(manifest_file), manifest_file.generic_string());
  auto root = manifest_file.parent_path();
  const desc::Node* outs = doc.get("outputs");
  if (!outs) return {"manifest has no outputs list"};
  std::size_t listed = 0;
  for (auto& o : outs->seq) {
    ++listed;
    auto* p = o.get("path");
    auto* h = o.get("sha256");
    if (!p || !h) { problems.push_back("entry without path or sha256"); continue; }
    auto f = root / p->str;
    if (!std::filesystem::exists(f)) { problems.push_back("missing file: " + p->str); continue; }
    if (sha256_file(f) != h->str) problems.push_back("checksum mismatch: " + p->str);
  }
  std::size_t present = 0;
  for (auto& e : std::filesystem::recursive_directory_iterator(root))
    if (e.is_regular_file() && e.path().filename() != "manifest.json") ++present;
  if (present != listed) problems.push_back("folder holds " + std::to_string(present) + " files, manifest lists " + std::to_string(listed));
  return problems;
}

}  // namespace stratonix
