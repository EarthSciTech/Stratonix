/// @file document.hpp
/// @brief Generic document tree read from YAML or JSON, a JSON writer, and a JSON Schema validator.
#pragma once
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace stratonix::desc {

/// Key/value pair of a mapping (a plain template so Node can contain it while incomplete).
template <class T>
struct KeyValue {
  std::string first;
  T second;
  KeyValue() = default;
  KeyValue(std::string k, T v) : first(std::move(k)), second(std::move(v)) {}
};

/// Node of a YAML/JSON document. Maps keep their order.
struct Node {
  enum class Type { Null, Bool, Number, String, Seq, Map };
  Type type = Type::Null;
  bool boolean = false;
  double number = 0;
  bool integer = false;
  std::string str;  ///< string value, or the literal text of a scalar
  std::vector<Node> seq;
  std::vector<KeyValue<Node>> map;
  int line = 0, col = 0;

  bool is_null() const { return type == Type::Null; }
  bool is_map() const { return type == Type::Map; }
  bool is_seq() const { return type == Type::Seq; }
  bool is_string() const { return type == Type::String; }
  bool is_number() const { return type == Type::Number; }
  const Node* get(std::string_view key) const;
  Node* get(std::string_view key);
  Node& set(const std::string& key, Node v);
  static Node make_string(std::string s);
  static Node make_number(double v);
  static Node make_bool(bool b);
  static Node make_map();
  static Node make_seq();
};

/// Parse YAML (1.1 block or flow style) or JSON text. Throws Error(Parse) with line and column.
Node parse_document(std::string_view text, const std::string& origin = "<text>");

/// Serialise to JSON with stable formatting (2-space indent, keys in document order).
std::string to_json(const Node& n, int indent = 2);

/// One validation finding.
struct Issue {
  std::string path;     ///< e.g. faults[1].throw.max
  std::string rule;     ///< rule identifier, e.g. "E1", "S2", "schema:required"
  std::string message;  ///< plain-language explanation
  bool error = true;    ///< false for warnings
  int line = 0;
};

/// Hook called when a `$ref` to a named definition fails; may return a more specific issue.
using RefHook = std::function<std::optional<Issue>(const std::string& def_name, const Node& value,
                                                   const std::string& path)>;

/// Validate a document against a JSON Schema (the subset used by Stratonix: type, enum, const,
/// properties, required, additionalProperties, items, minItems, minProperties, pattern,
/// minLength, minimum, oneOf, anyOf, local $ref).
std::vector<Issue> validate_schema(const Node& doc, const Node& schema, const RefHook& hook = {});

}  // namespace stratonix::desc
