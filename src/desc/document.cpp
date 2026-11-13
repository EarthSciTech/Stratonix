// Document tree: YAML/JSON reading through libyaml, JSON writing, JSON Schema validation.
#include "stratonix/desc/document.hpp"

#include <yaml.h>

#include <cmath>
#include <map>
#include <regex>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"

namespace stratonix::desc {

const Node* Node::get(std::string_view key) const {
  if (type != Type::Map) return nullptr;
  for (auto& [k, v] : map)
    if (k == key) return &v;
  return nullptr;
}
Node* Node::get(std::string_view key) {
  if (type != Type::Map) return nullptr;
  for (auto& [k, v] : map)
    if (k == key) return &v;
  return nullptr;
}
Node& Node::set(const std::string& key, Node v) {
  if (type != Type::Map) { type = Type::Map; map.clear(); }
  for (auto& [k, x] : map)
    if (k == key) { x = std::move(v); return x; }
  map.emplace_back(key, std::move(v));
  return map.back().second;
}
Node Node::make_string(std::string s) { Node n; n.type = Type::String; n.str = std::move(s); return n; }
Node Node::make_number(double v) {
  Node n; n.type = Type::Number; n.number = v; n.integer = std::floor(v) == v && std::fabs(v) < 1e15;
  n.str = n.integer ? std::to_string(static_cast<long long>(v)) : fmt_shortest(v);
  return n;
}
Node Node::make_bool(bool b) { Node n; n.type = Type::Bool; n.boolean = b; n.str = b ? "true" : "false"; return n; }
Node Node::make_map() { Node n; n.type = Type::Map; return n; }
Node Node::make_seq() { Node n; n.type = Type::Seq; return n; }

// ---- YAML reading --------------------------------------------------------------------------
namespace {

Node resolve_plain(const std::string& s) {
  Node n;
  n.str = s;
  if (s.empty() || s == "~" || s == "null" || s == "Null" || s == "NULL") { n.type = Node::Type::Null; return n; }
  if (s == "true" || s == "True" || s == "TRUE") { n.type = Node::Type::Bool; n.boolean = true; return n; }
  if (s == "false" || s == "False" || s == "FALSE") { n.type = Node::Type::Bool; n.boolean = false; return n; }
  static const std::regex num(R"(^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$)");
  if (std::regex_match(s, num)) {
    n.type = Node::Type::Number;
    n.number = parse_double(s);
    n.integer = s.find_first_of(".eE") == std::string::npos;
    return n;
  }
  n.type = Node::Type::String;
  return n;
}

struct YamlParser {
  yaml_parser_t p{};
  std::string origin;
  YamlParser(std::string_view text, std::string o) : origin(std::move(o)) {
    yaml_parser_initialize(&p);
    yaml_parser_set_input_string(&p, reinterpret_cast<const unsigned char*>(text.data()), text.size());
  }
  ~YamlParser() { yaml_parser_delete(&p); }
  void next(yaml_event_t& e) {
    if (!yaml_parser_parse(&p, &e)) {
      std::string msg = p.problem ? p.problem : "YAML syntax error";
      throw Error(ErrorCode::Parse, origin + ":" + std::to_string(p.problem_mark.line + 1) + ":" +
                                        std::to_string(p.problem_mark.column + 1) + ": " + msg);
    }
  }
  Node parse_node(yaml_event_t& e) {
    Node n;
    n.line = static_cast<int>(e.start_mark.line) + 1;
    n.col = static_cast<int>(e.start_mark.column) + 1;
    switch (e.type) {
      case YAML_SCALAR_EVENT: {
        std::string v(reinterpret_cast<char*>(e.data.scalar.value), e.data.scalar.length);
        if (e.data.scalar.style == YAML_PLAIN_SCALAR_STYLE) {
          int l = n.line, c = n.col;
          n = resolve_plain(v);
          n.line = l; n.col = c;
        } else {
          n.type = Node::Type::String;
          n.str = v;
        }
        yaml_event_delete(&e);
        return n;
      }
      case YAML_SEQUENCE_START_EVENT: {
        yaml_event_delete(&e);
        n.type = Node::Type::Seq;
        for (;;) {
          yaml_event_t c;
          next(c);
          if (c.type == YAML_SEQUENCE_END_EVENT) { yaml_event_delete(&c); break; }
          n.seq.push_back(parse_node(c));
        }
        return n;
      }
      case YAML_MAPPING_START_EVENT: {
        yaml_event_delete(&e);
        n.type = Node::Type::Map;
        for (;;) {
          yaml_event_t k;
          next(k);
          if (k.type == YAML_MAPPING_END_EVENT) { yaml_event_delete(&k); break; }
          if (k.type != YAML_SCALAR_EVENT)
            throw Error(ErrorCode::Parse, origin + ":" + std::to_string(k.start_mark.line + 1) +
                                              ": mapping keys must be plain scalars");
          std::string key(reinterpret_cast<char*>(k.data.scalar.value), k.data.scalar.length);
          int kl = static_cast<int>(k.start_mark.line) + 1;
          yaml_event_delete(&k);
          if (n.get(key))
            throw Error(ErrorCode::Parse, origin + ":" + std::to_string(kl) + ": duplicate key '" + key + "'");
          yaml_event_t v;
          next(v);
          n.map.emplace_back(key, parse_node(v));
        }
        return n;
      }
      case YAML_ALIAS_EVENT:
        yaml_event_delete(&e);
        throw Error(ErrorCode::Parse, origin + ":" + std::to_string(n.line) +
                                          ": YAML aliases are not supported in model descriptions");
      default:
        yaml_event_delete(&e);
        throw Error(ErrorCode::Parse, origin + ": unexpected YAML event");
    }
  }
  Node parse() {
    yaml_event_t e;
    next(e);  // stream start
    yaml_event_delete(&e);
    next(e);
    if (e.type == YAML_STREAM_END_EVENT) { yaml_event_delete(&e); return Node{}; }
    yaml_event_delete(&e);  // document start
    yaml_event_t r;
    next(r);
    Node root = parse_node(r);
    next(e);  // document end
    yaml_event_delete(&e);
    next(e);
    bool more = e.type == YAML_DOCUMENT_START_EVENT;
    yaml_event_delete(&e);
    if (more) throw Error(ErrorCode::Parse, origin + ": only one YAML document per file is allowed");
    return root;
  }
};
}  // namespace

Node parse_document(std::string_view text, const std::string& origin) {
  // Strip a UTF-8 byte-order mark if present.
  if (text.size() >= 3 && static_cast<unsigned char>(text[0]) == 0xEF &&
      static_cast<unsigned char>(text[1]) == 0xBB && static_cast<unsigned char>(text[2]) == 0xBF)
    text.remove_prefix(3);
  YamlParser p(text, origin);
  return p.parse();
}

// ---- JSON writing --------------------------------------------------------------------------
namespace {
void write_json(const Node& n, std::string& out, int indent, int level) {
  std::string pad(static_cast<std::size_t>(indent * (level + 1)), ' ');
  std::string pad0(static_cast<std::size_t>(indent * level), ' ');
  switch (n.type) {
    case Node::Type::Null: out += "null"; break;
    case Node::Type::Bool: out += n.boolean ? "true" : "false"; break;
    case Node::Type::Number:
      out += n.integer ? std::to_string(static_cast<long long>(n.number)) : fmt_shortest(n.number);
      break;
    case Node::Type::String: out += "\"" + json_escape(n.str) + "\""; break;
    case Node::Type::Seq: {
      if (n.seq.empty()) { out += "[]"; break; }
      bool scalars = true;
      for (auto& c : n.seq)
        if (c.is_map() || c.is_seq()) scalars = false;
      if (scalars) {
        out += "[";
        for (std::size_t i = 0; i < n.seq.size(); ++i) {
          if (i) out += ", ";
          write_json(n.seq[i], out, indent, level + 1);
        }
        out += "]";
        break;
      }
      out += "[\n";
      for (std::size_t i = 0; i < n.seq.size(); ++i) {
        out += pad;
        write_json(n.seq[i], out, indent, level + 1);
        out += i + 1 < n.seq.size() ? ",\n" : "\n";
      }
      out += pad0 + "]";
      break;
    }
    case Node::Type::Map: {
      if (n.map.empty()) { out += "{}"; break; }
      out += "{\n";
      for (std::size_t i = 0; i < n.map.size(); ++i) {
        out += pad + "\"" + json_escape(n.map[i].first) + "\": ";
        write_json(n.map[i].second, out, indent, level + 1);
        out += i + 1 < n.map.size() ? ",\n" : "\n";
      }
      out += pad0 + "}";
      break;
    }
  }
}
}  // namespace

std::string to_json(const Node& n, int indent) {
  std::string out;
  write_json(n, out, indent, 0);
  out += "\n";
  return out;
}

// ---- JSON Schema validation ----------------------------------------------------------------
namespace {
const char* type_name(Node::Type t) {
  switch (t) {
    case Node::Type::Null: return "null";
    case Node::Type::Bool: return "boolean";
    case Node::Type::Number: return "number";
    case Node::Type::String: return "string";
    case Node::Type::Seq: return "array";
    default: return "object";
  }
}
bool type_matches(const Node& v, const std::string& t) {
  if (t == "object") return v.is_map();
  if (t == "array") return v.is_seq();
  if (t == "string") return v.is_string();
  if (t == "number") return v.is_number();
  if (t == "integer") return v.is_number() && std::floor(v.number) == v.number;
  if (t == "boolean") return v.type == Node::Type::Bool;
  if (t == "null") return v.is_null();
  return false;
}
bool equal(const Node& a, const Node& b) {
  if (a.type != b.type) return false;
  switch (a.type) {
    case Node::Type::Null: return true;
    case Node::Type::Bool: return a.boolean == b.boolean;
    case Node::Type::Number: return a.number == b.number;
    case Node::Type::String: return a.str == b.str;
    case Node::Type::Seq:
      if (a.seq.size() != b.seq.size()) return false;
      for (std::size_t i = 0; i < a.seq.size(); ++i)
        if (!equal(a.seq[i], b.seq[i])) return false;
      return true;
    case Node::Type::Map:
      if (a.map.size() != b.map.size()) return false;
      for (auto& [k, v] : a.map) {
        auto* o = b.get(k);
        if (!o || !equal(v, *o)) return false;
      }
      return true;
  }
  return false;
}
std::string scalar_text(const Node& n) {
  if (n.is_string()) return "'" + n.str + "'";
  if (n.is_map()) return "a mapping";
  if (n.is_seq()) return "a list";
  return n.str;
}

struct Validator {
  const Node& root;
  const RefHook& hook;
  std::map<std::string, std::regex> regex_cache;

  const std::regex& rx(const std::string& p) {
    auto it = regex_cache.find(p);
    if (it == regex_cache.end()) it = regex_cache.emplace(p, std::regex(p, std::regex::ECMAScript)).first;
    return it->second;
  }
  const Node* deref(const std::string& ref) {
    if (!starts_with(ref, "#/")) return nullptr;
    const Node* n = &root;
    for (auto& part : split(ref.substr(2), '/')) {
      n = n->get(part);
      if (!n) return nullptr;
    }
    return n;
  }
  void run(const Node& v, const Node& s, const std::string& path, std::vector<Issue>& out) {
    if (!s.is_map()) return;
    auto add = [&](std::string rule, std::string msg) {
      out.push_back({path.empty() ? "(root)" : path, std::move(rule), std::move(msg), true, v.line});
    };
    if (auto* r = s.get("$ref")) {
      const Node* t = deref(r->str);
      if (!t) { add("schema:ref", "schema reference not found: " + r->str); return; }
      std::vector<Issue> sub;
      run(v, *t, path, sub);
      if (!sub.empty() && hook) {
        std::string name = r->str.substr(r->str.rfind('/') + 1);
        if (auto iss = hook(name, v, path)) { out.push_back(*iss); return; }
      }
      out.insert(out.end(), sub.begin(), sub.end());
      return;
    }
    if (auto* t = s.get("type")) {
      bool ok = false;
      if (t->is_string()) ok = type_matches(v, t->str);
      else for (auto& x : t->seq) ok = ok || type_matches(v, x.str);
      if (!ok) {
        add("schema:type", std::string("expected ") + (t->is_string() ? t->str : "one of several types") +
                               ", found " + type_name(v.type));
        return;
      }
    }
    if (auto* c = s.get("const"))
      if (!equal(v, *c)) add("schema:const", "must be " + scalar_text(*c) + ", found " + scalar_text(v));
    if (auto* e = s.get("enum")) {
      bool ok = false;
      std::string opts;
      for (auto& x : e->seq) {
        ok = ok || equal(v, x);
        opts += (opts.empty() ? "" : ", ") + x.str;
      }
      if (!ok) add("schema:enum", "must be one of [" + opts + "], found " + scalar_text(v));
    }
    if (v.is_string()) {
      if (auto* p = s.get("pattern"))
        if (!std::regex_search(v.str, rx(p->str))) add("schema:pattern", "text '" + v.str + "' does not have the required form");
      if (auto* m = s.get("minLength"))
        if (v.str.size() < static_cast<std::size_t>(m->number)) add("schema:minLength", "must not be empty");
    }
    if (v.is_number())
      if (auto* m = s.get("minimum"))
        if (v.number < m->number) add("schema:minimum", "must be at least " + m->str);
    if (v.is_map()) {
      if (auto* req = s.get("required"))
        for (auto& r : req->seq)
          if (!v.get(r.str))
            out.push_back({(path.empty() ? "" : path + ".") + r.str, "schema:required",
                           "required field is missing", true, v.line});
      if (auto* mp = s.get("minProperties"))
        if (v.map.size() < static_cast<std::size_t>(mp->number)) add("schema:minProperties", "must have at least " + mp->str + " entries");
      const Node* props = s.get("properties");
      const Node* addl = s.get("additionalProperties");
      for (auto& [k, cv] : v.map) {
        std::string cp = path.empty() ? k : path + "." + k;
        const Node* ps = props ? props->get(k) : nullptr;
        if (ps) run(cv, *ps, cp, out);
        else if (addl) {
          if (addl->type == Node::Type::Bool && !addl->boolean)
            out.push_back({cp, "schema:additionalProperties", "unknown field (not allowed here)", true, cv.line});
          else if (addl->is_map()) run(cv, *addl, cp, out);
        }
      }
    }
    if (v.is_seq()) {
      if (auto* mi = s.get("minItems"))
        if (v.seq.size() < static_cast<std::size_t>(mi->number)) add("schema:minItems", "must have at least " + mi->str + " item(s)");
      if (auto* it = s.get("items"))
        for (std::size_t i = 0; i < v.seq.size(); ++i) run(v.seq[i], *it, path + "[" + std::to_string(i) + "]", out);
    }
    if (auto* one = s.get("oneOf")) {
      int matches = 0;
      std::vector<Issue> best;
      for (auto& alt : one->seq) {
        std::vector<Issue> sub;
        run(v, alt, path, sub);
        if (sub.empty()) ++matches;
        else if (best.empty() || sub.size() < best.size()) best = sub;
      }
      if (matches != 1) {
        if (matches == 0) out.insert(out.end(), best.begin(), best.end());
        else add("schema:oneOf", "matches more than one allowed form");
      }
    }
    if (auto* any = s.get("anyOf")) {
      bool ok = false;
      for (auto& alt : any->seq) {
        std::vector<Issue> sub;
        run(v, alt, path, sub);
        ok = ok || sub.empty();
      }
      if (!ok) add("schema:anyOf", "does not match any allowed form");
    }
  }
};
}  // namespace

std::vector<Issue> validate_schema(const Node& doc, const Node& schema, const RefHook& hook) {
  Validator v{schema, hook, {}};
  std::vector<Issue> out;
  v.run(doc, schema, "", out);
  return out;
}

}  // namespace stratonix::desc
