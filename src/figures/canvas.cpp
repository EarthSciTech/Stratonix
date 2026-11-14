// Display list, text metrics, SVG writer, base64.
#include <cmath>
#include <cstdio>
#include <memory>
#include <sstream>

#include "stratonix/core/util.hpp"
#include "stratonix/figures/canvas.hpp"

namespace stratonix::fig {

Color Color::hex(const std::string& s) {
  if (s.size() != 7 || s[0] != '#') return rgb(128, 128, 128);
  auto h = [&](int i) { return static_cast<std::uint8_t>(std::stoi(s.substr(static_cast<std::size_t>(i), 2), nullptr, 16)); };
  return {h(1), h(3), h(5), 255};
}
std::string Color::svg() const {
  char buf[8];
  std::snprintf(buf, sizeof buf, "#%02x%02x%02x", r, g, b);
  return buf;
}

// ---- figure --------------------------------------------------------------------------------
void Figure::path(const std::vector<Pt>& p, bool closed, Color fill, Color stroke, double width, std::vector<double> dash) {
  if (p.size() < 2) return;
  Op o;
  o.kind = Op::Kind::Path;
  o.rings = {p};
  o.closed = closed;
  o.fill = fill;
  o.stroke = stroke;
  o.width = width;
  o.dash = std::move(dash);
  ops_.push_back(std::move(o));
}
void Figure::rings(const std::vector<std::vector<Pt>>& r, Color fill, Color stroke, double width) {
  Op o;
  o.kind = Op::Kind::Path;
  o.rings = r;
  o.closed = true;
  o.fill = fill;
  o.stroke = stroke;
  o.width = width;
  ops_.push_back(std::move(o));
}
void Figure::line(Pt a, Pt b, Color c, double width, std::vector<double> dash) { path({a, b}, false, none(), c, width, std::move(dash)); }
void Figure::rect(double x, double y, double w, double h, Color fill, Color stroke, double width) {
  path({{x, y}, {x + w, y}, {x + w, y + h}, {x, y + h}}, true, fill, stroke, width);
}
void Figure::text(Pt at, const std::string& s, double size, Anchor a, bool bold, Color c, double rotate, bool italic) {
  Op o;
  o.kind = Op::Kind::Text;
  o.at = at;
  o.text = s;
  o.size = size;
  o.anchor = a;
  o.bold = bold;
  o.italic = italic;
  o.color = c;
  o.rotate = rotate;
  ops_.push_back(std::move(o));
}
void Figure::image(double x, double y, double w, double h, const Image* img) {
  Op o;
  o.kind = Op::Kind::Image;
  o.x = x; o.y = y; o.w = w; o.h = h;
  o.image = img;
  ops_.push_back(std::move(o));
}
void Figure::clip(double x, double y, double w, double h) {
  Op o;
  o.kind = Op::Kind::ClipPush;
  o.x = x; o.y = y; o.w = w; o.h = h;
  ops_.push_back(std::move(o));
}
void Figure::unclip() {
  Op o;
  o.kind = Op::Kind::ClipPop;
  ops_.push_back(std::move(o));
}
const Image* Figure::own(Image img) {
  images_.push_back(std::make_unique<Image>(std::move(img)));
  return images_.back().get();
}

// ---- text metrics (Helvetica-compatible advance widths, 1/1000 em) --------------------------
double text_width(const std::string& s, double size, bool bold) {
  static const short w[95] = {278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278,  // ' '..'/'
                              556, 556, 556, 556, 556, 556, 556, 556, 556, 556,                                // 0-9
                              278, 278, 584, 584, 584, 556, 1015,                                              // :..@
                              667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778, 667,  // A-P
                              778, 722, 667, 611, 722, 667, 944, 667, 667, 611,                                // Q-Z
                              278, 278, 278, 469, 556, 333,                                                    // [..`
                              556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556, 556,  // a-p
                              556, 333, 500, 278, 556, 500, 722, 500, 500, 500,                                // q-z
                              334, 260, 334, 584};
  double t = 0;
  for (unsigned char c : s) t += (c >= 32 && c < 127) ? w[c - 32] : 556;
  return t / 1000.0 * size * (bold ? 1.07 : 1.0);
}

namespace {
static const short kDejaVu[95] = {318, 401, 460, 838, 636, 950, 780, 275, 390, 390, 500, 838, 318, 361, 318, 337, 636, 636, 636, 636, 636, 636, 636, 636, 636, 636, 337, 337, 838, 838, 838, 531, 1000, 684, 686, 698, 770, 632, 575, 775, 752, 295, 295, 656, 557, 863, 748, 787, 603, 787, 695, 635, 611, 732, 684, 989, 685, 611, 685, 390, 337, 390, 838, 500, 500, 613, 635, 550, 635, 615, 352, 635, 634, 278, 278, 579, 278, 974, 634, 612, 635, 635, 411, 521, 392, 634, 592, 818, 592, 592, 525, 636, 337, 636, 838};
static const short kDejaVuBold[95] = {348, 456, 521, 838, 696, 1002, 872, 306, 457, 457, 523, 838, 380, 415, 380, 365, 696, 696, 696, 696, 696, 696, 696, 696, 696, 696, 400, 400, 838, 838, 838, 580, 1000, 774, 762, 734, 830, 683, 683, 821, 837, 372, 372, 775, 637, 995, 837, 850, 733, 850, 770, 720, 682, 812, 774, 1103, 771, 724, 725, 457, 365, 457, 838, 500, 500, 675, 716, 593, 716, 678, 435, 716, 712, 343, 343, 665, 343, 1042, 712, 687, 716, 716, 493, 595, 478, 712, 652, 924, 645, 652, 582, 712, 365, 712, 838};
static const struct { unsigned cp; short w, wb; } kDejaVuExtra[] = {{229, 613, 675}, {248, 612, 687}, {230, 982, 1048}, {197, 684, 774}, {216, 787, 850}, {198, 974, 1085}, {233, 615, 678}, {232, 615, 678}, {252, 634, 712}, {246, 612, 687}, {228, 613, 675}, {214, 787, 850}, {196, 684, 774}, {220, 732, 812}, {8226, 590, 639}, {8211, 500, 500}, {8212, 1000, 1000}, {8722, 838, 838}, {8594, 838, 838}, {215, 838, 838}, {176, 500, 500}, {177, 838, 838}, {178, 401, 438}, {179, 401, 438}, {181, 636, 736}, {8217, 318, 380}, {8216, 318, 380}, {8220, 518, 657}, {8221, 518, 657}, {8230, 1000, 1000}, {183, 318, 380}, {8804, 838, 838}, {8805, 838, 838}};

/// Next code point of a UTF-8 string (invalid bytes are taken one by one).
unsigned next_cp(const std::string& s, std::size_t& i) {
  unsigned char c = static_cast<unsigned char>(s[i++]);
  if (c < 0x80) return c;
  int n = c >= 0xF0 ? 3 : c >= 0xE0 ? 2 : c >= 0xC0 ? 1 : 0;
  unsigned cp = n == 3 ? (c & 0x07u) : n == 2 ? (c & 0x0Fu) : (c & 0x1Fu);
  for (int k = 0; k < n && i < s.size(); ++k) cp = (cp << 6) | (static_cast<unsigned char>(s[i++]) & 0x3Fu);
  return cp;
}
}  // namespace

double text_width(const std::string& s, double size, bool bold, const std::string& family) {
  if (family.empty()) return text_width(s, size, bold);
  double t = 0;
  for (std::size_t i = 0; i < s.size();) {
    unsigned cp = next_cp(s, i);
    double w = 600;
    if (cp >= 32 && cp < 127) w = bold ? kDejaVuBold[cp - 32] : kDejaVu[cp - 32];
    else
      for (auto& e : kDejaVuExtra)
        if (e.cp == cp) { w = bold ? e.wb : e.w; break; }
    t += w;
  }
  return t / 1000.0 * size;
}

// ---- base64 --------------------------------------------------------------------------------
std::string base64(const std::string& d) {
  static const char* t = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  std::string o;
  o.reserve((d.size() + 2) / 3 * 4);
  std::size_t i = 0;
  for (; i + 2 < d.size(); i += 3) {
    unsigned v = (static_cast<unsigned char>(d[i]) << 16) | (static_cast<unsigned char>(d[i + 1]) << 8) | static_cast<unsigned char>(d[i + 2]);
    o += t[(v >> 18) & 63]; o += t[(v >> 12) & 63]; o += t[(v >> 6) & 63]; o += t[v & 63];
  }
  if (i < d.size()) {
    unsigned v = static_cast<unsigned char>(d[i]) << 16;
    if (i + 1 < d.size()) v |= static_cast<unsigned char>(d[i + 1]) << 8;
    o += t[(v >> 18) & 63]; o += t[(v >> 12) & 63];
    o += i + 1 < d.size() ? t[(v >> 6) & 63] : '=';
    o += '=';
  }
  return o;
}

// ---- SVG -----------------------------------------------------------------------------------
static std::string num(double v) { return fmt_fixed(v, 2); }

std::string to_svg(const Figure& f, const std::string& title) {
  std::ostringstream s;
  s << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" version=\"1.1\" width=\""
    << num(f.width()) << "pt\" height=\"" << num(f.height()) << "pt\" viewBox=\"0 0 " << num(f.width()) << " " << num(f.height())
    << "\">\n";
  if (!title.empty()) s << "<title>" << xml_escape(title) << "</title>\n";
  s << "<rect x=\"0\" y=\"0\" width=\"" << num(f.width()) << "\" height=\"" << num(f.height()) << "\" fill=\"#ffffff\"/>\n";
  int clip_id = 0, depth = 0;
  for (auto& o : f.ops()) {
    switch (o.kind) {
      case Op::Kind::Path: {
        s << "<path d=\"";
        for (auto& r : o.rings) {
          for (std::size_t i = 0; i < r.size(); ++i) s << (i ? " L" : "M") << num(r[i].x) << " " << num(r[i].y);
          if (o.closed) s << " Z ";
        }
        s << "\" fill=\"" << (o.fill.none() ? "none" : o.fill.svg()) << "\"";
        if (!o.fill.none() && o.fill.a < 255) s << " fill-opacity=\"" << fmt_fixed(o.fill.a / 255.0, 3) << "\"";
        if (o.rings.size() > 1) s << " fill-rule=\"evenodd\"";
        if (!o.stroke.none()) {
          s << " stroke=\"" << o.stroke.svg() << "\" stroke-width=\"" << num(o.width) << "\" stroke-linejoin=\"round\"";
          if (o.stroke.a < 255) s << " stroke-opacity=\"" << fmt_fixed(o.stroke.a / 255.0, 3) << "\"";
          if (!o.dash.empty()) {
            s << " stroke-dasharray=\"";
            for (std::size_t i = 0; i < o.dash.size(); ++i) s << (i ? "," : "") << num(o.dash[i]);
            s << "\"";
          }
        } else {
          s << " stroke=\"none\"";
        }
        s << "/>\n";
        break;
      }
      case Op::Kind::Text: {
        s << "<text x=\"" << num(o.at.x) << "\" y=\"" << num(o.at.y) << "\" font-family=\""
          << (f.font_family().empty() ? std::string("Arial, 'Liberation Sans', Helvetica, sans-serif")
                                      : "'" + f.font_family() + "', 'Bitstream Vera Sans', Verdana, sans-serif")
          << "\" font-size=\""
          << num(o.size) << "\"";
        if (o.bold) s << " font-weight=\"bold\"";
        if (o.italic) s << " font-style=\"italic\"";
        if (o.anchor != Anchor::Start) s << " text-anchor=\"" << (o.anchor == Anchor::Middle ? "middle" : "end") << "\"";
        if (o.rotate != 0) s << " transform=\"rotate(" << num(o.rotate) << " " << num(o.at.x) << " " << num(o.at.y) << ")\"";
        s << " fill=\"" << o.color.svg() << "\">" << xml_escape(o.text) << "</text>\n";
        break;
      }
      case Op::Kind::Image:
        s << "<image x=\"" << num(o.x) << "\" y=\"" << num(o.y) << "\" width=\"" << num(o.w) << "\" height=\"" << num(o.h)
          << "\" preserveAspectRatio=\"none\" xlink:href=\"data:image/png;base64," << base64(encode_png_rgba(*o.image)) << "\"/>\n";
        break;
      case Op::Kind::ClipPush:
        ++clip_id;
        ++depth;
        s << "<clipPath id=\"c" << clip_id << "\"><rect x=\"" << num(o.x) << "\" y=\"" << num(o.y) << "\" width=\"" << num(o.w)
          << "\" height=\"" << num(o.h) << "\"/></clipPath>\n<g clip-path=\"url(#c" << clip_id << ")\">\n";
        break;
      case Op::Kind::ClipPop:
        if (depth > 0) { s << "</g>\n"; --depth; }
        break;
    }
  }
  while (depth-- > 0) s << "</g>\n";
  s << "</svg>\n";
  return s.str();
}

}  // namespace stratonix::fig
