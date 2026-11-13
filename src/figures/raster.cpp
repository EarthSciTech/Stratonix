// Rasterisation (Cairo), PNG and TIFF (LZW) encoders.
#include <zlib.h>

#include <cmath>
#include <cstring>

#include "stratonix/core/error.hpp"
#include "stratonix/figures/canvas.hpp"

#ifdef STRATONIX_HAVE_CAIRO
#include <cairo.h>
#endif

namespace stratonix::fig {

// ---- PNG -----------------------------------------------------------------------------------
namespace {
void be32(std::string& s, std::uint32_t v) {
  s += static_cast<char>(v >> 24); s += static_cast<char>(v >> 16); s += static_cast<char>(v >> 8); s += static_cast<char>(v);
}
void chunk(std::string& out, const char* type, const std::string& data) {
  be32(out, static_cast<std::uint32_t>(data.size()));
  std::string td = std::string(type, 4) + data;
  out += td;
  be32(out, static_cast<std::uint32_t>(crc32(0, reinterpret_cast<const Bytef*>(td.data()), static_cast<uInt>(td.size()))));
}
std::string zcompress(const std::string& raw) {
  uLongf n = compressBound(static_cast<uLong>(raw.size()));
  std::string o(n, '\0');
  if (compress2(reinterpret_cast<Bytef*>(o.data()), &n, reinterpret_cast<const Bytef*>(raw.data()), static_cast<uLong>(raw.size()), 9) != Z_OK)
    throw Error(ErrorCode::Internal, "zlib compression failed");
  o.resize(n);
  return o;
}
std::string png(int w, int h, int channels, const std::uint8_t* px, double dpi) {
  std::string out = "\x89PNG\r\n\x1a\n";
  std::string ih;
  be32(ih, static_cast<std::uint32_t>(w));
  be32(ih, static_cast<std::uint32_t>(h));
  ih += static_cast<char>(8);
  ih += static_cast<char>(channels == 4 ? 6 : 2);
  ih += std::string(3, '\0');
  chunk(out, "IHDR", ih);
  if (dpi > 0) {
    std::string ph;
    auto ppm = static_cast<std::uint32_t>(std::lround(dpi / 0.0254));
    be32(ph, ppm); be32(ph, ppm);
    ph += static_cast<char>(1);
    chunk(out, "pHYs", ph);
  }
  std::string raw;
  raw.reserve(static_cast<std::size_t>(h) * (static_cast<std::size_t>(w) * channels + 1));
  for (int y = 0; y < h; ++y) {
    raw += '\0';
    raw.append(reinterpret_cast<const char*>(px + static_cast<std::size_t>(y) * w * channels), static_cast<std::size_t>(w) * channels);
  }
  chunk(out, "IDAT", zcompress(raw));
  chunk(out, "IEND", "");
  return out;
}
}  // namespace

std::string encode_png(const Raster& r, double dpi) { return png(r.w, r.h, 3, r.rgb.data(), dpi); }
std::string encode_png_rgba(const Image& im) { return png(im.w, im.h, 4, im.rgba.data(), 0); }

// ---- TIFF (baseline RGB, LZW) --------------------------------------------------------------
namespace {
class LzwWriter {
 public:
  std::string out;
  void put(unsigned code, int width) {
    acc_ = (acc_ << width) | code;
    nbits_ += width;
    while (nbits_ >= 8) {
      out += static_cast<char>((acc_ >> (nbits_ - 8)) & 0xFF);
      nbits_ -= 8;
    }
  }
  void flush() {
    if (nbits_ > 0) out += static_cast<char>((acc_ << (8 - nbits_)) & 0xFF);
    nbits_ = 0;
    acc_ = 0;
  }

 private:
  std::uint64_t acc_ = 0;
  int nbits_ = 0;
};

std::string lzw_encode(const std::uint8_t* d, std::size_t n) {
  // TIFF LZW: codes 256 = Clear, 257 = EOI; widths 9..12 bits with "early change".
  LzwWriter wr;
  thread_local std::vector<int> child(4096 * 256, -1);  // dense trie: child[code * 256 + byte]
  std::vector<std::size_t> used;                         // entries to clear on reset
  int next = 258, width = 9;
  auto reset = [&] {
    for (auto u : used) child[u] = -1;
    used.clear();
    next = 258;
    width = 9;
  };
  wr.put(256, width);
  if (n == 0) { wr.put(257, width); wr.flush(); return wr.out; }
  int w = d[0];
  for (std::size_t i = 1; i < n; ++i) {
    int c = d[i];
    int k = child[static_cast<std::size_t>(w) * 256 + c];
    if (k >= 0) { w = k; continue; }
    wr.put(static_cast<unsigned>(w), width);
    child[static_cast<std::size_t>(w) * 256 + c] = next++;
    used.push_back(static_cast<std::size_t>(w) * 256 + c);
    if (next + 1 > (1 << width) && width < 12) ++width;  // early change
    if (next >= 4094) {
      wr.put(256, width);
      reset();
    }
    w = c;
  }
  wr.put(static_cast<unsigned>(w), width);
  ++next;
  if (next + 1 > (1 << width) && width < 12) ++width;
  wr.put(257, width);
  wr.flush();
  reset();
  return wr.out;
}

void le16(std::string& s, std::uint16_t v) { s += static_cast<char>(v & 0xFF); s += static_cast<char>(v >> 8); }
void le32(std::string& s, std::uint32_t v) { for (int i = 0; i < 4; ++i) s += static_cast<char>((v >> (8 * i)) & 0xFF); }
}  // namespace

std::string encode_tiff_lzw(const Raster& r, double dpi, const std::string& software) {
  const int rows_per_strip = 64;
  const int nstrips = (r.h + rows_per_strip - 1) / rows_per_strip;
  std::vector<std::string> strips;
  for (int s = 0; s < nstrips; ++s) {
    int y0 = s * rows_per_strip, y1 = std::min(r.h, y0 + rows_per_strip);
    strips.push_back(lzw_encode(r.rgb.data() + static_cast<std::size_t>(y0) * r.w * 3, static_cast<std::size_t>(y1 - y0) * r.w * 3));
  }
  // Layout: header (8) | strips | bits-per-sample (6) | resolutions (16) | software | offsets | counts | IFD
  std::string o = "II";
  le16(o, 42);
  le32(o, 0);  // IFD offset, patched below
  std::vector<std::uint32_t> offs, counts;
  for (auto& st : strips) {
    offs.push_back(static_cast<std::uint32_t>(o.size()));
    counts.push_back(static_cast<std::uint32_t>(st.size()));
    o += st;
    if (o.size() % 2) o += '\0';
  }
  std::uint32_t bps_off = static_cast<std::uint32_t>(o.size());
  le16(o, 8); le16(o, 8); le16(o, 8);
  std::uint32_t xres_off = static_cast<std::uint32_t>(o.size());
  le32(o, static_cast<std::uint32_t>(std::lround(dpi))); le32(o, 1);
  std::uint32_t yres_off = static_cast<std::uint32_t>(o.size());
  le32(o, static_cast<std::uint32_t>(std::lround(dpi))); le32(o, 1);
  std::uint32_t sw_off = static_cast<std::uint32_t>(o.size());
  o += software;
  o += '\0';
  if (o.size() % 2) o += '\0';
  std::uint32_t so_off = static_cast<std::uint32_t>(o.size());
  for (auto v : offs) le32(o, v);
  std::uint32_t sc_off = static_cast<std::uint32_t>(o.size());
  for (auto v : counts) le32(o, v);
  std::uint32_t ifd = static_cast<std::uint32_t>(o.size());
  std::memcpy(o.data() + 4, &ifd, 4);  // little-endian host assumed below; fix explicitly
  o[4] = static_cast<char>(ifd & 0xFF); o[5] = static_cast<char>((ifd >> 8) & 0xFF);
  o[6] = static_cast<char>((ifd >> 16) & 0xFF); o[7] = static_cast<char>((ifd >> 24) & 0xFF);
  struct Tag { std::uint16_t tag, type; std::uint32_t count, value; };
  std::vector<Tag> tags = {
      {256, 4, 1, static_cast<std::uint32_t>(r.w)},            // ImageWidth
      {257, 4, 1, static_cast<std::uint32_t>(r.h)},            // ImageLength
      {258, 3, 3, bps_off},                                    // BitsPerSample
      {259, 3, 1, 5},                                          // Compression = LZW
      {262, 3, 1, 2},                                          // Photometric = RGB
      {273, 4, static_cast<std::uint32_t>(nstrips), nstrips == 1 ? offs[0] : so_off},  // StripOffsets
      {277, 3, 1, 3},                                          // SamplesPerPixel
      {278, 4, 1, static_cast<std::uint32_t>(rows_per_strip)}, // RowsPerStrip
      {279, 4, static_cast<std::uint32_t>(nstrips), nstrips == 1 ? counts[0] : sc_off},  // StripByteCounts
      {282, 5, 1, xres_off},                                   // XResolution
      {283, 5, 1, yres_off},                                   // YResolution
      {284, 3, 1, 1},                                          // PlanarConfiguration = chunky
      {296, 3, 1, 2},                                          // ResolutionUnit = inch
      {305, 2, static_cast<std::uint32_t>(software.size() + 1), sw_off},  // Software
      {317, 3, 1, 1},                                          // Predictor = none
  };
  le16(o, static_cast<std::uint16_t>(tags.size()));
  for (auto& t : tags) {
    le16(o, t.tag);
    le16(o, t.type);
    le32(o, t.count);
    if (t.type == 3 && t.count == 1) { le16(o, static_cast<std::uint16_t>(t.value)); le16(o, 0); }
    else le32(o, t.value);
  }
  le32(o, 0);
  return o;
}

// ---- rasterisation -------------------------------------------------------------------------
bool raster_available() {
#ifdef STRATONIX_HAVE_CAIRO
  return true;
#else
  return false;
#endif
}

Raster rasterise(const Figure& f, double dpi) {
#ifndef STRATONIX_HAVE_CAIRO
  (void)f; (void)dpi;
  throw Error(ErrorCode::Unsupported, "raster figures need Cairo (build with STRATONIX_WITH_CAIRO)");
#else
  const double sc = dpi / 72.0;
  int W = static_cast<int>(std::lround(f.width() * sc)), H = static_cast<int>(std::lround(f.height() * sc));
  cairo_surface_t* surf = cairo_image_surface_create(CAIRO_FORMAT_RGB24, W, H);
  cairo_t* cr = cairo_create(surf);
  cairo_scale(cr, sc, sc);
  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_paint(cr);
  auto setc = [&](Color c) { cairo_set_source_rgba(cr, c.r / 255.0, c.g / 255.0, c.b / 255.0, c.a / 255.0); };
#if defined(_WIN32)
  const char* family = "Arial";
#else
  const char* family = "Liberation Sans";
#endif
  for (auto& o : f.ops()) {
    switch (o.kind) {
      case Op::Kind::Path: {
        cairo_new_path(cr);
        for (auto& r : o.rings) {
          for (std::size_t i = 0; i < r.size(); ++i) (i ? cairo_line_to : cairo_move_to)(cr, r[i].x, r[i].y);
          if (o.closed) cairo_close_path(cr);
        }
        cairo_set_fill_rule(cr, o.rings.size() > 1 ? CAIRO_FILL_RULE_EVEN_ODD : CAIRO_FILL_RULE_WINDING);
        if (!o.fill.none()) {
          setc(o.fill);
          if (o.stroke.none()) cairo_fill(cr); else cairo_fill_preserve(cr);
        }
        if (!o.stroke.none()) {
          setc(o.stroke);
          cairo_set_line_width(cr, o.width);
          cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
          cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT);
          cairo_set_dash(cr, o.dash.empty() ? nullptr : o.dash.data(), static_cast<int>(o.dash.size()), 0);
          cairo_stroke(cr);
        }
        cairo_new_path(cr);
        break;
      }
      case Op::Kind::Text: {
        cairo_save(cr);
        cairo_select_font_face(cr, family, o.italic ? CAIRO_FONT_SLANT_ITALIC : CAIRO_FONT_SLANT_NORMAL,
                               o.bold ? CAIRO_FONT_WEIGHT_BOLD : CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, o.size);
        cairo_translate(cr, o.at.x, o.at.y);
        if (o.rotate != 0) cairo_rotate(cr, o.rotate * 3.14159265358979323846 / 180.0);
        cairo_text_extents_t ext;
        cairo_text_extents(cr, o.text.c_str(), &ext);
        double dx = o.anchor == Anchor::Start ? 0 : o.anchor == Anchor::Middle ? -ext.x_advance / 2 : -ext.x_advance;
        cairo_move_to(cr, dx, 0);
        setc(o.color);
        cairo_show_text(cr, o.text.c_str());
        cairo_restore(cr);
        break;
      }
      case Op::Kind::Image: {
        const Image& im = *o.image;
        cairo_surface_t* is = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, im.w, im.h);
        cairo_surface_flush(is);
        unsigned char* data = cairo_image_surface_get_data(is);
        int stride = cairo_image_surface_get_stride(is);
        for (int y = 0; y < im.h; ++y)
          for (int x = 0; x < im.w; ++x) {
            const std::uint8_t* p = &im.rgba[(static_cast<std::size_t>(y) * im.w + x) * 4];
            std::uint32_t a = p[3];
            std::uint32_t v = (a << 24) | ((p[0] * a / 255) << 16) | ((p[1] * a / 255) << 8) | (p[2] * a / 255);
            std::memcpy(data + static_cast<std::size_t>(y) * stride + x * 4, &v, 4);
          }
        cairo_surface_mark_dirty(is);
        cairo_save(cr);
        cairo_translate(cr, o.x, o.y);
        cairo_scale(cr, o.w / im.w, o.h / im.h);
        cairo_set_source_surface(cr, is, 0, 0);
        cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_GOOD);
        cairo_paint(cr);
        cairo_restore(cr);
        cairo_surface_destroy(is);
        break;
      }
      case Op::Kind::ClipPush:
        cairo_save(cr);
        cairo_rectangle(cr, o.x, o.y, o.w, o.h);
        cairo_clip(cr);
        break;
      case Op::Kind::ClipPop:
        cairo_restore(cr);
        break;
    }
  }
  cairo_surface_flush(surf);
  Raster r;
  r.w = W;
  r.h = H;
  r.rgb.resize(static_cast<std::size_t>(W) * H * 3);
  unsigned char* data = cairo_image_surface_get_data(surf);
  int stride = cairo_image_surface_get_stride(surf);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      std::uint32_t v;
      std::memcpy(&v, data + static_cast<std::size_t>(y) * stride + x * 4, 4);
      std::size_t k = (static_cast<std::size_t>(y) * W + x) * 3;
      r.rgb[k] = static_cast<std::uint8_t>((v >> 16) & 0xFF);
      r.rgb[k + 1] = static_cast<std::uint8_t>((v >> 8) & 0xFF);
      r.rgb[k + 2] = static_cast<std::uint8_t>(v & 0xFF);
    }
  cairo_destroy(cr);
  cairo_surface_destroy(surf);
  return r;
#endif
}

// ---- TIFF decoder for tests (our own files only: RGB, LZW, chunky) ---------------------------
Raster decode_tiff_lzw_for_test(const std::string& t) {
  auto u16 = [&](std::size_t o) { return static_cast<unsigned>(static_cast<unsigned char>(t[o]) | (static_cast<unsigned char>(t[o + 1]) << 8)); };
  auto u32 = [&](std::size_t o) { return u16(o) | (u16(o + 2) << 16); };
  std::size_t ifd = u32(4);
  unsigned n = u16(ifd);
  unsigned W = 0, H = 0, rps = 0, ns = 0;
  std::size_t so = 0, sc = 0;
  for (unsigned i = 0; i < n; ++i) {
    std::size_t e = ifd + 2 + 12 * i;
    unsigned tag = u16(e), type = u16(e + 2), cnt = u32(e + 4);
    unsigned val = type == 3 ? u16(e + 8) : u32(e + 8);
    if (tag == 256) W = val;
    if (tag == 257) H = val;
    if (tag == 278) rps = val;
    if (tag == 273) { ns = cnt; so = cnt == 1 ? e + 8 : val; }
    if (tag == 279) sc = cnt == 1 ? e + 8 : val;
  }
  Raster r;
  r.w = static_cast<int>(W);
  r.h = static_cast<int>(H);
  for (unsigned s = 0; s < ns; ++s) {
    std::size_t off = u32(so + 4 * s), len = u32(sc + 4 * s);
    // LZW decode
    std::vector<std::string> dict(4096);
    for (int i = 0; i < 256; ++i) dict[static_cast<std::size_t>(i)] = std::string(1, static_cast<char>(i));
    int next = 258, width = 9;
    std::uint64_t acc = 0;
    int nb = 0;
    std::size_t pos = off;
    std::string prev, out;
    bool first = true;
    for (;;) {
      while (nb < width && pos < off + len) { acc = (acc << 8) | static_cast<unsigned char>(t[pos++]); nb += 8; }
      if (nb < width) break;
      unsigned code = static_cast<unsigned>((acc >> (nb - width)) & ((1u << width) - 1));
      nb -= width;
      if (code == 257) break;
      if (code == 256) { next = 258; width = 9; first = true; prev.clear(); continue; }
      std::string entry;
      if (code < static_cast<unsigned>(next) && !dict[code].empty()) entry = dict[code];
      else entry = prev + prev[0];
      out += entry;
      if (!first) {
        dict[static_cast<std::size_t>(next++)] = prev + entry[0];
        if (next + 1 >= (1 << width) && width < 12) ++width;  // TIFF early change
      }
      first = false;
      prev = entry;
    }
    r.rgb.insert(r.rgb.end(), out.begin(), out.end());
  }
  (void)rps;
  return r;
}

}  // namespace stratonix::fig
