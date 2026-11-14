/// @file canvas.hpp
/// @brief Device-independent figures: a display list of paths, text and images in points
/// (1/72 inch), written to SVG (text kept as text), rasterised (PNG, TIFF) and placed in PPTX.
/// The same display list feeds every format, so content is identical across formats and platforms.
#pragma once
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace stratonix::fig {

struct Color {
  std::uint8_t r = 0, g = 0, b = 0, a = 255;
  static Color hex(const std::string& s);  ///< "#rrggbb"
  std::string svg() const;                  ///< "#rrggbb"
  bool none() const { return a == 0; }
};
inline Color rgb(int r, int g, int b) { return {static_cast<std::uint8_t>(r), static_cast<std::uint8_t>(g), static_cast<std::uint8_t>(b), 255}; }
inline Color none() { return {0, 0, 0, 0}; }

struct Pt { double x = 0, y = 0; };

/// RGBA image (row-major, top row first, 4 bytes per pixel, straight alpha).
struct Image {
  int w = 0, h = 0;
  std::vector<std::uint8_t> rgba;
};

enum class Anchor { Start, Middle, End };

struct Op {
  enum class Kind { Path, Text, Image, ClipPush, ClipPop };
  Kind kind = Kind::Path;
  // path
  std::vector<std::vector<Pt>> rings;  ///< subpaths
  bool closed = false;
  Color fill = none(), stroke = none();
  double width = 0.5;
  std::vector<double> dash;
  // text
  Pt at;
  std::string text;
  double size = 8;
  bool bold = false, italic = false;
  Anchor anchor = Anchor::Start;
  double rotate = 0;  ///< degrees, clockwise
  Color color = rgb(0, 0, 0);
  // image / clip rectangle
  double x = 0, y = 0, w = 0, h = 0;
  const Image* image = nullptr;
};

/// A figure: size in points, origin top-left, y downwards.
class Figure {
 public:
  Figure(double width_pt, double height_pt) : w_(width_pt), h_(height_pt) {}
  double width() const { return w_; }
  double height() const { return h_; }
  const std::vector<Op>& ops() const { return ops_; }

  void path(const std::vector<Pt>& p, bool closed, Color fill, Color stroke, double width = 0.5,
            std::vector<double> dash = {});
  void rings(const std::vector<std::vector<Pt>>& r, Color fill, Color stroke, double width = 0.5);
  void line(Pt a, Pt b, Color c, double width = 0.5, std::vector<double> dash = {});
  void rect(double x, double y, double w, double h, Color fill, Color stroke, double width = 0.5);
  void text(Pt at, const std::string& s, double size, Anchor a = Anchor::Start, bool bold = false, Color c = rgb(0, 0, 0),
            double rotate = 0, bool italic = false);
  void image(double x, double y, double w, double h, const Image* img);
  void clip(double x, double y, double w, double h);
  void unclip();
  /// Keep images alive with the figure.
  const Image* own(Image img);
  /// Font family for all text of this figure: empty (Arial / Liberation Sans, the default) or
  /// "DejaVu Sans". Text widths for layout use the matching built-in metrics.
  void set_font_family(std::string f) { font_ = std::move(f); }
  const std::string& font_family() const { return font_; }

 private:
  double w_, h_;
  std::string font_;
  std::vector<Op> ops_;
  std::vector<std::unique_ptr<Image>> images_;
};

/// Text width in points from built-in sans-serif metrics (Helvetica/Arial-compatible advance
/// widths), so layout is identical on every platform.
double text_width(const std::string& s, double size, bool bold = false);
/// Text width for a font family ("" = the default metrics above, "DejaVu Sans" = DejaVu Sans
/// advance widths). Text is UTF-8; characters without a table entry count as 0.6 em.
double text_width(const std::string& s, double size, bool bold, const std::string& family);

/// SVG document (text kept as <text>, images embedded as PNG).
std::string to_svg(const Figure& f, const std::string& title = {});

/// Rasterise at `dpi` on a white background (RGB, 3 bytes per pixel). Requires Cairo.
struct Raster {
  int w = 0, h = 0;
  std::vector<std::uint8_t> rgb;
};
bool raster_available();
Raster rasterise(const Figure& f, double dpi);

/// PNG (RGB, 8 bits, pHYs set to dpi) and TIFF (RGB, LZW, 300 dpi, no alpha) encoders.
std::string encode_png(const Raster& r, double dpi);
std::string encode_png_rgba(const Image& im);
std::string encode_tiff_lzw(const Raster& r, double dpi, const std::string& software);
/// Decoders used by tests and by readers.
Raster decode_tiff_lzw_for_test(const std::string& tiff);

/// One 16:9 slide holding a picture, with title, caption and speaker notes.
struct Slide {
  std::string title, caption, notes;
  std::string png;       ///< PNG bytes of the figure
  int px_w = 0, px_h = 0;
};
std::string make_pptx(const std::vector<Slide>& slides, const std::string& deck_title);

/// Minimal deterministic zip writer (deflate, fixed timestamps).
class ZipWriter {
 public:
  void add(const std::string& name, const std::string& data, bool compress = true);
  std::string finish();

 private:
  struct Entry { std::string name; std::uint32_t crc, csize, usize, offset; std::uint16_t method; };
  std::string out_;
  std::vector<Entry> entries_;
};

std::string base64(const std::string& data);

}  // namespace stratonix::fig
