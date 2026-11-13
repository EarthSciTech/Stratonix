// Deterministic zip writer and a minimal PowerPoint (PPTX) package: one 16:9 slide per figure
// with title, picture, caption and speaker notes.
#include <zlib.h>

#include <cstring>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/figures/canvas.hpp"

namespace stratonix::fig {

// ---- zip -----------------------------------------------------------------------------------
namespace {
void le16(std::string& s, std::uint32_t v) { s += static_cast<char>(v & 0xFF); s += static_cast<char>((v >> 8) & 0xFF); }
void le32(std::string& s, std::uint32_t v) { for (int i = 0; i < 4; ++i) s += static_cast<char>((v >> (8 * i)) & 0xFF); }
std::string raw_deflate(const std::string& in) {
  z_stream z{};
  if (deflateInit2(&z, 9, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) throw Error(ErrorCode::Internal, "deflateInit2 failed");
  std::string out(deflateBound(&z, static_cast<uLong>(in.size())) + 16, '\0');
  z.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(in.data()));
  z.avail_in = static_cast<uInt>(in.size());
  z.next_out = reinterpret_cast<Bytef*>(out.data());
  z.avail_out = static_cast<uInt>(out.size());
  int rc = deflate(&z, Z_FINISH);
  if (rc != Z_STREAM_END) { deflateEnd(&z); throw Error(ErrorCode::Internal, "deflate failed"); }
  out.resize(z.total_out);
  deflateEnd(&z);
  return out;
}
// Fixed DOS date/time (1980-01-01 00:00) so archives are byte-identical between runs.
constexpr std::uint16_t kDosTime = 0, kDosDate = (0 << 9) | (1 << 5) | 1;
}  // namespace

void ZipWriter::add(const std::string& name, const std::string& data, bool compress) {
  Entry e;
  e.name = name;
  e.crc = static_cast<std::uint32_t>(crc32(0, reinterpret_cast<const Bytef*>(data.data()), static_cast<uInt>(data.size())));
  e.usize = static_cast<std::uint32_t>(data.size());
  std::string body = compress ? raw_deflate(data) : data;
  e.method = compress ? 8 : 0;
  e.csize = static_cast<std::uint32_t>(body.size());
  e.offset = static_cast<std::uint32_t>(out_.size());
  le32(out_, 0x04034b50);
  le16(out_, 20); le16(out_, 0); le16(out_, e.method); le16(out_, kDosTime); le16(out_, kDosDate);
  le32(out_, e.crc); le32(out_, e.csize); le32(out_, e.usize);
  le16(out_, static_cast<std::uint32_t>(name.size())); le16(out_, 0);
  out_ += name;
  out_ += body;
  entries_.push_back(e);
}

std::string ZipWriter::finish() {
  std::uint32_t cd = static_cast<std::uint32_t>(out_.size());
  for (auto& e : entries_) {
    le32(out_, 0x02014b50);
    le16(out_, 20); le16(out_, 20); le16(out_, 0); le16(out_, e.method); le16(out_, kDosTime); le16(out_, kDosDate);
    le32(out_, e.crc); le32(out_, e.csize); le32(out_, e.usize);
    le16(out_, static_cast<std::uint32_t>(e.name.size())); le16(out_, 0); le16(out_, 0); le16(out_, 0); le16(out_, 0);
    le32(out_, 0); le32(out_, e.offset);
    out_ += e.name;
  }
  std::uint32_t cdsize = static_cast<std::uint32_t>(out_.size()) - cd;
  le32(out_, 0x06054b50);
  le16(out_, 0); le16(out_, 0);
  le16(out_, static_cast<std::uint32_t>(entries_.size())); le16(out_, static_cast<std::uint32_t>(entries_.size()));
  le32(out_, cdsize); le32(out_, cd); le16(out_, 0);
  return std::move(out_);
}

// ---- PPTX ----------------------------------------------------------------------------------
namespace {
const char* kXml = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
const char* kNsP = "xmlns:a=\"http://schemas.openxmlformats.org/drawingml/2006/main\" "
                   "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\" "
                   "xmlns:p=\"http://schemas.openxmlformats.org/presentationml/2006/main\"";
const char* kRelNs = "http://schemas.openxmlformats.org/officeDocument/2006/relationships";

std::string theme(const std::string& name) {
  return std::string(kXml) +
         "<a:theme xmlns:a=\"http://schemas.openxmlformats.org/drawingml/2006/main\" name=\"" + name + "\"><a:themeElements>"
         "<a:clrScheme name=\"Stratonix\"><a:dk1><a:srgbClr val=\"000000\"/></a:dk1><a:lt1><a:srgbClr val=\"FFFFFF\"/></a:lt1>"
         "<a:dk2><a:srgbClr val=\"1F3A5F\"/></a:dk2><a:lt2><a:srgbClr val=\"E8EEF4\"/></a:lt2>"
         "<a:accent1><a:srgbClr val=\"1F3A5F\"/></a:accent1><a:accent2><a:srgbClr val=\"C0504D\"/></a:accent2>"
         "<a:accent3><a:srgbClr val=\"9BBB59\"/></a:accent3><a:accent4><a:srgbClr val=\"8064A2\"/></a:accent4>"
         "<a:accent5><a:srgbClr val=\"4BACC6\"/></a:accent5><a:accent6><a:srgbClr val=\"F79646\"/></a:accent6>"
         "<a:hlink><a:srgbClr val=\"0000FF\"/></a:hlink><a:folHlink><a:srgbClr val=\"800080\"/></a:folHlink></a:clrScheme>"
         "<a:fontScheme name=\"Stratonix\"><a:majorFont><a:latin typeface=\"Arial\"/><a:ea typeface=\"\"/><a:cs typeface=\"\"/></a:majorFont>"
         "<a:minorFont><a:latin typeface=\"Arial\"/><a:ea typeface=\"\"/><a:cs typeface=\"\"/></a:minorFont></a:fontScheme>"
         "<a:fmtScheme name=\"Stratonix\"><a:fillStyleLst><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill>"
         "<a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill></a:fillStyleLst>"
         "<a:lnStyleLst><a:ln w=\"9525\"><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill></a:ln>"
         "<a:ln w=\"25400\"><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill></a:ln>"
         "<a:ln w=\"38100\"><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill></a:ln></a:lnStyleLst>"
         "<a:effectStyleLst><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle>"
         "<a:effectStyle><a:effectLst/></a:effectStyle></a:effectStyleLst>"
         "<a:bgFillStyleLst><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill>"
         "<a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill></a:bgFillStyleLst></a:fmtScheme></a:themeElements>"
         "<a:objectDefaults/><a:extraClrSchemeLst/></a:theme>";
}

std::string text_box(int id, const std::string& name, long x, long y, long cx, long cy, const std::string& text, int size_hpt,
                     bool bold) {
  std::string runs;
  for (auto& line : split(text, '\n'))
    runs += "<a:p><a:r><a:rPr lang=\"en-GB\" sz=\"" + std::to_string(size_hpt) + "\"" + (bold ? " b=\"1\"" : "") +
            " dirty=\"0\"><a:latin typeface=\"Arial\"/></a:rPr><a:t>" + xml_escape(line) + "</a:t></a:r></a:p>";
  return "<p:sp><p:nvSpPr><p:cNvPr id=\"" + std::to_string(id) + "\" name=\"" + name + "\"/><p:cNvSpPr txBox=\"1\"/><p:nvPr/></p:nvSpPr>"
         "<p:spPr><a:xfrm><a:off x=\"" + std::to_string(x) + "\" y=\"" + std::to_string(y) + "\"/><a:ext cx=\"" + std::to_string(cx) +
         "\" cy=\"" + std::to_string(cy) + "\"/></a:xfrm><a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom><a:noFill/></p:spPr>"
         "<p:txBody><a:bodyPr wrap=\"square\" rtlCol=\"0\"><a:normAutofit/></a:bodyPr><a:lstStyle/>" + runs + "</p:txBody></p:sp>";
}
const char* kEmptyTree =
    "<p:nvGrpSpPr><p:cNvPr id=\"1\" name=\"\"/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x=\"0\" y=\"0\"/>"
    "<a:ext cx=\"0\" cy=\"0\"/><a:chOff x=\"0\" y=\"0\"/><a:chExt cx=\"0\" cy=\"0\"/></a:xfrm></p:grpSpPr>";
}  // namespace

std::string make_pptx(const std::vector<Slide>& slides, const std::string& deck_title) {
  const long SW = 12192000, SH = 6858000;  // 16:9 in EMU
  ZipWriter z;
  std::string ct = std::string(kXml) +
                   "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
                   "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
                   "<Default Extension=\"xml\" ContentType=\"application/xml\"/><Default Extension=\"png\" ContentType=\"image/png\"/>"
                   "<Override PartName=\"/ppt/presentation.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml\"/>"
                   "<Override PartName=\"/ppt/slideMasters/slideMaster1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml\"/>"
                   "<Override PartName=\"/ppt/slideLayouts/slideLayout1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml\"/>"
                   "<Override PartName=\"/ppt/notesMasters/notesMaster1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.notesMaster+xml\"/>"
                   "<Override PartName=\"/ppt/theme/theme1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.theme+xml\"/>"
                   "<Override PartName=\"/ppt/theme/theme2.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.theme+xml\"/>"
                   "<Override PartName=\"/ppt/presProps.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.presProps+xml\"/>"
                   "<Override PartName=\"/ppt/viewProps.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.viewProps+xml\"/>"
                   "<Override PartName=\"/ppt/tableStyles.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.tableStyles+xml\"/>"
                   "<Override PartName=\"/docProps/core.xml\" ContentType=\"application/vnd.openxmlformats-package.core-properties+xml\"/>"
                   "<Override PartName=\"/docProps/app.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.extended-properties+xml\"/>";
  for (std::size_t i = 1; i <= slides.size(); ++i) {
    ct += "<Override PartName=\"/ppt/slides/slide" + std::to_string(i) + ".xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.slide+xml\"/>";
    ct += "<Override PartName=\"/ppt/notesSlides/notesSlide" + std::to_string(i) + ".xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.notesSlide+xml\"/>";
  }
  ct += "</Types>";
  z.add("[Content_Types].xml", ct);
  z.add("_rels/.rels", std::string(kXml) + "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                                           "<Relationship Id=\"rId1\" Type=\"" + kRelNs + "/officeDocument\" Target=\"ppt/presentation.xml\"/>"
                                           "<Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties\" Target=\"docProps/core.xml\"/>"
                                           "<Relationship Id=\"rId3\" Type=\"" + kRelNs + "/extended-properties\" Target=\"docProps/app.xml\"/></Relationships>");
  z.add("docProps/core.xml", std::string(kXml) +
                                 "<cp:coreProperties xmlns:cp=\"http://schemas.openxmlformats.org/package/2006/metadata/core-properties\" "
                                 "xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:dcterms=\"http://purl.org/dc/terms/\" "
                                 "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"><dc:title>" + xml_escape(deck_title) +
                                 "</dc:title><dc:creator>Stratonix</dc:creator></cp:coreProperties>");
  z.add("docProps/app.xml", std::string(kXml) +
                                "<Properties xmlns=\"http://schemas.openxmlformats.org/officeDocument/2006/extended-properties\">"
                                "<Application>Stratonix " + build_info().version + "</Application><Slides>" + std::to_string(slides.size()) +
                                "</Slides></Properties>");
  std::string pres = std::string(kXml) + "<p:presentation " + kNsP + " saveSubsetFonts=\"1\"><p:sldMasterIdLst><p:sldMasterId id=\"2147483648\" r:id=\"rId1\"/></p:sldMasterIdLst>"
                     "<p:notesMasterIdLst><p:notesMasterId r:id=\"rId2\"/></p:notesMasterIdLst><p:sldIdLst>";
  std::string prels = std::string(kXml) + "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                      "<Relationship Id=\"rId1\" Type=\"" + kRelNs + "/slideMaster\" Target=\"slideMasters/slideMaster1.xml\"/>"
                      "<Relationship Id=\"rId2\" Type=\"" + kRelNs + "/notesMaster\" Target=\"notesMasters/notesMaster1.xml\"/>"
                      "<Relationship Id=\"rId3\" Type=\"" + kRelNs + "/theme\" Target=\"theme/theme1.xml\"/>"
                      "<Relationship Id=\"rId4\" Type=\"" + kRelNs + "/presProps\" Target=\"presProps.xml\"/>"
                      "<Relationship Id=\"rId5\" Type=\"" + kRelNs + "/viewProps\" Target=\"viewProps.xml\"/>"
                      "<Relationship Id=\"rId6\" Type=\"" + kRelNs + "/tableStyles\" Target=\"tableStyles.xml\"/>";
  for (std::size_t i = 1; i <= slides.size(); ++i) {
    pres += "<p:sldId id=\"" + std::to_string(255 + i) + "\" r:id=\"rId" + std::to_string(10 + i) + "\"/>";
    prels += "<Relationship Id=\"rId" + std::to_string(10 + i) + "\" Type=\"" + kRelNs + "/slide\" Target=\"slides/slide" + std::to_string(i) + ".xml\"/>";
  }
  pres += "</p:sldIdLst><p:sldSz cx=\"" + std::to_string(SW) + "\" cy=\"" + std::to_string(SH) + "\"/><p:notesSz cx=\"6858000\" cy=\"9144000\"/>"
          "</p:presentation>";
  prels += "</Relationships>";
  z.add("ppt/presentation.xml", pres);
  z.add("ppt/_rels/presentation.xml.rels", prels);
  z.add("ppt/presProps.xml", std::string(kXml) + "<p:presentationPr " + kNsP + "/>");
  z.add("ppt/viewProps.xml", std::string(kXml) + "<p:viewPr " + kNsP + "><p:normalViewPr><p:restoredLeft sz=\"15620\"/><p:restoredTop sz=\"94660\"/></p:normalViewPr></p:viewPr>");
  z.add("ppt/tableStyles.xml", std::string(kXml) + "<a:tblStyleLst xmlns:a=\"http://schemas.openxmlformats.org/drawingml/2006/main\" def=\"{5C22544A-7EE6-4342-B048-85BDC9FD1C3A}\"/>");
  z.add("ppt/theme/theme1.xml", theme("Stratonix"));
  z.add("ppt/theme/theme2.xml", theme("Stratonix notes"));
  const std::string clrmap = "<p:clrMap bg1=\"lt1\" tx1=\"dk1\" bg2=\"lt2\" tx2=\"dk2\" accent1=\"accent1\" accent2=\"accent2\" accent3=\"accent3\" "
                             "accent4=\"accent4\" accent5=\"accent5\" accent6=\"accent6\" hlink=\"hlink\" folHlink=\"folHlink\"/>";
  z.add("ppt/slideMasters/slideMaster1.xml", std::string(kXml) + "<p:sldMaster " + kNsP + "><p:cSld><p:spTree>" + kEmptyTree +
                                                 "</p:spTree></p:cSld>" + clrmap +
                                                 "<p:sldLayoutIdLst><p:sldLayoutId id=\"2147483649\" r:id=\"rId1\"/></p:sldLayoutIdLst></p:sldMaster>");
  z.add("ppt/slideMasters/_rels/slideMaster1.xml.rels", std::string(kXml) +
                                                            "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                                                            "<Relationship Id=\"rId1\" Type=\"" + kRelNs + "/slideLayout\" Target=\"../slideLayouts/slideLayout1.xml\"/>"
                                                            "<Relationship Id=\"rId2\" Type=\"" + kRelNs + "/theme\" Target=\"../theme/theme1.xml\"/></Relationships>");
  z.add("ppt/slideLayouts/slideLayout1.xml", std::string(kXml) + "<p:sldLayout " + kNsP + " type=\"blank\" preserve=\"1\"><p:cSld name=\"Blank\"><p:spTree>" +
                                                 kEmptyTree + "</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sldLayout>");
  z.add("ppt/slideLayouts/_rels/slideLayout1.xml.rels", std::string(kXml) +
                                                            "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                                                            "<Relationship Id=\"rId1\" Type=\"" + kRelNs + "/slideMaster\" Target=\"../slideMasters/slideMaster1.xml\"/></Relationships>");
  z.add("ppt/notesMasters/notesMaster1.xml", std::string(kXml) + "<p:notesMaster " + kNsP + "><p:cSld><p:spTree>" + kEmptyTree +
                                                 "</p:spTree></p:cSld>" + clrmap + "</p:notesMaster>");
  z.add("ppt/notesMasters/_rels/notesMaster1.xml.rels", std::string(kXml) +
                                                            "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                                                            "<Relationship Id=\"rId1\" Type=\"" + kRelNs + "/theme\" Target=\"../theme/theme2.xml\"/></Relationships>");
  for (std::size_t i = 0; i < slides.size(); ++i) {
    const Slide& s = slides[i];
    std::string n = std::to_string(i + 1);
    // picture placement: below the title, above the caption, aspect kept
    long top = 640000, bottom_space = 700000, margin = 300000;
    long availw = SW - 2 * margin, availh = SH - top - bottom_space;
    double ar = s.px_h > 0 ? static_cast<double>(s.px_w) / s.px_h : 1.6;
    long pw = availw, ph = static_cast<long>(availw / ar);
    if (ph > availh) { ph = availh; pw = static_cast<long>(availh * ar); }
    long px = (SW - pw) / 2, py = top + (availh - ph) / 2;
    std::string slide = std::string(kXml) + "<p:sld " + kNsP + "><p:cSld><p:spTree>" + kEmptyTree +
                        text_box(2, "Title", margin, 120000, SW - 2 * margin, 480000, s.title, 2200, true) +
                        "<p:pic><p:nvPicPr><p:cNvPr id=\"3\" name=\"Figure\" descr=\"" + xml_escape(s.title) + "\"/><p:cNvPicPr><a:picLocks noChangeAspect=\"1\"/></p:cNvPicPr><p:nvPr/></p:nvPicPr>"
                        "<p:blipFill><a:blip r:embed=\"rId2\"/><a:stretch><a:fillRect/></a:stretch></p:blipFill>"
                        "<p:spPr><a:xfrm><a:off x=\"" + std::to_string(px) + "\" y=\"" + std::to_string(py) + "\"/><a:ext cx=\"" +
                        std::to_string(pw) + "\" cy=\"" + std::to_string(ph) + "\"/></a:xfrm><a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom></p:spPr></p:pic>" +
                        text_box(4, "Caption", margin, SH - bottom_space + 60000, SW - 2 * margin, bottom_space - 120000, s.caption, 1100, false) +
                        "</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>";
    z.add("ppt/slides/slide" + n + ".xml", slide);
    z.add("ppt/slides/_rels/slide" + n + ".xml.rels",
          std::string(kXml) + "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                              "<Relationship Id=\"rId1\" Type=\"" + kRelNs + "/slideLayout\" Target=\"../slideLayouts/slideLayout1.xml\"/>"
                              "<Relationship Id=\"rId2\" Type=\"" + kRelNs + "/image\" Target=\"../media/image" + n + ".png\"/>"
                              "<Relationship Id=\"rId3\" Type=\"" + kRelNs + "/notesSlide\" Target=\"../notesSlides/notesSlide" + n + ".xml\"/></Relationships>");
    z.add("ppt/media/image" + n + ".png", s.png, false);
    std::string notes = std::string(kXml) + "<p:notes " + kNsP + "><p:cSld><p:spTree>" + kEmptyTree +
                        "<p:sp><p:nvSpPr><p:cNvPr id=\"2\" name=\"Notes\"/><p:cNvSpPr><a:spLocks noGrp=\"1\"/></p:cNvSpPr><p:nvPr><p:ph type=\"body\" idx=\"1\"/></p:nvPr></p:nvSpPr>"
                        "<p:spPr><a:xfrm><a:off x=\"685800\" y=\"4400550\"/><a:ext cx=\"5486400\" cy=\"4000500\"/></a:xfrm></p:spPr><p:txBody><a:bodyPr/><a:lstStyle/>";
    for (auto& line : split(s.notes, '\n')) notes += "<a:p><a:r><a:rPr lang=\"en-GB\" dirty=\"0\"/><a:t>" + xml_escape(line) + "</a:t></a:r></a:p>";
    notes += "</p:txBody></p:sp></p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:notes>";
    z.add("ppt/notesSlides/notesSlide" + n + ".xml", notes);
    z.add("ppt/notesSlides/_rels/notesSlide" + n + ".xml.rels",
          std::string(kXml) + "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                              "<Relationship Id=\"rId1\" Type=\"" + kRelNs + "/notesMaster\" Target=\"../notesMasters/notesMaster1.xml\"/>"
                              "<Relationship Id=\"rId2\" Type=\"" + kRelNs + "/slide\" Target=\"../slides/slide" + n + ".xml\"/></Relationships>");
  }
  return z.finish();
}

}  // namespace stratonix::fig
