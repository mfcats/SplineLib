#ifndef SRC_IO_IO_CONVERTER_H_
#define SRC_IO_IO_CONVERTER_H_

#include "iges_reader.h"
#include "xml_writer.h"

namespace io {
template<int DIM>
class IOConverter {
 public:
  IOConverter() = default;

  void ConvertIGESFileToXMLFile(const char *input_filename, const char *output_filename) {
    io::IGESReader iges_reader;
    std::vector<std::any> input_splines = iges_reader.ReadIGESFile(input_filename);
    std::vector<std::any> output_splines;
    for (const auto &spline : input_splines) {
      if (GetSplineType(spline) == DIM) {
        output_splines.emplace_back(spline);
      }
    }
    io::XMLWriter<DIM> xml_writer;
    xml_writer.WriteXMLFile(output_splines, output_filename);
  }

 private:
  int GetSplineType(std::any spline) {
    try {
      std::any_cast<std::shared_ptr<spl::BSpline<1>>>(spline);
      return 1;
    } catch (std::bad_any_cast &msg) {
      try {
        std::any_cast<std::shared_ptr<spl::NURBS<1>>>(spline);
        return 1;
      } catch (std::bad_any_cast &msg) {
        try {
          std::any_cast<std::shared_ptr<spl::BSpline<2>>>(spline);
          return 2;
        } catch (std::bad_any_cast &msg) {
          try {
            std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline);
            return 2;
          } catch (std::bad_any_cast &msg) {
            return 0;
          }
        }
      }
    }
  }
};
}  // namespace io

#endif  // SRC_IO_IO_CONVERTER_H_
