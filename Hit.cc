#include "Hit.h"

void MCHitInfo::reset()
{
}

void print(std::string label, const MeasurementState& s)
{
  std::cout << label << std::endl;
  std::cout << "x: "  << s.parameters()[0]
            << " y: " << s.parameters()[1]
            << " z: " << s.parameters()[2] << std::endl
            << "errors: " << std::endl;
  dumpMatrix(s.errors());
  std::cout << std::endl;
}
