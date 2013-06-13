#include "Tools/interface/FileExists.h"

bool fexists(const std::string filename, bool fail)
{
  std::ifstream ifile(filename.c_str());
  if (fail && !ifile) {
    std::cout << "File does not exist: " << filename << std::endl;
    assert(ifile);
  }
  return ifile;
}

