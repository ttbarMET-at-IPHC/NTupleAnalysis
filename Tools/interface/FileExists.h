#ifndef FileExists_h
#define FileExists_h
#include <iostream>
#include <fstream>
#include <assert.h>

  /**
   * Will test whether a file exists
   * if fail is set to true, executable will exit
   */

bool fexists(const std::string filename, bool fail);

#endif
