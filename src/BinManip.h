// created 09/14/2017 by Brennan Young
// last modified 09/14/2017 by Brennan Young

// General-use structures and functions for binary file manipulation.

#ifndef YOUNG_STDLIB_BINMANIP_20170914
#define YOUNG_STDLIB_BINMANIP_20170914

#include <fstream>  // std::ifstream, std::ofstream
#include <string>   // std::string

namespace bystd { // Brennan Young standard namespace

// created 07/10/2017 by Brennan Young
// modified 09/14/2017 by Brennan Young
//   * moved to BinManip library, under bystd namespace
std::string readStrBinary(std::ifstream & file){
    size_t len;                                 // read length of string
    file.read((char *) &len, sizeof(size_t));
    char * temp = new char[len+1];              // read string into char array
    file.read(temp, len);
    temp[len] = '\0';                           // append null character
    std::string s = temp;
    delete [] temp;                             // clean up
    return s;
}
void writeStrBinary(std::ofstream & file, const std::string & s){
    size_t len = s.length();                    // length of string
    file.write((char*) &len, sizeof(size_t));   // write length of string
    file.write(s.c_str(), len);                 // write string
}

} // namespace bystd

#endif // YOUNG_STDLIB_BINMANIP_20170914