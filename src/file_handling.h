#ifndef FILEHANDLING_H
#define FILEHANDLING_H

#include <fstream>
#include <string>
#include <TString.h>
#include <array>

std::array<float, 2> getELifetime(const std::string& cintyp);
TString getOutputNameChannelPerTheta(const std::string& cintyp);
TString getOutputNameChannelByAngle(const std::string& cintyp);
TString getOutputNameChannel(const std::string& cintyp);
TString getOutputNameAng(const std::string& cintyp);
TString getOutputNameYZ(const std::string& cintyp);
TString getOutputNameX(const std::string& cintyp);
TString getOutputNamedQdx(const std::string& cintyp);
TString getOutputNameOfflinePitch(const std::string& cintyp);
std::ifstream openInputFile(const std::string& cintyp);

#endif
