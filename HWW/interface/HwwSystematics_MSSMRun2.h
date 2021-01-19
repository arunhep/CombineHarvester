#ifndef HWW_HwwSystematics_MSSMRun2_h
#define HWW_HwwSystematics_MSSMRun2_h
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {
// Run2 MSSM analysis systematics
//void AddMSSMFullRun2Systematics(CombineHarvester& cb, int highmass, std::string mssmorindep, std::string model);
void AddMSSMFullRun2Systematics(CombineHarvester& cb, int highmass, bool DNNdiscrim, std::string model);
}
#endif
