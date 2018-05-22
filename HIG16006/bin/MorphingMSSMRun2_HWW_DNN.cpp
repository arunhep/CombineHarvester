#include <string>
#include <map>
#include <set>
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib>
#include "boost/algorithm/string/predicate.hpp"
#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/CardWriter.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"
#include "CombineHarvester/CombineTools/interface/AutoRebin.h"
#include "CombineHarvester/CombinePdfs/interface/MorphFunctions.h"
#include "CombineHarvester/HIG16006/interface/HwwSystematics_MSSMRun2.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TH2.h"
#include "TF1.h"

using namespace std;
using boost::starts_with;
namespace po = boost::program_options;

template <typename T>
void To1Bin(T* proc)
{
    std::unique_ptr<TH1> originalHist = proc->ClonedScaledShape();
    TH1F *hist = new TH1F("hist","hist",1,0,1);
    double err = 0;
    double rate =
        originalHist->IntegralAndError(0, originalHist->GetNbinsX() + 1, err);
    hist->SetDirectory(0);
    hist->SetBinContent(1, rate);
    hist->SetBinError(1, err);
    proc->set_shape(*hist, true);  // True means adjust the process rate to the
                                   // integral of the hist
}

bool BinIsControlRegion(ch::Object const* obj)
{
    return ((obj->channel() == std::string("dytt")) || (obj->channel() == std::string("top")) || (obj->channel() == std::string("eedytt")) || (obj->channel() == std::string("eetop")) || (obj->channel() == std::string("mmdytt")) || (obj->channel() == std::string("mmtop"))    );
}

// Useful to have the inverse sometimes too
bool BinIsNotControlRegion(ch::Object const* obj)
{
    return !BinIsControlRegion(obj);
}



int main(int argc, char** argv) {
  // First define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  string SM125= "";
  string mass = "";
  string output_folder = "mssm_run2";
  // TODO: option to pick up cards from different dirs depending on channel?
  // ^ Something like line 90?
  string input_folder_em="RWTH/";
  string postfix="-mti";
  bool auto_rebin = false;
  bool manual_rebin = false;
  bool real_data = true;
  int control_region = 0;
  bool check_neg_bins = false;
  bool poisson_bbb = false;
  bool do_w_weighting = false;
  bool include_VBF = true;
  bool include_SF = false;
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
    ("mass,m", po::value<string>(&mass)->default_value(mass))
    ("input_folder_em", po::value<string>(&input_folder_em)->default_value("RWTH"))
    ("postfix", po::value<string>(&postfix)->default_value(""))
    ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(false))
    ("real_data", po::value<bool>(&real_data)->default_value(true))
    ("manual_rebin", po::value<bool>(&manual_rebin)->default_value(false))
    ("output_folder", po::value<string>(&output_folder)->default_value("mssm_run2"))
    ("SM125,h", po::value<string>(&SM125)->default_value(SM125))
    ("control_region", po::value<int>(&control_region)->default_value(0))
    ("check_neg_bins", po::value<bool>(&check_neg_bins)->default_value(false))
    ("poisson_bbb", po::value<bool>(&poisson_bbb)->default_value(false))
    ("w_weighting", po::value<bool>(&do_w_weighting)->default_value(false));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  typedef vector<string> VString;
  typedef vector<pair<int, string>> Categories;
  std::map<string, string> input_dir;
  input_dir["em"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HIG16006/plots_HWWhighMass2016_OF_peter.root";
  input_dir["dytt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HIG16006/plots_HWWhighMass2016_OF_peter.root";
  input_dir["top"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HIG16006/plots_HWWhighMass2016_OF_peter.root";
  if (include_SF){
  input_dir["ee"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HIG16006/plots_HWWhighMass2016_SF_embeddedWeights_ee_VBF_V2.root";
  input_dir["eedytt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HIG16006/plots_HWWhighMass2016_SF_embeddedWeights_ee_VBF_V2.root";
  input_dir["eetop"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HIG16006/plots_HWWhighMass2016_SF_embeddedWeights_ee_VBF_V2.root";
  input_dir["mm"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HIG16006/plots_HWWhighMass2016_SF_embeddedWeights_mm_VBF_V2.root";
  input_dir["mmdytt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HIG16006/plots_HWWhighMass2016_SF_embeddedWeights_mm_VBF_V2.root";
  input_dir["mmtop"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HIG16006/plots_HWWhighMass2016_SF_embeddedWeights_mm_VBF_V2.root";
  }

  VString chns =
      {"em","dytt","top"};
  if (include_SF){
  chns = {"em","dytt","top", "ee","eedytt","eetop", "mm","mmdytt","mmtop"};
  }
  /*if (control_region==1){
  chns =
      {"dytt"};
  }
  if (control_region==2){
  chns =
      {"top"};
  }*/

  //RooRealVar mA(mass.c_str(), mass.c_str(), 90., 3200.);
  //mA.setConstant(true);
  string massvar="";
  if (mass=="indep"){
    massvar="MH";
  }else{
    massvar="mH";
  }
  RooRealVar mH(massvar.c_str(), massvar.c_str(), 90., 3200.);
  if (mass=="indep"){
    mH.setConstant(true);
  }
  //RooRealVar mh("mh", "mh", 90., 3200.);

  map<string, VString> bkg_procs;
  bkg_procs["em"] = {"WW", "DY", "VVV", "top", "Fake", "ggWW", "VZ", "Vg", "VgS", "WW2J", "qqWWqq"};//{"ggWW", "Vg", "VVV", "top", "Fake", "WW", "DY", "VZ", "VgS"};
  bkg_procs["dytt"] = {"WW", "DY", "VVV", "top", "Fake", "ggWW", "VZ", "Vg", "VgS", "WW2J", "qqWWqq"};
  bkg_procs["top"] = {"WW", "DY", "VVV", "top", "Fake", "ggWW", "VZ", "Vg", "VgS", "WW2J", "qqWWqq"};
  if (include_SF){
  bkg_procs["ee"] = {"WW", "DY", "VVV", "top", "Fake", "ggWW", "VZ", "Vg", "VgS", "WW2J", "qqWWqq"};
  bkg_procs["eedytt"] = {"WW", "DY", "VVV", "top", "Fake", "ggWW", "VZ", "Vg", "VgS", "WW2J", "qqWWqq"};
  bkg_procs["eetop"] = {"WW", "DY", "VVV", "top", "Fake", "ggWW", "VZ", "Vg", "VgS", "WW2J", "qqWWqq"};
  bkg_procs["mm"] = {"WW", "DY", "VVV", "top", "Fake", "ggWW", "VZ", "Vg", "VgS", "WW2J", "qqWWqq"};
  bkg_procs["mmdytt"] = {"WW", "DY", "VVV", "top", "Fake", "ggWW", "VZ", "Vg", "VgS", "WW2J", "qqWWqq"};
  bkg_procs["mmtop"] = {"WW", "DY", "VVV", "top", "Fake", "ggWW", "VZ", "Vg", "VgS", "WW2J", "qqWWqq"};
  }

  VString SM_procs = {"WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt"};//{"ggH_hww_SM125", "qqH_hww_SM125", "ZH_hww_SM125", "ggZH_hww_SM125","WH_hww_SM125","bbH_hww_SM125"};
  //VString SM_tauprocs = {"H_htt"};//{"ggH_htt_SM125", "qqH_htt_SM125", "ZH_htt_SM125", "WH_htt_SM125"};

  //Example - could fill this map with hardcoded binning for different
  //categories if manual_rebin is turned on
  map<string, vector<double> > binning;
  binning["em_of0j"] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900,4000};
  binning["em_of1j"] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900,4000};
  binning["em_of2j"] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900,4000};
  if (include_VBF){
  binning["em_of2j_vbf"] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900,4000};
  }

  // Create an empty CombineHarvester instance that will hold all of the
  // datacard configuration and histograms etc.
  ch::CombineHarvester cb;
  // Uncomment this next line to see a *lot* of debug information
  // cb.SetVerbosity(3);

  // Here we will just define two categories for an 8TeV analysis. Each entry in
  // the vector below specifies a bin name and corresponding bin_id.
  //
  map<string,Categories> cats;
  if (include_VBF){
  cats["em_13TeV"] = {
    {8, "of0j"},
    {9, "of1j"},
    {10, "of2j"},
    {11, "ofVBF"}//"of2j_vbf"
    };
  cats["dytt_13TeV"] = {
    {8, "dy_of0j"},
    {9, "dy_of1j"},
    {10, "dy_of2j"},
    {11, "dy_ofVBF"}//"of2j_vbf"
    };
  cats["top_13TeV"] = {
    {8, "top_of0j"},
    {9, "top_of1j"},
    {10, "top_of2j"},
    {11, "top_ofVBF"}//"of2j_vbf"
    };
  }else{
  cats["em_13TeV"] = {
    {8, "of_0j"},
    {9, "of_1j"},
    {10, "of2j"}
    };
  cats["dytt_13TeV"] = {
    {8, "dytt_of0j"},
    {9, "dytt_of1j"},
    {10, "dytt_of2j"}
    };
  cats["top_13TeV"] = {
    {8, "top_of0j"},
    {9, "top_of1j"},
    {10, "top_of2j"}
    };
  }
  if (include_SF){
  cats["ee_13TeV"] = {
    {11, "e_e_2j_VBF"}
    };
  cats["eedytt_13TeV"] = {
    {11, "dy_e_e_2j_VBF"}
    };
  cats["eetop_13TeV"] = {
    {11, "top_e_e_2j_VBF"}
    };
  cats["mm_13TeV"] = {
    {11, "mu_mu_2j_VBF"}
    };
  cats["mmdytt_13TeV"] = {
    {11, "dy_mu_mu_2j_VBF"}
    };
  cats["mmtop_13TeV"] = {
    {11, "top_mu_mu_2j_VBF"}
    };
  }

  vector<string> masses = {"300", "350", "400", "450", "500", "550", "650", "700", "750", "800", "900", "1000", "1500", "2000", "2500", "3000"};//{"120","124","126","130","135","140","145","150","155","160","165","170","175","180","190","200","210","230", "250","270","300","350","400","450","500","550","650","700","750","800","900","1000","1500","2000","2500","3000"};//{"90","100","110","120","130","140","160","180", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900","1000","1200","1400","1500","1600","1800","2000","2300","2600","2900","3200"};
  vector<string> massesadd = {"200", "210", "230", "250", "270"};

  map<string, VString> signal_types = {
    {"ggH", {"ggH_HWW"}},//,
    {"qqH", {"qqH_HWW"}},
    {"ggH_SBI", {"ggH_SBIHWW"}},
    {"qqH_SBI", {"qqH_SBIHWW"}}
    //{"bbH", {"bbh_htautau", "bbH_Htautau", "bbA_Atautau"}}
  };
  if (!include_VBF){
signal_types = {
    {"ggH", {"ggH_HWW"}},//,
    {"ggH_SBI", {"ggH_SBIHWW"}}
    //{"bbH", {"bbh_htautau", "bbH_Htautau", "bbA_Atautau"}}
  };
  }
/*  if(mass=="MH"){
    signal_types = {
      {"ggH", {"ggH"}}//,
      //{"bbH", {"bbH"}}
    };
  }*/
    vector<string> sig_procs = {"ggH","qqH","ggH_SBI","qqH_SBI"};//,"bbH"
  if (!include_VBF){
    sig_procs = {"ggH","ggH_SBI"};
  }
  for(auto chn : chns){
    cb.AddObservations({"*"}, {"hww"}, {"13TeV"}, {chn}, cats[chn+"_13TeV"]);

    cb.AddProcesses({"*"}, {"hww"}, {"13TeV"}, {chn}, bkg_procs[chn], cats[chn+"_13TeV"], false);

    cb.AddProcesses(masses, {"hww"}, {"13TeV"}, {chn}, signal_types["ggH"], cats[chn+"_13TeV"], true);
    if (include_VBF){
    cb.AddProcesses(masses, {"hww"}, {"13TeV"}, {chn}, signal_types["qqH"], cats[chn+"_13TeV"], true);
    }
    cb.AddProcesses(masses, {"hww"}, {"13TeV"}, {chn}, signal_types["ggH_SBI"], cats[chn+"_13TeV"], true);
    if (include_VBF){
    cb.AddProcesses(masses, {"hww"}, {"13TeV"}, {chn}, signal_types["qqH_SBI"], cats[chn+"_13TeV"], true);
    }
    //cb.AddProcesses(masses, {"hww"}, {"13TeV"}, {chn}, signal_types["bbH"], cats[chn+"_13TeV"], true);
    //cb.AddProcesses({"*"}, {"hww"}, {"13TeV"}, {chn}, SM_tauprocs, cats[chn+"_13TeV"], false);
    if (chn == "em" or chn == "dytt" or chn == "top"){
    cb.AddProcesses({"*"}, {"hww"}, {"13TeV"}, {chn}, SM_procs, cats[chn+"_13TeV"], false);
    }
    }


  //Remove signal from ControlRegion
  cb.FilterAll([](ch::Object const* obj) {
              return (BinIsControlRegion(obj) && obj->signal());
              });

  ch::AddMSSMRun2Systematics_HWW_lat(cb,control_region);
  //! [part7]
  for (string chn:chns){
  if (chn == "eedytt" or chn == "eetop" or chn == "mmdytt" or chn == "mmtop"){

    std::cout << "Test1 \n ";
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn],
        "hww2l2v_13TeV_$BIN/events/histo_$PROCESS",
        "hww2l2v_13TeV_$BIN/events/histo_$PROCESS_$SYSTEMATIC");
    std::cout << "Test2 \n ";
/*    cb.cp().channel({chn}).process(SM_procs).ExtractShapes(
         input_dir[chn],
         "hwwhm_13TeV_$BIN/mTi/histo_$PROCESS",
         "hwwhm_13TeV_$BIN/mTi/histo_$PROCESS_$SYSTEMATIC");*/
    std::cout << "Test3 \n ";
    cb.cp().channel({chn}).process(signal_types["ggH"]).ExtractShapes(
        input_dir[chn],
        "hww2l2v_13TeV_$BIN/events/histo_ggH_hww_$MASS_c10brn00",
        "hww2l2v_13TeV_$BIN/events/histo_ggH_hww_$MASS_c10brn00_$SYSTEMATIC");
    if (include_VBF){
    cb.cp().channel({chn}).process(signal_types["qqH"]).ExtractShapes(
        input_dir[chn],
        "hww2l2v_13TeV_$BIN/events/histo_qqH_hww_$MASS_c10brn00",
        "hww2l2v_13TeV_$BIN/events/histo_qqH_hww_$MASS_c10brn00_$SYSTEMATIC");
    }
    cb.cp().channel({chn}).process(signal_types["ggH_SBI"]).ExtractShapes(
        input_dir[chn],
        "hww2l2v_13TeV_$BIN/events/histo_ggH_hww_SBI$MASS_c10brn00",
        "hww2l2v_13TeV_$BIN/events/histo_ggH_hww_SBI$MASS_c10brn00_$SYSTEMATIC");
    if (include_VBF){
    cb.cp().channel({chn}).process(signal_types["qqH_SBI"]).ExtractShapes(
        input_dir[chn],
        "hww2l2v_13TeV_$BIN/events/histo_qqH_hww_SBI$MASS_c10brn00",
        "hww2l2v_13TeV_$BIN/events/histo_qqH_hww_SBI$MASS_c10brn00_$SYSTEMATIC");
    }
    std::cout << "Test4 \n ";
  }else if (chn == "ee" or chn == "mm"){
    std::cout << "Test5 \n ";
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn],
        "hwwhm_13TeV_$BIN/ml_max_score/histo_$PROCESS",
        "hwwhm_13TeV_$BIN/ml_max_score/histo_$PROCESS_$SYSTEMATIC");
    std::cout << "Test6 \n ";
/*    cb.cp().channel({chn}).process(SM_procs).ExtractShapes(
         input_dir[chn],
         "hwwhm_13TeV_$BIN/mTi/histo_$PROCESS",
         "hwwhm_13TeV_$BIN/mTi/histo_$PROCESS_$SYSTEMATIC");*/
    std::cout << "Test7 \n ";
    cb.cp().channel({chn}).process(signal_types["ggH"]).ExtractShapes(
        input_dir[chn],
        "hwwhm_13TeV_$BIN/ml_max_score/histo_ggH_hww_$MASS_c10brn00",
        "hwwhm_13TeV_$BIN/ml_max_score/histo_ggH_hww_$MASS_c10brn00_$SYSTEMATIC");
    if (include_VBF){
    cb.cp().channel({chn}).process(signal_types["qqH"]).ExtractShapes(
        input_dir[chn],
        "hwwhm_13TeV_$BIN/ml_max_score/histo_qqH_hww_$MASS_c10brn00",
        "hwwhm_13TeV_$BIN/ml_max_score/histo_qqH_hww_$MASS_c10brn00_$SYSTEMATIC");
    }
    cb.cp().channel({chn}).process(signal_types["ggH_SBI"]).ExtractShapes(
        input_dir[chn],
        "hwwhm_13TeV_$BIN/ml_max_score/histo_ggH_hww_SBI$MASS_c10brn00",
        "hwwhm_13TeV_$BIN/ml_max_score/histo_ggH_hww_SBI$MASS_c10brn00_$SYSTEMATIC");
    if (include_VBF){
    cb.cp().channel({chn}).process(signal_types["qqH_SBI"]).ExtractShapes(
        input_dir[chn],
        "hwwhm_13TeV_$BIN/ml_max_score/histo_qqH_hww_SBI$MASS_c10brn00",
        "hwwhm_13TeV_$BIN/ml_max_score/histo_qqH_hww_SBI$MASS_c10brn00_$SYSTEMATIC");
    }
    std::cout << "Test8 \n ";
   // cb.cp().channel({chn}).process(signal_types["bbH"]).ExtractShapes(
   //     input_dir[chn] + "hww_"+chn+".inputs-mssm-13TeV"+postfix+".root",
   //     "$BIN/bbH$MASS",
   //     "$BIN/bbH$MASS_$SYSTEMATIC");
  }
  }










  for(auto chn : chns){
    cb.AddProcesses(massesadd, {"hww"}, {"13TeV"}, {chn}, signal_types["ggH"], cats[chn+"_13TeV"], true);
    if (include_VBF){
    cb.AddProcesses(massesadd, {"hww"}, {"13TeV"}, {chn}, signal_types["qqH"], cats[chn+"_13TeV"], true);
    }
    cb.AddProcesses(massesadd, {"hww"}, {"13TeV"}, {chn}, signal_types["ggH_SBI"], cats[chn+"_13TeV"], true);
    if (include_VBF){
    cb.AddProcesses(massesadd, {"hww"}, {"13TeV"}, {chn}, signal_types["qqH_SBI"], cats[chn+"_13TeV"], true);
    }
    }

  //Remove signal from ControlRegion
  cb.FilterAll([](ch::Object const* obj) {
              return (BinIsControlRegion(obj) && obj->signal());
              });

  //! [part7]
  for (string chn:chns){
  if (chn == "dytt" or chn == "top"){

    std::cout << "Test1 \n ";
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn],
        "hww2l2v_13TeV_$BIN/events/histo_$PROCESS",
        "hww2l2v_13TeV_$BIN/events/histo_$PROCESS_$SYSTEMATIC");
    std::cout << "Test2 \n ";
/*    cb.cp().channel({chn}).process(SM_procs).ExtractShapes(
         input_dir[chn],
         "hwwhm_13TeV_$BIN/mTi/histo_$PROCESS",
         "hwwhm_13TeV_$BIN/mTi/histo_$PROCESS_$SYSTEMATIC");*/
    std::cout << "Test3 \n ";
    cb.cp().channel({chn}).process(signal_types["ggH"]).ExtractShapes(
        input_dir[chn],
        "hww2l2v_13TeV_$BIN/events/histo_ggH_hww_$MASS_c10brn00",
        "hww2l2v_13TeV_$BIN/events/histo_ggH_hww_$MASS_c10brn00_$SYSTEMATIC");
    if (include_VBF){
    cb.cp().channel({chn}).process(signal_types["qqH"]).ExtractShapes(
        input_dir[chn],
        "hww2l2v_13TeV_$BIN/events/histo_qqH_hww_$MASS_c10brn00",
        "hww2l2v_13TeV_$BIN/events/histo_qqH_hww_$MASS_c10brn00_$SYSTEMATIC");
    }
    cb.cp().channel({chn}).process(signal_types["ggH_SBI"]).ExtractShapes(
        input_dir[chn],
        "hww2l2v_13TeV_$BIN/events/histo_ggH_hww_SBI$MASS_c10brn00",
        "hww2l2v_13TeV_$BIN/events/histo_ggH_hww_SBI$MASS_c10brn00_$SYSTEMATIC");
    if (include_VBF){
    cb.cp().channel({chn}).process(signal_types["qqH_SBI"]).ExtractShapes(
        input_dir[chn],
        "hww2l2v_13TeV_$BIN/events/histo_qqH_hww_SBI$MASS_c10brn00",
        "hww2l2v_13TeV_$BIN/events/histo_qqH_hww_SBI$MASS_c10brn00_$SYSTEMATIC");
    }
    std::cout << "Test4 \n ";
  }else if (chn == "em"){
    std::cout << "Test5 \n ";
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn],
        "hwwhm_13TeV_$BIN/ml_max_score/histo_$PROCESS",
        "hwwhm_13TeV_$BIN/ml_max_score/histo_$PROCESS_$SYSTEMATIC");
    std::cout << "Test6 \n ";
/*    cb.cp().channel({chn}).process(SM_procs).ExtractShapes(
         input_dir[chn],
         "hwwhm_13TeV_$BIN/mTi/histo_$PROCESS",
         "hwwhm_13TeV_$BIN/mTi/histo_$PROCESS_$SYSTEMATIC");*/
    std::cout << "Test7 \n ";
    cb.cp().channel({chn}).process(signal_types["ggH"]).ExtractShapes(
        input_dir[chn],
        "hwwhm_13TeV_$BIN/ml_max_score/histo_ggH_hww_$MASS_c10brn00",
        "hwwhm_13TeV_$BIN/ml_max_score/histo_ggH_hww_$MASS_c10brn00_$SYSTEMATIC");
    if (include_VBF){
    cb.cp().channel({chn}).process(signal_types["qqH"]).ExtractShapes(
        input_dir[chn],
        "hwwhm_13TeV_$BIN/ml_max_score/histo_qqH_hww_$MASS_c10brn00",
        "hwwhm_13TeV_$BIN/ml_max_score/histo_qqH_hww_$MASS_c10brn00_$SYSTEMATIC");
    }
    cb.cp().channel({chn}).process(signal_types["ggH_SBI"]).ExtractShapes(
        input_dir[chn],
        "hwwhm_13TeV_$BIN/ml_max_score/histo_ggH_hww_SBI$MASS_c10brn00",
        "hwwhm_13TeV_$BIN/ml_max_score/histo_ggH_hww_SBI$MASS_c10brn00_$SYSTEMATIC");
    if (include_VBF){
    cb.cp().channel({chn}).process(signal_types["qqH_SBI"]).ExtractShapes(
        input_dir[chn],
        "hwwhm_13TeV_$BIN/ml_max_score/histo_qqH_hww_SBI$MASS_c10brn00",
        "hwwhm_13TeV_$BIN/ml_max_score/histo_qqH_hww_SBI$MASS_c10brn00_$SYSTEMATIC");
     }
    std::cout << "Test8 \n ";
   // cb.cp().channel({chn}).process(signal_types["bbH"]).ExtractShapes(
   //     input_dir[chn] + "hww_"+chn+".inputs-mssm-13TeV"+postfix+".root",
   //     "$BIN/bbH$MASS",
   //     "$BIN/bbH$MASS_$SYSTEMATIC");
  }
  }










 //Now delete processes with 0 yield
 cb.FilterProcs([&](ch::Process *p) {
  bool null_yield = !(p->rate() > 0.);// || BinIsControlRegion(p)
  if (null_yield){
     std::cout << "[Null yield] Removing process with null yield: \n ";
     std::cout << ch::Process::PrintHeader << *p << "\n"; 
     cb.FilterSysts([&](ch::Systematic *s){
       bool remove_syst = (MatchingProcess(*p,*s));
       return remove_syst;
    });
  }
  return null_yield;
 });



    /* cb.cp().process(SM_tauprocs).ForEachProc([&](ch::Process * proc) {
       proc->set_rate(proc->rate()*0.0627);
      });
     cb.cp().process({"ggH_htt_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*44.14);
     });
     cb.cp().process({"qqH_htt_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*3.782);
     });
     cb.cp().process({"WH_htt_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*0.6863);
     });
     cb.cp().process({"ZH_htt_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*0.8839);
     });*/
/*     cb.cp().process({"TTH_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*0.5071);
     });*/


//  // And convert any shapes in the CRs to lnN: ??????????
//  // Convert all shapes to lnN at this stage
//  cb.cp().FilterSysts(BinIsNotControlRegion).syst_type({"shape"}).ForEachSyst([](ch::Systematic *sys) {
//    sys->set_type("lnN");
//  });

   //Replacing observation with the sum of the backgrounds (asimov) - nice to ensure blinding
    auto bins = cb.cp().bin_set();
    // For control region bins data (should) = sum of bkgs already
    // useful to be able to check this, so don't do the replacement
    // for these
  if(!real_data){
      for (auto b : cb.cp().FilterAll(BinIsControlRegion).bin_set()) {
          std::cout << " - Replacing data with asimov in bin " << b << "\n";
          cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
            obs->set_shape(cb.cp().bin({b}).backgrounds().GetShape(), true);
          });
        }
  }


  // Merge to one bin for control region bins
  //cb.cp().FilterAll(BinIsNotControlRegion).ForEachProc(To1Bin<ch::Process>);
  //cb.cp().FilterAll(BinIsNotControlRegion).ForEachObs(To1Bin<ch::Observation>);

  // Rebinning
  // --------------------
  // Keep track of shapes before and after rebinning for comparison
  // and for checking the effect of negative bin contents
  std::map<std::string, TH1F> before_rebin;
  std::map<std::string, TH1F> after_rebin;
  std::map<std::string, TH1F> after_rebin_neg;
  if (check_neg_bins) {
    for (auto b : bins) {
      before_rebin[b] = cb.cp().bin({b}).backgrounds().GetShape();
    }
  }


  auto rebin = ch::AutoRebin()
    .SetBinThreshold(0.)
    // .SetBinUncertFraction(0.5)
    .SetRebinMode(1)
    .SetPerformRebin(true)
    .SetVerbosity(1);
  if(auto_rebin) rebin.Rebin(cb, cb);

  if(manual_rebin) {
    for(auto b : bins) {
      std::cout << "Rebinning by hand for bin: " << b <<  std::endl;
      cb.cp().bin({b}).VariableRebin(binning[b]);
    }
  }

  if (check_neg_bins) {
    for (auto b : bins) {
      after_rebin[b] = cb.cp().bin({b}).backgrounds().GetShape();
      // std::cout << "Bin: " << b << " (before)\n";
      // before_rebin[b].Print("range");
      // std::cout << "Bin: " << b << " (after)\n";
      // after_rebin[b].Print("range");
      // Build a sum-of-bkgs TH1 that doesn't truncate the negative yields
      // like the CH GetShape does
      for (auto p : cb.cp().bin({b}).backgrounds().process_set()) {
        TH1F proc_hist;
        cb.cp().bin({b}).process({p}).ForEachProc([&](ch::Process *proc) {
          proc_hist = proc->ShapeAsTH1F();
          proc_hist.Scale(proc->no_norm_rate());
          for (int i = 1; i <= proc_hist.GetNbinsX(); ++i) {
            if (proc_hist.GetBinContent(i) < 0.) {
              std::cout << p << " bin " << i << ": " << proc_hist.GetBinContent(i) << "\n";
            }
          }
        });
        if (after_rebin_neg.count(b)) {
          after_rebin_neg[b].Add(&proc_hist);
        } else {
          after_rebin_neg[b] = proc_hist;
        }
      }
      std::cout << "Bin: " << b << "\n";
      for (int i = 1; i <= after_rebin[b].GetNbinsX(); ++i) {
        double offset = after_rebin[b].GetBinContent(i) - after_rebin_neg[b].GetBinContent(i);
        double offset_by_yield = offset / after_rebin[b].GetBinContent(i);
        double offset_by_err = offset / after_rebin[b].GetBinError(i);
        printf("%-2i offset %-10.4f tot %-10.4f err %-10.4f off/tot %-10.4f off/err %-10.4f\n", i , offset, after_rebin[b].GetBinContent(i), after_rebin[b].GetBinError(i), offset_by_yield, offset_by_err);
      }
    }
  }

  // Uncomment this to inject 1 obs event in the last bin of every signal-region
  // category
  // if(!real_data){
  //     for (auto b : cb.cp().FilterAll(BinIsControlRegion).bin_set()) {
  //       std::cout << " - Adjusting data in bin " << b << "\n";
  //         cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
  //           TH1F new_obs = cb.cp().bin({b}).GetObservedShape();
  //           new_obs.SetBinContent(new_obs.GetNbinsX(), 1.);
  //           new_obs.Print("range");
  //           obs->set_shape(new_obs, true);
  //         });
  //       }
  // }

  // At this point we can fix the negative bins
  cb.ForEachProc([](ch::Process *p) {
    if (ch::HasNegativeBins(p->shape())) {
      std::cout << "[Negative bins] Fixing negative bins for " << p->bin()
                << "," << p->process() << "\n";
      // std::cout << "[Negative bins] Before:\n";
      // p->shape()->Print("range");
      auto newhist = p->ClonedShape();
      ch::ZeroNegativeBins(newhist.get());
      // Set the new shape but do not change the rate, we want the rate to still
      // reflect the total integral of the events
      p->set_shape(std::move(newhist), false);
      // std::cout << "[Negative bins] After:\n";
      // p->shape()->Print("range");
    }
  });

  cb.ForEachSyst([](ch::Systematic *s) {
    if (s->type().find("shape") == std::string::npos) return;
    if (ch::HasNegativeBins(s->shape_u()) || ch::HasNegativeBins(s->shape_d())) {
      std::cout << "[Negative bins] Fixing negative bins for syst" << s->bin()
                << "," << s->process() << "," << s->name() << "\n";
       std::cout << "[Negative bins] Before:\n";
       s->shape_u()->Print("range");
       s->shape_d()->Print("range");
      auto newhist_u = s->ClonedShapeU();
      auto newhist_d = s->ClonedShapeD();
      ch::ZeroNegativeBins(newhist_u.get());
      ch::ZeroNegativeBins(newhist_d.get());
      // Set the new shape but do not change the rate, we want the rate to still
      // reflect the total integral of the events
      s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
       std::cout << "[Negative bins] After:\n";
       s->shape_u()->Print("range");
       s->shape_d()->Print("range");
    }
  });

  cout << "Generating bbb uncertainties...";
  auto bbb = ch::BinByBinFactory()
    .SetPattern("CMS_$ANALYSIS_$BIN_$ERA_$PROCESS_bin_$#")
    .SetAddThreshold(0.)
    .SetMergeThreshold(0.4)
    .SetFixNorm(true)
    // .SetMergeZeroBins(false)
    .SetPoissonErrors(poisson_bbb);
  for (auto chn : chns) {
    std::cout << " - Doing bbb for channel " << chn << "\n";
    bbb.MergeAndAdd(cb.cp().channel({chn}).process({"WW", "DY", "VVV", "top", "FakeOF", "FakeSF", "ggWW", "VZ", "Vg", "VgS", "WW2J", "qqWWqq"}).FilterAll([](ch::Object const* obj) {
                return BinIsControlRegion(obj);
                }), cb);
  }
  // And now do bbb for the control region with a slightly different config:
  auto bbb_ctl = ch::BinByBinFactory()
    .SetPattern("CMS_$ANALYSIS_$BIN_$ERA_$PROCESS_bin_$#")
    .SetAddThreshold(0.)
    .SetMergeThreshold(0.4)
    .SetFixNorm(false)  // contrary to signal region, bbb *should* change yield here
    .SetVerbosity(1);
  // Will merge but only for non W and QCD processes, to be on the safe side
  bbb_ctl.MergeBinErrors(cb.cp().process({"DY", "top"}, false).FilterProcs(BinIsNotControlRegion));
  bbb_ctl.AddBinByBin(cb.cp().process({"DY", "top"}, false).FilterProcs(BinIsNotControlRegion), cb);
  cout << " done\n";

  //Switch JES over to lnN:
  //cb.cp().syst_name({"CMS_scale_j_13TeV"}).ForEachSyst([](ch::Systematic *sys) { sys->set_type("lnN");});

  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  // which is commonly used in the hww analyses
  ch::SetStandardBinNames(cb);
  //! [part8]


  //! [part9]
  // First we generate a set of bin names:
  RooWorkspace ws("hww", "hww");

  TFile demo("hww_mssm_demo.root", "RECREATE");

  bool do_morphing = true;
  map<string, RooAbsReal *> mass_var = {
    {"ggH_HWW", &mH}, {"qqH_HWW", &mH}, {"ggH_SBIHWW", &mH}, {"qqH_SBIHWW", &mH}//,
    //{"bbh_htautau", &mh}, {"bbH_Htautau", &mH}, {"bbA_Atautau", &mA}
  };
  if (!include_VBF){
mass_var = {
    {"ggH_HWW", &mH},{"ggH_SBIHWW", &mH}  };
  }
/*  if(mass=="MH"){
    mass_var = {
      {"ggH", &mA}//,
      //{"bbH", &mA}
    };
  }*/
  if (do_morphing) {
    auto bins = cb.bin_set();
    for (auto b : bins) {
      std::cout << "Bin:" << b << "\n";
      auto procs = cb.cp().bin({b}).process(ch::JoinStr({signal_types["ggH"], signal_types["qqH"], signal_types["ggH_SBI"], signal_types["qqH_SBI"]})).process_set();//, signal_types["bbH"]
    if (!include_VBF){
      procs = cb.cp().bin({b}).process(ch::JoinStr({signal_types["ggH"], signal_types["ggH_SBI"]})).process_set();
    }
      for (auto p : procs) {
        std::cout << "Process:" << p << "\n";
        ch::BuildRooMorphing(ws, cb, b, p, *(mass_var[p]),
                             "norm", true, false, false, &demo);
      }
    }
  }
  demo.Close();
  cb.AddWorkspace(ws);
  if (do_morphing) {
  if (include_VBF){
  cb.cp().process(ch::JoinStr({signal_types["ggH"], signal_types["qqH"], signal_types["ggH_SBI"], signal_types["qqH_SBI"]})).ExtractPdfs(cb, "hww", "$BIN_$PROCESS_morph");//, signal_types["bbH"]
  }else{
  cb.cp().process(ch::JoinStr({signal_types["ggH"], signal_types["ggH_SBI"]})).ExtractPdfs(cb, "hww", "$BIN_$PROCESS_morph");//, signal_types["bbH"]
  }
  }
  cb.PrintAll();


 //Write out datacards. Naming convention important for rest of workflow. We
 //make one directory per chn-cat, one per chn and cmb. In this code we only
 //store the individual datacards for each directory to be combined later, but
 //note that it's also possible to write out the full combined card with CH
  ch::CardWriter writer("output/" + output_folder + "/$TAG/$BIN.txt",
                        "output/" + output_folder + "/$TAG/$BIN_input.root");
  // We're not using mass as an identifier - which we need to tell the CardWriter
  // otherwise it will see "*" as the mass value for every object and skip it
  writer.SetWildcardMasses({});
  writer.SetVerbosity(1);

  writer.WriteCards("cmb", cb);
  for (auto chn : chns) {
    // per-channel
    writer.WriteCards(chn, cb.cp().channel({chn}));
    // And per-channel-category
    //writer.WriteCards("hww_"+chn+"_8_13TeV", cb.cp().channel({chn}).bin_id({8, 10, 11, 12}));
    //writer.WriteCards("hww_"+chn+"_9_13TeV", cb.cp().channel({chn}).bin_id({9, 13, 14, 15}));
  }
  // For btag/nobtag areas want to include control regions. This will
  // work even if the extra categories aren't there.
  writer.WriteCards("hww_cmb_0jet_13TeV", cb.cp().bin_id({8}));
  writer.WriteCards("hww_cmb_1jet_13TeV", cb.cp().bin_id({9}));
  writer.WriteCards("hww_cmb_2jet_13TeV", cb.cp().bin_id({10}));
  writer.WriteCards("hww_cmb_vbf_13TeV", cb.cp().bin_id({11}));
  //writer.WriteCards("hww_cmb_9_13TeV", cb.cp().bin_id({9, 13, 14, 15}));

  cb.PrintAll();
  cout << " done\n";
}
