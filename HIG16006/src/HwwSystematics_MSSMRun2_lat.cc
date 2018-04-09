#include "CombineHarvester/HIG16006/interface/HwwSystematics_MSSMRun2.h"
#include <vector>
#include <string>
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"

namespace ch {

using ch::syst::SystMap;
using ch::syst::SystMapAsymm;
using ch::syst::era;
using ch::syst::channel;
using ch::syst::bin_id;
using ch::syst::process;
using ch::syst::bin;
using ch::syst::mass;
using ch::JoinStr;

void AddMSSMRun2Systematics_HWW_lat(CombineHarvester & cb, int control_region = 0) {
  // Create a CombineHarvester clone that only contains the signal
  // categories
  CombineHarvester cb_sig = cb.cp();

  //std::vector<std::string> ggH = {"ggH", "ggH_HWW"};
  //std::vector<std::string> bbH = {"bbH", "bbH_Htautau", "bbA_Atautau", "bbh_htautau"};
/*
  if (control_region == 1){
    // we only want to cosider systematic uncertainties in the signal region.
    // limit to only the btag and nobtag categories
    cb_sig.bin_id({8,9});
  }
*/
  //std::vector<std::string> SM_procsone = {"ggH_hww_SM125", "qqH_hww_SM125"};
  //std::vector<std::string> SM_procstwo = {"ZH_hww_SM125", "ggZH_hww_SM125","WH_hww_SM125","bbH_hww_SM125"};
  //std::vector<std::string> SM_tauprocs = {"ggH_htt_SM125", "qqH_htt_SM125", "ZH_htt_SM125", "WH_htt_SM125"};

  //auto signal = Set2Vec(cb.cp().signals().SetFromProcs(
  //    std::mem_fn(&Process::process)));
  //std::vector<std::string> signal = {"ggH", "qqH"};

  //signal = JoinStr({signal,SM_procs});

  cb.cp().process({"top"}).AddSyst(cb,
    "TopPtRew", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggH_hww", "H_htt"}).AddSyst(cb, "QCDscale_ggH", "lnN",
    SystMapAsymm<process>::init
    ({"ggH_hww"},0.919,1.076)
    ({"H_htt"},0.919,1.076));

  cb.cp().process({"qqH", "qqH_hww"}).AddSyst(cb, "QCDscale_qqH", "lnN",
    SystMapAsymm<process, mass>::init
    ({"qqH_hww"},{"*"},0.997,1.004)
    ({"qqH"},{"200"},0.998,1.003)
    ({"qqH"},{"210"},0.998,1.003)
    ({"qqH"},{"230"},0.999,1.003)
    ({"qqH"},{"250"},0.999,1.003)
    ({"qqH"},{"270"},0.999,1.003)
    ({"qqH"},{"300"},0.999,1.003)
    ({"qqH"},{"350"},0.999,1.003)
    ({"qqH"},{"400"},0.999,1.003)
    ({"qqH"},{"450"},0.999,1.003)
    ({"qqH"},{"500"},0.999,1.003)
    ({"qqH"},{"550"},0.998,1.003)
    ({"qqH"},{"650"},0.997,1.003)
    ({"qqH"},{"700"},0.996,1.003)
    ({"qqH"},{"750"},0.996,1.003)
    ({"qqH"},{"800"},0.995,1.003)
    ({"qqH"},{"900"},0.994,1.003)
    ({"qqH"},{"1000"},0.993,1.003)
    ({"qqH"},{"1500"},0.989,1.004)
    ({"qqH"},{"2000"},0.984,1.007)
    ({"qqH"},{"2500"},0.979,1.010)
    ({"qqH"},{"3000"},0.973,1.012));

  cb.cp().process({"WH_hww", "ZH_hww"}).AddSyst(cb, "QCDscale_VH", "lnN",
    SystMapAsymm<process>::init
    ({"WH_hww"},0.993,1.005)
    ({"ZH_hww"},0.969,1.038));

  cb.cp().process({"ggZH_hww"}).AddSyst(cb, "QCDscale_ggZH", "lnN",
    SystMapAsymm<>::init(0.919,1.076));

  cb.cp().process({"qqH_hww", "WH_hww", "ZH_hww", "VZ", "qqH", "qqH_SBI"}).AddSyst(cb,
    "QCDscale_qqbar_accept", "lnN", SystMap<process, mass>::init
    ({"qqH_hww"},{"*"},1.007)
    ({"WH_hww"},{"*"},1.05)
    ({"ZH_hww"},{"*"},1.04)
    ({"VZ"},{"*"},1.029)
    ({"qqH", "qqH_SBI"},{"200"},1.003)
    ({"qqH", "qqH_SBI"},{"210"},1.003)
    ({"qqH", "qqH_SBI"},{"230"},1.003)
    ({"qqH", "qqH_SBI"},{"250"},1.003)
    ({"qqH", "qqH_SBI"},{"270"},1.003)
    ({"qqH", "qqH_SBI"},{"300"},1.003)
    ({"qqH", "qqH_SBI"},{"350"},1.003)
    ({"qqH", "qqH_SBI"},{"400"},1.003)
    ({"qqH", "qqH_SBI"},{"450"},1.004)
    ({"qqH", "qqH_SBI"},{"500"},1.007)
    ({"qqH", "qqH_SBI"},{"550"},1.009)
    ({"qqH", "qqH_SBI"},{"650"},1.014)
    ({"qqH", "qqH_SBI"},{"700"},1.016)
    ({"qqH", "qqH_SBI"},{"750"},1.018)
    ({"qqH", "qqH_SBI"},{"800"},1.020)
    ({"qqH", "qqH_SBI"},{"900"},1.023)
    ({"qqH", "qqH_SBI"},{"1000"},1.027)
    ({"qqH", "qqH_SBI"},{"1500"},1.037)
    ({"qqH", "qqH_SBI"},{"2000"},1.042)
    ({"qqH", "qqH_SBI"},{"2500"},1.048)
    ({"qqH", "qqH_SBI"},{"3000"},1.060));

  cb.cp().process({"ggWW", "ggH_hww", "H_htt", "ggZH_hww", "ggH", "ggH_SBI"}).AddSyst(cb,
    "QCDscale_gg_accept", "lnN", SystMap<process, mass>::init
    ({"ggWW"},{"*"},1.027)
    ({"ggH_hww"},{"*"},1.027)
    ({"H_htt"},{"*"},1.027)
    ({"ggZH_hww"},{"*"},1.027)
    ({"ggH", "ggH_SBI"},{"200"},1.054)
    ({"ggH", "ggH_SBI"},{"210"},1.055)
    ({"ggH", "ggH_SBI"},{"230"},1.055)
    ({"ggH", "ggH_SBI"},{"250"},1.056)
    ({"ggH", "ggH_SBI"},{"270"},1.056)
    ({"ggH", "ggH_SBI"},{"300"},1.057)
    ({"ggH", "ggH_SBI"},{"350"},1.058)
    ({"ggH", "ggH_SBI"},{"400"},1.059)
    ({"ggH", "ggH_SBI"},{"450"},1.059)
    ({"ggH", "ggH_SBI"},{"500"},1.060)
    ({"ggH", "ggH_SBI"},{"550"},1.061)
    ({"ggH", "ggH_SBI"},{"650"},1.062)
    ({"ggH", "ggH_SBI"},{"700"},1.062)
    ({"ggH", "ggH_SBI"},{"750"},1.062)
    ({"ggH", "ggH_SBI"},{"800"},1.062)
    ({"ggH", "ggH_SBI"},{"900"},1.063)
    ({"ggH", "ggH_SBI"},{"1000"},1.062)
    ({"ggH", "ggH_SBI"},{"1500"},1.058)
    ({"ggH", "ggH_SBI"},{"2000"},1.050)
    ({"ggH", "ggH_SBI"},{"2500"},1.041)
    ({"ggH", "ggH_SBI"},{"3000"},1.035));

  cb.cp().process({"ggH_hww", "ggH"}).bin({"8","9","10","11"}).AddSyst(cb,
    "QCDscale", "shape", SystMap<>::init(1.00));
  cb.cp().process({"ggH_hww", "ggH"}).bin({"8","9","10","11"}).AddSyst(cb,
    "QCDscale1in", "shape", SystMap<>::init(1.00));
  cb.cp().process({"ggH_hww", "ggH"}).bin({"8","9","10","11"}).AddSyst(cb,
    "QCDscale2in", "shape", SystMap<>::init(1.00));
  cb.cp().process({"ggH_hww", "ggH"}).bin({"8","9","10","11"}).AddSyst(cb,
    "QCDscale3in", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggH_hww", "H_htt", "ggZH_hww", "ggH", "ggH_SBI"}).AddSyst(cb, "pdf_gg", "lnN",
    SystMap<process, mass>::init
    ({"ggH_hww", "H_htt"},{"*"},1.018)
    ({"ggZH_hww"},{"*"},1.013)
    ({"ggH", "ggH_SBI"},{"200"},1.017)
    ({"ggH", "ggH_SBI"},{"210"},1.017)
    ({"ggH", "ggH_SBI"},{"230"},1.017)
    ({"ggH", "ggH_SBI"},{"250"},1.018)
    ({"ggH", "ggH_SBI"},{"270"},1.018)
    ({"ggH", "ggH_SBI"},{"300"},1.018)
    ({"ggH", "ggH_SBI"},{"350"},1.018)
    ({"ggH", "ggH_SBI"},{"400"},1.019)
    ({"ggH", "ggH_SBI"},{"450"},1.020)
    ({"ggH", "ggH_SBI"},{"500"},1.022)
    ({"ggH", "ggH_SBI"},{"550"},1.023)
    ({"ggH", "ggH_SBI"},{"650"},1.026)
    ({"ggH", "ggH_SBI"},{"700"},1.028)
    ({"ggH", "ggH_SBI"},{"750"},1.030)
    ({"ggH", "ggH_SBI"},{"800"},1.032)
    ({"ggH", "ggH_SBI"},{"900"},1.036)
    ({"ggH", "ggH_SBI"},{"1000"},1.040)
    ({"ggH", "ggH_SBI"},{"1500"},1.063)
    ({"ggH", "ggH_SBI"},{"2000"},1.087)
    ({"ggH", "ggH_SBI"},{"2500"},1.116)
    ({"ggH", "ggH_SBI"},{"3000"},1.150));

  cb.cp().process({"qqH_hww", "WH_hww", "ZH_hww", "qqH", "qqH_SBI"}).AddSyst(cb, "pdf_qqbar", "lnN",
    SystMap<process, mass>::init
    ({"qqH_hww"},{"*"},1.021)
    ({"WH_hww"},{"*"},1.017)
    ({"ZH_hww"},{"*"},1.013)
    ({"qqH", "qqH_SBI"},{"200"},1.018)
    ({"qqH", "qqH_SBI"},{"210"},1.018)
    ({"qqH", "qqH_SBI"},{"230"},1.018)
    ({"qqH", "qqH_SBI"},{"250"},1.018)
    ({"qqH", "qqH_SBI"},{"270"},1.019)
    ({"qqH", "qqH_SBI"},{"300"},1.019)
    ({"qqH", "qqH_SBI"},{"350"},1.019)
    ({"qqH", "qqH_SBI"},{"400"},1.019)
    ({"qqH", "qqH_SBI"},{"450"},1.020)
    ({"qqH", "qqH_SBI"},{"500"},1.020)
    ({"qqH", "qqH_SBI"},{"550"},1.021)
    ({"qqH", "qqH_SBI"},{"650"},1.021)
    ({"qqH", "qqH_SBI"},{"700"},1.022)
    ({"qqH", "qqH_SBI"},{"750"},1.022)
    ({"qqH", "qqH_SBI"},{"800"},1.023)
    ({"qqH", "qqH_SBI"},{"900"},1.024)
    ({"qqH", "qqH_SBI"},{"1000"},1.025)
    ({"qqH", "qqH_SBI"},{"1500"},1.030)
    ({"qqH", "qqH_SBI"},{"2000"},1.036)
    ({"qqH", "qqH_SBI"},{"2500"},1.045)
    ({"qqH", "qqH_SBI"},{"3000"},1.057));

  cb.cp().process({"ggWW", "ggH_hww", "H_htt", "ggZH_hww", "ggH", "ggH_SBI"}).AddSyst(cb, "pdf_gg_accept","lnN",
    SystMap<process, mass>::init
    ({"ggWW", "ggH_hww", "H_htt", "ggZH_hww"},{"*"},1.005)
    ({"ggH"},{"200"},1.007)
    ({"ggH"},{"210"},1.007)
    ({"ggH"},{"230"},1.007)
    ({"ggH"},{"250"},1.007)
    ({"ggH"},{"270"},1.007)
    ({"ggH"},{"300"},1.007)
    ({"ggH"},{"350"},1.007)
    ({"ggH"},{"400"},1.007)
    ({"ggH"},{"450"},1.007)
    ({"ggH"},{"500"},1.007)
    ({"ggH"},{"550"},1.007)
    ({"ggH"},{"650"},1.007)
    ({"ggH"},{"700"},1.007)
    ({"ggH"},{"750"},1.007)
    ({"ggH"},{"800"},1.007)
    ({"ggH"},{"900"},1.007)
    ({"ggH"},{"1000"},1.007)
    ({"ggH"},{"1500"},1.012)
    ({"ggH"},{"2000"},1.012)
    ({"ggH"},{"2500"},1.012)
    ({"ggH"},{"3000"},1.012)
    ({"ggH_SBI"},{"200"},1.010)
    ({"ggH_SBI"},{"210"},1.010)
    ({"ggH_SBI"},{"230"},1.010)
    ({"ggH_SBI"},{"250"},1.010)
    ({"ggH_SBI"},{"270"},1.010)
    ({"ggH_SBI"},{"300"},1.010)
    ({"ggH_SBI"},{"350"},1.010)
    ({"ggH_SBI"},{"400"},1.010)
    ({"ggH_SBI"},{"450"},1.010)
    ({"ggH_SBI"},{"500"},1.010)
    ({"ggH_SBI"},{"550"},1.010)
    ({"ggH_SBI"},{"650"},1.010)
    ({"ggH_SBI"},{"700"},1.010)
    ({"ggH_SBI"},{"750"},1.010)
    ({"ggH_SBI"},{"800"},1.010)
    ({"ggH_SBI"},{"900"},1.010)
    ({"ggH_SBI"},{"1000"},1.010)
    ({"ggH_SBI"},{"1500"},1.035)
    ({"ggH_SBI"},{"2000"},1.035)
    ({"ggH_SBI"},{"2500"},1.035)
    ({"ggH_SBI"},{"3000"},1.035));

  cb.cp().process({"qqH_hww", "WH_hww", "ZH_hww", "VZ", "qqH", "qqH_SBI"}).AddSyst(cb,
    "pdf_qqbar_accept","lnN", SystMap<process, mass>::init
    ({"qqH_hww"},{"*"},1.011)
    ({"WH_hww"},{"*"},1.007)
    ({"ZH_hww"},{"*"},1.012)
    ({"VZ"},{"*"},1.005)
    ({"qqH", "qqH_SBI"},{"200"},1.005)
    ({"qqH", "qqH_SBI"},{"210"},1.005)
    ({"qqH", "qqH_SBI"},{"230"},1.005)
    ({"qqH", "qqH_SBI"},{"250"},1.005)
    ({"qqH", "qqH_SBI"},{"270"},1.005)
    ({"qqH", "qqH_SBI"},{"300"},1.005)
    ({"qqH", "qqH_SBI"},{"350"},1.005)
    ({"qqH", "qqH_SBI"},{"400"},1.005)
    ({"qqH", "qqH_SBI"},{"450"},1.005)
    ({"qqH", "qqH_SBI"},{"500"},1.005)
    ({"qqH", "qqH_SBI"},{"550"},1.005)
    ({"qqH", "qqH_SBI"},{"650"},1.005)
    ({"qqH", "qqH_SBI"},{"700"},1.005)
    ({"qqH", "qqH_SBI"},{"750"},1.005)
    ({"qqH", "qqH_SBI"},{"800"},1.005)
    ({"qqH", "qqH_SBI"},{"900"},1.005)
    ({"qqH", "qqH_SBI"},{"1000"},1.015)
    ({"qqH", "qqH_SBI"},{"1500"},1.015)
    ({"qqH", "qqH_SBI"},{"2000"},1.015)
    ({"qqH", "qqH_SBI"},{"2500"},1.015)
    ({"qqH", "qqH_SBI"},{"3000"},1.015));

  cb.cp().process({"ggWW"}).AddSyst(cb,
    "QCDscale_ggWW","lnN", SystMap<>::init(1.15));

  cb.cp().process({"WW"}).AddSyst(cb,
    "WWresum0j", "shape", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWresum1j", "shape", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"WW"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb, //"em","dytt","top",
    "WWresum2j", "shape", SystMap<bin_id>::init({10, 11}, 1.00));

  cb.cp().process({"WW"}).AddSyst(cb,
    "WWqscale0j", "shape", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWqscale1j", "shape", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"WW"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb, //"em","dytt","top",
    "WWqscale2j", "shape", SystMap<bin_id>::init({10, 11}, 1.00));

  cb.cp().process({"WW"}).AddSyst(cb,
    "PS", "shape", SystMap<>::init(1.00));

  cb.cp().process({"WW"}).AddSyst(cb,
    "UE", "shape", SystMap<>::init(1.00));

  cb.cp().process({"WgS", "VgS"}).AddSyst(cb,
    "WgStarScale","lnN", SystMap<>::init(1.25));

  cb.cp().process({"DY"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_DYttnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"DY"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_DYttnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"DY"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_DYttnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"DY"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_DYttnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"WW"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_WWofnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_WWnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"WW"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_WWnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"WW"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_WWnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"top"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_Topnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"top"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_Topnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"top"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_Topnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"top"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_Topnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"top"}).channel({"em","dytt","top"}).AddSyst(cb,
    "singleTopToTTbar", "shape", SystMap<>::init(1.00));

  cb.cp().process({"top"}).AddSyst(cb,
    "TopPS", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggWW", "Vg", "VgS", "VZ", "VVV", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "lumi_13TeV", "lnN", SystMap<>::init(1.025));

  cb.cp().process({"Fake"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwof_fake_syst", "lnN", SystMap<>::init(1.30));

  cb.cp().process({"Fake"}).AddSyst(cb,
    "fake_ele_hww", "shape", SystMap<>::init(1.00));

  cb.cp().process({"Fake"}).AddSyst(cb,
    "fake_ele_stat_hww", "shape", SystMap<>::init(1.00));

  cb.cp().process({"Fake"}).AddSyst(cb,
    "fake_mu_hww", "shape", SystMap<>::init(1.00));

  cb.cp().process({"Fake"}).AddSyst(cb,
    "fake_mu_stat_hww", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "WW", "ggWW", "VVV", "VZ", "top", "Vg", "VgS", "WH_hww", "ZH_hww", "ggH_hww", "qqH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "btag_heavy", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "VVV", "VZ", "WW", "ggWW", "top", "Vg", "VgS", "WH_hww", "ZH_hww", "ggH_hww", "qqH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "btag_light", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "trigger", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "eff_e", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "scale_e", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "hww_elePtCor", "shape", SystMap<>::init(1.00));//"ggH_hww", "qqH_hww", 

  cb.cp().process({"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "hww_eleEtaCor", "shape", SystMap<>::init(1.00));//"ggH_hww", "qqH_hww", 

  cb.cp().process({"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "eff_m", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggWW", "WW", "DY", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "scale_m", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggWW", "WW", "DY", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).AddSyst(cb,
    "scale_j", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggWW", "WW", "DY", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "bbH_hww", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).channel({"em","dytt","top"}).AddSyst(cb,
    "scale_met", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggWW", "WW", "DY", "top", "VZ", "VVV", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_SBI", "qqH", "qqH_SBI"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
    "scale_met", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY"}).channel({"em","dytt","top"}).AddSyst(cb,
    "QCDscale_V", "shape", SystMap<>::init(1.00));

  cb.cp().process({"qqWWqq", "WW2J"}).AddSyst(cb,
    "WW2J_QCDscale_V", "shape", SystMap<>::init(1.00));



  cb.cp().process({"DY"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
    "DYMetRew", "shape", SystMap<>::init(1.00));

  cb.cp().process({"VW"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
    "QCDscale_VW","lnN", SystMap<>::init(1.03));

  cb.cp().process({"ggH_hww"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).bin({"8"}).AddSyst(cb,
    "QCDscale_ggH0j","lnN", SystMap<>::init(1.08539));

  cb.cp().process({"ggWW"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
    "kfactggww","lnN", SystMap<>::init(1.15));

  cb.cp().process({"DY"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_DYnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"DY"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_DYnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"DY"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_DYnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"DY"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_DYnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"DY"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_DYnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"DY"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_DYnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"DY"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_DYnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"DY"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_DYnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"WW"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_WWnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_WWnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"WW"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_WWnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"WW"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_WWnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"WW"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_WWnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_WWnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"WW"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_WWnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"WW"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_WWnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"top"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_Topnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"top"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_Topnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"top"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_Topnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"top"}).channel({"ee","eedytt","eetop"}).AddSyst(cb,
    "CMS_hwwhmee_Topnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"top"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_Topnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"top"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_Topnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"top"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_Topnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"top"}).channel({"mm","mmdytt","mmtop"}).AddSyst(cb,
    "CMS_hwwhmmm_Topnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"Fake"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
    "fake_syst", "lnN", SystMap<>::init(1.30));

  cb.cp().process({"top"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
    "tttwTh", "shape", SystMap<>::init(1.00));

  //cb.cp().process({"DY"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
  //  "DYptRew", "shape", SystMap<>::init(1.00));

  }

}
