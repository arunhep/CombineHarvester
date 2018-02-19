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

  cb.cp().process({"DY"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
    "DYMetRew", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggWW", "Vg", "VgS", "VZ", "VVV", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "lumi_13TeV", "lnN", SystMap<>::init(1.025));

  cb.cp().process({"ggH_hww"}).bin({"8"}).AddSyst(cb,
    "QCDscale_ggH0j","lnN", SystMap<>::init(1.08539));

  cb.cp().process({"H_htt"}).AddSyst(cb, "QCDscale_ggH", "lnN",
    SystMapAsymm<>::init(0.919,1.076));

  cb.cp().process({"qqH", "qqH_hww"}).AddSyst(cb, "QCDscale_qqH", "lnN",
    SystMapAsymm<process>::init
    ({"qqH_hww"},0.997,1.004)
    ({"qqH"},0.973,1.012));

  cb.cp().process({"WH_hww", "ZH_hww"}).AddSyst(cb, "QCDscale_VH", "lnN",
    SystMapAsymm<process>::init
    ({"WH_hww"},0.993,1.005)
    ({"ZH_hww"},0.969,1.038));

  cb.cp().process({"ggZH_hww"}).AddSyst(cb, "QCDscale_ggZH", "lnN",
    SystMapAsymm<>::init(0.919,1.076));

  cb.cp().process({"qqH_hww", "WH_hww", "ZH_hww", "qqH"}).AddSyst(cb,
    "QCDscale_qqbar_accept", "lnN", SystMap<process>::init
    ({"qqH_hww"},1.007)
    ({"WH_hww"},1.05)
    ({"ZH_hww"},1.04)
    ({"qqH"},1.02));

  cb.cp().process({"ggWW", "ggH_hww", "H_htt", "ggZH_hww", "ggH", "ggH_INT"}).AddSyst(cb,
    "QCDscale_gg_accept", "lnN", SystMap<>::init(1.027));

  cb.cp().process({"ggH_hww"}).bin({"8","9","10","11"}).AddSyst(cb,
    "QCDscale", "shape", SystMap<>::init(1.00));
  cb.cp().process({"ggH_hww"}).bin({"8","9","10","11"}).AddSyst(cb,
    "QCDscale1in", "shape", SystMap<>::init(1.00));
  cb.cp().process({"ggH_hww"}).bin({"8","9","10","11"}).AddSyst(cb,
    "QCDscale2in", "shape", SystMap<>::init(1.00));
  cb.cp().process({"ggH_hww"}).bin({"8","9","10","11"}).AddSyst(cb,
    "QCDscale3in", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggH_hww", "H_htt", "ggZH_hww", "ggH", "ggH_INT"}).AddSyst(cb, "pdf_gg", "lnN",
    SystMapAsymm<process>::init
    ({"ggH_hww", "H_htt", "ggZH_hww"},0.919,1.076)
    ({"ggH", "ggH_INT"},0.945,1.034));

  cb.cp().process({"qqH_hww", "WH_hww", "ZH_hww", "qqH"}).AddSyst(cb, "pdf_qqbar", "lnN",
    SystMapAsymm<process>::init
    ({"qqH_hww"},0.997,1.004)
    ({"WH_hww"},0.993,1.005)
    ({"ZH_hww"},0.969,1.038)
    ({"qqH"},0.973,1.012));

  cb.cp().process({"ggWW", "ggH_hww", "H_htt", "ggZH_hww", "ggH", "ggH_INT"}).AddSyst(cb,
    "pdf_gg_accept","lnN", SystMap<>::init(1.005));

  cb.cp().process({"qqH_hww", "WH_hww", "ZH_hww", "qqH"}).AddSyst(cb,
    "pdf_qqbar_accept","lnN", SystMap<process>::init
    ({"qqH_hww"},1.011)
    ({"WH_hww"},1.007)
    ({"ZH_hww"},1.012)
    ({"qqH"},1.011));

  cb.cp().process({"ggWW"}).AddSyst(cb,
    "kfactggww","lnN", SystMap<>::init(1.15));

  cb.cp().process({"WW"}).AddSyst(cb,
    "WWresum0j", "shape", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWresum1j", "shape", SystMap<bin_id>::init({9}, 1.00));
  //cb.cp().process({"WW"}).AddSyst(cb,
  //  "WWresum2j", "shape", SystMap<bin_id>::init({10, 11}, 1.00));

  cb.cp().process({"WW"}).AddSyst(cb,
    "WWqscale0j", "shape", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWqscale1j", "shape", SystMap<bin_id>::init({9}, 1.00));
  //cb.cp().process({"WW"}).AddSyst(cb,
  //  "WWqscale2j", "shape", SystMap<bin_id>::init({10, 11}, 1.00));

  cb.cp().process({"WW"}).AddSyst(cb,
    "UE", "shape", SystMap<>::init(1.00));//, "ggH_hww", "qqH_hww"

  cb.cp().process({"WW"}).AddSyst(cb,
    "PS", "shape", SystMap<>::init(1.00));

  cb.cp().process({"VgS"}).AddSyst(cb,
    "WgStarScale","lnN", SystMap<>::init(1.25));

  cb.cp().process({"DY"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_DYttnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"DY"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_DYttnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"DY"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_DYttnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"DY"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_DYttnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

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

  cb.cp().process({"WW"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_WWofnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_WWnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"WW"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_WWnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"WW"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_WWnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

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

  cb.cp().process({"top"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_Topnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"top"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_Topnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"top"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_Topnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"top"}).channel({"em","dytt","top"}).AddSyst(cb,
    "CMS_hwwhmof_Topnorm2jVBF", "rateParam", SystMap<bin_id>::init({11}, 1.00));

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

  cb.cp().process({"top"}).channel({"em","dytt","top"}).AddSyst(cb,
    "singleTopToTTbar", "shape", SystMap<>::init(1.00));
  cb.cp().process({"top"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
    "tttwTh", "shape", SystMap<>::init(1.00));

  cb.cp().process({"top"}).AddSyst(cb,
    "TopPS", "shape", SystMap<>::init(1.00));

  //cb.cp().process({"DY"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
  //  "DYptRew", "shape", SystMap<>::init(1.00));

  cb.cp().process({"Fake"}).AddSyst(cb,
    "fake_syst", "lnN", SystMap<>::init(1.30));

  cb.cp().process({"Fake"}).AddSyst(cb,
    "fake_ele_hww", "shape", SystMap<>::init(1.00));

  cb.cp().process({"Fake"}).AddSyst(cb,
    "fake_ele_stat_hww", "shape", SystMap<>::init(1.00));

  cb.cp().process({"Fake"}).AddSyst(cb,
    "fake_mu_hww", "shape", SystMap<>::init(1.00));

  cb.cp().process({"Fake"}).AddSyst(cb,
    "fake_mu_stat_hww", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "WW", "ggWW", "VVV", "VZ", "top", "Vg", "VgS", "WH_hww", "ZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "btag_heavy", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "VVV", "VZ", "WW", "ggWW", "top", "Vg", "VgS", "WH_hww", "ZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "btag_light", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "trigger", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "eff_e", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "scale_e", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "hww_elePtCor", "shape", SystMap<>::init(1.00));//"ggH_hww", "qqH_hww", 

  cb.cp().process({"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "hww_eleEtaCor", "shape", SystMap<>::init(1.00));//"ggH_hww", "qqH_hww", 

  cb.cp().process({"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "eff_m", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggWW", "WW", "DY", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "scale_m", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggWW", "WW", "DY", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).AddSyst(cb,
    "scale_j", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggWW", "WW", "DY", "top", "VZ", "VVV", "Vg", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).channel({"em","dytt","top"}).AddSyst(cb,
    "scale_met", "shape", SystMap<>::init(1.00));
  cb.cp().process({"ggWW", "WW", "DY", "top", "VZ", "VVV", "VgS", "WH_hww", "ZH_hww", "ggZH_hww", "ggH_hww", "qqH_hww", "H_htt", "ggH", "ggH_INT", "qqH"}).channel({"ee","eedytt","eetop","mm","mmdytt","mmtop"}).AddSyst(cb,
    "scale_met", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY"}).channel({"em","dytt","top"}).AddSyst(cb,
    "QCDscale_V", "shape", SystMap<>::init(1.00));











/*

  cb.cp().process(JoinStr({signal, {"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS"}})).AddSyst(cb,
    "elePtCor", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal, {"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS"}})).AddSyst(cb,
    "eleEtaCor", "shape", SystMap<>::init(1.00));





  cb.cp().process({"DY"}).AddSyst(cb,
    "DYQCDscale", "shape", SystMap<>::init(1.00));

  cb.cp().process(JoinStr({signal})).AddSyst(cb,
    "THU_ggH_Mu", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal})).AddSyst(cb,
    "THU_ggH_Res", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal})).AddSyst(cb,
    "THU_ggH_Mig01", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal})).AddSyst(cb,
    "THU_ggH_Mig12", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal})).AddSyst(cb,
    "THU_ggH_PT60", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal})).AddSyst(cb,
    "THU_ggH_PT120", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal})).AddSyst(cb,
    "THU_ggH_VBF2j", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal})).AddSyst(cb,
    "THU_ggH_VBF3j", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal})).AddSyst(cb,
    "THU_ggH_qmtop", "shape", SystMap<>::init(1.00));

  cb.cp().process({"DY"}).AddSyst(cb,
    "QCDscale_CRSR_accept_dytt", "lnN", SystMap<>::init(1.02));

  cb.cp().process({"top"}).AddSyst(cb,
    "QCDscale_CRSR_accept_top", "lnN", SystMap<>::init(1.01));



  cb.cp().process({"DY"}).AddSyst(cb,
    "DYQCDscale", "shape", SystMap<>::init(1.00));

  cb.cp().process({"top"}).AddSyst(cb,
    "TopPtRew", "shape", SystMap<>::init(1.00));


  cb.cp().process({"WW"}).AddSyst(cb,
    "WWresumvbf", "shape", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"WW"}).AddSyst(cb,
    "WWqscalevbf", "shape", SystMap<bin_id>::init({11}, 1.00));
*/



/*  cb.cp().process({"TTH_SM125"}).AddSyst(cb,"QCDscale_ttH","lnN", 
    SystMapAsymm<>::init(1.058,0.908));*/
/*
  cb.cp().process({"ggH_htt_SM125"}).AddSyst(cb, "pdf_Higgs_gg","lnN",
    SystMap<>::init(1.031));
*/
/*  cb.cp().process({"TTH_SM125"}).AddSyst(cb, "pdf_Higgs_ttH","lnN",
    SystMap<>::init(1.036));*/
/*
  cb.cp().process({"ZH_htt_SM125","WH_htt_SM125","qqH_htt_SM125"}).AddSyst(cb, "pdf_Higgs_qqbar","lnN",
    SystMap<process>::init
      ({"ZH_htt_SM125"},1.016)
      ({"WH_htt_SM125"},1.019)
      ({"qqH_htt_SM125"},1.021)); 
    */
  }

}
