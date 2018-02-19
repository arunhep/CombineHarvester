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

void AddMSSMRun2Systematics_HWW(CombineHarvester & cb, int control_region = 0) {
  // Create a CombineHarvester clone that only contains the signal
  // categories
  CombineHarvester cb_sig = cb.cp();

  std::vector<std::string> ggH = {"ggH", "ggH_HWW"};
  //std::vector<std::string> bbH = {"bbH", "bbH_Htautau", "bbA_Atautau", "bbh_htautau"};

  if (control_region == 1){
    // we only want to cosider systematic uncertainties in the signal region.
    // limit to only the btag and nobtag categories
    cb_sig.bin_id({8,9});
  }

  std::vector<std::string> SM_procsone = {"ggH_hww_SM125", "qqH_hww_SM125"};
  std::vector<std::string> SM_procstwo = {"ZH_hww_SM125", "ggZH_hww_SM125","WH_hww_SM125","bbH_hww_SM125"};
  std::vector<std::string> SM_tauprocs = {"ggH_htt_SM125", "qqH_htt_SM125", "ZH_htt_SM125", "WH_htt_SM125"};

  //auto signal = Set2Vec(cb.cp().signals().SetFromProcs(
  //    std::mem_fn(&Process::process)));
  std::vector<std::string> signal = {"ggH"};

  //signal = JoinStr({signal,SM_procs});


  cb.cp().process(JoinStr({signal, {"DY", "top", "WW", "ggWW", "Vg", "VgS", "VZ", "VVV"}})).AddSyst(cb,
    "lumi_13TeV", "lnN", SystMap<>::init(1.025));

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

  cb.cp().process(JoinStr({signal, {"DY", "WW", "ggWW", "VVV", "VZ", "top", "Vg", "VgS"}})).AddSyst(cb,
    "Full2016_btag_bc", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal, {"DY", "VVV", "VZ", "WW", "ggWW", "top", "Vg", "VgS"}})).AddSyst(cb,
    "Full2016_btag_udsg", "shape", SystMap<>::init(1.00));

  cb.cp().process(JoinStr({signal, {"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS"}})).AddSyst(cb,
    "trigger", "shape", SystMap<>::init(1.00));

  cb.cp().process(JoinStr({signal, {"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS"}})).AddSyst(cb,
    "idiso_ele", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal, {"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS"}})).AddSyst(cb,
    "scale_e", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal, {"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS"}})).AddSyst(cb,
    "elePtCor", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal, {"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS"}})).AddSyst(cb,
    "eleEtaCor", "shape", SystMap<>::init(1.00));

  cb.cp().process(JoinStr({signal, {"DY", "VVV", "VZ", "ggWW", "WW", "top", "Vg", "VgS"}})).AddSyst(cb,
    "idiso_mu", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal, {"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS"}})).AddSyst(cb,
    "scale_m", "shape", SystMap<>::init(1.00));

  cb.cp().process(JoinStr({signal, {"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS"}})).AddSyst(cb,
    "scale_j", "shape", SystMap<>::init(1.00));

  cb.cp().process(JoinStr({signal, {"DY", "ggWW", "WW", "top", "VZ", "VVV", "Vg", "VgS"}})).AddSyst(cb,
    "scale_met", "shape", SystMap<>::init(1.00));

  cb.cp().process(JoinStr({signal, {"WW"}})).AddSyst(cb,
    "PS", "shape", SystMap<>::init(1.00));
  cb.cp().process(JoinStr({signal, {"WW"}})).AddSyst(cb,
    "UE", "shape", SystMap<>::init(1.00));

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

  cb.cp().process(JoinStr({signal, {"ggWW"}})).AddSyst(cb,
    "QCDscale_gg_accept", "lnN", SystMap<>::init(1.027));

  cb.cp().process({"ggH_SM125"}).AddSyst(cb,
    "pdf_gg","lnN", SystMap<>::init(1.031));

  cb.cp().process({"ggH_SM125", "ggWW"}).AddSyst(cb,
    "pdf_gg_accept","lnN", SystMap<>::init(1.005));

  cb.cp().process({"DY"}).AddSyst(cb,
    "DYQCDscale", "shape", SystMap<>::init(1.00));

  cb.cp().process({"top"}).AddSyst(cb,
    "TopPtRew", "shape", SystMap<>::init(1.00));

  cb.cp().process({"ggWW"}).AddSyst(cb,
    "kfactggww","lnN", SystMap<>::init(1.15));

  cb.cp().process({"WW"}).AddSyst(cb,
    "WWresum0j", "shape", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWresum1j", "shape", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWresum2j", "shape", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWresumvbf", "shape", SystMap<bin_id>::init({11}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWqscale0j", "shape", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWqscale1j", "shape", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWqscale2j", "shape", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWqscalevbf", "shape", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"VgS"}).AddSyst(cb,
    "WgStarScale","lnN", SystMap<>::init(1.25));

  cb.cp().process({"DY"}).AddSyst(cb,
    "DYttnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"DY"}).AddSyst(cb,
    "DYttnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"DY"}).AddSyst(cb,
    "DYttnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"DY"}).AddSyst(cb,
    "DYttnormvbf", "rateParam", SystMap<bin_id>::init({11}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"WW"}).AddSyst(cb,
    "WWnormvbf", "rateParam", SystMap<bin_id>::init({11}, 1.00));
  cb.cp().process({"top"}).AddSyst(cb,
    "Topnorm0j", "rateParam", SystMap<bin_id>::init({8}, 1.00));
  cb.cp().process({"top"}).AddSyst(cb,
    "Topnorm1j", "rateParam", SystMap<bin_id>::init({9}, 1.00));
  cb.cp().process({"top"}).AddSyst(cb,
    "Topnorm2j", "rateParam", SystMap<bin_id>::init({10}, 1.00));
  cb.cp().process({"top"}).AddSyst(cb,
    "Topnormvbf", "rateParam", SystMap<bin_id>::init({11}, 1.00));

  cb.cp().process({"top"}).AddSyst(cb,
    "tttwTh", "shape", SystMap<>::init(1.00));
  cb.cp().process({"top"}).AddSyst(cb,
    "TopPtRew", "shape", SystMap<>::init(1.00));

  //SM theory uncertainties
  cb.cp().process({"ggH_htt_SM125"}).AddSyst(cb, "QCDscale_ggH", "lnN",
    SystMapAsymm<>::init(1.076,0.919));

  cb.cp().process({"qqH_htt_SM125"}).AddSyst(cb, "QCDscale_qqH", "lnN",
    SystMapAsymm<>::init(1.004,0.997));


  cb.cp().process({"ZH_htt_SM125","WH_htt_SM125"}).AddSyst(cb, "QCDscale_VH", "lnN",
    SystMapAsymm<process>::init
     ({"ZH_htt_SM125"}, 1.038, 0.969)
     ({"WH_htt_SM125"},1.005,0.993));

/*  cb.cp().process({"TTH_SM125"}).AddSyst(cb,"QCDscale_ttH","lnN", 
    SystMapAsymm<>::init(1.058,0.908));*/

  cb.cp().process({"ggH_htt_SM125"}).AddSyst(cb, "pdf_Higgs_gg","lnN",
    SystMap<>::init(1.031));

/*  cb.cp().process({"TTH_SM125"}).AddSyst(cb, "pdf_Higgs_ttH","lnN",
    SystMap<>::init(1.036));*/

  cb.cp().process({"ZH_htt_SM125","WH_htt_SM125","qqH_htt_SM125"}).AddSyst(cb, "pdf_Higgs_qqbar","lnN",
    SystMap<process>::init
      ({"ZH_htt_SM125"},1.016)
      ({"WH_htt_SM125"},1.019)
      ({"qqH_htt_SM125"},1.021)); 
    
  }
}
