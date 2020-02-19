#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
    
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileCollection.h"
#include "TH1.h"

#include "TMath.h"
#include "TComplex.h"
#include "TH2.h"

#include "../include/initAnalysisBinning_pt_v2.h" ///IMPORTANT: here I define analysis binnning
///#include "initAnalysisBinning_pt.h" //for wider bins - forward

///Histograms for BDT selection: use optimized cuts for BDT, done using mass scan in (centrality, y, pT)
TFile *file_bdtcut_cent0to30 = TFile::Open("../HistogramsForBDTselection/BDTCuts_finePtBins_PrompD0_PbPb_cent0to30_trkPt1GeV.root");
TH1D *hist_bdtcut_cent0to30 = (TH1D*)file_bdtcut_cent0to30->Get("hist_bdtcut");
TFile *file_bdtcut_cent30to50 = TFile::Open("../HistogramsForBDTselection/BDTCuts_finePtBins_PrompD0_PbPb_cent30to50_trkPt1GeV.root");
TH1D *hist_bdtcut_cent30to50 = (TH1D*)file_bdtcut_cent30to50->Get("hist_bdtcut");
TFile *file_bdtcut_cent50to80 = TFile::Open("../HistogramsForBDTselection/BDTCuts_finePtBins_PrompD0_PbPb_cent50to80_trkPt1GeV.root");
TH1D *hist_bdtcut_cent50to80 = (TH1D*)file_bdtcut_cent50to80->Get("hist_bdtcut");

///Function to get optimized BDT cuts
Double_t GetMVACut(Double_t y, Double_t pt, Int_t centrality){

Double_t mvacut = -1.0;
if(fabs(y)>=2.4) return 999.9;

if(!hist_bdtcut_cent0to30)return mvacut;
if(!hist_bdtcut_cent30to50)return mvacut;
if(!hist_bdtcut_cent50to80)return mvacut;

if(centrality<60){
  mvacut = hist_bdtcut_cent0to30->GetBinContent(hist_bdtcut_cent0to30->GetXaxis()->FindBin(y),hist_bdtcut_cent0to30->GetYaxis()->FindBin(pt));
  if(pt>=8.0) mvacut = hist_bdtcut_cent0to30->GetBinContent(hist_bdtcut_cent0to30->GetXaxis()->FindBin(y),hist_bdtcut_cent0to30->GetYaxis()->FindBin(7.99));
  if(pt<=1.5) mvacut = hist_bdtcut_cent0to30->GetBinContent(hist_bdtcut_cent0to30->GetXaxis()->FindBin(y),hist_bdtcut_cent0to30->GetYaxis()->FindBin(1.51));
}
else if(60<=centrality && centrality<100){
  mvacut = hist_bdtcut_cent30to50->GetBinContent(hist_bdtcut_cent30to50->GetXaxis()->FindBin(y),hist_bdtcut_cent30to50->GetYaxis()->FindBin(pt));
  if(pt>=8.0) mvacut = hist_bdtcut_cent30to50->GetBinContent(hist_bdtcut_cent30to50->GetXaxis()->FindBin(y),hist_bdtcut_cent30to50->GetYaxis()->FindBin(7.99));
  if(pt<=1.5) mvacut = hist_bdtcut_cent30to50->GetBinContent(hist_bdtcut_cent30to50->GetXaxis()->FindBin(y),hist_bdtcut_cent30to50->GetYaxis()->FindBin(1.51));
}else{
  mvacut = hist_bdtcut_cent50to80->GetBinContent(hist_bdtcut_cent50to80->GetXaxis()->FindBin(y),hist_bdtcut_cent50to80->GetYaxis()->FindBin(pt));
  if(pt>=8.0) mvacut = hist_bdtcut_cent50to80->GetBinContent(hist_bdtcut_cent50to80->GetXaxis()->FindBin(y),hist_bdtcut_cent50to80->GetYaxis()->FindBin(7.99));
  if(pt<=1.5) mvacut = hist_bdtcut_cent50to80->GetBinContent(hist_bdtcut_cent50to80->GetXaxis()->FindBin(y),hist_bdtcut_cent50to80->GetYaxis()->FindBin(1.51));
}

return mvacut;

}


///Function to calculate vn and mass of D0, D0bar, D0+D0bar, D0-D0bar
void build_vnHistograms_version8_pt(std::string datalist = "", std::string outputfile = ""){

cout << datalist << endl;
cout << outputfile << endl;

TH1::StatOverflows(kTRUE);

TChain *tree = new TChain("d0ana/VertexCompositeNtuple");
TFileCollection *fcData = new TFileCollection(datalist.c_str(), "", datalist.c_str());
tree->AddFileInfoList(fcData->GetList());
std::cout << "d0 ready" << endl;


Int_t MAXD0CandSize=100;
tree->SetBranchStatus("*",0);     // disable all branches
tree->SetBranchStatus("candSize",1);
tree->SetBranchStatus("centrality",1);
tree->SetBranchStatus("bestvtxZ",1);
tree->SetBranchStatus("ephfpAngle",1);
tree->SetBranchStatus("ephfmAngle",1);
tree->SetBranchStatus("ephfpQ",1);
tree->SetBranchStatus("ephfmQ",1);
tree->SetBranchStatus("pT",1);
tree->SetBranchStatus("phi",1);
tree->SetBranchStatus("mva",1);
tree->SetBranchStatus("y",1);
tree->SetBranchStatus("flavor",1);
tree->SetBranchStatus("mass",1);
tree->SetBranchStatus("dca",1);
tree->SetBranchStatus("eptkAngle",1);
tree->SetBranchStatus("eptkQ",1);

// SetBranchAddress
Int_t candSize;
Int_t centrality;
Float_t bestvtxZ;
Float_t ephfpAngle[3];
Float_t ephfmAngle[3];
Float_t ephfpQ[3];
Float_t ephfmQ[3];
Float_t pT[MAXD0CandSize];
Float_t phi[MAXD0CandSize];
Float_t mva[MAXD0CandSize];
Float_t y[MAXD0CandSize];
Float_t flavor[MAXD0CandSize];
Float_t mass[MAXD0CandSize];
Float_t dca[MAXD0CandSize];
Float_t eptkAngle[2];
Float_t eptkQ[2];
tree->SetBranchAddress("candSize", &candSize);
tree->SetBranchAddress("centrality", &centrality);
tree->SetBranchAddress("bestvtxZ", &bestvtxZ);
tree->SetBranchAddress("ephfpAngle", &ephfpAngle);
tree->SetBranchAddress("ephfmAngle", &ephfmAngle);
tree->SetBranchAddress("ephfpQ", &ephfpQ);
tree->SetBranchAddress("ephfmQ", &ephfmQ);
tree->SetBranchAddress("pT", &pT);
tree->SetBranchAddress("phi", &phi);
tree->SetBranchAddress("mva", &mva);
tree->SetBranchAddress("y", &y);
tree->SetBranchAddress("flavor", &flavor);
tree->SetBranchAddress("mass", &mass);
tree->SetBranchAddress("dca", &dca);
tree->SetBranchAddress("eptkAngle", &eptkAngle);
tree->SetBranchAddress("eptkQ", &eptkQ);


// array for vn denominator for each centrality classification

TH1D *hist_Q2Q2_HFmHFp_Re_centbin[N_CENTBIN];//centrality
TH1D *hist_Q2Q2_HFmHFp_Im_centbin[N_CENTBIN];//centrality
TH1D *aux_hist_Q2Q2_HFmHFp = new TH1D("aux_hist_Q2Q2_HFmHFp","",60000,-30000,30000);

TH1D *hist_Q2Q2_HFmTrk_Re_centbin[N_CENTBIN];//centrality
TH1D *hist_Q2Q2_HFmTrk_Im_centbin[N_CENTBIN];//centrality
TH1D *aux_hist_Q2Q2_HFmTrk = new TH1D("aux_hist_Q2Q2_HFmTrk","",60000,-30000,30000);

TH1D *hist_Q2Q2_HFpTrk_Re_centbin[N_CENTBIN];//centrality
TH1D *hist_Q2Q2_HFpTrk_Im_centbin[N_CENTBIN];//centrality
TH1D *aux_hist_Q2Q2_HFpTrk = new TH1D("aux_hist_Q2Q2_HFpTrk","",60000,-30000,30000);

TH1D *hist_Q3Q3_HFmHFp_Re_centbin[N_CENTBIN];//centrality
TH1D *hist_Q3Q3_HFmHFp_Im_centbin[N_CENTBIN];//centrality
TH1D *aux_hist_Q3Q3_HFmHFp = new TH1D("aux_hist_Q3Q3_HFmHFp","",30000,-15000,15000);

TH1D *hist_Q3Q3_HFmTrk_Re_centbin[N_CENTBIN];//centrality
TH1D *hist_Q3Q3_HFmTrk_Im_centbin[N_CENTBIN];//centrality
TH1D *aux_hist_Q3Q3_HFmTrk = new TH1D("aux_hist_Q3Q3_HFmTrk","",30000,-15000,15000);

TH1D *hist_Q3Q3_HFpTrk_Re_centbin[N_CENTBIN];//centrality
TH1D *hist_Q3Q3_HFpTrk_Im_centbin[N_CENTBIN];//centrality
TH1D *aux_hist_Q3Q3_HFpTrk = new TH1D("aux_hist_Q3Q3_HFpTrk","",30000,-15000,15000);

TH1D *hist_Q1Q1_HFmHFp_Re_centbin[N_CENTBIN];//centrality
TH1D *hist_Q1Q1_HFmHFp_Im_centbin[N_CENTBIN];//centrality
TH1D *aux_hist_Q1Q1_HFmHFp = new TH1D("aux_hist_Q1Q1_HFmHFp","",60000,-30000,30000);

TH1D *hist_Q1Q1Q2Q2_HFmHFpHFmHFp_Re_centbin[N_CENTBIN];//centrality
TH1D *hist_Q1Q1Q2Q2_HFmHFpHFmHFp_Im_centbin[N_CENTBIN];//centrality
TH1D *aux_hist_Q1Q1Q2Q2_HFmHFpHFmHFp = new TH1D("aux_hist_Q1Q1Q2Q2_HFmHFpHFmHFp","",20000,-1000000,1000000);

TH1D *hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Re_centbin[N_CENTBIN];//centrality
TH1D *hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Im_centbin[N_CENTBIN];//centrality
TH1D *aux_hist_Q2Q2Q2Q2_HFmHFmHFpHFp = new TH1D("aux_hist_Q2Q2Q2Q2_HFmHFmHFpHFp","",20000,-1000000,1000000);

for(Int_t i=0; i<N_CENTBIN; i++){
  hist_Q2Q2_HFmHFp_Re_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFmHFp->Clone("Q2Q2_HFmHFp_Re_"+label_centbin[i]);
  hist_Q2Q2_HFmHFp_Im_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFmHFp->Clone("Q2Q2_HFmHFp_Im_"+label_centbin[i]);
  hist_Q2Q2_HFmTrk_Re_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFmTrk->Clone("Q2Q2_HFmTrk_Re_"+label_centbin[i]);
  hist_Q2Q2_HFmTrk_Im_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFmTrk->Clone("Q2Q2_HFmTrk_Im_"+label_centbin[i]);
  hist_Q2Q2_HFpTrk_Re_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFpTrk->Clone("Q2Q2_HFpTrk_Re_"+label_centbin[i]);
  hist_Q2Q2_HFpTrk_Im_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFpTrk->Clone("Q2Q2_HFpTrk_Im_"+label_centbin[i]);

  hist_Q3Q3_HFmHFp_Re_centbin[i] = (TH1D *)aux_hist_Q3Q3_HFmHFp->Clone("Q3Q3_HFmHFp_Re_"+label_centbin[i]);
  hist_Q3Q3_HFmHFp_Im_centbin[i] = (TH1D *)aux_hist_Q3Q3_HFmHFp->Clone("Q3Q3_HFmHFp_Im_"+label_centbin[i]);
  hist_Q3Q3_HFmTrk_Re_centbin[i] = (TH1D *)aux_hist_Q3Q3_HFmTrk->Clone("Q3Q3_HFmTrk_Re_"+label_centbin[i]);
  hist_Q3Q3_HFmTrk_Im_centbin[i] = (TH1D *)aux_hist_Q3Q3_HFmTrk->Clone("Q3Q3_HFmTrk_Im_"+label_centbin[i]);
  hist_Q3Q3_HFpTrk_Re_centbin[i] = (TH1D *)aux_hist_Q3Q3_HFpTrk->Clone("Q3Q3_HFpTrk_Re_"+label_centbin[i]);
  hist_Q3Q3_HFpTrk_Im_centbin[i] = (TH1D *)aux_hist_Q3Q3_HFpTrk->Clone("Q3Q3_HFpTrk_Im_"+label_centbin[i]);

  hist_Q1Q1_HFmHFp_Re_centbin[i] = (TH1D *)aux_hist_Q1Q1_HFmHFp->Clone("Q1Q1_HFmHFp_Re_"+label_centbin[i]);
  hist_Q1Q1_HFmHFp_Im_centbin[i] = (TH1D *)aux_hist_Q1Q1_HFmHFp->Clone("Q1Q1_HFmHFp_Im_"+label_centbin[i]);

  hist_Q1Q1Q2Q2_HFmHFpHFmHFp_Re_centbin[i] = (TH1D *)aux_hist_Q1Q1Q2Q2_HFmHFpHFmHFp->Clone("Q1Q1Q2Q2_HFmHFpHFmHFp_Re_"+label_centbin[i]);
  hist_Q1Q1Q2Q2_HFmHFpHFmHFp_Im_centbin[i] = (TH1D *)aux_hist_Q1Q1Q2Q2_HFmHFpHFmHFp->Clone("Q1Q1Q2Q2_HFmHFpHFmHFp_Im_"+label_centbin[i]);  

  hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Re_centbin[i] = (TH1D *)aux_hist_Q2Q2Q2Q2_HFmHFmHFpHFp->Clone("Q2Q2Q2Q2_HFmHFmHFpHFp_Re_"+label_centbin[i]);
  hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Im_centbin[i] = (TH1D *)aux_hist_Q2Q2Q2Q2_HFmHFmHFpHFp->Clone("Q2Q2Q2Q2_HFmHFmHFpHFp_Im_"+label_centbin[i]);

}


///Prepare arrays for vn numerator

TH1D *h_D0_mass[N_CENTBIN][N_PTBIN][N_YBIN];//pT,y,centrality
TH1D *h_D0bar_mass[N_CENTBIN][N_PTBIN][N_YBIN];//pT,y,centrality
TH1D *h_D0plusD0bar_mass[N_CENTBIN][N_PTBIN][N_YBIN];//pT,y,centrality
TH2D *h_D0plusD0bar_mass_dca[N_CENTBIN][N_PTBIN][N_YBIN];//pT,y,centrality
TH1D *hist_D0_mass_aux = new TH1D("D0_mass_aux","",100,1.74,2.0);
TH1D *hist_D0bar_mass_aux = new TH1D("D0bar_mass_aux","",100,1.74,2.0);
TH1D *hist_D0plusD0bar_mass_aux = new TH1D("D0plusD0bar_mass_aux","",100,1.74,2.0);
TH2D *hist_D0plusD0bar_mass_dca_aux = new TH2D("D0plusD0bar_mass_dca_aux","",100,1.74,2.0,100,0.0,0.1);

TH2D *hist_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin = new TH2D("aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin","",N_MASSBIN,massbinning,400,-200.0,200.0);

TH2D *hist_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q3Q3_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q3Q3_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q3Q3_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q3Q3_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *aux_Q3Q3_D0D0barHF_centbin_pTbin_ybin_massbin = new TH2D("aux_Q3Q3_D0D0barHF_centbin_pTbin_ybin_massbin","",N_MASSBIN,massbinning,400,-200.0,200.0);

TH2D *hist_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *aux_Q1Q1_D0D0barHF_centbin_pTbin_ybin_massbin = new TH2D("aux_Q1Q1_D0D0barHF_centbin_pTbin_ybin_massbin","",N_MASSBIN,massbinning,400,-200.0,200.0);

TH2D *hist_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1Q2_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1Q2_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1Q2_D0minusD0barHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q1Q1Q2_D0minusD0barHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *aux_Q1Q1Q2_D0D0barHFHF_centbin_pTbin_ybin_massbin = new TH2D("aux_Q1Q1Q2_D0D0barHFHF_centbin_pTbin_ybin_massbin","",N_MASSBIN,massbinning,400,-15000.0,15000.0);

TH2D *hist_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];//centrality,pT,y
TH2D *aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin = new TH2D("aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin","",N_MASSBIN,massbinning,400,-100000.0,100000.0);




for(Int_t i_centbin=0; i_centbin<N_CENTBIN; i_centbin++){
  for(Int_t i_pTbin=0; i_pTbin<N_PTBIN; i_pTbin++){
    for(Int_t i_ybin=0; i_ybin<N_YBIN; i_ybin++){

        h_D0_mass[i_centbin][i_pTbin][i_ybin] = (TH1D*)hist_D0_mass_aux->Clone("hist_D0_mass_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        h_D0bar_mass[i_centbin][i_pTbin][i_ybin] = (TH1D*)hist_D0bar_mass_aux->Clone("hist_D0bar_mass_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        h_D0plusD0bar_mass[i_centbin][i_pTbin][i_ybin] = (TH1D*)hist_D0plusD0bar_mass_aux->Clone("hist_D0plusD0bar_mass_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);      
        h_D0plusD0bar_mass_dca[i_centbin][i_pTbin][i_ybin] = (TH2D*)hist_D0plusD0bar_mass_dca_aux->Clone("hist_D0plusD0bar_mass_dca_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0HF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0HF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
 
        hist_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);        
        hist_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0plusD0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0plusD0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
         
        hist_Q2Q2_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0minusD0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0minusD0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

 
        hist_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q3Q3_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q3Q3_D0HF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q3Q3_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q3Q3_D0HF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
 
        hist_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q3Q3_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q3Q3_D0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);        
        hist_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q3Q3_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q3Q3_D0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q3Q3_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q3Q3_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q3Q3_D0plusD0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q3Q3_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q3Q3_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q3Q3_D0plusD0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        
        hist_Q3Q3_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q3Q3_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q3Q3_D0minusD0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q3Q3_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q3Q3_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q3Q3_D0minusD0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
       
       
        hist_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1_D0HF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1_D0HF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
 
        hist_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1_D0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);        
        hist_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1_D0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q1Q1_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1_D0plusD0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q1Q1_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1_D0plusD0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
         
        hist_Q1Q1_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1_D0minusD0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q1Q1_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1_D0minusD0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

      
        hist_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1Q2_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1Q2_D0HFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1Q2_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1Q2_D0HFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
 
        hist_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1Q2_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1Q2_D0barHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);        
        hist_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1Q2_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1Q2_D0barHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q1Q1Q2_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1Q2_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1Q2_D0plusD0barHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q1Q1Q2_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1Q2_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1Q2_D0plusD0barHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
         
        hist_Q1Q1Q2_D0minusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1Q2_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1Q2_D0minusD0barHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q1Q1Q2_D0minusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q1Q1Q2_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q1Q1Q2_D0minusD0barHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

       
        hist_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0HFHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0HFHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
 
        hist_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0barHFHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);        
        hist_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0barHFHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
         
        hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
 
 
    }  
  }
}


//Long64_t n_entries = tree->GetEntries();
Long64_t n_entries = 1000;
//Long64_t n_entries = tree->GetEntriesFast();
for(Long64_t ii=0; ii<n_entries; ii++){//loop in events

  ///Long64_t ientry = tree->LoadTree(ii);
  ///if (ientry < 0) break;

  tree->GetEntry(ii);
  if(ii % 500000 == 0)
    printf("current entry = %lld out of %lld : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);

  ///PV Z-position condition
  if(fabs(bestvtxZ) >= 15)continue; 

  ///For vn denominator - EP resolution
  TComplex aux_Q2_HFm(ephfmQ[1]*TMath::Cos(2.0*ephfmAngle[1]),ephfmQ[1]*TMath::Sin(2.0*ephfmAngle[1]),0);
  TComplex aux_Q2_HFp(ephfpQ[1]*TMath::Cos(2.0*ephfpAngle[1]),ephfpQ[1]*TMath::Sin(2.0*ephfpAngle[1]),0); 
  TComplex aux_Q2_Trk(eptkQ[0]*TMath::Cos(2.0*eptkAngle[0]),eptkQ[0]*TMath::Sin(2.0*eptkAngle[0]),0);

  TComplex aux_Q3_HFm(ephfmQ[2]*TMath::Cos(3.0*ephfmAngle[2]),ephfmQ[2]*TMath::Sin(3.0*ephfmAngle[2]),0);
  TComplex aux_Q3_HFp(ephfpQ[2]*TMath::Cos(3.0*ephfpAngle[2]),ephfpQ[2]*TMath::Sin(3.0*ephfpAngle[2]),0);
  TComplex aux_Q3_Trk(eptkQ[1]*TMath::Cos(3.0*eptkAngle[1]),eptkQ[1]*TMath::Sin(3.0*eptkAngle[1]),0); 

  TComplex aux_Q1_HFm(ephfmQ[0]*TMath::Cos(1.0*ephfmAngle[0]),ephfmQ[0]*TMath::Sin(1.0*ephfmAngle[0]),0);
  TComplex aux_Q1_HFp(ephfpQ[0]*TMath::Cos(1.0*ephfpAngle[0]),ephfpQ[0]*TMath::Sin(1.0*ephfpAngle[0]),0);

  for(Int_t i=0; i<N_CENTBIN; i++){
    if(min_centbin[i]<=centrality && centrality<max_centbin[i]){

       Double_t aux_Q2Q2_HFmHFp_Re = (aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)).Re();
       Double_t aux_Q2Q2_HFmHFp_Im = (aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)).Im();
       Double_t aux_Q2Q2_HFmTrk_Re = (aux_Q2_HFm*TComplex::Conjugate(aux_Q2_Trk)).Re();
       Double_t aux_Q2Q2_HFmTrk_Im = (aux_Q2_HFm*TComplex::Conjugate(aux_Q2_Trk)).Im();
       Double_t aux_Q2Q2_HFpTrk_Re = (aux_Q2_HFp*TComplex::Conjugate(aux_Q2_Trk)).Re();
       Double_t aux_Q2Q2_HFpTrk_Im = (aux_Q2_HFp*TComplex::Conjugate(aux_Q2_Trk)).Im();
       hist_Q2Q2_HFmHFp_Re_centbin[i]->Fill(aux_Q2Q2_HFmHFp_Re);
       hist_Q2Q2_HFmHFp_Im_centbin[i]->Fill(aux_Q2Q2_HFmHFp_Im);
       hist_Q2Q2_HFmTrk_Re_centbin[i]->Fill(aux_Q2Q2_HFmTrk_Re);
       hist_Q2Q2_HFmTrk_Im_centbin[i]->Fill(aux_Q2Q2_HFmTrk_Im);
       hist_Q2Q2_HFpTrk_Re_centbin[i]->Fill(aux_Q2Q2_HFpTrk_Re);
       hist_Q2Q2_HFpTrk_Im_centbin[i]->Fill(aux_Q2Q2_HFpTrk_Im);

       Double_t aux_Q3Q3_HFmHFp_Re = (aux_Q3_HFm*TComplex::Conjugate(aux_Q3_HFp)).Re();
       Double_t aux_Q3Q3_HFmHFp_Im = (aux_Q3_HFm*TComplex::Conjugate(aux_Q3_HFp)).Im();
       Double_t aux_Q3Q3_HFmTrk_Re = (aux_Q3_HFm*TComplex::Conjugate(aux_Q3_Trk)).Re();
       Double_t aux_Q3Q3_HFmTrk_Im = (aux_Q3_HFm*TComplex::Conjugate(aux_Q3_Trk)).Im();
       Double_t aux_Q3Q3_HFpTrk_Re = (aux_Q3_HFp*TComplex::Conjugate(aux_Q3_Trk)).Re();
       Double_t aux_Q3Q3_HFpTrk_Im = (aux_Q3_HFp*TComplex::Conjugate(aux_Q3_Trk)).Im();
       hist_Q3Q3_HFmHFp_Re_centbin[i]->Fill(aux_Q3Q3_HFmHFp_Re);
       hist_Q3Q3_HFmHFp_Im_centbin[i]->Fill(aux_Q3Q3_HFmHFp_Im);
       hist_Q3Q3_HFmTrk_Re_centbin[i]->Fill(aux_Q3Q3_HFmTrk_Re);
       hist_Q3Q3_HFmTrk_Im_centbin[i]->Fill(aux_Q3Q3_HFmTrk_Im);
       hist_Q3Q3_HFpTrk_Re_centbin[i]->Fill(aux_Q3Q3_HFpTrk_Re);
       hist_Q3Q3_HFpTrk_Im_centbin[i]->Fill(aux_Q3Q3_HFpTrk_Im);

       Double_t aux_Q1Q1_HFmHFp_Re = (aux_Q1_HFm*TComplex::Conjugate(aux_Q1_HFp)).Re();
       Double_t aux_Q1Q1_HFmHFp_Im = (aux_Q1_HFm*TComplex::Conjugate(aux_Q1_HFp)).Im();
       hist_Q1Q1_HFmHFp_Re_centbin[i]->Fill(aux_Q1Q1_HFmHFp_Re);
       hist_Q1Q1_HFmHFp_Im_centbin[i]->Fill(aux_Q1Q1_HFmHFp_Im);

       Double_t aux_Q1Q1Q2Q2_HFmHFpHFmHFp_Re = (aux_Q1_HFm*TComplex::Conjugate(aux_Q1_HFp)*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)).Re();
       Double_t aux_Q1Q1Q2Q2_HFmHFpHFmHFp_Im = (aux_Q1_HFm*TComplex::Conjugate(aux_Q1_HFp)*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)).Im();
       hist_Q1Q1Q2Q2_HFmHFpHFmHFp_Re_centbin[i]->Fill(aux_Q1Q1Q2Q2_HFmHFpHFmHFp_Re);
       hist_Q1Q1Q2Q2_HFmHFpHFmHFp_Im_centbin[i]->Fill(aux_Q1Q1Q2Q2_HFmHFpHFmHFp_Im);

       Double_t aux_Q2Q2Q2Q2_HFmHFmHFpHFp_Re = (aux_Q2_HFm*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Re();
       Double_t aux_Q2Q2Q2Q2_HFmHFmHFpHFp_Im = (aux_Q2_HFm*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Im();
       hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Re_centbin[i]->Fill(aux_Q2Q2Q2Q2_HFmHFmHFpHFp_Re);
       hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Im_centbin[i]->Fill(aux_Q2Q2Q2Q2_HFmHFmHFpHFp_Im);

    }     
  }  

  ///For vn numerator
  for(Int_t i_centbin=0; i_centbin<N_CENTBIN; i_centbin++){
 
    for(Int_t i_pTbin=0; i_pTbin<N_PTBIN; i_pTbin++){

      for(Int_t i_ybin=0; i_ybin<N_YBIN; i_ybin++){   
    
          for(Int_t i=0; i<candSize; i++){//loop in D0 + D0bar Candidates

            //apply optimized BDT cuts
            if(mva[i]<GetMVACut(y[i],pT[i],centrality))continue; //default
            ///if(mva[i]<(GetMVACut(y[i],pT[i],centrality)+0.05))continue; //sys: tight DBT by 0.05
            ///if(mva[i]<(GetMVACut(y[i],pT[i],centrality)-0.05))continue; //sys: loose DBT by 0.05
          
            if( (min_centbin[i_centbin]<=centrality && centrality<max_centbin[i_centbin]) && 
                (min_pTbin[i_pTbin]<=pT[i] && pT[i]<max_pTbin[i_pTbin]) &&
                (min_ybin[i_ybin]<=y[i] && y[i]<max_ybin[i_ybin])
              ){
              if(flavor[i]>0){
                h_D0_mass[i_centbin][i_pTbin][i_ybin]->Fill(mass[i]);//Fill mass histogram for D0  
                h_D0plusD0bar_mass[i_centbin][i_pTbin][i_ybin]->Fill(mass[i]);//Fill mass histogram for D0plusD0bar 
                h_D0plusD0bar_mass_dca[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],dca[i]);//Fill mass vs dca histogram for D0plusD0bar 
                if(y[i]>0){

                  TComplex aux_Q2_D0(TMath::Cos(2.0*phi[i]),TMath::Sin(2.0*phi[i]),0);
                  Double_t aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_HFm)).Re();
                  Double_t aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_HFm)).Im();
                  hist_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin);   
                  hist_Q2Q2_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin);

                  TComplex aux_Q3_D0(TMath::Cos(3.0*phi[i]),TMath::Sin(3.0*phi[i]),0);
                  Double_t aux_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin = (aux_Q3_D0*TComplex::Conjugate(aux_Q3_HFm)).Re();
                  Double_t aux_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin = (aux_Q3_D0*TComplex::Conjugate(aux_Q3_HFm)).Im();
                  hist_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin);

                  TComplex aux_Q1_D0(TMath::Cos(1.0*phi[i]),TMath::Sin(1.0*phi[i]),0);
                  //Double_t aux_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFm)).Re();
                  //Double_t aux_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFm)).Im();
                  Double_t aux_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFm) + aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFp)).Re();//for v1_odd
                  Double_t aux_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFm) + aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFp)).Im();//for v1_odd
                  hist_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin);   
                  hist_Q1Q1_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin);

                  Double_t aux_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0*aux_Q1_HFp*TComplex::Conjugate(aux_Q2_HFm)).Re();
                  Double_t aux_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0*aux_Q1_HFp*TComplex::Conjugate(aux_Q2_HFm)).Im();
                  hist_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin);   
                  hist_Q1Q1Q2_D0minusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0minusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin); 

                  Double_t aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFp*TComplex::Conjugate(aux_Q2_HFm)*TComplex::Conjugate(aux_Q2_HFm)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFp*TComplex::Conjugate(aux_Q2_HFm)*TComplex::Conjugate(aux_Q2_HFm)).Im();
                  hist_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin);   
                  hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin); 

                }
                else{

                  TComplex aux_Q2_D0(TMath::Cos(2.0*phi[i]),TMath::Sin(2.0*phi[i]),0);
                  Double_t aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_HFp)).Re();
                  Double_t aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_HFp)).Im();
                  hist_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin);
                 
                  TComplex aux_Q3_D0(TMath::Cos(3.0*phi[i]),TMath::Sin(3.0*phi[i]),0);
                  Double_t aux_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin = (aux_Q3_D0*TComplex::Conjugate(aux_Q3_HFp)).Re();
                  Double_t aux_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin = (aux_Q3_D0*TComplex::Conjugate(aux_Q3_HFp)).Im();
                  hist_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin);

                  TComplex aux_Q1_D0(TMath::Cos(1.0*phi[i]),TMath::Sin(1.0*phi[i]),0);
                  //Double_t aux_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFp)).Re();
                  //Double_t aux_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFp)).Im();
                  Double_t aux_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFp) + aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFm)).Re();//for v1_odd
                  Double_t aux_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFp) + aux_Q1_D0*TComplex::Conjugate(aux_Q1_HFm)).Im();//for v1_odd
                  hist_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin);
 
                  Double_t aux_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0*aux_Q1_HFm*TComplex::Conjugate(aux_Q2_HFp)).Re();
                  Double_t aux_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0*aux_Q1_HFm*TComplex::Conjugate(aux_Q2_HFp)).Im();
                  hist_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin);   
                  hist_Q1Q1Q2_D0minusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0minusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin); 

                  Double_t aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Im();
                  hist_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin);   
                  hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin); 

                }
              }//end if D0
              else{
                h_D0bar_mass[i_centbin][i_pTbin][i_ybin]->Fill(mass[i]);//Fill mass histogram for D0bar 
                h_D0plusD0bar_mass[i_centbin][i_pTbin][i_ybin]->Fill(mass[i]);//Fill mass histogram for D0plusD0bar 
                h_D0plusD0bar_mass_dca[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],dca[i]);//Fill mass vs dca histogram for D0plusD0bar 
                if(y[i]>0){

                  TComplex aux_Q2_D0bar(TMath::Cos(2.0*phi[i]),TMath::Sin(2.0*phi[i]),0); 
                  Double_t aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_HFm)).Re();
                  Double_t aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_HFm)).Im();
                  hist_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin);
                 
                  TComplex aux_Q3_D0bar(TMath::Cos(3.0*phi[i]),TMath::Sin(3.0*phi[i]),0);
                  Double_t aux_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin = (aux_Q3_D0bar*TComplex::Conjugate(aux_Q3_HFm)).Re();
                  Double_t aux_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin = (aux_Q3_D0bar*TComplex::Conjugate(aux_Q3_HFm)).Im();
                  hist_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin); 
                  hist_Q3Q3_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin);

                  TComplex aux_Q1_D0bar(TMath::Cos(1.0*phi[i]),TMath::Sin(1.0*phi[i]),0); 
                  //Double_t aux_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFm)).Re();
                  //Double_t aux_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFm)).Im();
                  Double_t aux_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFm) + aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFp)).Re();//for v1_odd
                  Double_t aux_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFm) + aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFp)).Im();//for v1_odd
                  hist_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin);

                  Double_t aux_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*aux_Q1_HFp*TComplex::Conjugate(aux_Q2_HFm)).Re();
                  Double_t aux_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*aux_Q1_HFp*TComplex::Conjugate(aux_Q2_HFm)).Im();
                  hist_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin);   
                  hist_Q1Q1Q2_D0minusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0minusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin); 

                  Double_t aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFp*TComplex::Conjugate(aux_Q2_HFm)*TComplex::Conjugate(aux_Q2_HFm)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFp*TComplex::Conjugate(aux_Q2_HFm)*TComplex::Conjugate(aux_Q2_HFm)).Im();
                  hist_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin);   
                  hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin); 

                }
                else{

                  TComplex aux_Q2_D0bar(TMath::Cos(2.0*phi[i]),TMath::Sin(2.0*phi[i]),0);
                  Double_t aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_HFp)).Re();
                  Double_t aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_HFp)).Im();
                  hist_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin);       
                  hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin);

                  TComplex aux_Q3_D0bar(TMath::Cos(3.0*phi[i]),TMath::Sin(3.0*phi[i]),0);
                  Double_t aux_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin = (aux_Q3_D0bar*TComplex::Conjugate(aux_Q3_HFp)).Re();
                  Double_t aux_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin = (aux_Q3_D0bar*TComplex::Conjugate(aux_Q3_HFp)).Im();
                  hist_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q3Q3_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin);

                  TComplex aux_Q1_D0bar(TMath::Cos(1.0*phi[i]),TMath::Sin(1.0*phi[i]),0);
                  //Double_t aux_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFp)).Re();
                  //Double_t aux_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFp)).Im();
                  Double_t aux_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFp) + aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFm)).Re();//for v1_odd
                  Double_t aux_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFp) + aux_Q1_D0bar*TComplex::Conjugate(aux_Q1_HFm)).Im();//for v1_odd
                  hist_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin);       
                  hist_Q1Q1_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin);

                  Double_t aux_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*aux_Q1_HFm*TComplex::Conjugate(aux_Q2_HFp)).Re();
                  Double_t aux_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q1_D0bar*aux_Q1_HFm*TComplex::Conjugate(aux_Q2_HFp)).Im();
                  hist_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin);   
                  hist_Q1Q1Q2_D0minusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q1Q1Q2_D0minusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin); 

                  Double_t aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Im();
                  hist_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin);   
                  hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin);
                  hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],-1.0*aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin); 

                }

              }//end if D0bar 

            }//if condition for all bins

          }//D0 + D0bar Candidates

      }//y bins
    
    }//pT bins
  
  }//end loop in centrality bins 
  
}//end loop in events


///Save histograms
//TString outstr = outputfile; 
TFile* fout = new TFile(outputfile.c_str(),"RECREATE");
fout->cd();


for(Int_t i_centbin=0; i_centbin<N_CENTBIN; i_centbin++){

  ///v2 - denominator
  hist_Q2Q2_HFmHFp_Re_centbin[i_centbin]->Write(); 
  hist_Q2Q2_HFmHFp_Im_centbin[i_centbin]->Write();
  hist_Q2Q2_HFmTrk_Re_centbin[i_centbin]->Write();
  hist_Q2Q2_HFmTrk_Im_centbin[i_centbin]->Write();
  hist_Q2Q2_HFpTrk_Re_centbin[i_centbin]->Write();
  hist_Q2Q2_HFpTrk_Im_centbin[i_centbin]->Write();

  ///v3 - denominator
  hist_Q3Q3_HFmHFp_Re_centbin[i_centbin]->Write();
  hist_Q3Q3_HFmHFp_Im_centbin[i_centbin]->Write();
  hist_Q3Q3_HFmTrk_Re_centbin[i_centbin]->Write();
  hist_Q3Q3_HFmTrk_Im_centbin[i_centbin]->Write();
  hist_Q3Q3_HFpTrk_Re_centbin[i_centbin]->Write();
  hist_Q3Q3_HFpTrk_Im_centbin[i_centbin]->Write();

  ///v1 - denominator
  hist_Q1Q1_HFmHFp_Re_centbin[i_centbin]->Write();
  hist_Q1Q1_HFmHFp_Im_centbin[i_centbin]->Write();

  ///v1 (mixed harmonics) - denominator
  hist_Q1Q1Q2Q2_HFmHFpHFmHFp_Re_centbin[i_centbin]->Write();
  hist_Q1Q1Q2Q2_HFmHFpHFmHFp_Im_centbin[i_centbin]->Write();

  ///v2_4 - denominator (other terms are computed in other parts above)
  hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Re_centbin[i_centbin]->Write();
  hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Im_centbin[i_centbin]->Write(); 

  for(Int_t i_pTbin=0; i_pTbin<N_PTBIN; i_pTbin++){
    for(Int_t i_ybin=0; i_ybin<N_YBIN; i_ybin++){
      h_D0_mass[i_centbin][i_pTbin][i_ybin]->Write(); 
      h_D0bar_mass[i_centbin][i_pTbin][i_ybin]->Write(); 
      h_D0plusD0bar_mass[i_centbin][i_pTbin][i_ybin]->Write();
      h_D0plusD0bar_mass_dca[i_centbin][i_pTbin][i_ybin]->Write();

      ///v2 - numerator
      hist_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();

      ///v3 - numerator
      hist_Q3Q3_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q3Q3_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q3Q3_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q3Q3_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q3Q3_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q3Q3_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q3Q3_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q3Q3_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();

      ///v1 - numerator 
      hist_Q1Q1_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1_D0minusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1_D0minusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();

      ///v1 (mixed harmonics) - numerator
      hist_Q1Q1Q2_D0HFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1Q2_D0HFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1Q2_D0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1Q2_D0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1Q2_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1Q2_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1Q2_D0minusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q1Q1Q2_D0minusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();

      ///v2_4 - numerator (other terms are computed in other parts above)
      hist_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0minusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();

    }
  }
}



fout->Close();


std::cout<<" \n Done ... "<<std::endl;

}
