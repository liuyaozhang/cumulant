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

#include "../include/initAnalysisBinning_pt_v2_build15_ExactCumulant_HFsumET.h" ///IMPORTANT: here I define analysis binnning
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
//this function to test the exact cumulant formula. 
void build_vnHistograms_version15_pt_ExactCumulant(std::string datalist = "", std::string outputfile = ""){

cout << datalist << endl;
cout << outputfile << endl;

TH1::StatOverflows(kTRUE);

TChain *tree = new TChain("d0ana/VertexCompositeNtuple");
TFileCollection *fcData = new TFileCollection(datalist.c_str(), "", datalist.c_str());
tree->AddFileInfoList(fcData->GetList());
std::cout << "d0 ready" << endl;


Int_t MAXD0CandSize=100;
tree->SetBranchStatus("*",0);     // disable all branches
tree->SetBranchStatus("Ntrkoffline",1);
tree->SetBranchStatus("HFsumETPlus",1); 
tree->SetBranchStatus("HFsumETMinus",1); 
tree->SetBranchStatus("candSize",1);
tree->SetBranchStatus("centrality",1);
tree->SetBranchStatus("bestvtxZ",1);
tree->SetBranchStatus("ephfpAngle",1);
tree->SetBranchStatus("ephfmAngle",1);
tree->SetBranchStatus("ephfpQ",1);
tree->SetBranchStatus("ephfmQ",1);
tree->SetBranchStatus("ephfpM",1);
tree->SetBranchStatus("ephfmM",1);

tree->SetBranchStatus("eptkAngle_2p4",1);
tree->SetBranchStatus("eptkAngle_2p4m",1);
tree->SetBranchStatus("eptkAngle_2p4p",1);
tree->SetBranchStatus("eptkQ_2p4",1);
tree->SetBranchStatus("eptkQ_2p4m",1);
tree->SetBranchStatus("eptkQ_2p4p",1);
tree->SetBranchStatus("eptkM_2p4",1);
tree->SetBranchStatus("eptkM_2p4m",1);
tree->SetBranchStatus("eptkM_2p4p",1);

tree->SetBranchStatus("pT",1);
tree->SetBranchStatus("phi",1);
tree->SetBranchStatus("mva",1);
tree->SetBranchStatus("y",1);
tree->SetBranchStatus("flavor",1);
tree->SetBranchStatus("mass",1);
tree->SetBranchStatus("dca",1);

// SetBranchAddress
Int_t Ntrkoffline;
Float_t HFsumETPlus;
Float_t HFsumETMinus;
Int_t candSize;
Int_t centrality;
Float_t bestvtxZ;
Float_t ephfpAngle[4];
Float_t ephfmAngle[4];
Float_t ephfpQ[4];
Float_t ephfmQ[4];
Float_t ephfpM;
Float_t ephfmM;

Float_t eptkAngle_2p4[3];
Float_t eptkAngle_2p4m[3];
Float_t eptkAngle_2p4p[3];
Float_t eptkQ_2p4[3];
Float_t eptkQ_2p4m[3];
Float_t eptkQ_2p4p[3];
Float_t eptkM_2p4;
Float_t eptkM_2p4m;
Float_t eptkM_2p4p;

Float_t pT[MAXD0CandSize];
Float_t phi[MAXD0CandSize];
Float_t mva[MAXD0CandSize];
Float_t y[MAXD0CandSize];
Float_t flavor[MAXD0CandSize];
Float_t mass[MAXD0CandSize];
Float_t dca[MAXD0CandSize];

tree->SetBranchAddress("Ntrkoffline", &Ntrkoffline);
tree->SetBranchAddress("HFsumETPlus", &HFsumETPlus);
tree->SetBranchAddress("HFsumETMinus", &HFsumETMinus);
tree->SetBranchAddress("candSize", &candSize);
tree->SetBranchAddress("centrality", &centrality);
tree->SetBranchAddress("bestvtxZ", &bestvtxZ);
tree->SetBranchAddress("ephfpAngle", &ephfpAngle);
tree->SetBranchAddress("ephfmAngle", &ephfmAngle);
tree->SetBranchAddress("ephfpQ", &ephfpQ);
tree->SetBranchAddress("ephfmQ", &ephfmQ);
tree->SetBranchAddress("ephfpM", &ephfpM);
tree->SetBranchAddress("ephfmM", &ephfmM);

tree->SetBranchAddress("eptkAngle_2p4", &eptkAngle_2p4);
tree->SetBranchAddress("eptkAngle_2p4m", &eptkAngle_2p4m);
tree->SetBranchAddress("eptkAngle_2p4p", &eptkAngle_2p4p);
tree->SetBranchAddress("eptkQ_2p4", &eptkQ_2p4);
tree->SetBranchAddress("eptkQ_2p4m", &eptkQ_2p4m);
tree->SetBranchAddress("eptkQ_2p4p", &eptkQ_2p4p);
tree->SetBranchAddress("eptkM_2p4", &eptkM_2p4);
tree->SetBranchAddress("eptkM_2p4m", &eptkM_2p4m);
tree->SetBranchAddress("eptkM_2p4p", &eptkM_2p4p);

tree->SetBranchAddress("pT", &pT);
tree->SetBranchAddress("phi", &phi);
tree->SetBranchAddress("mva", &mva);
tree->SetBranchAddress("y", &y);
tree->SetBranchAddress("flavor", &flavor);
tree->SetBranchAddress("mass", &mass);
tree->SetBranchAddress("dca", &dca);


// array for vn denominator for each centrality classification
//v2{4} trakcer: 
TH1D* hist_Q2Q2_TrkmTrkp_Re_centbin[N_CENTBIN];
TH1D* hist_Q2Q2_TrkmTrkp_Im_centbin[N_CENTBIN];
TH1D* aux_hist_Q2Q2_TrkmTrkp = new TH1D("aux_hist_Q2Q2_TrkmTrkp","",500,-0.06,0.06);

TH1D* hist_Q4Q4_TrkmTrkp_Re_centbin[N_CENTBIN];
TH1D* hist_Q4Q4_TrkmTrkp_Im_centbin[N_CENTBIN];
TH1D* aux_hist_Q4Q4_TrkmTrkp = new TH1D("aux_hist_Q4Q4_TrkmTrkp","",500,-0.06,0.06);

TH1D* hist_Q4Q2Q2_TrkmTrkpTrkp_Re_centbin[N_CENTBIN];
TH1D* hist_Q4Q2Q2_TrkmTrkpTrkp_Im_centbin[N_CENTBIN];
TH1D* aux_hist_Q4Q2Q2_TrkmTrkpTrkp = new TH1D("aux_hist_Q4Q2Q2_TrkmTrkpTrkp","",500,-0.06,0.06);

TH1D* hist_Q2Q2Q4_TrkmTrkmTrkp_Re_centbin[N_CENTBIN];
TH1D* hist_Q2Q2Q4_TrkmTrkmTrkp_Im_centbin[N_CENTBIN];
TH1D* aux_hist_Q2Q2Q4_TrkmTrkmTrkp = new TH1D("aux_hist_Q2Q2Q4_TrkmTrkmTrkp","",500,-0.06,0.06);

TH1D* hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Re_centbin[N_CENTBIN];
TH1D* hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Im_centbin[N_CENTBIN];
TH1D* aux_hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp = new TH1D("aux_hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp","",500,-0.06,0.06);

//v2{2} HF&trakcer
TH1D* hist_Q2Q2_HFmHFp_Re_centbin[N_CENTBIN];//centrality
TH1D* hist_Q2Q2_HFmHFp_Im_centbin[N_CENTBIN];//centrality
TH1D* aux_hist_Q2Q2_HFmHFp = new TH1D("aux_hist_Q2Q2_HFmHFp","",200,-0.06,0.06);//60000,-30000,30000
//v2{4} 
TH1D* hist_Q2Q2_HFmTrk_Re_centbin[N_CENTBIN];//centrality
TH1D* hist_Q2Q2_HFmTrk_Im_centbin[N_CENTBIN];//centrality
TH1D* aux_hist_Q2Q2_HFmTrk = new TH1D("aux_hist_Q2Q2_HFmTrk","",200,-0.06,0.06);

TH1D* hist_Q2Q2_HFpTrk_Re_centbin[N_CENTBIN];//centrality
TH1D* hist_Q2Q2_HFpTrk_Im_centbin[N_CENTBIN];//centrality
TH1D* aux_hist_Q2Q2_HFpTrk = new TH1D("aux_hist_Q2Q2_HFpTrk","",200,-0.06,0.06);

TH1D* hist_Q4Q4_HFmHFp_Re_centbin[N_CENTBIN]; 
TH1D* hist_Q4Q4_HFmHFp_Im_centbin[N_CENTBIN];
TH1D* aux_hist_Q4Q4_HFmHFp = new TH1D("aux_hist_Q4Q4_HFmHFp","",200,-0.06,0.06);

TH1D* hist_Q4Q2Q2_HFmHFpHFp_Re_centbin[N_CENTBIN];
TH1D* hist_Q4Q2Q2_HFmHFpHFp_Im_centbin[N_CENTBIN];
TH1D* aux_hist_Q4Q2Q2_HFmHFpHFp = new TH1D("aux_hist_Q4Q2Q2_HFmHFpHFp","",200,-0.06,0.06);

TH1D* hist_Q2Q2Q4_HFmHFmHFp_Re_centbin[N_CENTBIN];
TH1D* hist_Q2Q2Q4_HFmHFmHFp_Im_centbin[N_CENTBIN];
TH1D* aux_hist_Q2Q2Q4_HFmHFmHFp = new TH1D("aux_hist_Q2Q2Q4_HFmHFmHFp","",200,-0.06,0.06);

TH1D* hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Re_centbin[N_CENTBIN];//centrality
TH1D* hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Im_centbin[N_CENTBIN];//centrality
TH1D* aux_hist_Q2Q2Q2Q2_HFmHFmHFpHFp = new TH1D("aux_hist_Q2Q2Q2Q2_HFmHFmHFpHFp","",200,-0.06,0.06);

for(Int_t i=0; i<N_CENTBIN; i++){
  hist_Q2Q2_TrkmTrkp_Re_centbin[i] = (TH1D*)aux_hist_Q2Q2_TrkmTrkp->Clone("Q2Q2_TrkmTrkp_Re_"+label_centbin[i]);
  hist_Q2Q2_TrkmTrkp_Im_centbin[i] = (TH1D*)aux_hist_Q2Q2_TrkmTrkp->Clone("Q2Q2_TrkmTrkp_Im_"+label_centbin[i]);
  
  hist_Q4Q4_TrkmTrkp_Re_centbin[i] = (TH1D*)aux_hist_Q4Q4_TrkmTrkp->Clone("Q4Q4_TrkmTrkp_Re_"+label_centbin[i]);
  hist_Q4Q4_TrkmTrkp_Im_centbin[i] = (TH1D*)aux_hist_Q4Q4_TrkmTrkp->Clone("Q4Q4_TrkmTrkp_Im_"+label_centbin[i]);
  
  hist_Q4Q2Q2_TrkmTrkpTrkp_Re_centbin[i] = (TH1D*)aux_hist_Q4Q2Q2_TrkmTrkpTrkp->Clone("Q4Q2Q2_TrkmTrkpTrkp_Re_"+label_centbin[i]);
  hist_Q4Q2Q2_TrkmTrkpTrkp_Im_centbin[i] = (TH1D*)aux_hist_Q4Q2Q2_TrkmTrkpTrkp->Clone("Q4Q2Q2_TrkmTrkpTrkp_Im_"+label_centbin[i]);
  
  hist_Q2Q2Q4_TrkmTrkmTrkp_Re_centbin[i] = (TH1D*)aux_hist_Q2Q2Q4_TrkmTrkmTrkp->Clone("Q2Q2Q4_TrkmTrkmTrkp_Re_"+label_centbin[i]);
  hist_Q2Q2Q4_TrkmTrkmTrkp_Im_centbin[i] = (TH1D*)aux_hist_Q2Q2Q4_TrkmTrkmTrkp->Clone("Q2Q2Q4_TrkmTrkmTrkp_Im_"+label_centbin[i]);
  
  hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Re_centbin[i] = (TH1D*)aux_hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp->Clone("Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Re_"+label_centbin[i]);
  hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Im_centbin[i] = (TH1D*)aux_hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp->Clone("Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Im_"+label_centbin[i]);
   
  //HF
  hist_Q2Q2_HFmHFp_Re_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFmHFp->Clone("Q2Q2_HFmHFp_Re_"+label_centbin[i]);
  hist_Q2Q2_HFmHFp_Im_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFmHFp->Clone("Q2Q2_HFmHFp_Im_"+label_centbin[i]);
  hist_Q2Q2_HFmTrk_Re_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFmTrk->Clone("Q2Q2_HFmTrk_Re_"+label_centbin[i]);
  hist_Q2Q2_HFmTrk_Im_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFmTrk->Clone("Q2Q2_HFmTrk_Im_"+label_centbin[i]);
  hist_Q2Q2_HFpTrk_Re_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFpTrk->Clone("Q2Q2_HFpTrk_Re_"+label_centbin[i]);
  hist_Q2Q2_HFpTrk_Im_centbin[i] = (TH1D *)aux_hist_Q2Q2_HFpTrk->Clone("Q2Q2_HFpTrk_Im_"+label_centbin[i]);
  
  hist_Q4Q4_HFmHFp_Re_centbin[i] = (TH1D *)aux_hist_Q4Q4_HFmHFp->Clone("Q4Q4_HFmHFp_Re_"+label_centbin[i]);  
  hist_Q4Q4_HFmHFp_Im_centbin[i] = (TH1D *)aux_hist_Q4Q4_HFmHFp->Clone("Q4Q4_HFmHFp_Im_"+label_centbin[i]);
  
  hist_Q4Q2Q2_HFmHFpHFp_Re_centbin[i] = (TH1D *)aux_hist_Q4Q2Q2_HFmHFpHFp->Clone("Q4Q2Q2_HFmHFpHFp_Re_"+label_centbin[i]);
  hist_Q4Q2Q2_HFmHFpHFp_Im_centbin[i] = (TH1D *)aux_hist_Q4Q2Q2_HFmHFpHFp->Clone("Q4Q2Q2_HFmHFpHFp_Im_"+label_centbin[i]);
  
  hist_Q2Q2Q4_HFmHFmHFp_Re_centbin[i] = (TH1D *)aux_hist_Q2Q2Q4_HFmHFmHFp->Clone("Q2Q2Q4_HFmHFmHFp_Re_"+label_centbin[i]);
  hist_Q2Q2Q4_HFmHFmHFp_Im_centbin[i] = (TH1D *)aux_hist_Q2Q2Q4_HFmHFmHFp->Clone("Q2Q2Q4_HFmHFmHFp_Im_"+label_centbin[i]);
 
  hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Re_centbin[i] = (TH1D *)aux_hist_Q2Q2Q2Q2_HFmHFmHFpHFp->Clone("Q2Q2Q2Q2_HFmHFmHFpHFp_Re_"+label_centbin[i]);
  hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Im_centbin[i] = (TH1D *)aux_hist_Q2Q2Q2Q2_HFmHFmHFpHFp->Clone("Q2Q2Q2Q2_HFmHFmHFpHFp_Im_"+label_centbin[i]);

}


///Prepare arrays for vn numerator
TH1D* h_D0_mass[N_CENTBIN][N_PTBIN][N_YBIN];//pT,y,centrality
TH1D* h_D0bar_mass[N_CENTBIN][N_PTBIN][N_YBIN];//pT,y,centrality
TH1D* h_D0plusD0bar_mass[N_CENTBIN][N_PTBIN][N_YBIN];//pT,y,centrality
TH2D* h_D0plusD0bar_mass_dca[N_CENTBIN][N_PTBIN][N_YBIN];//pT,y,centrality
TH1D* hist_D0_mass_aux = new TH1D("D0_mass_aux","",100,1.74,2.0);
TH1D* hist_D0bar_mass_aux = new TH1D("D0bar_mass_aux","",100,1.74,2.0);
TH1D* hist_D0plusD0bar_mass_aux = new TH1D("D0plusD0bar_mass_aux","",100,1.74,2.0);
TH2D* hist_D0plusD0bar_mass_dca_aux = new TH2D("D0plusD0bar_mass_dca_aux","",100,1.74,2.0,100,0.0,0.1);

//from Tracker
TH2D* hist_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2_D0plusD0barTrk_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2_D0plusD0barTrk_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* aux_Q2Q2_D0D0barTrk_centbin_pTbin_ybin_massbin = new TH2D("aux_Q2Q2_D0D0barTrk_centbin_pTbin_ybin_massbin","",N_MASSBIN,massbinning,20000,-0.4,0.4);

TH2D* hist_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q4_D0plusD0barTrkTrk_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q4_D0plusD0barTrkTrk_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* aux_Q2Q2Q4_D0D0barTrkTrk_centbin_pTbin_ybin_massbin = new TH2D("aux_Q2Q2Q4_D0D0barTrkTrk_centbin_pTbin_ybin_massbin","",N_MASSBIN,massbinning,20000,-0.2,0.2);

TH2D* hist_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* aux_Q2Q2Q2Q2_D0D0barTrkTrkTrk_centbin_pTbin_ybin_massbin = new TH2D("aux_Q2Q2Q2Q2_D0D0barTrkTrkTrk_centbin_pTbin_ybin_massbin","",N_MASSBIN,massbinning,20000,-0.2,0.2);

//from HF
TH2D* hist_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin = new TH2D("aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin","",N_MASSBIN,massbinning,20000,-0.2,0.2);//400,-200,200

TH2D* hist_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q4_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q4_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* aux_Q2Q2Q4_D0D0barHFHF_centbin_pTbin_ybin_massbin = new TH2D("aux_Q2Q2Q4_D0D0barHFHF_centbin_pTbin_ybin_massbin","",N_MASSBIN, massbinning, 20000, -0.05, 0.05);

TH2D* hist_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[N_CENTBIN][N_PTBIN][N_YBIN];
TH2D* aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin = new TH2D("aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin","",N_MASSBIN,massbinning,20000,-0.05,0.05);

for(Int_t i_centbin=0; i_centbin<N_CENTBIN; i_centbin++){
  for(Int_t i_pTbin=0; i_pTbin<N_PTBIN; i_pTbin++){
    for(Int_t i_ybin=0; i_ybin<N_YBIN; i_ybin++){
        h_D0_mass[i_centbin][i_pTbin][i_ybin] = (TH1D*)hist_D0_mass_aux->Clone("hist_D0_mass_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        h_D0bar_mass[i_centbin][i_pTbin][i_ybin] = (TH1D*)hist_D0bar_mass_aux->Clone("hist_D0bar_mass_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        h_D0plusD0bar_mass[i_centbin][i_pTbin][i_ybin] = (TH1D*)hist_D0plusD0bar_mass_aux->Clone("hist_D0plusD0bar_mass_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);      
        h_D0plusD0bar_mass_dca[i_centbin][i_pTbin][i_ybin] = (TH2D*)hist_D0plusD0bar_mass_dca_aux->Clone("hist_D0plusD0bar_mass_dca_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        //tracker
        hist_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0Trk_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0Trk_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0barTrk_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0barTrk_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2_D0plusD0barTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0plusD0barTrk_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2_D0plusD0barTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0plusD0barTrk_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        
        hist_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0TrkTrk_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0TrkTrk_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0barTrkTrk_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0barTrkTrk_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q4_D0plusD0barTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0plusD0barTrkTrk_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q4_D0plusD0barTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0plusD0barTrkTrk_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
       
        hist_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barTrkTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0TrkTrkTrk_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barTrkTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0TrkTrkTrk_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barTrkTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barTrkTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barTrkTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barTrkTrkTrk_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        //HF
        hist_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0HF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0HF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
 
        hist_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);        
        hist_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0plusD0barHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2_D0D0barHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2_D0plusD0barHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
         
        //Q2Q2Q2_D0HFHF
        hist_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0HFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0HFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0barHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0barHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q2Q2Q4_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0plusD0barHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q4_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q4_D0D0barHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q4_D0plusD0barHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0HFHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0HFHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
 
        hist_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0barHFHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);        
        hist_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0barHFHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);

        hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
        hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin] = (TH2D*)aux_Q2Q2Q2Q2_D0D0barHFHFHF_centbin_pTbin_ybin_massbin->Clone("hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_"+label_centbin[i_centbin]+"_"+label_pTbin[i_pTbin]+"_"+label_ybin[i_ybin]);
         
    }  
  }
}


//Long64_t n_entries = tree->GetEntries();
Long64_t n_entries = 100000;
//Long64_t n_entries = tree->GetEntriesFast();
for(Long64_t ii=0; ii<n_entries; ii++){//loop in events

  ///Long64_t ientry = tree->LoadTree(ii);
  ///if (ientry < 0) break;

  tree->GetEntry(ii);
  if(ii % 50000 == 0)
    printf("current entry = %lld out of %lld : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);
  ///PV Z-position condition
  if(fabs(bestvtxZ) >= 15)continue; 

  ///For vn denominator - EP resolution
  TComplex aux_Q2_HFm(ephfmQ[1]*TMath::Cos(2.0*ephfmAngle[1]),ephfmQ[1]*TMath::Sin(2.0*ephfmAngle[1]),0);
  TComplex aux_Q2_HFp(ephfpQ[1]*TMath::Cos(2.0*ephfpAngle[1]),ephfpQ[1]*TMath::Sin(2.0*ephfpAngle[1]),0); 
  TComplex aux_Q4_HFm(ephfmQ[3]*TMath::Cos(4.0*ephfmAngle[3]),ephfmQ[3]*TMath::Sin(4.0*ephfmAngle[3]),0);
  TComplex aux_Q4_HFp(ephfpQ[3]*TMath::Cos(4.0*ephfpAngle[3]),ephfpQ[3]*TMath::Sin(4.0*ephfpAngle[3]),0); 

  TComplex aux_Q2_Trk(eptkQ_2p4[0]*TMath::Cos(2.0*eptkAngle_2p4[0]),eptkQ_2p4[0]*TMath::Sin(2.0*eptkAngle_2p4[0]),0);
  TComplex aux_Q2_Trkm(eptkQ_2p4m[0]*TMath::Cos(2.0*eptkAngle_2p4m[0]),eptkQ_2p4m[0]*TMath::Sin(2.0*eptkAngle_2p4m[0]),0);
  TComplex aux_Q2_Trkp(eptkQ_2p4p[0]*TMath::Cos(2.0*eptkAngle_2p4p[0]),eptkQ_2p4p[0]*TMath::Sin(2.0*eptkAngle_2p4p[0]),0);
  
  TComplex aux_Q4_Trk(eptkQ_2p4[2]*TMath::Cos(4.0*eptkAngle_2p4[2]),eptkQ_2p4[2]*TMath::Sin(4.0*eptkAngle_2p4[2]),0);
  TComplex aux_Q4_Trkm(eptkQ_2p4m[2]*TMath::Cos(4.0*eptkAngle_2p4m[2]),eptkQ_2p4m[2]*TMath::Sin(4.0*eptkAngle_2p4m[2]),0);
  TComplex aux_Q4_Trkp(eptkQ_2p4p[2]*TMath::Cos(4.0*eptkAngle_2p4p[2]),eptkQ_2p4p[2]*TMath::Sin(4.0*eptkAngle_2p4p[2]),0);
  
  //cout <<"========================================="<<endl;
  //cout <<"ephfmQ[1]"<< ":" <<ephfmQ[1] << " "<< "ephfmQ[3]" << ":" << ephfmQ[3] << " " << "ephfpQ[1]" << ":" << ephfpQ[1] << " " << "ephfpQ[3]" << ":" << ephfpQ[3]<< endl; 
  //cout <<"ephfmAngle[1]"<< "=" << ephfmAngle[1] << " " << "ephfmAngle[3]" <<  " " << ephfmAngle[3] << " " << "ephfpAngle[1]" << " " << ephfpAngle[1] << " " <<"ephfpAngle[3]" << " " << ephfpAngle[3]<< endl;  
  //cout << "eptkQ_2p4m[0]" << ":" << eptkQ_2p4m[0] << " " << "eptkQ_2p4m[2]" << ":" << eptkQ_2p4m[2] << " " <<"eptkQ_2p4p[0]" << ":" << eptkQ_2p4p[0] << " " << "eptkQ_2p4p[2]"  << ":" << eptkQ_2p4p[2] <<endl;  
  ////cout << "eptkAngle_2p4m[0]" << "=" << eptkAngle_2p4m[0] << "eptkAngle_2p4m[2]"<<" " << eptkAngle_2p4m[2] <<"eptkAngle_2p4p[0]" << "=" << eptkAngle_2p4p[0] << "eptkAngle_2p4p[2]"<<" " << eptkAngle_2p4p[2] << endl;  
  //cout <<"HFsumETPlus"<< "="<< HFsumETPlus <<" " << "HFsumETMinus"<<"="<< HFsumETMinus<<endl; 
  //cout <<"ephfpM" << ":" << ephfpM << " " << "ephfmM" << ":" << ephfmM << endl; 
  //cout <<"eptkM_2p4m"<< ":" <<eptkM_2p4m << " " << "eptkM_2p4p" << ":" << eptkM_2p4p << endl; 
  //cout <<"----------------------------------------"<<endl; 
  for(Int_t i=0; i<N_CENTBIN; i++){
    if(min_centbin[i]<=centrality && centrality<max_centbin[i]){

       Double_t aux_Q2Q2_TrkmTrkp_Re = (aux_Q2_Trkm*TComplex::Conjugate(aux_Q2_Trkp)).Re(); 
       Double_t aux_Q2Q2_TrkmTrkp_Im = (aux_Q2_Trkm*TComplex::Conjugate(aux_Q2_Trkp)).Im(); 
        
       if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q2Q2_TrkmTrkp_Re_centbin[i]->Fill(aux_Q2Q2_TrkmTrkp_Re * 1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
       if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q2Q2_TrkmTrkp_Im_centbin[i]->Fill(aux_Q2Q2_TrkmTrkp_Im * 1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
       
       Double_t aux_Q4Q4_TrkmTrkp_Re = (aux_Q4_Trkm*TComplex::Conjugate(aux_Q4_Trkp)).Re(); 
       Double_t aux_Q4Q4_TrkmTrkp_Im = (aux_Q4_Trkm*TComplex::Conjugate(aux_Q4_Trkp)).Im(); 
       if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q4Q4_TrkmTrkp_Re_centbin[i]->Fill(aux_Q4Q4_TrkmTrkp_Re * 1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
       if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q4Q4_TrkmTrkp_Im_centbin[i]->Fill(aux_Q4Q4_TrkmTrkp_Im * 1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
       
       Double_t aux_Q4Q2Q2_TrkmTrkpTrkp_Re = (aux_Q4_Trkm*TComplex::Conjugate(aux_Q2_Trkp)*TComplex::Conjugate(aux_Q2_Trkp)).Re(); 
       Double_t aux_Q4Q2Q2_TrkmTrkpTrkp_Im = (aux_Q4_Trkm*TComplex::Conjugate(aux_Q2_Trkp)*TComplex::Conjugate(aux_Q2_Trkp)).Im(); 
       if(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1) != 0) hist_Q4Q2Q2_TrkmTrkpTrkp_Re_centbin[i]->Fill(aux_Q4Q2Q2_TrkmTrkpTrkp_Re * 1./(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1));
       if(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1) != 0) hist_Q4Q2Q2_TrkmTrkpTrkp_Im_centbin[i]->Fill(aux_Q4Q2Q2_TrkmTrkpTrkp_Im * 1./(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1));
       
       Double_t aux_Q2Q2Q4_TrkmTrkmTrkp_Re = (aux_Q2_Trkm*aux_Q2_Trkm*TComplex::Conjugate(aux_Q4_Trkp)).Re(); 
       Double_t aux_Q2Q2Q4_TrkmTrkmTrkp_Im = (aux_Q2_Trkm*aux_Q2_Trkm*TComplex::Conjugate(aux_Q4_Trkp)).Im(); 
       if(eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p != 0) hist_Q2Q2Q4_TrkmTrkmTrkp_Re_centbin[i]->Fill(aux_Q2Q2Q4_TrkmTrkmTrkp_Re * 1./(eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p), eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p);
       if(eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p != 0) hist_Q2Q2Q4_TrkmTrkmTrkp_Im_centbin[i]->Fill(aux_Q2Q2Q4_TrkmTrkmTrkp_Im * 1./(eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p), eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p);
       
       Double_t aux_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Re = (aux_Q2_Trkm*aux_Q2_Trkm*TComplex::Conjugate(aux_Q2_Trkp)*TComplex::Conjugate(aux_Q2_Trkp)).Re(); 
       Double_t aux_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Im = (aux_Q2_Trkm*aux_Q2_Trkm*TComplex::Conjugate(aux_Q2_Trkp)*TComplex::Conjugate(aux_Q2_Trkp)).Im();    
       if(eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p*(eptkM_2p4p-1) != 0) hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Re_centbin[i]->Fill(aux_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Re * 1./(eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p*(eptkM_2p4p-1));
       if(eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p*(eptkM_2p4p-1) != 0) hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Im_centbin[i]->Fill(aux_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Im * 1./(eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*(eptkM_2p4m-1)*eptkM_2p4p*(eptkM_2p4p-1));
      
       //HF
       Double_t aux_Q2Q2_HFmHFp_Re = (aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)).Re();
       Double_t aux_Q2Q2_HFmHFp_Im = (aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)).Im();
       Double_t aux_Q2Q2_HFmTrk_Re = (aux_Q2_HFm*TComplex::Conjugate(aux_Q2_Trk)).Re();
       Double_t aux_Q2Q2_HFmTrk_Im = (aux_Q2_HFm*TComplex::Conjugate(aux_Q2_Trk)).Im();
       Double_t aux_Q2Q2_HFpTrk_Re = (aux_Q2_HFp*TComplex::Conjugate(aux_Q2_Trk)).Re();
       Double_t aux_Q2Q2_HFpTrk_Im = (aux_Q2_HFp*TComplex::Conjugate(aux_Q2_Trk)).Im();
       if(HFsumETMinus*HFsumETPlus != 0) hist_Q2Q2_HFmHFp_Re_centbin[i]->Fill(aux_Q2Q2_HFmHFp_Re * 1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
       if(HFsumETMinus*HFsumETPlus != 0) hist_Q2Q2_HFmHFp_Im_centbin[i]->Fill(aux_Q2Q2_HFmHFp_Im * 1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
       if(HFsumETMinus*Ntrkoffline != 0) hist_Q2Q2_HFmTrk_Re_centbin[i]->Fill(aux_Q2Q2_HFmTrk_Re * 1./(HFsumETMinus*Ntrkoffline), HFsumETMinus*Ntrkoffline);
       if(HFsumETMinus*Ntrkoffline != 0) hist_Q2Q2_HFmTrk_Im_centbin[i]->Fill(aux_Q2Q2_HFmTrk_Im * 1./(HFsumETMinus*Ntrkoffline), HFsumETMinus*Ntrkoffline);
       if(HFsumETPlus*Ntrkoffline  != 0) hist_Q2Q2_HFpTrk_Re_centbin[i]->Fill(aux_Q2Q2_HFpTrk_Re * 1./(HFsumETPlus*Ntrkoffline), HFsumETPlus*Ntrkoffline);
       if(HFsumETPlus*Ntrkoffline  != 0) hist_Q2Q2_HFpTrk_Im_centbin[i]->Fill(aux_Q2Q2_HFpTrk_Im * 1./(HFsumETPlus*Ntrkoffline), HFsumETPlus*Ntrkoffline);
        
       Double_t aux_Q4Q4_HFmHFp_Re = (aux_Q4_HFm*TComplex::Conjugate(aux_Q4_HFp)).Re();
       Double_t aux_Q4Q4_HFmHFp_Im = (aux_Q4_HFm*TComplex::Conjugate(aux_Q4_HFp)).Im();
       if(HFsumETMinus*HFsumETPlus != 0) hist_Q4Q4_HFmHFp_Re_centbin[i]->Fill(aux_Q4Q4_HFmHFp_Re * 1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
       if(HFsumETMinus*HFsumETPlus != 0) hist_Q4Q4_HFmHFp_Im_centbin[i]->Fill(aux_Q4Q4_HFmHFp_Im * 1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
       
       Double_t aux_Q4Q2Q2_HFmHFpHFp_Re = (aux_Q4_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Re();
       Double_t aux_Q4Q2Q2_HFmHFpHFp_Im = (aux_Q4_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Im();
       if(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1) != 0) hist_Q4Q2Q2_HFmHFpHFp_Re_centbin[i]->Fill(aux_Q4Q2Q2_HFmHFpHFp_Re * 1./(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)), (HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)));
       if(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1) != 0) hist_Q4Q2Q2_HFmHFpHFp_Im_centbin[i]->Fill(aux_Q4Q2Q2_HFmHFpHFp_Im * 1./(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)), (HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)));
       
       Double_t aux_Q2Q2Q4_HFmHFmHFp_Re = (aux_Q2_HFm*aux_Q2_HFm*TComplex::Conjugate(aux_Q4_HFp)).Re();
       Double_t aux_Q2Q2Q4_HFmHFmHFp_Im = (aux_Q2_HFm*aux_Q2_HFm*TComplex::Conjugate(aux_Q4_HFp)).Im();
       if(HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus != 0) hist_Q2Q2Q4_HFmHFmHFp_Re_centbin[i]->Fill(aux_Q2Q2Q4_HFmHFmHFp_Re * 1./(HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus), HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus);
       if(HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus != 0) hist_Q2Q2Q4_HFmHFmHFp_Im_centbin[i]->Fill(aux_Q2Q2Q4_HFmHFmHFp_Im * 1./(HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus), HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus);
       Double_t aux_Q2Q2Q2Q2_HFmHFmHFpHFp_Re = (aux_Q2_HFm*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Re();
       Double_t aux_Q2Q2Q2Q2_HFmHFmHFpHFp_Im = (aux_Q2_HFm*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Im();
       if(HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus*(HFsumETPlus-1) != 0) hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Re_centbin[i]->Fill(aux_Q2Q2Q2Q2_HFmHFmHFpHFp_Re * 1./(HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus*(HFsumETPlus-1)), (HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus*(HFsumETPlus-1)));
       if(HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus*(HFsumETPlus-1) != 0) hist_Q2Q2Q2Q2_HFmHFmHFpHFp_Im_centbin[i]->Fill(aux_Q2Q2Q2Q2_HFmHFmHFpHFp_Im * 1./(HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus*(HFsumETPlus-1)), (HFsumETMinus*(HFsumETMinus-1)*HFsumETPlus*(HFsumETPlus-1)));

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
                  //Tracker
                  Double_t aux_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_Trkm)).Re();
                  Double_t aux_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_Trkm)).Im();
                
                  if(eptkM_2p4m !=0 ) hist_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin*1./eptkM_2p4m, eptkM_2p4m); 
                  if(eptkM_2p4m !=0 ) hist_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin*1./eptkM_2p4m, eptkM_2p4m); 
                  if(eptkM_2p4m !=0 ) hist_Q2Q2_D0plusD0barTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin*1./eptkM_2p4m, eptkM_2p4m);; 
                  if(eptkM_2p4m !=0 ) hist_Q2Q2_D0plusD0barTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin*1./eptkM_2p4m, eptkM_2p4m);; 
                  
                  Double_t aux_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_Trkp*TComplex::Conjugate(aux_Q4_Trkm)).Re();
                  Double_t aux_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_Trkp*TComplex::Conjugate(aux_Q4_Trkm)).Im();
                     
                  if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p); 
                  if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p); 
                  if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q2Q2Q4_D0plusD0barTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p); 
                  if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q2Q2Q4_D0plusD0barTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p); 
              
                  Double_t aux_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_Trkp*TComplex::Conjugate(aux_Q2_Trkm)*TComplex::Conjugate(aux_Q2_Trkm)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_Trkp*TComplex::Conjugate(aux_Q2_Trkm)*TComplex::Conjugate(aux_Q2_Trkm)).Im();
                  
                  if(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1) != 0) hist_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)), eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)); 
                  if(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1) != 0) hist_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)), eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)); 
                  if(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1) != 0) hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)), eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)); 
                  if(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1) != 0) hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)), eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)); 

                  //HF
                  Double_t aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_HFm)).Re();
                  Double_t aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_HFm)).Im();
                  if(HFsumETMinus != 0) hist_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin*1./HFsumETMinus, HFsumETMinus);
                  if(HFsumETMinus != 0) hist_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin*1./HFsumETMinus, HFsumETMinus);
                  if(HFsumETMinus != 0) hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin*1./HFsumETMinus, HFsumETMinus);
                  if(HFsumETMinus != 0) hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin*1./HFsumETMinus, HFsumETMinus);   
                  Double_t aux_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFp*TComplex::Conjugate(aux_Q4_HFm)).Re();
                  Double_t aux_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFp*TComplex::Conjugate(aux_Q4_HFm)).Im();
                  if(HFsumETMinus*HFsumETPlus != 0) hist_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus != 0) hist_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus != 0) hist_Q2Q2Q4_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus != 0) hist_Q2Q2Q4_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  
                  Double_t aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFp*TComplex::Conjugate(aux_Q2_HFm)*TComplex::Conjugate(aux_Q2_HFm)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFp*TComplex::Conjugate(aux_Q2_HFm)*TComplex::Conjugate(aux_Q2_HFm)).Im();
                  if(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1) != 0) hist_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1)), HFsumETPlus*HFsumETMinus*(HFsumETMinus-1));
                  if(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1) != 0) hist_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1)), HFsumETPlus*HFsumETMinus*(HFsumETMinus-1));
                  if(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1) != 0) hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1)), HFsumETPlus*HFsumETMinus*(HFsumETMinus-1));
                  if(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1) != 0) hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1)), HFsumETPlus*HFsumETMinus*(HFsumETMinus-1));   

                }
                else{
                  TComplex aux_Q2_D0(TMath::Cos(2.0*phi[i]),TMath::Sin(2.0*phi[i]),0);
                  //Tracker
                  Double_t aux_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_Trkp)).Re();
                  Double_t aux_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_Trkp)).Im();

                  if(eptkM_2p4p != 0) hist_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin*1./eptkM_2p4p, eptkM_2p4p);
                  if(eptkM_2p4p != 0) hist_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin*1./eptkM_2p4p, eptkM_2p4p);
                  if(eptkM_2p4p != 0) hist_Q2Q2_D0plusD0barTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin*1./eptkM_2p4p, eptkM_2p4p);;
                  if(eptkM_2p4p != 0) hist_Q2Q2_D0plusD0barTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin*1./eptkM_2p4p, eptkM_2p4p);;

                  Double_t aux_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_Trkm*TComplex::Conjugate(aux_Q4_Trkp)).Re();
                  Double_t aux_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_Trkm*TComplex::Conjugate(aux_Q4_Trkp)).Im();

                  if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
                  if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
                  if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q2Q2Q4_D0plusD0barTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
                  if(eptkM_2p4m*eptkM_2p4p != 0) hist_Q2Q2Q4_D0plusD0barTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);

                  Double_t aux_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_Trkm*TComplex::Conjugate(aux_Q2_Trkp)*TComplex::Conjugate(aux_Q2_Trkp)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_Trkm*TComplex::Conjugate(aux_Q2_Trkp)*TComplex::Conjugate(aux_Q2_Trkp)).Im();

                  if(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1) != 0) hist_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1));
                  if(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1) != 0) hist_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1));
                  if(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1) != 0) hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1));
                  if(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1) != 0) hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1));
                 

                  //HF
                  Double_t aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_HFp)).Re();
                  Double_t aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*TComplex::Conjugate(aux_Q2_HFp)).Im();
                  if(HFsumETPlus !=0) hist_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin*1./HFsumETPlus, HFsumETPlus);
                  if(HFsumETPlus !=0) hist_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin*1./HFsumETPlus, HFsumETPlus);
                  if(HFsumETPlus !=0) hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin*1./HFsumETPlus, HFsumETPlus);
                  if(HFsumETPlus !=0) hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin*1./HFsumETPlus, HFsumETPlus);
                   
                  Double_t aux_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFm*TComplex::Conjugate(aux_Q4_HFp)).Re();
                  Double_t aux_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFm*TComplex::Conjugate(aux_Q4_HFp)).Im();
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);

                  Double_t aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Im();
                  if(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1) !=0) hist_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)), HFsumETMinus*HFsumETPlus*(HFsumETPlus-1));
                  if(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1) !=0) hist_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)), HFsumETMinus*HFsumETPlus*(HFsumETPlus-1));
                  if(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1) !=0) hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)), HFsumETMinus*HFsumETPlus*(HFsumETPlus-1));
                  if(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1) !=0) hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)), HFsumETMinus*HFsumETPlus*(HFsumETPlus-1));   

                }
              }//end if D0
              else{
                h_D0bar_mass[i_centbin][i_pTbin][i_ybin]->Fill(mass[i]);//Fill mass histogram for D0bar 
                h_D0plusD0bar_mass[i_centbin][i_pTbin][i_ybin]->Fill(mass[i]);//Fill mass histogram for D0plusD0bar 
                h_D0plusD0bar_mass_dca[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],dca[i]);//Fill mass vs dca histogram for D0plusD0bar 
                if(y[i]>0){
                  TComplex aux_Q2_D0bar(TMath::Cos(2.0*phi[i]),TMath::Sin(2.0*phi[i]),0); 
                  Double_t aux_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_Trkm)).Re();
                  Double_t aux_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_Trkm)).Im();

                  if(eptkM_2p4m !=0) hist_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin*1./eptkM_2p4m, eptkM_2p4m);
                  if(eptkM_2p4m !=0) hist_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin*1./eptkM_2p4m, eptkM_2p4m);
                  if(eptkM_2p4m !=0) hist_Q2Q2_D0plusD0barTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin*1./eptkM_2p4m, eptkM_2p4m);;
                  if(eptkM_2p4m !=0) hist_Q2Q2_D0plusD0barTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin*1./eptkM_2p4m, eptkM_2p4m);;

                  Double_t aux_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_Trkp*TComplex::Conjugate(aux_Q4_Trkm)).Re();
                  Double_t aux_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_Trkp*TComplex::Conjugate(aux_Q4_Trkm)).Im();

                  if(eptkM_2p4m*eptkM_2p4p !=0) hist_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
                  if(eptkM_2p4m*eptkM_2p4p !=0) hist_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
                  if(eptkM_2p4m*eptkM_2p4p !=0) hist_Q2Q2Q4_D0plusD0barTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
                  if(eptkM_2p4m*eptkM_2p4p !=0) hist_Q2Q2Q4_D0plusD0barTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);

                  Double_t aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_Trkp*TComplex::Conjugate(aux_Q2_Trkm)*TComplex::Conjugate(aux_Q2_Trkm)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_Trkp*TComplex::Conjugate(aux_Q2_Trkm)*TComplex::Conjugate(aux_Q2_Trkm)).Im();

                  if(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1) !=0) hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)), eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1));
                  if(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1) !=0) hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)), eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1));
                  if(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1) !=0) hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)), eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1));
                  if(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1) !=0) hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1)), eptkM_2p4p*eptkM_2p4m*(eptkM_2p4m-1));

                  //HF
                  Double_t aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_HFm)).Re();
                  Double_t aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_HFm)).Im();
                  if(HFsumETMinus !=0) hist_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin*1./HFsumETMinus, HFsumETMinus);
                  if(HFsumETMinus !=0) hist_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin*1./HFsumETMinus, HFsumETMinus);
                  if(HFsumETMinus !=0) hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin*1./HFsumETMinus, HFsumETMinus);
                  if(HFsumETMinus !=0) hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin*1./HFsumETMinus, HFsumETMinus);
                 
                  Double_t aux_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFp*TComplex::Conjugate(aux_Q4_HFm)).Re();
                  Double_t aux_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFp*TComplex::Conjugate(aux_Q4_HFm)).Im();
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
  
                  Double_t aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFp*TComplex::Conjugate(aux_Q2_HFm)*TComplex::Conjugate(aux_Q2_HFm)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFp*TComplex::Conjugate(aux_Q2_HFm)*TComplex::Conjugate(aux_Q2_HFm)).Im();
                  if(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1) !=0) hist_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1)), HFsumETPlus*HFsumETMinus*(HFsumETMinus-1));
                  if(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1) !=0) hist_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1)), HFsumETPlus*HFsumETMinus*(HFsumETMinus-1));
                  if(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1) !=0) hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1)), HFsumETPlus*HFsumETMinus*(HFsumETMinus-1));
                  if(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1) !=0) hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETPlus*HFsumETMinus*(HFsumETMinus-1)), HFsumETPlus*HFsumETMinus*(HFsumETMinus-1));   

                }
                else{
                  TComplex aux_Q2_D0bar(TMath::Cos(2.0*phi[i]),TMath::Sin(2.0*phi[i]),0);
                  Double_t aux_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_Trkp)).Re();
                  Double_t aux_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_Trkp)).Im();

                  if(eptkM_2p4p !=0) hist_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin*1./eptkM_2p4p, eptkM_2p4p);
                  if(eptkM_2p4p !=0) hist_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin*1./eptkM_2p4p, eptkM_2p4p);
                  if(eptkM_2p4p !=0) hist_Q2Q2_D0plusD0barTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin*1./eptkM_2p4p, eptkM_2p4p);;
                  if(eptkM_2p4p !=0) hist_Q2Q2_D0plusD0barTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin*1./eptkM_2p4p, eptkM_2p4p);;

                  Double_t aux_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_Trkm*TComplex::Conjugate(aux_Q4_Trkp)).Re();
                  Double_t aux_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_Trkm*TComplex::Conjugate(aux_Q4_Trkp)).Im();

                  if(eptkM_2p4m*eptkM_2p4p !=0) hist_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
                  if(eptkM_2p4m*eptkM_2p4p !=0) hist_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
                  if(eptkM_2p4m*eptkM_2p4p !=0) hist_Q2Q2Q4_D0plusD0barTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);
                  if(eptkM_2p4m*eptkM_2p4p !=0) hist_Q2Q2Q4_D0plusD0barTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p), eptkM_2p4m*eptkM_2p4p);

                  Double_t aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_Trkm*TComplex::Conjugate(aux_Q2_Trkp)*TComplex::Conjugate(aux_Q2_Trkp)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_Trkm*TComplex::Conjugate(aux_Q2_Trkp)*TComplex::Conjugate(aux_Q2_Trkp)).Im();

                  if(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1) !=0) hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1));
                  if(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1) !=0) hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1));
                  if(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1) !=0) hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1));
                  if(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1) !=0) hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i], aux_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin*1./(eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1)), eptkM_2p4m*eptkM_2p4p*(eptkM_2p4p-1));
                  
                  //HF
                  Double_t aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_HFp)).Re();
                  Double_t aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*TComplex::Conjugate(aux_Q2_HFp)).Im();
                  if(HFsumETPlus !=0) hist_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin*1./HFsumETPlus, HFsumETPlus);
                  if(HFsumETPlus !=0) hist_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin*1./HFsumETPlus, HFsumETPlus);
                  if(HFsumETPlus !=0) hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin*1./HFsumETPlus, HFsumETPlus);       
                  if(HFsumETPlus !=0) hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin*1./HFsumETPlus, HFsumETPlus);

                  Double_t aux_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFm*TComplex::Conjugate(aux_Q4_HFp)).Re();
                  Double_t aux_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFm*TComplex::Conjugate(aux_Q4_HFp)).Im();
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  if(HFsumETMinus*HFsumETPlus !=0) hist_Q2Q2Q4_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus), HFsumETMinus*HFsumETPlus);
                  
                  Double_t aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Re();
                  Double_t aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin = (aux_Q2_D0bar*aux_Q2_HFm*TComplex::Conjugate(aux_Q2_HFp)*TComplex::Conjugate(aux_Q2_HFp)).Im();
                  if(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1) !=0) hist_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)), HFsumETMinus*HFsumETPlus*(HFsumETPlus-1));
                  if(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1) !=0) hist_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)), HFsumETMinus*HFsumETPlus*(HFsumETPlus-1));
                  if(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1) !=0) hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)), HFsumETMinus*HFsumETPlus*(HFsumETPlus-1));
                  if(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1) !=0) hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Fill(mass[i],aux_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin*1./(HFsumETMinus*HFsumETPlus*(HFsumETPlus-1)), HFsumETMinus*HFsumETPlus*(HFsumETPlus-1));   

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
  //tracker
  hist_Q2Q2_TrkmTrkp_Re_centbin[i_centbin]->Write();
  hist_Q2Q2_TrkmTrkp_Im_centbin[i_centbin]->Write();

  hist_Q4Q4_TrkmTrkp_Re_centbin[i_centbin]->Write();
  hist_Q4Q4_TrkmTrkp_Im_centbin[i_centbin]->Write();

  hist_Q4Q2Q2_TrkmTrkpTrkp_Re_centbin[i_centbin]->Write();
  hist_Q4Q2Q2_TrkmTrkpTrkp_Im_centbin[i_centbin]->Write();

  hist_Q2Q2Q4_TrkmTrkmTrkp_Re_centbin[i_centbin]->Write();
  hist_Q2Q2Q4_TrkmTrkmTrkp_Im_centbin[i_centbin]->Write();

  hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Re_centbin[i_centbin]->Write();
  hist_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Im_centbin[i_centbin]->Write();
  
  //HF
  hist_Q2Q2_HFmHFp_Re_centbin[i_centbin]->Write(); 
  hist_Q2Q2_HFmHFp_Im_centbin[i_centbin]->Write();
  hist_Q2Q2_HFmTrk_Re_centbin[i_centbin]->Write();
  hist_Q2Q2_HFmTrk_Im_centbin[i_centbin]->Write();
  hist_Q2Q2_HFpTrk_Re_centbin[i_centbin]->Write();
  hist_Q2Q2_HFpTrk_Im_centbin[i_centbin]->Write();

  hist_Q4Q4_HFmHFp_Re_centbin[i_centbin]->Write();
  hist_Q4Q4_HFmHFp_Im_centbin[i_centbin]->Write();
  
  hist_Q4Q2Q2_HFmHFpHFp_Re_centbin[i_centbin]->Write();
  hist_Q4Q2Q2_HFmHFpHFp_Im_centbin[i_centbin]->Write();

  hist_Q2Q2Q4_HFmHFmHFp_Re_centbin[i_centbin]->Write();
  hist_Q2Q2Q4_HFmHFmHFp_Im_centbin[i_centbin]->Write();

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
      //Tracker
      hist_Q2Q2_D0Trk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0Trk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0barTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0barTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0plusD0barTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0plusD0barTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      
      hist_Q2Q2Q4_D0TrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q4_D0TrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q4_D0barTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q4_D0barTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q4_D0plusD0barTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q4_D0plusD0barTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      
      hist_Q2Q2Q2Q2_D0TrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0TrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();

      //HF
      hist_Q2Q2_D0HF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0HF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0plusD0barHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2_D0plusD0barHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      
      //
      hist_Q2Q2Q4_D0HFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q4_D0HFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q4_D0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q4_D0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q4_D0plusD0barHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();       
      hist_Q2Q2Q4_D0plusD0barHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();

      ///v2_4 - numerator (other terms are computed in other parts above)
      hist_Q2Q2Q2Q2_D0HFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0HFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();
      hist_Q2Q2Q2Q2_D0plusD0barHFHFHF_Im_centbin_pTbin_ybin_massbin[i_centbin][i_pTbin][i_ybin]->Write();

    }
  }
}



fout->Close();


std::cout<<" \n Done ... "<<std::endl;

}
