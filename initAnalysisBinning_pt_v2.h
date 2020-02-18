#ifndef INITANALYSISBINNING_PT_V2_H
#define INITANALYSISBINNING_PT_V2_H

/// centrality ///
#define N_CENTBIN 3
Int_t min_centbin[N_CENTBIN]={0 , 20 ,60};//the total N_CENTBIIN=200. So 0-20 represents 0-10% in centrality
Int_t max_centbin[N_CENTBIN]={20, 60 ,100};
TString label_centbin[N_CENTBIN]={"cent0to10","cent10to30","cent30to50"};
TString centbin[N_CENTBIN]={"0 <= cent < 10 %","10 <= cent < 30 %","30 <= cent < 50 %"};
//Int_t min_centbin[N_CENTBIN]={0 };//the total N_CENTBIIN=200. So 0-20 represents 0-10% in centrality
//Int_t max_centbin[N_CENTBIN]={100};
//TString label_centbin[N_CENTBIN]={"cent0to100"};
//TString centbin[N_CENTBIN]={"0 <= cent < 100 %"};


/// rapidity ///
#define N_YBIN 3
Float_t min_ybin[N_YBIN]={-1.0,-2.0,1.0};
Float_t max_ybin[N_YBIN]={1.0,-1.0,2.0};
TString label_ybin[N_YBIN]={"y-1.0to1.0","y-2.0to-1.0","y1.0to2.0"};
TString ybin[N_YBIN]={"-1.0 < y < 1.0","-2.0 < y < -1.0","1.0 < y < 2.0"};
//Float_t min_ybin[N_YBIN]={-2.0};
//Float_t max_ybin[N_YBIN]={2.0};
//TString label_ybin[N_YBIN]={"y-2.0to2.0"};
//TString ybin[N_YBIN]={"-2.0 < y < 2.0"};


//#define N_PTBIN 20
//Float_t min_pTbin[N_PTBIN]={1.0,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,10.0,12.5,15.0,20.0,30.0,40.0,60.0};
//Float_t max_pTbin[N_PTBIN]={2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,10.0,12.5,15.0,20.0,30.0,40.0,60.0,100.0};
//TString label_pTbin[N_PTBIN]={"pT1.0to2.0","pT2.0to2.5","pT2.5to3.0","pT3.0to3.5","pT3.5to4.0","pT4.0to4.5","pT4.5to5.0","pT5.0to5.5","pT5.5to6.0","pT6.0to6.5","pT6.5to7.0","pT7.0to8.0","pT8.0to10.0","pT10.0to12.5","pT12.5to15.0","pT15.0to20.0","pT20.0to30.0","pT30.0to40.0","pT40.0to60.0","pT60.0to100.0"};
//TString pTbin[N_PTBIN]={"1.0 < pT < 2.0 GeV","2.0 < pT < 2.5 GeV","2.5 < pT < 3.0 GeV","3.0 < pT < 3.5 GeV","3.5 < pT < 4.0 GeV","4.0 < pT < 4.5 GeV","4.5 < pT < 5.0 GeV","5.0 < pT < 5.5 GeV","5.5 < pT < 6.0 GeV","6.0 < pT < 6.5 GeV","6.5 < pT < 7.0 GeV","7.0 < pT < 8.0 GeV","8.0 < pT < 10.0 GeV","10.0 < pT < 12.5 GeV","12.5 < pT < 15.0 GeV","15.0 < pT < 20.0 GeV","20.0 < pT < 30.0 GeV","30.0 < pT < 40.0 GeV","40.0 < pT < 60.0 GeV","60.0 < pT < 100.0 GeV"};
#define N_PTBIN  10
Float_t min_pTbin[N_PTBIN]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0, 15.0, 30.0};
Float_t max_pTbin[N_PTBIN]={2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0,15.0, 30.0, 100.0};
TString label_pTbin[N_PTBIN]={"pT1.0to2.0","pT2.0to3.0","pT3.0to4.0","pT4.0to5.0","pT5.0to6.0","pT6.0to7.0","pT7.0to10.0","pT10.0to15.0","pT15.0to30.0","pT30.0to100.0"};
TString pTbin[N_PTBIN]={"1.0 < pT < 2.0 GeV", "2.0 < pT < 3.0 GeV", "3.0 < pT < 4.0 GeV","4.0 < pT < 5.0 GeV","5.0 < pT < 6.0 GeV","6.0 < pT < 7.0 GeV","7.0 < pT < 10.0 GeV", "10.0 < pT < 15.0","15.0 < pT < 30.0 GeV","30.0 < pT < 100.0 GeV"};


/// mass ///
//mass bins for v2 vs mass
#define N_MASSBIN 13
Float_t min_massbin[N_MASSBIN]={1.74,1.78,1.80,1.82,1.84,1.85,1.86 ,1.865,1.87,1.88,1.90,1.92,1.96};
Float_t max_massbin[N_MASSBIN]={1.78,1.80,1.82,1.84,1.85,1.86,1.865,1.87,1.88,1.90,1.92,1.96,2.00};
Double_t massbinning[N_MASSBIN+1]={1.74,1.78,1.80,1.82,1.84,1.85,1.86,1.865,1.87,1.88,1.90,1.92,1.96,2.00};
Float_t median_massbin[N_MASSBIN]={1.76,1.79,1.81,1.83,1.845,1.855,1.8625,1.8675,1.875,1.885,1.91,1.94,1.98};//for massVsV2 histogram filling
TString label_massbin[N_MASSBIN]={"m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12","m13"};

#endif
