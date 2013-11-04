//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------- Script to create the MC Pileup distribution -------------------------------------------- 
//------------------------------------------------------------------------------------------------------------------------------------
#include "TChain.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>

const int numTrigger = 8;

TH1D* readTH1(const TString &fileName, const TString &histName, const TString &newHistName, bool useCurrentStyle); 
TH1D* createTH1(TH1D *histo,const TString &HistoName,const TString &HistoTitle,int nBin, double bound_low, double bound_up,const TString &xTitle);
void plotTH1(TH1* histo,const TString &filename);

void calcPUMCdist(){
  
  // Create Histogram for PU distribution
  TH1D *hPUMCDist = new TH1D("hPUMCDist","PU Distribution in MC",600,0,60);
  hPUMCDist -> Sumw2();
  hPUMCDist -> SetXTitle("# PU");
  
  // Add Chain to MC Sample
  TString DataFilename = "/scratch/hh/lustre/cms/user/telenz/mc/PhotonJetTuple2012/";
  DataFilename += "ak5FastPF_*.root";  
  TChain* chain = new TChain("GammaJetTree");
  chain -> Add(DataFilename); 
  
  // Set branch address to interesting variables
  float weight;
  float PUMCNumTruth; 
  chain->SetBranchAddress("Weight",&weight);
  chain->SetBranchAddress("PUMCNumTruth",&PUMCNumTruth);
  
  // Go through alle events and fill the MC Pileup-distribution with the true Pileup number
  cout<<"There are "<<chain->GetEntries()<<" entries in the Sample!!"<<endl<<endl;
  for(int n = 0; n < chain->GetEntries(); n++) {
    chain->GetEntry(n);
    if( (n+1)%500000 == 0 ) std::cout << "Event " << (n+1) << std::endl;
    hPUMCDist -> Fill(PUMCNumTruth);
  }
  
  // Normalize the histogram to one
  hPUMCDist -> Scale(1.0/hPUMCDist->Integral());     
  
  // Save and plot the PU distribution histogram
  TCanvas* canvas = new TCanvas("PUMC","PUMCDist",0,0,500,500);
  canvas -> SetLeftMargin(0.17);  
  hPUMCDist ->GetYaxis()->SetTitleOffset(1.5); 
  canvas -> cd();
  hPUMCDist -> Draw(); 
  canvas -> SaveAs("PUMCdist.pdf");
  TFile *f = new TFile("PUMCdist.root","RECREATE");
  f -> WriteTObject(hPUMCDist);
  f ->Close();
  delete f;

  // Read Data PU Distribution and calculate weight histograms for each trigger
  
  // Read Pileup-distribution from Data (seperately for different Triggers)
  TH1D* hPUDataDist[numTrigger] = {0};
  TH1D* hPUWeight[numTrigger]   = {0};
  TString histoName;
  char tname[100];
  for(int i=0; i<numTrigger; i++){
    sprintf(tname,"PU-Weights for %i. pT Bin",i);
    histoName.Form("PUWeight%i",i);
    hPUWeight[i] = createTH1(hPUWeight[i],histoName,tname,600,0,60,"# Pileup");
  }
  
  TString PUVersion = "cmssw5_2_5_Tag04-01-01PixelLumiCalcFineBinning/";
  hPUDataDist[0] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon20_CaloIdVL_IsoL.root","pileup;1","pileup0", 0);
  hPUDataDist[1] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon30_CaloIdVL_IsoL.root","pileup;1","pileup1", 0);
  hPUDataDist[2] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon50_CaloIdVL_IsoL.root","pileup;1","pileup2", 0);
  hPUDataDist[3] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon75_CaloIdVL_IsoL.root","pileup;1","pileup3", 0);
  hPUDataDist[4] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon90_CaloIdVL_IsoL.root","pileup;1","pileup4", 0);
  hPUDataDist[5] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon135.root","pileup;1","pileup5", 0);
  hPUDataDist[6] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon150.root","pileup;1","pileup6", 0);
  hPUDataDist[7] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon160.root","pileup;1","pileup7", 0);
  
  // Normalize PUData Distribution to 1.0
  for(int i = 0; i<numTrigger; i++ ) hPUDataDist[i] -> Scale(1.0/hPUDataDist[i]->Integral());
  
  // Calculate weights and normalize weight histogram to 1
  for(int i=0;i<numTrigger;i++) {
    hPUWeight[i] -> Divide(hPUDataDist[i],hPUMCDist);   
    hPUWeight[i] -> Scale(1.0/hPUWeight[i]->Integral());          
  }
    
  // Plot the weights for all different pTGamma bins
  
  for(int i=0; i<numTrigger; i++){
    sprintf(tname,"PUWeights_%i",i);
    plotTH1(hPUWeight[i],tname);
  }
   
}

//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------- Script to create the MC Pileup distribution -------------------------------------------- 
//------------------------------------------------------------------------------------------------------------------------------------


TH1D* readTH1(const TString &fileName, const TString &histName, const TString &newHistName, bool useCurrentStyle) {
 
  TFile f(fileName,"READ");
  TH1D *h = 0;
  f.GetObject(histName,h);
  if( h ) {
    h->SetDirectory(0);
    if( useCurrentStyle ) h->UseCurrentStyle();
    if( newHistName.Length() ) h->SetName(newHistName);
  } else {
    std::cerr << "ERROR in FileOps::readTH1: Histogram with name '" << histName << "' does not exist in file '" << fileName << "'\n.";
    f.Close();
    exit(-1);
  }
  f.Close();
  
  return h;
}


void plotTH1(TH1* histo,const TString &filename){

  TString tot_filename;
  TFile* f;

  TCanvas* canvas = new TCanvas(filename,filename,0,0,500,500);
  
  canvas -> SetLeftMargin(0.17);
 
  histo ->GetYaxis()->SetTitleOffset(1.5); 
  canvas -> cd();
  histo -> Draw(); 
  
  
  tot_filename = filename + ".pdf";
  canvas -> SaveAs(tot_filename);
  
  tot_filename = filename + ".root";
  f = new TFile(tot_filename,"RECREATE");
  f -> WriteTObject(histo);
  f ->Close();
  delete f;
}


TH1D* createTH1(TH1D *histo,const TString &HistoName,const TString &HistoTitle,int nBin, double bound_low, double bound_up,const TString &xTitle){
  
  histo = new TH1D(HistoName,HistoTitle,nBin,bound_low,bound_up);
  histo -> Sumw2();
  histo -> SetXTitle(xTitle);
    
  return histo;
}
