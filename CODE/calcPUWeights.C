#include "myDeclarations.h"
#include "myFunctions.h"
#include <iostream>

void calcPUWeights(){
  
  
  TString PUVersion;
  if(date == 2012) PUVersion = "PUData_2012/cmssw5_2_5_Tag04-01-01PixelLumiCalcFineBinning/";
  else if(date == 2011) PUVersion = "PUData_2011/newVersion/";
  
  // Read Pileup-distribution from Data (seperately for different Triggers)
  hPUDataDist[0] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon20_CaloIdVL_IsoL.root","pileup","pileup0", 0);
  hPUDataDist[1] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon30_CaloIdVL_IsoL.root","pileup","pileup1", 0);
  hPUDataDist[2] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon50_CaloIdVL_IsoL.root","pileup","pileup2", 0);
  hPUDataDist[3] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon75_CaloIdVL_IsoL.root","pileup","pileup3", 0);
  hPUDataDist[4] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon90_CaloIdVL_IsoL.root","pileup","pileup4", 0);
  
  if(date==2012){
    hPUDataDist[5] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon135.root","pileup","pileup5", 0);
    hPUDataDist[6] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon150.root","pileup","pileup6", 0);
    hPUDataDist[7] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon160.root","pileup","pileup7", 0);
  }
  
  /*hPUDataDist[0] =  readTH1(PUVersion + "VtxN0_PF_data.root","13","pileup0", 0);
    hPUDataDist[1] =  readTH1(PUVersion + "VtxN1_PF_data.root","14","pileup1", 0);
    hPUDataDist[2] =  readTH1(PUVersion + "VtxN2_PF_data.root","15","pileup2", 0);
    hPUDataDist[3] =  readTH1(PUVersion + "VtxN3_PF_data.root","16","pileup3", 0);
    hPUDataDist[4] =  readTH1(PUVersion + "VtxN4_PF_data.root","17","pileup4", 0);
    hPUDataDist[5] =  readTH1(PUVersion + "VtxN5_PF_data.root","18","pileup5", 0);
    hPUDataDist[6] =  readTH1(PUVersion + "VtxN6_PF_data.root","19","pileup6", 0);
    hPUDataDist[7] =  readTH1(PUVersion + "VtxN7_PF_data.root","20","pileup7", 0);
  */
  
  char tname[100];
  
  // 1d PUMC Distribution
  for(int i=0; i<numTrigger; i++){
    if(i<numTrigger-1) sprintf(tname,"MC PU-Distribution in MC for %4.1f GeV < p_{T} < %4.1f GeV",bd[i],bd[i+1]);
    else               sprintf(tname,"MC PU-Distribution in MC for %4.1f GeV < p_{T}",bd[i]);
  }
  
  hPUMCDist = createTH1(hPUMCDist,"hPUMCDist",tname,600,0,60,"#PU");
  
  cout<<"numTrigger = "<<numTrigger<<endl;
  
  // Go through alle events and fill the MC Pileup-distribution with the true Pileup number
  for(int n = 0; n < 1000000; n++) {
  //for(int n = 0; n < chain->GetEntries(); n++) {
    if( n%(chain->GetEntries())/2. == 0 ) std::cout << "Event " << n << std::endl;
    chain->GetEntry(n);

    hPUMCDist->Fill(PUMCNumTruth,weight);
    //cout<<"hPUMCDist"<<PUMCNumTruth<<endl;
    //hPUMCDist->Fill(vtxN);
  }
  
  // Normalize PUData Distribution to PUMC entries
  for(int i = 0; i<numTrigger; i++ ) hPUDataDist[i] -> Scale(hPUMCDist->Integral()/hPUDataDist[i]->Integral());

  // Calculate weights and normalize weight histogram to 1
  for(int i=0;i<numTrigger;i++) {
    hPUWeight[i] -> Divide(hPUDataDist[i],hPUMCDist);
    hPUWeight[i] -> Scale(1.0/hPUWeight[i]->Integral());      
  }
    
  // Plot the weights for all different pTGamma bins
  for(int i=0; i<numTrigger; i++){
    sprintf(tname,"PUWeights_%i",i);
    plotTH1(hPUWeight[i],tname,0);
  }
  cout<<endl;
}
