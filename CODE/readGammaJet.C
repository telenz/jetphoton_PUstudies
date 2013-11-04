#include "myDeclarations.h"
#include "myClasses.h"

#include <iostream>
#include "TChain.h"
#include "TString.h"


void readGammaJet(int nEvents) {

  chain  = new TChain("GammaJetTree");

  TString DataFilename;
  TString PUVersion;
  TString set = 2;  // 1=AB; 2= ABC; 3=AB rereco; 4=C
  
  // Read Ntuple
  if(date == 2012){
    if(type == 1){
      if(set == 1){
	DataFilename = "/scratch/hh/lustre/cms/user/telenz/data/PhotonJetTuple2012AB/";      
	if(jetType==1)      DataFilename += "ak5FastPF_*.root";
	else if(jetType==2) DataFilename += "ak5PFCHS_*.root"; 
	else cout<<"No input File available!"<<endl;
      }
      else if(set == 2){
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/data/PhotonJetTuple2012ABC/OnlyTightPhotons/";      
	if(jetType==1)      DataFilename += "ak5FastPF*.root";
	else if(jetType==2) DataFilename += "ak5PFCHS*.root";     
	else cout<<"No input File available!"<<endl;
      }
      else if(set == 3){
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/data/PhotonJetTuple2012ABC/otherSamples/OnlyTightPhotons/";      
	if(jetType==1)      DataFilename += "ak5FastPFAB.root";
	else if(jetType==2) DataFilename += "ak5PFCHSAB.root";     
	else cout<<"No input File available!"<<endl;
      }
      else if(set == 4){
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/data/PhotonJetTuple2012ABC/OnlyTightPhotons/";      
	if(jetType==1)      DataFilename += "ak5FastPFC*.root";
	else if(jetType==2) DataFilename += "ak5PFCHSC*.root";     
	else cout<<"No input File available!"<<endl;
      }
    }
    else if (type == 2){
      if(set == 1){
	DataFilename = "/scratch/hh/lustre/cms/user/telenz/mc/PhotonJetTuple2012/allTogether/OnlyTightPhotons/";
	if(jetType==1)      DataFilename += "ak5FastPF_*.root";  
	else if(jetType==2) DataFilename += "ak5PFCHS_*.root";  
	else cout<<"No input File available!"<<endl;
      }
      else if(set == 2 || set == 3 || set == 4){
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/mc/pythia_flat_53/MCPhoton2012/OnlyTightPhotons/";      
	//DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/mc/pythia_flat_535/MCPhoton2012/OnlyTightPhotons/";      
	if(jetType==1)      DataFilename += "ak5FastPF.root";
	else if(jetType==2) DataFilename += "ak5PFCHS.root";    
	else cout<<"No input File available!"<<endl;
      }

    }
    else cout<<"No such \"type\" available. Please choose a different number for variable type."<<endl;
  }
  else if (date == 2011){
    if(type == 1){
      DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/data/PhotonJetTuple2011/OnlyTightPhotons/";
      if(jetType==1) DataFilename += "ak5FastJet.root";
      else cout<<"No input File available!"<<endl;
    }
    else if(type == 2){
      DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/mc/mc_2011/OnlyTightPhotons/";
      if(jetType==1) DataFilename += "ak5FastJet.root";
      else cout<<"No input File available!"<<endl;
    }
  }
  else  cout<<"No input File available! (no none 2012 data available)"<<endl;
  cout<<"filename = "<<DataFilename<<endl; 
  chain->Add(DataFilename); 

  // Set branch addresses
  chain->SetBranchAddress("NobjPhoton",&nobjPhoton);
  chain->SetBranchAddress("NobjGenPhoton",&nobjGenPhoton);
  chain->SetBranchAddress("PhotonPt",photonPt);
  chain->SetBranchAddress("PhotonE",photonE);
  chain->SetBranchAddress("PhotonEta",photonEta);
  chain->SetBranchAddress("PhotonPhi",photonPhi);
  chain->SetBranchAddress("GenPhotonPt",genPhotonPt);
  chain->SetBranchAddress("GenPhotonEta",genPhotonEta);
  chain->SetBranchAddress("GenPhotonPhi",genPhotonPhi);
  chain->SetBranchAddress("PhotonIDTight",tight);

  chain->SetBranchAddress("NobjJet",&nobjJet);
  chain->SetBranchAddress("JetPt",jetPt);
  chain->SetBranchAddress("JetE",jetE);
  chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  chain->SetBranchAddress("JetCorrL1",jetCorrL1);
  if(date == 2012) chain->SetBranchAddress("JetCorrUncert",jetUncert);
  chain->SetBranchAddress("JetPhi",jetPhi);
  chain->SetBranchAddress("JetEta",jetEta);
  chain->SetBranchAddress("JetIDTight",JetTight);
  chain->SetBranchAddress("VtxN",&vtxN);

  chain->SetBranchAddress("PassesECALDeadCellBEFilter",&EcalBEFilter);
  chain->SetBranchAddress("PassesECALDeadCellTPFilter",&EcalTPFilter);

  chain->SetBranchAddress("RunNumber",&runNum);
  chain->SetBranchAddress("LumiBlockNumber",&lumiNum);
  chain->SetBranchAddress("EventNumber",&eventNum);  

  chain->SetBranchAddress("PhotonIsoECAL04",photonIsoEcal);  
  chain->SetBranchAddress("PhotonIsoHCAL04",photonIsoHcal);
  chain->SetBranchAddress("PhotonIsoTrk04",photonIsoTrk);

  chain->SetBranchAddress("JetFRBX",jetFRBX);
  chain->SetBranchAddress("JetFHPD",jetFHPD);
  chain->SetBranchAddress("JetEMF",jetEMF);
  
  chain->SetBranchAddress("Weight",&weight);
  chain->SetBranchAddress("Rho",&rho);
  chain->SetBranchAddress("PUMCNumTruth",&PUMCNumTruth);

  chain->SetBranchAddress("GenJetColPt",genJetPt);
  //chain->SetBranchAddress("GenJetPt",genJetPt);
  chain->SetBranchAddress("GenJetColEta",genJetEta);
  chain->SetBranchAddress("GenJetColPhi",genJetPhi);
  chain->SetBranchAddress("GenPartId_algo",genJetID);
  chain->SetBranchAddress("GenJetColJetIdx",genJetColJetIdx);
  chain->SetBranchAddress("NobjGenJet",&nobjGenJet);

  chain->SetBranchAddress("HltPhoton20",&hltPhoton[0]);
  chain->SetBranchAddress("HltPhoton30",&hltPhoton[1]);
  chain->SetBranchAddress("HltPhoton50",&hltPhoton[2]);
  chain->SetBranchAddress("HltPhoton75",&hltPhoton[3]);
  chain->SetBranchAddress("HltPhoton90",&hltPhoton[4]);  

  if(date == 2012){
    chain->SetBranchAddress("HltPhoton135",&hltPhoton[5]);
    chain->SetBranchAddress("HltPhoton150",&hltPhoton[6]);
    chain->SetBranchAddress("HltPhoton160",&hltPhoton[7]);
  }

  if(set == 1 && date == 2012)      PUVersion = "PUData_2012/weightsCMSSW5_2_5/";
  else if(set == 2 && date == 2012) PUVersion = "PUData_2012/weightsCMSSW5_3_3/";
  else if(set == 3 && date == 2012) PUVersion = "PUData_2012/weightsCMSSW5_3_3_AB/";
  else if(set == 4 && date == 2012) PUVersion = "PUData_2012/weightsCMSSW5_3_3_C/";
  else if(date == 2011) PUVersion = "PUData_2011/newVersion/";
  else cout<<"No PU data histograms available !!! (look in CODE/readGammaJets.C)"<<endl;
 
  
  // Read PU Weight-distribution  (seperately for different Triggers)
  hPUWeight[0] =  readTH1(PUVersion + "PUWeights_0.root","PUWeight0","PUWeight0", 0);
  hPUWeight[1] =  readTH1(PUVersion + "PUWeights_1.root","PUWeight1","PUWeight1", 0);
  hPUWeight[2] =  readTH1(PUVersion + "PUWeights_2.root","PUWeight2","PUWeight2", 0);
  hPUWeight[3] =  readTH1(PUVersion + "PUWeights_3.root","PUWeight3","PUWeight3", 0);
  hPUWeight[4] =  readTH1(PUVersion + "PUWeights_4.root","PUWeight4","PUWeight4", 0);
  
  if(date==2012){
    hPUWeight[5] =  readTH1(PUVersion + "PUWeights_5.root","PUWeight5","PUWeight5", 0);
    hPUWeight[6] =  readTH1(PUVersion + "PUWeights_6.root","PUWeight6","PUWeight6", 0);
    hPUWeight[7] =  readTH1(PUVersion + "PUWeights_7.root","PUWeight7","PUWeight7", 0);
  }

  // If tree contains less entries than nEvents use less events
  nMax = nEvents;
  if( chain->GetEntries() < nEvents ) nMax = chain->GetEntries();
  cout<<"There are "<<chain->GetEntries()<<" events in this sample!"<<endl;
  
}
