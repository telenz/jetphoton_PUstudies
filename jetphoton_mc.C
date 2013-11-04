//  Script to calculate JET Energy Response in Gamma Jet Events --
//  (see: CMS AN-2010/076 and CMS AN-2010/421)                  --
// ---------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLorentzVector.h"

const int jetType = 2;      // 1 = PF_L1Fast , 2 = PF_CHS Jets , 3 = Calo_L1FastJet 
const int date    = 2012;
const int type    = 2;      //1 = Data , 2 = MC
//class CsysUncerMC;

#include "CODE/myDeclarations.h"
#include "CODE/myClasses.h"
#include "CODE/myFunctions.h"
#include "CODE/myClasses.C"
#include "CODE/myFunctions.C"
#include "CODE/calcPath.C"
#include "CODE/bookHistos.C"
#include "CODE/readGammaJet.C"
#include "CODE/calcPUWeights.C"
#include "CODE/applyCuts.C"
#include "CODE/draw.C"

void calcSample();
void calcScale();

  
// ---- Function implementations ---------------
// ---------------------------------------------

// Run main script
// This is the actual command to be typed
// used in the ROOT session.
// ---------------------------------------------
int runGammaJet(int nEvts, int step) {
  time_t start,end;

  time (&start);
  
  // -------------------------------------
  if(date == 2011 && numTrigger!=5){
    cout<<"Please correct the number of Triggers used!"<<endl;  
    return 0;
  }
  if(date == 2012 && numTrigger!=8){
    cout<<"Please correct the number of Triggers used!"<<endl;  
    return 0;
  }

  //-------- print out of Information --------------
  cout<<endl;
  if(type == 1) cout<<"Data!!"<<endl<<endl;
  else if(type == 2) cout<<"MC!!"<<endl<<endl;

  if(date == 2012) cout<<"date = 2012"<<endl<<endl;
  else if(date == 2011) cout<<"date = 2011"<<endl<<endl;
  
  if(jetType == 1) cout<<"jetType = PF_L1Fast"<<endl<<endl<<endl;
  else if(jetType == 2) cout<<"jetType = PF_CHSJets"<<endl<<endl<<endl;
  
  cout<<"calcPath() is executed!"<<endl<<endl;
  calcPath(step); 
  cout<<"bookHistos() is executed!"<<endl<<endl;;
  bookHistos();
  cout<<"readGammaJet() is executed!"<<endl<<endl;
  readGammaJet(nEvts);
  if(step==1){
    cout<<"calcSample() is executed!"<<endl<<endl;
    calcSample();
  }
  if(step==2){
    cout<<"calcScale() is executed!"<<endl<<endl;
    calcScale();
  }
  //calcScale(JetResponsePUle10, "_PUle10");
  //calcScale(JetResponsePUgt10le15, "_PUgt10le15");
  //calcScale(JetResponsePUgt15le20, "_PUgt15le20");
  //calcScale(JetResponsePUgt20, "_PUgt20");
  //cout<<"draw() is executed!"<<endl<<endl;
  //draw();
  //----------------
  time (&end);
  
  cout<<"elapsed time:"<<difftime(end,start)/60<<" minutes"<<endl;

  return 0;
}


// Aplly all cuts and fill histograms
// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void calcSample() {

  noDefMatch0 = 0;
  noDefMatch2 = 0;
  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){
      for(int k=0; k<alpha_int; k++){

	JetResponseJetHemisphere[i][j][k]         = new CResponse(1); 	
	JetImbalanceJetHemisphere[i][j][k]        = new CResponse(3);
	JetResponsePhotonHemisphere[i][j][k]      = new CResponse(1); 
	JetImbalancePhotonHemisphere[i][j][k]     = new CResponse(3);
	JetIntrinsic[i][j][k]                     = new CResponse(2);
	JetIntrinsicle5[i][j][k]                  = new CResponse(2);
	JetIntrinsicle10gt5[i][j][k]              = new CResponse(2);
	JetIntrinsicle15gt10[i][j][k]             = new CResponse(2);
	JetIntrinsicle20gt15[i][j][k]             = new CResponse(2);
	JetIntrinsicgt20[i][j][k]                 = new CResponse(2);

	//JetIntrinsicGluon[i][j][k] = new CResponse;	
	//JetIntrinsicQuark[i][j][k] = new CResponse;
	//JetResponsePUle10[i][j][k] = new CResponse;
	//JetResponsePUgt10le15[i][j][k] = new CResponse;
	//JetResponsePUgt15le20[i][j][k] = new CResponse;
	//JetResponsePUgt20[i][j][k] = new CResponse;
	
      }
      /*JetScaleResAlpha[i][j]  = new CScaleResAlpha;
      JetIntrinsicAlpha[i][j] = new CScaleResAlpha;
      JetImbalanceAlpha[i][j] = new CScaleResAlpha;   
      JetIntrinsicAlphaGluon[i][j] = new CScaleResAlpha;  
      JetIntrinsicAlphaQuark[i][j] = new CScaleResAlpha;
      */
    }
    /*
    JetScaleRes[j]         = new CScaleRes;
    JetIntrinsicPt[j]      = new CScaleRes;
    JetImbalancePt[j]      = new CScaleRes;
    JetIntrinsicPtGluon[j] = new CScaleRes;    
    JetIntrinsicPtQuark[j] = new CScaleRes;
    */

    //sysUncer[j] = new CsysUncerMC;
  }
 

  //------------------------------------------------------------------------------------------------------------
  // Loop over nMax entries, read variables and fill histogram

  std::cout << "Processing events" << std::endl;
  
  
  for(int n = 0; n < nMax; n++) {
    
    if( (n+1)%1000000 == 0 ) std::cout << "Event " << (n+1) << std::endl;
    
    // Get this event i.e. entry n and fill values into variables
    chain->GetEntry(n);
    
    // Calculating PUWeight for this event    
    for(int i=0; i<numTrigger-1; i++){
      if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){
	
	PUWeight = hPUWeight[i]->GetBinContent(hPUWeight[i]->FindBin(PUMCNumTruth));
	hPUgenMC[i]     -> Fill(PUMCNumTruth,weight*PUWeight);
	break;
      }
      else if(photonPt[0]>=bd[numTrigger-1]){
	
	PUWeight = hPUWeight[numTrigger-1]->GetBinContent(hPUWeight[numTrigger-1]->FindBin(PUMCNumTruth));
	hPUgenMC[numTrigger-1]    -> Fill(PUMCNumTruth,weight*PUWeight);
	break;
      }
    }

   
    deltaRJets = 0;
    //PUWeight=1;
    //weight = 1;
      
    lead_jet = 0;
    jet_2    = 0;
    alpha    = 0;

    
    
    //---------------------------------------------------------------------------------------------
    // Sort Jets and apply all cuts with following function
    bool testVar = applyCuts();
    if(!testVar) continue;    
    //---------------------------------------------------------------------------------------------
    
    
    
    /*
      for(int i=0; i<numTrigger-1; i++){
      if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){
      hVtxPtbinned[i] -> Fill(vtxN,weight*PUWeight);
      //hPUgenMC[i]     -> Fill(PUMCNumTruth,weight*PUWeight);
      break;
      }
      else if(photonPt[0]>=bd[numTrigger-1]){
      hVtxPtbinned[numTrigger-1]-> Fill(vtxN,weight*PUWeight); 
      //hPUgenMC[numTrigger-1]    -> Fill(PUMCNumTruth,weight*PUWeight);
      break;
      }
      }      
      
      hPhotonIsoEcal -> Fill(photonIsoEcal[0],weight*PUWeight);
      hPhotonIsoHcal -> Fill(photonIsoHcal[0],weight*PUWeight);
      hPhotonIsoTrk  -> Fill(photonIsoTrk[0],weight*PUWeight);
    */
   
    

    // Fill Response functions for whole sample
    float RatioPhotonPt = 0;
    
    for(int k=0; k<alpha_int; k++){

      if(alpha >= alphaBin[k] && alpha < alphaBin[k+1]){ 
		
	for(int j=0; j<eta_int; j++){
	  
	  if(std::abs(jetEta[corrJets.idx(lead_jet)])>= etaBin[j] && std::abs(jetEta[corrJets.idx(lead_jet)]) < etaBin[j+1] ){
	    
	    
	    for(int i=0; i<pt_int; i++){ 
	      
	      //if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){
	      if(genJetPt[genJetidx] >= bd[i] && genJetPt[genJetidx] < bd[i+1]){

		//if(jetPt[corrJets.idx(lead_jet)] >= bd[i] && jetPt[corrJets.idx(lead_jet)] < bd[i+1]){
	      

		/*
		  if(i==pt_int-1){
		  hImbalanceAlpha->Fill(imbalance,alpha,weight*PUWeight);	
		  }	
		
		  if(i==pt_int-1 && k==alpha_int-1){

		  hJet1PtPhoton1Pt->Fill(jetPt[corrJets.idx(lead_jet)],photonPt[0],weight*PUWeight);
		  hImbalanceJetN->Fill(imbalance,nobjJet,weight*PUWeight);
		  hImbalancePhoton1Pt->Fill(imbalance,photonPt[0],weight*PUWeight);
		  hImbalanceVtx->Fill(imbalance,vtxN,weight*PUWeight);
		  hImbalanceJetPt->Fill(imbalance,jetPt[corrJets.idx(lead_jet)],weight*PUWeight);
		  hImbalanceJetEta->Fill(imbalance,std::abs(jetEta[corrJets.idx(lead_jet)]),weight*PUWeight);
		  hImbalancePhotonEta->Fill(imbalance,std::abs(photonEta[0]),weight*PUWeight);
		  hImbalanceWeights->Fill(imbalance,weight*PUWeight,1);
		  hImbalanceJet2Pt->Fill(imbalance,jetPt[corrJets.idx(jet_2)],weight*PUWeight);
		 
		  
		  
		  hImbalanceDeltaRPhoton1stJet->Fill(imbalance,deltaRPhoton1stJet,weight*PUWeight);
		  hImbalanceDeltaRPhoton2ndJet->Fill(imbalance,deltaRPhoton2ndJet,weight*PUWeight);
		  hImbalanceDeltaR1stJet2ndJet->Fill(imbalance,deltaR1stJet2ndJet,weight*PUWeight);
		  hImbalanceEta2ndJet->Fill(imbalance,std::abs(jetEta[corrJets.idx(jet_2)]),weight*PUWeight);
		  hImbalanceRatioPhotonPt12->Fill(imbalance,RatioPhotonPt,weight*PUWeight);
		  hImbalanceRatioPhotongenPhoton->Fill(imbalance,photonPt[0]/genPhotonPt[0],weight*PUWeight);
		  hImbalanceDeltaPhiPhoton2ndJet->Fill(imbalance,std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(jet_2)]-photonPhi[0])),weight*PUWeight);
		  if(nobjPhoton>1) hImbalancePhoton2Pt->Fill(imbalance,photonPt[1],weight*PUWeight);
		  hImbalanceDeltaEtaPhoton2ndJet->Fill(imbalance, std::abs(jetEta[corrJets.idx(jet_2)] - photonEta[0]),weight*PUWeight);
		  //hImbalanceDeltaEtaPhoton2ndJet->Fill(imbalance, std::abs(jetEta[corrJets.idx(jet_2)] - jetEta[corrJets.idx(lead_jet)]),weight*PUWeight);
		  hImbalanceNobjPhoton->Fill(imbalance,nobjPhoton,weight*PUWeight);
		  hImbalance2ndJetCorr->Fill(imbalance,jetCorrL1[corrJets.idx(jet_2)]*jetCorrL2L3[corrJets.idx(jet_2)],weight*PUWeight);
		  hImbalance1stGenJetID->Fill(imbalance,std::abs(genJetID[genJetidx]),weight*PUWeight);
		  hPhotonEta2ndJet->Fill(photonEta[0],jetEta[corrJets.idx(jet_2)],weight*PUWeight);
		  hEnergyBalance -> Fill(sumJet.Pt()/photon.Pt(),PUWeight*weight);
		  hImbalanceRatioPhotonPtJetColPt -> Fill(imbalance,photonPt[0]/jetPt[corrJets.idx(photonidx)],PUWeight*weight);
		  }

		  else if(i==pt_int-1 && k==0){

		  hJet1PtPhoton1Pt2->Fill(jetPt[corrJets.idx(lead_jet)],photonPt[0],weight*PUWeight);
		  hImbalanceJetN2->Fill(imbalance,nobjJet,weight*PUWeight);
		  hImbalancePhoton1Pt2->Fill(imbalance,photonPt[0],weight*PUWeight);
		  hImbalanceVtx2->Fill(imbalance,vtxN,weight*PUWeight);
		  hImbalanceJetPt2->Fill(imbalance,jetPt[corrJets.idx(lead_jet)],weight*PUWeight);
		  hImbalanceJetEta2->Fill(imbalance,std::abs(jetEta[corrJets.idx(lead_jet)]),weight*PUWeight);
		  hImbalancePhotonEta2->Fill(imbalance,std::abs(photonEta[0]),weight*PUWeight);
		  hImbalanceWeights2->Fill(imbalance,weight*PUWeight,1);
		  hImbalanceJet2Pt2->Fill(imbalance,jetPt[corrJets.idx(jet_2)],weight*PUWeight);
		  hImbalanceAlpha2->Fill(imbalance,alpha,weight*PUWeight);
		  hImbalanceDeltaRPhoton1stJet2->Fill(imbalance,deltaRPhoton1stJet,weight*PUWeight);
		  hImbalanceDeltaRPhoton2ndJet2->Fill(imbalance,deltaRPhoton2ndJet,weight*PUWeight);
		  hImbalanceDeltaR1stJet2ndJet2->Fill(imbalance,deltaR1stJet2ndJet,weight*PUWeight);
		  hImbalanceEta2ndJet2->Fill(imbalance,std::abs(jetEta[corrJets.idx(jet_2)]),weight*PUWeight);
		  hImbalanceRatioPhotonPt122->Fill(imbalance,RatioPhotonPt,weight*PUWeight);	
		  hImbalanceRatioPhotongenPhoton2->Fill(imbalance,photonPt[0]/genPhotonPt[0],weight*PUWeight);
		  hImbalanceDeltaPhiPhoton2ndJet2->Fill(imbalance,std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(jet_2)]-photonPhi[0])),weight*PUWeight);
		  if(nobjPhoton>1) hImbalancePhoton2Pt2->Fill(imbalance,photonPt[1],weight*PUWeight);
		  hImbalanceDeltaEtaPhoton2ndJet2->Fill(imbalance, std::abs(jetEta[corrJets.idx(jet_2)] - photonEta[0]),weight*PUWeight);
		  hImbalanceNobjPhoton2->Fill(imbalance,nobjPhoton,weight*PUWeight);
		  hImbalance2ndJetCorr2->Fill(imbalance,jetCorrL1[corrJets.idx(jet_2)]*jetCorrL2L3[corrJets.idx(jet_2)],weight*PUWeight);
		  hImbalance1stGenJetID2->Fill(imbalance,std::abs(genJetID[genJetidx]),weight*PUWeight);
		  hPhotonEta2ndJet2->Fill(photonEta[0],jetEta[corrJets.idx(jet_2)],weight*PUWeight);
		  hEnergyBalance2 -> Fill(sumJet.Pt()/photon.Pt(),PUWeight*weight);
		  hImbalanceRatioPhotonPtJetColPt2 -> Fill(imbalance,photonPt[0]/jetPt[corrJets.idx(photonidx)],PUWeight*weight);
		  }
		*/
		
		float deltaphi1stJetPhoton = 0;
		float deltaphi2ndJetPhoton = 0;
		if(jet_2 != -1){
		  if(std::abs(TVector2::Phi_mpi_pi((jetPhi[corrJets.idx(lead_jet)]+photonPhi[0])/2. - jetPhi[corrJets.idx(jet_2)])) < TMath::Pi()/2.){
		    
		    deltaphi1stJetPhoton = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(lead_jet)]-photonPhi[0]));
		  }
		  else deltaphi1stJetPhoton = TMath::Pi()*2. - std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(lead_jet)]-photonPhi[0]));
		  
		  
		  deltaphi2ndJetPhoton = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(jet_2)]-photonPhi[0])); 
		}
	       		
		  
		if(deltaphi2ndJetPhoton > deltaphi1stJetPhoton/2. || jet_2 == -1){

		  JetResponseJetHemisphere[i][j][k] -> hResponse      -> Fill(response,weight*PUWeight);
		  JetImbalanceJetHemisphere[i][j][k]-> hResponse      -> Fill(imbalance,weight*PUWeight);
		  JetResponseJetHemisphere[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);
		  JetResponseJetHemisphere[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);	
		  JetResponseJetHemisphere[i][j][k] -> hPtLeadingJet  -> Fill(jetPt[corrJets.idx(lead_jet)],weight*PUWeight);  

		  		  	
		}  
		else{
		  JetResponsePhotonHemisphere[i][j][k] -> hResponse -> Fill(response,weight*PUWeight);
		  JetImbalancePhotonHemisphere[i][j][k]-> hResponse -> Fill(imbalance,weight*PUWeight); 
		  JetResponsePhotonHemisphere[i][j][0] -> hPt       -> Fill(photonPt[0],weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][k] -> hAlpha    -> Fill(alpha,weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][k] -> hPtLeadingJet  -> Fill(jetPt[corrJets.idx(lead_jet)],weight*PUWeight);	
		}

		JetIntrinsic[i][j][k] -> hResponse      -> Fill(intrinsic,weight*PUWeight); 
		JetIntrinsic[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);	
		JetIntrinsic[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);
		JetIntrinsic[i][j][k] -> hPtLeadingJet  -> Fill(jetPt[corrJets.idx(lead_jet)],weight*PUWeight);	
		
		
		if(PUMCNumTruth <= 8){
		  JetIntrinsicle5[i][j][k] -> hResponse      -> Fill(intrinsic,weight*PUWeight); 
		  JetIntrinsicle5[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);	
		  JetIntrinsicle5[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);
		  JetIntrinsicle5[i][j][k] -> hPtLeadingJet  -> Fill(jetPt[corrJets.idx(lead_jet)],weight*PUWeight);	
		  hVtxN1[i][j]->Fill(PUMCNumTruth,weight*PUWeight); 
		}
		else if(PUMCNumTruth > 8 && PUMCNumTruth <= 16){
		  JetIntrinsicle10gt5[i][j][k] -> hResponse      -> Fill(intrinsic,weight*PUWeight); 
		  JetIntrinsicle10gt5[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);	
		  JetIntrinsicle10gt5[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);
		  JetIntrinsicle10gt5[i][j][k] -> hPtLeadingJet  -> Fill(jetPt[corrJets.idx(lead_jet)],weight*PUWeight);	
		  hVtxN2[i][j]->Fill(PUMCNumTruth,weight*PUWeight); 
		}
		else if(PUMCNumTruth > 16 && PUMCNumTruth <= 24){
		  JetIntrinsicle15gt10[i][j][k] -> hResponse      -> Fill(intrinsic,weight*PUWeight); 
		  JetIntrinsicle15gt10[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);	
		  JetIntrinsicle15gt10[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);
		  JetIntrinsicle15gt10[i][j][k] -> hPtLeadingJet  -> Fill(jetPt[corrJets.idx(lead_jet)],weight*PUWeight);	
		  hVtxN3[i][j]->Fill(PUMCNumTruth,weight*PUWeight); 
		}
		else if(PUMCNumTruth > 24 && PUMCNumTruth <= 32){
		  JetIntrinsicle20gt15[i][j][k] -> hResponse      -> Fill(intrinsic,weight*PUWeight); 
		  JetIntrinsicle20gt15[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);	
		  JetIntrinsicle20gt15[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);
		  JetIntrinsicle20gt15[i][j][k] -> hPtLeadingJet  -> Fill(jetPt[corrJets.idx(lead_jet)],weight*PUWeight);	
		  hVtxN4[i][j]->Fill(PUMCNumTruth,weight*PUWeight); 
		}
		else{
		  JetIntrinsicgt20[i][j][k] -> hResponse      -> Fill(intrinsic,weight*PUWeight); 
		  JetIntrinsicgt20[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);	
		  JetIntrinsicgt20[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);
		  JetIntrinsicgt20[i][j][k] -> hPtLeadingJet  -> Fill(jetPt[corrJets.idx(lead_jet)],weight*PUWeight);	
		  hVtxN5[i][j]->Fill(PUMCNumTruth,weight*PUWeight); 
		}

		hNPU -> Fill(PUMCNumTruth, weight*PUWeight);
		hNPV -> Fill(vtxN, weight*PUWeight);
  		

		if(jet_2 != -1 && i==pt_int-1 && k==alpha_int-1){
		  hDeltaPhi1st2ndJet ->Fill(TVector2::Phi_0_2pi(jetPhi[corrJets.idx(lead_jet)]-jetPhi[corrJets.idx(jet_2)]),weight*PUWeight); 
		  hDeltaPhi1st2ndJetDeltaPt ->Fill(TVector2::Phi_0_2pi(jetPhi[corrJets.idx(lead_jet)]-jetPhi[corrJets.idx(jet_2)]),std::abs(photonPt[0]-corrJets.pt(lead_jet)),weight*PUWeight); 
		}

		
		
		/*
		  if(PUMCNumTruth<= 10){
		  JetResponsePUle10[i][j][k] -> hResponse -> Fill(response,weight*PUWeight); 
		  JetResponsePUle10[i][j][0] -> hPt       -> Fill(photonPt[0],weight*PUWeight);
		  JetResponsePUle10[i][j][k] -> hAlpha    -> Fill(alpha,weight*PUWeight);		  
		  hPUsysY[0][i][j] -> Fill(PUMCNumTruth,weight*PUWeight);
		  }
		  else if(PUMCNumTruth> 10 && PUMCNumTruth<=15){
		  JetResponsePUgt10le15[i][j][k] -> hResponse -> Fill(response,weight*PUWeight); 
		  JetResponsePUgt10le15[i][j][0] -> hPt       -> Fill(photonPt[0],weight*PUWeight);
		  JetResponsePUgt10le15[i][j][k] -> hAlpha    -> Fill(alpha,weight*PUWeight);	
		  hPUsysY[1][i][j] -> Fill(PUMCNumTruth,weight*PUWeight);	  
		  }
		  else if(PUMCNumTruth> 15 && PUMCNumTruth<=20){
		  JetResponsePUgt15le20[i][j][k] -> hResponse -> Fill(response,weight*PUWeight); 
		  JetResponsePUgt15le20[i][j][0] -> hPt       -> Fill(photonPt[0],weight*PUWeight);
		  JetResponsePUgt15le20[i][j][k] -> hAlpha    -> Fill(alpha,weight*PUWeight);	
		  hPUsysY[2][i][j] -> Fill(PUMCNumTruth,weight*PUWeight);	  
		  }
		  else if(PUMCNumTruth> 20 ){
		  JetResponsePUgt20[i][j][k] -> hResponse -> Fill(response,weight*PUWeight); 
		  JetResponsePUgt20[i][j][0] -> hPt       -> Fill(photonPt[0],weight*PUWeight);
		  JetResponsePUgt20[i][j][k] -> hAlpha    -> Fill(alpha,weight*PUWeight);		  
		  hPUsysY[3][i][j] -> Fill(PUMCNumTruth,weight*PUWeight);
		  }
		*/
		break;
	      }
	      else{
		//cout<<"i = "<<i<<endl;
                
                //cout<<"photonPt = "<<photonPt[0]<<endl;
		//cout<<"Something wrong with Trigger and Pt Bin bounds (MC)"<<endl<<endl;
	      }
	    }
	    break;
	  }
	}
	break;
      }
    }
    

   
    /*
    // Fill Response functions for subsamples gluon and light quarks
    for(int k=0; k<alpha_int; k++){
    if(alpha >= alphaBin[k] && alpha < alphaBin[k+1]){ 
	
    for(int j=0; j<eta_int; j++){
	  
    if(std::abs(jetEta[corrJets.idx(lead_jet)])>= etaBin[j] && std::abs(jetEta[corrJets.idx(lead_jet)]) < etaBin[j+1] ){
	    
    for(int i=0; i<pt_int; i++){ 
	      
    if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){
		
    if(std::abs(genJetID[corrJets.idx(lead_jet)])<=3){
		  
    JetIntrinsicQuark[i][j][k] -> hResponse -> Fill(intrinsic,weight*PUWeight); 		  
    JetIntrinsicQuark[i][j][0] -> hPt       -> Fill(photonPt[0],weight*PUWeight);
    JetIntrinsicQuark[i][j][k] -> hAlpha    -> Fill(alpha,weight*PUWeight);
    //cout<<"Quark!"<<endl;
                  
    }
		
    else if(std::abs(genJetID[corrJets.idx(lead_jet)])==21){
 
    JetIntrinsicGluon[i][j][k] -> hResponse -> Fill(intrinsic,weight*PUWeight); 	 
    JetIntrinsicGluon[i][j][0] -> hPt       -> Fill(photonPt[0],weight*PUWeight);
    JetIntrinsicGluon[i][j][k] -> hAlpha    -> Fill(alpha,weight*PUWeight);
    //cout<<"Gluon!"<<endl;
    }
    break;
    }
    else{
    //cout<<"i = "<<i<<endl;
    //cout<<"photonPt = "<<photonPt[0]<<endl;
    //cout<<"Something wrong with Trigger and Pt Bin bounds (MC)"<<endl<<endl;
    }
    }
    break;
    }
    }
    break;
    }
    }
    */
    
    // SMALL RESPONSE DIAGNOSTIC::Plots just for small response region
    /*if( corrJets.pt(lead_jet)/photonPt[0] <= 0.3 ){
      hPhotonEta_low_resp_reg -> Fill(photonEta[0],weight*PUWeight);
      
      // Write files for strange events in 2011 data
      if(!filestr.is_open()) cout<<"wrong"<<endl;
      filestr<<runNum<<":"<<lumiNum<<":"<<eventNum<<endl;	   
      if(!EventVariables.is_open()) cout<<"EventVariable.txt is not open"<<endl;
      EventVariables<<photonPt[0]<<":"<<corrJets.pt(lead_jet)<<":"<<photonEta[0]<<":"<<jetEta[corrJets.idx(lead_jet)]<<":"<<photonPhi[0]<<":"<<jetPhi[corrJets.idx(lead_jet)]<<":"<<photonIsoEcal[0]<<":"<<photonIsoHcal[0]<<":"<<photonIsoTrk[0]<<endl;
      }
    
      if(corrJets.pt(lead_jet)/photonPt[0]<=1.0 && photonPt[0]>130.){
      hSmallJetPtResponse -> Fill(corrJets.pt(lead_jet),corrJets.pt(lead_jet)/photonPt[0],weight*PUWeight);
      hPhotonEta_high_pt_reg_Response -> Fill(photonEta[0],corrJets.pt(lead_jet)/photonPt[0],weight*PUWeight);	 
      //cout<<"eventNum="<<eventNum<<endl;         
      }*/
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
    
    /*// Fill 2d Diagram Response L1 correctio
      hResponseL1 -> Fill(corrJets.pt(lead_jet)/photonPt[0],jetCorrL1[corrJets.idx(lead_jet)],weight*PUWeight);
      hResponseL2L3 -> Fill(corrJets.pt(lead_jet)/photonPt[0],jetCorrL2L3[corrJets.idx(lead_jet)],weight*PUWeight);
      
      //Fill 2d histogram Response agaist Photon Eta
      hRespPhotonEta -> Fill(corrJets.pt(lead_jet)/photonPt[0],photonEta[0],weight*PUWeight);
      
      // Fill 2d Diagram Response against JetEMF JetFHPD JetFRBX
      hResponseJetEMF -> Fill(response,jetEMF[corrJets.idx(lead_jet)],weight*PUWeight);
      hResponseJetFHPD -> Fill(response,jetFHPD[corrJets.idx(lead_jet)],weight*PUWeight);
      hResponseJetFRBX -> Fill(response,jetFRBX[corrJets.idx(lead_jet)],weight*PUWeight);
      
      // Fill 2d Diagram Response against PhotonEMF PhotonFHPD PhotonFRBX
      if(photonidx!=-1){
      hResponsePhotonEMF -> Fill(response,jetEMF[corrJets.idx(photonidx)],weight*PUWeight);
      hResponsePhotonFHPD -> Fill(response,jetFHPD[corrJets.idx(photonidx)],weight*PUWeight);
      hResponsePhotonFRBX -> Fill(response,jetFRBX[corrJets.idx(photonidx)],weight*PUWeight);
      if(response<=1 && corrJets.pt(photonidx)/photonPt[0]<=2)    hResponsePhotonPtRatio -> Fill(response,corrJets.pt(photonidx)/photonPt[0],weight*PUWeight);
      }*/
    
    // Fill Histogram for first Photon pT 
    hPhoton1Pt -> Fill(photonPt[0],weight*PUWeight);
    hJet2Pt    -> Fill(jetPt2nd,weight*PUWeight);
    hJet1Pt    -> Fill(jetPt1stJet,weight*PUWeight);
    
    // Fill 2d histogram with first and second Photon of one event
    hPhoton12Pt->Fill(photonPt[0],photonPt[1],weight*PUWeight);
    
    nocut = nocut +1; 
   
  } // End of loop over entries


  TFile *f;

  TString TotFilename;  

  // Scale Response functions to MC statistics
  for(int i=0;i<pt_int;i++){
    for(int j=0; j<eta_int; j++){

      JetImbalanceJetHemisphere[i][j][0]    -> hPt    = JetResponseJetHemisphere[i][j][0]    -> hPt;  
      JetImbalancePhotonHemisphere[i][j][0] -> hPt    = JetResponsePhotonHemisphere[i][j][0] -> hPt;       
          
      TotFilename = RootPath + "hPt_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";
      f = new TFile(TotFilename,"RECREATE");
      f -> WriteTObject(JetResponseJetHemisphere[i][j][0]->hPt,"histo");
      f->Close();
      delete f;
      
      TotFilename = RootPath + "hPt_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";  
      f = new TFile(TotFilename,"RECREATE");
      f -> WriteTObject(JetResponsePhotonHemisphere[i][j][0]->hPt,"histo");
      f->Close();
      delete f;
  
      TotFilename = RootPath + "hPt_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";
      f = new TFile(TotFilename,"RECREATE");
      f -> WriteTObject(JetImbalanceJetHemisphere[i][j][0]->hPt,"histo");
      f->Close();
      delete f;
 
      TotFilename = RootPath + "hPt_imbalance_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";      
      f = new TFile(TotFilename,"RECREATE");
      f -> WriteTObject(JetImbalancePhotonHemisphere[i][j][0]->hPt,"histo");
      f->Close();
      delete f;

      TotFilename = RootPath + "hPt_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";          
      f = new TFile(TotFilename,"RECREATE");
      f -> WriteTObject(JetIntrinsic[i][j][0]->hPt,"histo");
      f->Close();
      delete f;
      TotFilename = RootPath + "hPt_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + "_le5.root";          
      f = new TFile(TotFilename,"RECREATE");
      f -> WriteTObject(JetIntrinsicle5[i][j][0]->hPt,"histo");
      f->Close();
      delete f;
      TotFilename = RootPath + "hPt_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + "_le10gt5.root";          
      f = new TFile(TotFilename,"RECREATE");
      f -> WriteTObject(JetIntrinsicle10gt5[i][j][0]->hPt,"histo");
      f->Close();
      delete f;
      TotFilename = RootPath + "hPt_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + "_le15gt10.root";          
      f = new TFile(TotFilename,"RECREATE");
      f -> WriteTObject(JetIntrinsicle15gt10[i][j][0]->hPt,"histo");
      f->Close();
      delete f;
      TotFilename = RootPath + "hPt_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + "_le20gt15.root";          
      f = new TFile(TotFilename,"RECREATE");
      f -> WriteTObject(JetIntrinsicle20gt15[i][j][0]->hPt,"histo");
      f->Close();
      delete f;
      TotFilename = RootPath + "hPt_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + "_gt20.root";          
      f = new TFile(TotFilename,"RECREATE");
      f -> WriteTObject(JetIntrinsicgt20[i][j][0]->hPt,"histo");
      f->Close();
      delete f;


      for(int k=0; k<alpha_int; k++){	

	JetImbalanceJetHemisphere[i][j][k]    -> hAlpha = JetResponseJetHemisphere[i][j][k]    -> hAlpha;	                
	JetImbalancePhotonHemisphere[i][j][k] -> hAlpha = JetResponsePhotonHemisphere[i][j][k] -> hAlpha;
  
	TotFilename = RootPath + "hAlpha_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetResponseJetHemisphere[i][j][k]->hAlpha,"histo");
	f->Close();
	delete f;

	TotFilename = RootPath + "hAlpha_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";         
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetResponsePhotonHemisphere[i][j][k]->hAlpha,"histo");
	f->Close();
	delete f;

	TotFilename = RootPath + "hAlpha_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";     
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetImbalanceJetHemisphere[i][j][k]->hAlpha,"histo");
	f->Close();
	delete f;

	TotFilename = RootPath + "hAlpha_imbalance_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";  
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetImbalancePhotonHemisphere[i][j][k]->hAlpha,"histo");
	f->Close();
	delete f;

	TotFilename = RootPath + "hAlpha_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";          
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsic[i][j][k]->hAlpha,"histo");
	f->Close();
	delete f;
	TotFilename = RootPath + "hAlpha_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + "_le5.root";      
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsicle5[i][j][k]->hAlpha,"histo");
	f->Close();
	delete f;
	TotFilename = RootPath + "hAlpha_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + "_le10gt5.root"; 
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsicle10gt5[i][j][k]->hAlpha,"histo");
	f->Close();
	delete f;
	TotFilename = RootPath + "hAlpha_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + "_le15gt10.root";
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsicle15gt10[i][j][k]->hAlpha,"histo");
	f->Close();
	delete f;
	TotFilename = RootPath + "hAlpha_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + "_le20gt15.root";
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsicle20gt15[i][j][k]->hAlpha,"histo");
	f->Close();
	delete f;
	TotFilename = RootPath + "hAlpha_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + "_gt20.root"; 
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsicgt20[i][j][k]->hAlpha,"histo");
	f->Close();
	delete f;

	if(JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries() != 0 && JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries() != 0 ){
	  JetResponseJetHemisphere[i][j][k]->hResponse->Scale((JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries())/(JetResponseJetHemisphere[i][j][k]->hResponse->Integral()));
	  JetResponsePhotonHemisphere[i][j][k]->hResponse->Scale((JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries())/(JetResponsePhotonHemisphere[i][j][k]->hResponse->Integral()));
	  JetImbalancePhotonHemisphere[i][j][k]->hResponse->Scale((JetImbalancePhotonHemisphere[i][j][k]->hResponse->GetEntries())/(JetImbalancePhotonHemisphere[i][j][k]->hResponse->Integral()));
	  JetImbalanceJetHemisphere[i][j][k]->hResponse->Scale((JetImbalanceJetHemisphere[i][j][k]->hResponse->GetEntries())/(JetImbalanceJetHemisphere[i][j][k]->hResponse->Integral()));
	  JetIntrinsic[i][j][k]->hResponse->Scale((JetIntrinsic[i][j][k]->hResponse->GetEntries())/(JetIntrinsic[i][j][k]->hResponse->Integral())); 
	  JetIntrinsicle5[i][j][k]->hResponse->Scale((JetIntrinsicle5[i][j][k]->hResponse->GetEntries())/(JetIntrinsicle5[i][j][k]->hResponse->Integral())); 
	  JetIntrinsicle10gt5[i][j][k]->hResponse->Scale((JetIntrinsicle10gt5[i][j][k]->hResponse->GetEntries())/(JetIntrinsicle10gt5[i][j][k]->hResponse->Integral())); 
	  JetIntrinsicle15gt10[i][j][k]->hResponse->Scale((JetIntrinsicle15gt10[i][j][k]->hResponse->GetEntries())/(JetIntrinsicle15gt10[i][j][k]->hResponse->Integral())); 
	  JetIntrinsicle20gt15[i][j][k]->hResponse->Scale((JetIntrinsicle20gt15[i][j][k]->hResponse->GetEntries())/(JetIntrinsicle20gt15[i][j][k]->hResponse->Integral())); 
	  JetIntrinsicgt20[i][j][k]->hResponse->Scale((JetIntrinsicgt20[i][j][k]->hResponse->GetEntries())/(JetIntrinsicgt20[i][j][k]->hResponse->Integral())); 
	}


	TotFilename = RootPath + "response_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";         
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetResponseJetHemisphere[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;
		
	TotFilename = RootPath + "response_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";   
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetImbalanceJetHemisphere[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;
	
	TotFilename = RootPath + "response_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";           
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetResponsePhotonHemisphere[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;
	
	TotFilename = RootPath + "response_imbalance_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root"; 
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetImbalancePhotonHemisphere[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;
	
	TotFilename = RootPath + "response_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";       
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsic[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;
	TotFilename = RootPath + "response_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + "_le5.root"; 
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsicle5[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;
	TotFilename = RootPath + "response_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + "_le10gt5.root";
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsicle10gt5[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;
	TotFilename = RootPath + "response_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + "_le15gt10.root";
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsicle15gt10[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;
	TotFilename = RootPath + "response_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + "_le20gt15.root";       
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsicle20gt15[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;
	TotFilename = RootPath + "response_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + "_gt20.root";
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetIntrinsicgt20[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;

      }
    }
  }

  // Save some control plots


  TotFilename = RootPath + "hPhoton1Pt.root";
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hPhoton1Pt,"histo");
  f->Close();
  delete f;
  TotFilename = RootPath + "hJet1Pt.root";
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hJet1Pt,"histo");
  f->Close();
  delete f;
  TotFilename = RootPath + "hJet2Pt.root";
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hJet2Pt,"histo");
  f->Close();
  delete f;
  TotFilename = RootPath + "hDeltaPhi.root";
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hDeltaPhi,"histo");
  f->Close();
  delete f;
  TotFilename = RootPath + "hNPU.root";
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hNPU,"histo");
  f->Close();
  delete f;
  TotFilename = RootPath + "hNPV.root";
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hNPV,"histo");
  f->Close();
  delete f;

  

  for(int i=0;i<pt_int;i++){
    for(int j=0; j<eta_int; j++){
          
    
  TotFilename = RootPath + "hVtxN1_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hVtxN1[i][j],"histo");
  f->Close();
  delete f;

  TotFilename = RootPath + "hVtxN2_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hVtxN2[i][j],"histo");
  f->Close();
  delete f;

  TotFilename = RootPath + "hVtxN3_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";  
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hVtxN3[i][j],"histo");
  f->Close();
  delete f;
  TotFilename = RootPath + "hVtxN4_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";
  
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hVtxN4[i][j],"histo");
  f->Close();
  delete f;
  TotFilename = RootPath + "hVtxN5_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";
  
  f = new TFile(TotFilename,"RECREATE");
  f -> WriteTObject(hVtxN5[i][j],"histo");
  f->Close();
  delete f;

    }
  }


  // CUTflow
  cout<<endl;
  cout<<"CUTFLOW:"<<endl;
  cout<<"cut1 = "<<cut1<<endl;
  cout<<"cut2 = "<<cut2<<endl;
  cout<<"cut3 = "<<cut3<<endl;
  cout<<"cut4 = "<<cut4<<endl;
  cout<<"cut5 = "<<cut5<<endl;
  cout<<"cut6 = "<<cut6<<endl;
  cout<<"cut7 = "<<cut7<<endl;
  cout<<"cut8 = "<<cut8<<endl;
  cout<<"cut9 = "<<cut9<<endl;
  cout<<"cut10 = "<<cut10<<endl;
  cout<<"cut11 = "<<cut11<<endl;
  for(int i=0; i<numTrigger; i++) cout<<"cut12["<<i<<"] = "<<cut12[i]<<endl;
  for(int i=0; i<3; i++) cout<<"cut13["<<i<<"] = "<<cut13[i]<<endl;
  cout<<"nocut = "<<nocut<<endl;
  cout<<"number of all events = "<<cut1+cut2+cut3+cut4+cut5+cut6+cut7+cut8+cut9+cut10+cut11+cut12[0]+cut12[1]+cut12[2]+cut12[3]+cut12[4]+cut12[5]+cut12[6]+cut12[7]+cut13[0]+cut13[1]+cut13[2]+nocut<<endl<<endl;

  
  delete chain;
  filestr.close();   
  EventVariables.close();
  
}


// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// Now take the response functions and fit to the core region a gauss function 
// so that one gets the jet energy scale and jet energy resolution


// Calculate sigma and mean for all Response functions
void calcScale(){

  std::cout << "Read Response Histograms!" << std::endl<< std::endl;

  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){
      for(int k=0; k<alpha_int; k++){

	JetResponseJetHemisphere[i][j][k]         = new CResponse(1); 	
	JetImbalanceJetHemisphere[i][j][k]        = new CResponse(3);
	JetResponsePhotonHemisphere[i][j][k]      = new CResponse(1); 
	JetImbalancePhotonHemisphere[i][j][k]     = new CResponse(3);
	JetIntrinsic[i][j][k]                     = new CResponse(2);
	
      }
    }
  }

  TFile *file;
  TString TotFilename;  
  for(int i=0;i<pt_int;i++){
    for(int j=0; j<eta_int; j++){

        
      TotFilename = RootPath + "hPt_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";    
      file = new TFile(TotFilename);      
      JetResponseJetHemisphere[i][j][0]->hPt =  (TH1D*) gDirectory->Get("histo");
      JetResponseJetHemisphere[i][j][0]->hPt -> SetDirectory(0);
      delete file;  
      
      TotFilename = RootPath + "hPt_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";    
      file = new TFile(TotFilename);      
      JetResponsePhotonHemisphere[i][j][0]->hPt =  (TH1D*) gDirectory->Get("histo");
      JetResponsePhotonHemisphere[i][j][0]->hPt -> SetDirectory(0);
      delete file;       
      
      TotFilename = RootPath + "hPt_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";     
      file = new TFile(TotFilename);      
      JetImbalanceJetHemisphere[i][j][0]->hPt =  (TH1D*) gDirectory->Get("histo");
      JetImbalanceJetHemisphere[i][j][0]->hPt -> SetDirectory(0);
      delete file;           
       
      TotFilename = RootPath + "hPt_imbalance_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";    
      file = new TFile(TotFilename);      
      JetImbalancePhotonHemisphere[i][j][0]->hPt =  (TH1D*) gDirectory->Get("histo");
      JetImbalancePhotonHemisphere[i][j][0]->hPt -> SetDirectory(0);
      delete file;           
      TString FilenameIntrinsic = "_le5.root";
      TotFilename = RootPath + "hPt_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + FilenameIntrinsic;    
      file = new TFile(TotFilename);      
      JetIntrinsic[i][j][0]->hPt =  (TH1D*) gDirectory->Get("histo");
      JetIntrinsic[i][j][0]->hPt -> SetDirectory(0);
      delete file;                         
      

      for(int k=0; k<alpha_int; k++){
	
	TotFilename = RootPath + "response_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetResponseJetHemisphere[i][j][k]->hResponse =  (TH1D*) gDirectory->Get("histo");
	JetResponseJetHemisphere[i][j][k]->hResponse -> SetDirectory(0);
	delete file; 
		
	TotFilename = RootPath + "response_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetImbalanceJetHemisphere[i][j][k]->hResponse =  (TH1D*) gDirectory->Get("histo");
	JetImbalanceJetHemisphere[i][j][k]->hResponse -> SetDirectory(0);
	delete file;  
		
	TotFilename = RootPath + "response_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetResponsePhotonHemisphere[i][j][k]->hResponse =  (TH1D*) gDirectory->Get("histo");
	JetResponsePhotonHemisphere[i][j][k]->hResponse -> SetDirectory(0);
	delete file; 
	  	
	TotFilename = RootPath + "response_imbalance_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";
	file = new TFile(TotFilename);      
	JetImbalancePhotonHemisphere[i][j][k]->hResponse =  (TH1D*) gDirectory->Get("histo");
	JetImbalancePhotonHemisphere[i][j][k]->hResponse -> SetDirectory(0);
	delete file;   
		
	TotFilename = RootPath + "response_intrinsic_in_" + (long)(i+1) +"_Pt_bin_"+ (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + FilenameIntrinsic;
	file = new TFile(TotFilename);      
	JetIntrinsic[i][j][k]->hResponse =  (TH1D*) gDirectory->Get("histo");
	JetIntrinsic[i][j][k]->hResponse -> SetDirectory(0);
	delete file;    
	
	TotFilename = RootPath + "hAlpha_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetResponseJetHemisphere[i][j][k]->hAlpha =  (TH1D*) gDirectory->Get("histo");
	JetResponseJetHemisphere[i][j][k]->hAlpha -> SetDirectory(0);
	delete file; 
     
	TotFilename = RootPath + "hAlpha_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetResponsePhotonHemisphere[i][j][k]->hAlpha =  (TH1D*) gDirectory->Get("histo");
	JetResponsePhotonHemisphere[i][j][k]->hAlpha -> SetDirectory(0);
	delete file; 
	 
	TotFilename = RootPath + "hAlpha_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetImbalanceJetHemisphere[i][j][k]->hAlpha =  (TH1D*) gDirectory->Get("histo");
	JetImbalanceJetHemisphere[i][j][k]->hAlpha -> SetDirectory(0);
	delete file; 
	
	TotFilename = RootPath + "hAlpha_imbalance_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root"; 
	file = new TFile(TotFilename);      
	JetImbalancePhotonHemisphere[i][j][k]->hAlpha =  (TH1D*) gDirectory->Get("histo");
	JetImbalancePhotonHemisphere[i][j][k]->hAlpha -> SetDirectory(0);
	delete file;              
      
	TotFilename = RootPath + "hAlpha_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + FilenameIntrinsic;
	file = new TFile(TotFilename);      
	JetIntrinsic[i][j][k]->hAlpha =  (TH1D*) gDirectory->Get("histo");
	JetIntrinsic[i][j][k]->hAlpha -> SetDirectory(0);
	delete file; 
		
      }
    }
  }
  
  // Remove old object from root_files folder and pdf folder
  TString command;
  command = ".! rm -r " + RootPath + "jet_energy_*";
  gROOT->ProcessLine(command);
  command = ".! rm -r " + RootPath + "Resolution_*";
  gROOT->ProcessLine(command);
  command = ".! rm -r " + RootPath + "Scale_*";
  gROOT->ProcessLine(command);
  command = ".! rm -r " + PDFPath + "jet_energy_*";
  gROOT->ProcessLine(command); 
  command = ".! rm -r " + PDFPath + "Resolution_*";
  gROOT->ProcessLine(command);
  command = ".! rm -r " + PDFPath + "Scale_*";
  gROOT->ProcessLine(command);


  std::cout << "Calculation of Jet energy scale and resolution" << std::endl;

  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){
      
      JetScaleResAlpha[i][j]  = new CScaleResAlpha;
      JetIntrinsicAlpha[i][j] = new CScaleResAlpha;
      JetImbalanceAlpha[i][j] = new CScaleResAlpha;      
    }
    JetScaleRes[j]         = new CScaleRes;
    JetIntrinsicPt[j]      = new CScaleRes;
    JetImbalancePt[j]      = new CScaleRes;
    
  }
  
  // Fit a Gaussian to all Response Functions (conducted in function calculate()) 
  
  for(int i=0; i<pt_int; i++){
    for(int j=0; j<eta_int; j++){
      
      JetResponseJetHemisphere[i][j][0]          -> calculatePt(); 
      JetImbalanceJetHemisphere[i][j][0]         -> calculatePt();
      JetResponsePhotonHemisphere[i][j][0]       -> calculatePt(); 
      JetImbalancePhotonHemisphere[i][j][0]      -> calculatePt();
      JetIntrinsic[i][j][0]                      -> calculatePt(); 
      
      for(int k=0; k<alpha_int; k++){

	// Fit Gaussian Functions to Response Histogram  
	if(JetIntrinsic[i][j][k] ->hResponse -> GetEntries() >= 100 /*JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries()>100  && JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries()>100*/){
	  
	  JetResponseJetHemisphere[i][j][k]       -> calculate(); 
	  JetResponsePhotonHemisphere[i][j][k]    -> calculate();
	  JetImbalanceJetHemisphere[i][j][k]      -> calculate();
	  JetImbalancePhotonHemisphere[i][j][k]   -> calculate();
	  JetIntrinsic[i][j][k]                   -> calculate(); 	  
        }
	
      }
    }
  }
  
  TString etaRegion[eta_int];
  TString ptRegion[pt_int];
  TString filename;
  TString graphTitle;
  int idx,length;
  
  
  // Some description declarations
  for(int j=0;j<eta_int;j++){
    etaRegion[j].Form("%4.1f < |#eta| < %4.1f ",etaBin[j],etaBin[j+1]);  
  }

  for(int i=0;i<pt_int;i++){
    if(i!=pt_int-1)  ptRegion[i].Form("and %4.0f GeV < p_{T}^{#gamma} < %4.0f GeV",bd[i],bd[i+1]);
    else             ptRegion[i].Form("and %4.0f GeV < p_{T}^{#gamma} ",bd[i]);
  }


  // Write in other object all relevant numbers (only if they are different from zero)
  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){
  
      length = alpha_int;     
      idx = 0;      
      for(int k=0; k<alpha_int;k++){
	

        if(JetResponseJetHemisphere[i][j][k]->mean_array == 0 ){
	  length = length - 1;  
	  continue;
	}

	// Calculate the resulting Resolution from the two histograms in photon and jet hemissphere and write it into 

	// Calculate weighted mean  
	float weightPhoton = 1/TMath::Power(JetResponsePhotonHemisphere[i][j][k] -> sigma_error,2);
	//weightPhoton = 0 ;
	float weightJet = 1/TMath::Power(JetResponseJetHemisphere[i][j][k] -> sigma_error,2);	
	//weightJet = 0 ;
	float sigma = ((JetResponseJetHemisphere[i][j][k]->sigma_array)*weightJet + (JetResponsePhotonHemisphere[i][j][k]->sigma_array)*weightPhoton)/(weightJet+weightPhoton);
	float sigmaError = TMath::Sqrt(1./(weightJet+weightPhoton));
	weightPhoton = 1/TMath::Power(JetResponsePhotonHemisphere[i][j][k] -> alpha_error,2);
	//weightPhoton = 0 ;
	weightJet = 1/TMath::Power(JetResponseJetHemisphere[i][j][k] -> alpha_error,2);	
	//weightJet = 0 ;
	float alphaVal = ((JetResponseJetHemisphere[i][j][k]->alpha_array)*weightJet + (JetResponsePhotonHemisphere[i][j][k]->alpha_array)*weightPhoton)/(weightJet+weightPhoton);
	float alphaError = TMath::Sqrt(1./(weightJet+weightPhoton));
	
	JetScaleResAlpha[i][j]  -> alpha[idx]      = alphaVal;
        JetScaleResAlpha[i][j]  -> alphaError[idx] = alphaError;
	JetScaleResAlpha[i][j]  -> sigma[idx]      = sigma;         
        JetScaleResAlpha[i][j]  -> sigmaError[idx] = sigmaError;
	//JetScaleResAlpha[i][j]  -> alpha[idx]       = JetResponseJetHemisphere[i][j][k] -> alpha_array;         
        //JetScaleResAlpha[i][j]  -> alphaError[idx]  = JetResponseJetHemisphere[i][j][k] -> alpha_error;
	//JetScaleResAlpha[i][j]  -> sigma[idx]       = JetResponseJetHemisphere[i][j][k] -> sigma_array;         
        //JetScaleResAlpha[i][j]  -> sigmaError[idx]  = JetResponseJetHemisphere[i][j][k] -> sigma_error;

	//weightPhoton = 1/TMath::Power(JetResponsePhotonHemisphere[i][j][k] -> mean_error,2);
	//weightJet = 1/TMath::Power(JetResponseJetHemisphere[i][j][k] -> mean_error,2);		
	//float mean = ((JetResponseJetHemisphere[i][j][k]->mean_array)*weightJet + (JetResponsePhotonHemisphere[i][j][k]->mean_array)*weightPhoton)/(weightJet+weightPhoton);
	//float meanError = TMath::Sqrt(1./(weightJet+weightPhoton));
	//JetScaleResAlpha[i][j]  -> mean[idx]       = mean;  
        //JetScaleResAlpha[i][j]  -> meanError[idx]  = meanError;
	
	JetScaleResAlpha[i][j]  -> mean[idx]       = JetResponseJetHemisphere[i][j][k] -> mean_array;         
        JetScaleResAlpha[i][j]  -> meanError[idx]  = JetResponseJetHemisphere[i][j][k] -> mean_error;
	
        
	weightPhoton = 1/TMath::Power(JetImbalancePhotonHemisphere[i][j][k] -> sigma_error,2);
	//weightPhoton = 0 ;
	weightJet = 1/TMath::Power(JetImbalanceJetHemisphere[i][j][k] -> sigma_error,2);	
	//weightJet = 0 ;
	sigma = ((JetImbalanceJetHemisphere[i][j][k]->sigma_array)*weightJet + (JetImbalancePhotonHemisphere[i][j][k]->sigma_array)*weightPhoton)/(weightJet+weightPhoton);
	sigmaError = TMath::Sqrt(1./(weightJet+weightPhoton));
	weightPhoton = 1/TMath::Power(JetImbalancePhotonHemisphere[i][j][k] -> alpha_error,2);
	//weightPhoton = 0 ;
	weightJet    = 1/TMath::Power(JetImbalanceJetHemisphere[i][j][k] -> alpha_error,2);
	//weightJet = 0 ;	
	alphaVal = ((JetImbalanceJetHemisphere[i][j][k]->alpha_array)*weightJet + (JetImbalancePhotonHemisphere[i][j][k]->alpha_array)*weightPhoton)/(weightJet+weightPhoton);
	alphaError = TMath::Sqrt(1./(weightJet+weightPhoton));
	
	JetImbalanceAlpha[i][j] -> alpha[idx]      = alphaVal;  
        JetImbalanceAlpha[i][j] -> alphaError[idx] = alphaError; 
	JetImbalanceAlpha[i][j] -> sigma[idx]      = sigma;             
        JetImbalanceAlpha[i][j] -> sigmaError[idx] = sigmaError;   
	//JetImbalanceAlpha[i][j] -> alpha[idx]       = JetImbalanceJetHemisphere[i][j][k] -> alpha_array;         
	//JetImbalanceAlpha[i][j] -> alphaError[idx]  = JetImbalanceJetHemisphere[i][j][k] -> alpha_error;
	//JetImbalanceAlpha[i][j] -> sigma[idx]       = JetImbalanceJetHemisphere[i][j][k] -> sigma_array;         
        //JetImbalanceAlpha[i][j] -> sigmaError[idx]  = JetImbalanceJetHemisphere[i][j][k] -> sigma_error;

	//weightPhoton = 1/TMath::Power(JetImbalancePhotonHemisphere[i][j][k] -> mean_error,2);
	//weightJet = 1/TMath::Power(JetImbalanceJetHemisphere[i][j][k] -> mean_error,2);		
	//mean = ((JetImbalanceJetHemisphere[i][j][k]->mean_array)*weightJet + (JetImbalancePhotonHemisphere[i][j][k]->mean_array)*weightPhoton)/(weightJet+weightPhoton);
	//meanError = TMath::Sqrt(1./(weightJet+weightPhoton));
	//JetImbalanceAlpha[i][j] -> mean[idx]       = mean;         
	//JetImbalanceAlpha[i][j] -> meanError[idx]  = meanError;	
	JetImbalanceAlpha[i][j] -> mean[idx]       = JetImbalanceJetHemisphere[i][j][k] -> mean_array;         
        JetImbalanceAlpha[i][j] -> meanError[idx]  = JetImbalanceJetHemisphere[i][j][k] -> mean_error;


	JetIntrinsicAlpha[i][j] -> alpha[idx]      = JetIntrinsic[i][j][k] -> alpha_array;  
        JetIntrinsicAlpha[i][j] -> alphaError[idx] = JetIntrinsic[i][j][k] -> alpha_error;
        JetIntrinsicAlpha[i][j] -> mean[idx]       = JetIntrinsic[i][j][k] -> mean_array;         
        JetIntrinsicAlpha[i][j] -> meanError[idx]  = JetIntrinsic[i][j][k] -> mean_error;       
        JetIntrinsicAlpha[i][j] -> sigma[idx]      = JetIntrinsic[i][j][k] -> sigma_array;         
        JetIntrinsicAlpha[i][j] -> sigmaError[idx] = JetIntrinsic[i][j][k] -> sigma_error;
	    	
	idx = idx + 1;  
	
      }

      // Calculate weighted mean of Photon Pt
      float weightPhoton = 1/TMath::Power(JetResponsePhotonHemisphere[i][j][0] -> pT_error,2);
      //weightPhoton = 0 ;
      float weightJet = 1/TMath::Power(JetResponseJetHemisphere[i][j][0] -> pT_error,2);	
      float pT = ((JetResponseJetHemisphere[i][j][0]->pT_array)*weightJet + (JetResponsePhotonHemisphere[i][j][0]->pT_array)*weightPhoton)/(weightJet+weightPhoton);
      //weightJet = 0 ;
      float pTError = TMath::Sqrt(1./(weightJet+weightPhoton));

      /* OLD
      //JetScaleResAlpha[i][j]  -> pT         = JetResponseJetHemisphere[i][j][0] -> pT_array;
      //JetScaleResAlpha[i][j]  -> pTError    = JetResponseJetHemisphere[i][j][0] -> pT_error;
      //JetImbalanceAlpha[i][j] -> pT         = JetScaleResAlpha[i][j] -> pT;  
      //JetImbalanceAlpha[i][j] -> pTError    = JetScaleResAlpha[i][j] -> pTError;
      */

      JetScaleResAlpha[i][j]  -> pT        = pT;
      JetScaleResAlpha[i][j]  -> pTError   = pTError;
      JetImbalanceAlpha[i][j] -> pT        = pT;  
      JetImbalanceAlpha[i][j] -> pTError   = pTError;

      cout<<"pT = "<<pT<<endl;
      

      JetIntrinsicAlpha[i][j] -> pT         = JetIntrinsic[i][j][0] -> pT_array;  
      JetIntrinsicAlpha[i][j] -> pTError    = JetIntrinsic[i][j][0] -> pT_error;

      cout<<"JetIntrinsicAlpha[i][j] -> pT = "<<JetIntrinsicAlpha[i][j] -> pT<<endl;       
      

   	    
      cout<<"ptBoundLow = "<<bd[i]<<endl;
      cout<<"Eta bin = "<<j+1<<endl;

      JetIntrinsicAlpha[i][j] -> calculate(length, 2);
           
      JetImbalanceAlpha[i][j] -> calculate(length, 3);
   
      JetScaleResAlpha[i][j]  -> qprime      = JetImbalanceAlpha[i][j] -> qprime;
      JetScaleResAlpha[i][j]  -> qprimeError = JetImbalanceAlpha[i][j] -> qprimeError;
      JetScaleResAlpha[i][j]  -> mprime      = JetImbalanceAlpha[i][j] -> mprime;
      JetScaleResAlpha[i][j]  -> cprime      = JetIntrinsicAlpha[i][j] -> cprime;
      JetScaleResAlpha[i][j]  -> cprimeError = JetIntrinsicAlpha[i][j] -> cprimeError;
      JetScaleResAlpha[i][j]  -> q           = JetImbalanceAlpha[i][j] -> q;
      JetScaleResAlpha[i][j]  -> qError      = JetImbalanceAlpha[i][j] -> qError;

      JetScaleResAlpha[i][j]  -> calculate(length, 4);
            
      char legEntry[100];
      
      if(JetScaleResAlpha[i][j]->scale !=0){

	filename.Form("jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin",j+1,i+1);
        graphTitle = "Scale for " + etaRegion[j] + " " + ptRegion[i];

        plotTGraphErrors(JetScaleResAlpha[i][j] -> gJetScaleAlpha, filename, graphTitle, "alpha","JES",0,alphaBin[alpha_int],0.8,1.2,"",0);
        filename.Form("jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_intrinsic",j+1,i+1);
        graphTitle = "Scale for " + etaRegion[j] + " " + ptRegion[i] + " (intrinsic)";

	plotTGraphErrors(JetIntrinsicAlpha[i][j] -> gJetScaleAlpha, filename, graphTitle, "alpha","JES",0,alphaBin[alpha_int],0.8,1.2,"",0);
        filename.Form("jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_imbalance",j+1,i+1);
        graphTitle = "Scale for " + etaRegion[j] + " " + ptRegion[i] + " (imbalance)";

        plotTGraphErrors(JetImbalanceAlpha[i][j] -> gJetScaleAlpha, filename, graphTitle, "alpha","JES",0,alphaBin[alpha_int],0.8,1.2,"",0);  

        TFile *f;
        TF1* totalScale = new TF1("totalScale"," [0]*(1. - [1] - [2] *TMath::Power(x,2))",0,600); 
        totalScale -> SetParameters(JetIntrinsicAlpha[i][j]->c, JetImbalanceAlpha[i][j]->q,JetImbalanceAlpha[i][j]->m);
        filename.Form("jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_total",j+1,i+1);
        TString tot_filename = RootPath + filename + DataType + ".root";
        f = new TFile(tot_filename,"RECREATE");
        f -> WriteTObject(totalScale);
        f ->Close();
        delete f;
        
      
	filename.Form("jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin",j+1,i+1);
        graphTitle = "Resolution for " + etaRegion[j] + " " + ptRegion[i];
	sprintf(legEntry,"#splitline{Chi^2 = %4.2f}{ndof = %i}",JetScaleResAlpha[i][j] -> gJetResolutionAlpha->GetFunction("fResolutionAlpha")->GetChisquare(),JetScaleResAlpha[i][j] -> gJetResolutionAlpha->GetFunction("fResolutionAlpha") ->GetNDF());
        plotTGraphErrors(JetScaleResAlpha[i][j] -> gJetResolutionAlpha, filename, graphTitle, "alpha","JER",0,alphaBin[alpha_int],0.,0.4,legEntry,1);
        filename.Form("jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_intrinsic",j+1,i+1);
        graphTitle = "Resolution for " + etaRegion[j] + " " + ptRegion[i] + " (intrinsic)";
	sprintf(legEntry,"#splitline{Chi^2 = %4.2f}{ndof = %i}",JetIntrinsicAlpha[i][j]->gJetResolutionAlpha->GetFunction("fResolutionAlpha")->GetChisquare(),JetIntrinsicAlpha[i][j] -> gJetResolutionAlpha->GetFunction("fResolutionAlpha") ->GetNDF());
        plotTGraphErrors(JetIntrinsicAlpha[i][j] -> gJetResolutionAlpha, filename, graphTitle, "alpha","JER",0,alphaBin[alpha_int],0.,0.4,legEntry,1);
        filename.Form("jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_imbalance",j+1,i+1);
        graphTitle = "Resolution for " + etaRegion[j] + " " + ptRegion[i] + " (imbalance)";
	sprintf(legEntry,"#splitline{Chi^2 = %4.2f}{ndof = %i}",JetImbalanceAlpha[i][j]->gJetResolutionAlpha->GetFunction("fResolutionAlpha")->GetChisquare(),JetImbalanceAlpha[i][j]->gJetResolutionAlpha->GetFunction("fResolutionAlpha")->GetNDF());
        plotTGraphErrors(JetImbalanceAlpha[i][j] -> gJetResolutionAlpha, filename, graphTitle, "alpha","JER",0,alphaBin[alpha_int],0.,0.4,legEntry,1); 

        TF1* totalResolution = new TF1("totalResolution","TMath::Sqrt(TMath::Power([0],2) + TMath::Power([1],2) + 2*[1]*[2]*x + TMath::Power(([2]*x),2) )",0,600);
	totalResolution -> SetParameters(JetIntrinsicAlpha[i][j]->cprime, JetImbalanceAlpha[i][j]->qprime,JetImbalanceAlpha[i][j]->mprime);
        filename.Form("jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_total",j+1,i+1);  
        tot_filename = RootPath + filename + DataType + ".root";
        f = new TFile(tot_filename,"RECREATE");
        f -> WriteTObject(totalResolution);
        f ->Close();
        delete f; 
      }

      
    }
  }
   
  // DRAW LAST PLOTS
  char legEntry[100];
  double q[pt_int],qprime[pt_int],qError[pt_int],qprimeError[pt_int]; 
  double pT[pt_int],pTError[pt_int];
  
  
  for(int j=0; j<eta_int; j++){

    length = pt_int;
    idx = 0;
      

    for(int i=0; i<pt_int; i++){

      q[i]           = JetImbalanceAlpha[i][j] -> q;
      qError[i]      = JetImbalanceAlpha[i][j] -> qError;
      qprime[i]      = JetImbalanceAlpha[i][j] -> qprime;
      qprimeError[i] = JetImbalanceAlpha[i][j] -> qprimeError;
      pT[i]          = JetScaleResAlpha[i][j]  -> pT;
      pTError[i]     = JetScaleResAlpha[i][j]  -> pTError;

     
      if(JetScaleResAlpha[i][j] -> scale == 0 ){
        length = length - 1;       
        continue;
      }
     
      JetScaleRes[j] -> pT[idx]              = JetScaleResAlpha[i][j] -> pT;
      JetScaleRes[j] -> pTError[idx]         = JetScaleResAlpha[i][j] -> pTError;
      JetImbalancePt[j] -> pT[idx]           = JetScaleRes[j] -> pT[idx];
      JetImbalancePt[j] -> pTError[idx]      = JetScaleRes[j] -> pTError[idx];
      
      JetIntrinsicPt[j] -> pT[idx]           = JetIntrinsicAlpha[i][j] -> pT;
      JetIntrinsicPt[j] -> pTError[idx]      = JetIntrinsicAlpha[i][j] -> pTError;
      
     
      JetScaleRes[j] -> scale[idx]           = JetScaleResAlpha[i][j]  -> scale;
      JetScaleRes[j] -> scaleError[idx]      = JetScaleResAlpha[i][j]  -> scaleError;
      JetScaleRes[j] -> resolution[idx]      = JetScaleResAlpha[i][j]  -> resolution;
      JetScaleRes[j] -> resolutionError[idx] = JetScaleResAlpha[i][j]  -> resolutionError;
      JetScaleRes[j] -> q[idx]               = JetImbalanceAlpha[i][j] -> q;
      JetScaleRes[j] -> qError[idx]          = JetImbalanceAlpha[i][j] -> qError;
      JetScaleRes[j] -> qprime[idx]          = JetImbalanceAlpha[i][j] -> qprime;
      JetScaleRes[j] -> qprimeError[idx]     = JetImbalanceAlpha[i][j] -> qprimeError;
      

      JetIntrinsicPt[j] -> scale[idx]           = JetIntrinsicAlpha[i][j] -> scale;
      JetIntrinsicPt[j] -> scaleError[idx]      = JetIntrinsicAlpha[i][j] -> scaleError;
      JetIntrinsicPt[j] -> resolution[idx]      = JetIntrinsicAlpha[i][j] -> resolution;
      JetIntrinsicPt[j] -> resolutionError[idx] = JetIntrinsicAlpha[i][j] -> resolutionError;
      
      JetImbalancePt[j] -> scale[idx]           = JetImbalanceAlpha[i][j] -> scale;
      JetImbalancePt[j] -> scaleError[idx]      = JetImbalanceAlpha[i][j] -> scaleError;
      JetImbalancePt[j] -> resolution[idx]      = JetImbalanceAlpha[i][j] -> resolution;
      JetImbalancePt[j] -> resolutionError[idx] = JetImbalanceAlpha[i][j] -> resolutionError;  

      idx = idx +1;
      
    }
      
     
     
    TGraphErrors *gq          = new TGraphErrors(pt_int, pT , q, pTError, qError);
    TGraphErrors *gqprime     = new TGraphErrors(pt_int, pT , qprime, pTError, qprimeError);  

    // Calculate Fit variable for Resolution and Scale
    
    JetScaleRes[j]    -> calculate(length);
    JetIntrinsicPt[j] -> calculate(length);
    JetImbalancePt[j] -> calculate(length);
    
  
   
     
    // "Scale" Plots
    TString filenameSubsample;
    filename.Form("Scale_for_%i_eta_bin",j+1);
    filename += filenameSubsample;
    graphTitle = "Jet Energy Scale for " + etaRegion[j];
    sprintf(legEntry,"scale = %4.3f #pm %4.3f ",JetScaleRes[j] -> scalenumber,JetScaleRes[j] ->scalenumberError);
    plotTGraphErrors(JetScaleRes[j] -> gScale, filename, graphTitle, "p_{T}^{#gamma}","JES",0,600,0.7,1.1,legEntry,1); 

    filename.Form("Scale_for_%i_eta_bin_intrinsic",j+1);
    filename += filenameSubsample;
    graphTitle = "Jet Energy Scale for " + etaRegion[j] + " (intrinsic)";
    sprintf(legEntry,"scale = %4.3f #pm %4.3f ",JetIntrinsicPt[j] -> scalenumber,JetIntrinsicPt[j] ->scalenumberError);
    plotTGraphErrors(JetIntrinsicPt[j] -> gScale, filename, graphTitle, "p_{T}^{#gamma}","JES",0,600,0.7,1.1,legEntry,1); 

    filename.Form("Scale_for_%i_eta_bin_imbalance",j+1);
    filename += filenameSubsample;
    graphTitle = "Jet Energy Scale for " + etaRegion[j] + " (imbalance)";
    sprintf(legEntry,"scale = %4.3f #pm %4.3f ",JetImbalancePt[j] -> scalenumber,JetImbalancePt[j] ->scalenumberError);
    plotTGraphErrors(JetImbalancePt[j] -> gScale, filename, graphTitle, "p_{T}^{#gamma}","JES",0,600,0.7,1.1,legEntry,1); 

    filename.Form("Residual_imbalance_Scale_for_%i_eta_bin",j+1);
    filename += filenameSubsample;
    graphTitle = "Residual Imbalance for " + etaRegion[j] + " (Scale)";
    plotTGraphErrors(JetScaleRes[j] -> gq, filename, graphTitle, "p_{T}^{#gamma}","Residual Imbalance",0,600,-0.04,0.1,"a",0);

    filename.Form("Residual_imbalance_Scale_for_%i_eta_bin_all_pt_bins",j+1);
    filename += filenameSubsample;
    graphTitle = "Residual Imbalance for " + etaRegion[j] + " (Scale)";
    plotTGraphErrors(gq, filename, graphTitle, "p_{T}^{#gamma}","Residual Imbalance",0,600,-0.04,0.1,"a",0);
    
 

    // Resolution plots
    filename.Form("Resolution_for_%i_eta_bin",j+1);
    filename += filenameSubsample;
    graphTitle = "Jet Energy Resolution for " + etaRegion[j];
    sprintf(legEntry,"#splitline{N = %4.3f #pm %4.3f, S = %4.3f #pm %4.3f}{m = %4.3f #pm %4.3f, C = %4.3f #pm %4.3f}",JetScaleRes[j]->N,JetScaleRes[j]->NErr,JetScaleRes[j]->S,JetScaleRes[j]->SErr,JetScaleRes[j]->m,JetScaleRes[j]->mErr,JetScaleRes[j]->C,JetScaleRes[j]->CErr);
    plotTGraphErrors(JetScaleRes[j]->gResolution, filename, graphTitle, "p_{T}^{#gamma}","JER",0,600,-0.1,0.3,legEntry,1); 

    filename.Form("Resolution_for_%i_eta_bin_intrinsic",j+1);
    filename += filenameSubsample;
    graphTitle = "Jet Energy Resolution for " + etaRegion[j] + "(intrinsic)";
    sprintf(legEntry,"#splitline{N = %4.3f #pm %4.3f, S = %4.3f #pm %4.3f}{m = %4.3f #pm %4.3f, C = %4.3f #pm %4.3f}",JetIntrinsicPt[j]->N,JetIntrinsicPt[j]->NErr,JetIntrinsicPt[j]->S,JetIntrinsicPt[j]->SErr,JetIntrinsicPt[j]->m,JetIntrinsicPt[j]->mErr,JetIntrinsicPt[j]->C,JetIntrinsicPt[j]->CErr);
    plotTGraphErrors(JetIntrinsicPt[j]->gResolution, filename, graphTitle, "p_{T}^{#gamma}","JER",0,600,-0.1,0.3,legEntry,1); 

    filename.Form("Resolution_for_%i_eta_bin_imbalance",j+1);
    filename += filenameSubsample;
    graphTitle = "Jet Energy Resolution for " + etaRegion[j] + "(imbalance)";
    plotTGraphErrors(JetImbalancePt[j]->gResolution, filename, graphTitle, "p_{T}^{#gamma}","JER",0,600,-0.1,0.3,"a",0); 

    filename.Form("Residual_imbalance_Resolution_for_%i_eta_bin",j+1);
    filename += filenameSubsample;
    graphTitle = "Residual Imbalance for " + etaRegion[j]+ " (Resolution)";
    plotTGraphErrors(JetScaleRes[j] -> gqprime, filename, graphTitle, "p_{T}^{#gamma}","Residual Imbalance",0,600,-0.04,0.1,"a",0);

    filename.Form("Residual_imbalance_Resolution_for_%i_eta_bin_all_pt_bins",j+1);
    filename += filenameSubsample;
    graphTitle = "Residual Imbalance for " + etaRegion[j]+ " (Resolution)";
    plotTGraphErrors(gqprime, filename, graphTitle, "p_{T}^{#gamma}","Residual Imbalance",0,600,-0.04,0.1,"a",0);

    
  } // End of Eta loop
 
}





