// ----------------------------------------------------------------------------------------------
//  Script to calculate JET Energy Response in Gamma Jet Events
//
//  For questions on ROOT see also the manual at
//  http://root.cern.ch/drupal/content/users-guide
// ----------------------------------------------------------------------------------------------

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
#include "TLatex.h"

const int jetType = 2;      // 1 = PF_L1Fast , 2 = PF_CHS Jets , 3 = Calo_L1FastJet 
const int date    = 2012;
const int type    = 1;      //1 = Data , 2 = MC

#include "CODE/myDeclarations.h"
#include "CODE/myClasses.h"
#include "CODE/myFunctions.h"
#include "CODE/myClasses.C"
#include "CODE/myFunctions.C"
#include "CODE/calcPath.C"
#include "CODE/bookHistos.C"
#include "CODE/readGammaJet.C"
#include "CODE/applyCuts.C"
#include "CODE/draw.C"

// --------- Function and Class declarations -----------------
// -----------------------------------------------------------
void calcScale();
void calcSample();


// Run main script
// This is the actual command to be typed
// used in the ROOT session.
// ---------------------------------------------
int runGammaJet(int nEvts, int step) {
  
  time_t start,end;
  double dif;
  time (&start);

  // -------------------------------------------------------------------------------------
  if(date == 2011 && numTrigger!=5){
    cout<<"Please correct the number of Triggers used (in CODE/myDeclarations.h)!"<<endl;  
    return 0;
  }
  if(date == 2012 && numTrigger!=8){
    cout<<"Please correct the number of Triggers used (in CODE/myDeclarations.h)!"<<endl;  
    return 0;
  }
  //-------- print out of Information ----------------------------------------------------
  cout<<endl;
  
  if(type == 1)          cout<<"Data!!"<<endl<<endl;
  else if(type == 2)     cout<<"MC!!"<<endl<<endl;
  
  if(date == 2012)       cout<<"date = 2012"<<endl<<endl;
  else if(date == 2011)  cout<<"date = 2011"<<endl<<endl;
  
  if(jetType == 1)       cout<<"jetType = PF_L1Fast"<<endl<<endl;
  else if(jetType == 2)  cout<<"jetType = PF_CHSJets"<<endl<<endl;
  
  cout<<endl;
  //--------------------------------------------------------------------------------------
  cout<<"calcPath() is executed!"<<endl<<endl;
  calcPath(step);
  cout<<"bookHistos() is executed!"<<endl<<endl;
  bookHistos();
  cout<<"readGammaJet() is executed!"<<endl<<endl;
  readGammaJet(nEvts);
  if(step == 1){
    cout<<"calcSample() is executed!"<<endl<<endl;
    calcSample();
  }
  if(step == 2){
    cout<<"calcScale() is executed!"<<endl<<endl;
    calcScale();
  }
  //cout<<"draw() is executed!"<<endl<<endl;
  //draw();
  //--------------------------------------------------------------------------------------

  time (&end);
  dif = difftime(end,start);
  cout<<"elapsed time:"<<dif/60<<" minutes"<<endl<<endl;

  return 0;
}


// Read GammaJet from file and fill histograms. 
// -------------------------------------------------------------------------------
void calcSample() {

  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){
      for(int k=0; k<alpha_int; k++){
	
	JetResponseJetHemisphere[i][j][k]          = new CResponse(1); 	
	JetResponsePhotonHemisphere[i][j][k]       = new CResponse(1); 	      
      }   
    }
  }
  //-----------------------------------------------------------------------------
  
  // Loop over nMax entries, read variables and fill histogram
  
  std::cout << "Processing events" << std::endl;
  int testCalc = 0;
  
  for(int n = 0; n < nMax; n++) {
    
    if( (n+1)%1000000 == 0 ) std::cout << "Event " << (n+1) << std::endl;
    
    // Get this event i.e. entry n and fill values into variables
    chain->GetEntry(n);
    
    //Just MC has PUweights; in data you have to set them to One
    PUWeight= 1.;
    weight  = 1.;  
    
    /*
      float photonEt =0;
      if(photonIsoEcal[0]> 4.0 + 0.012*photonEt){
      //cout<<"because Ecal thrown out."<<endl;
      
      continue;
      }
      if(photonIsoHcal[0]> 5.5 + 0.005*photonEt){
      //cout<<"because Hcal thrown out."<<endl;
      
      continue;
      }
      if(photonIsoTrk[0]> 4.0 + 0.002*photonEt){
      //cout<<"because Trk thrown out."<<endl;
      
      continue;
      }
    */
    
    
    //--------------------------------------------------------------------------
    // Sort Jets and apply Cuts
    bool testVar = applyCuts();
    if(!testVar){
      testCalc += 1;
      continue;
    }
    //--------------------------------------------------------------------------
    
    /*
      Number of Vertices for different pt bins
      for(int i=0; i<numTrigger-1; i++){
      if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1] && hltPhoton[i]==1){
      hVtxPtbinned[i] ->Fill(vtxN,weight*PUWeight);
      break;
      }
      else if(photonPt[0]>=bd[numTrigger-1] && hltPhoton[numTrigger-1]==1){
      hVtxPtbinned[numTrigger-1]->Fill(vtxN,weight*PUWeight); 
      break;
      }
      }  
      
      hPhotonIsoEcal -> Fill(photonIsoEcal[0],weight*PUWeight);
      hPhotonIsoHcal -> Fill(photonIsoHcal[0],weight*PUWeight);
      hPhotonIsoTrk  -> Fill(photonIsoTrk[0],weight*PUWeight);  
    */
    
    // Fill Response Functions    
    for(int k=0; k<alpha_int; k++){
      if(alpha >= alphaBin[k] && alpha < alphaBin[k+1]){ 
	
	for(int j=0; j<eta_int; j++){
	  
	  if(std::abs(jetEta[corrJets.idx(lead_jet)])>= etaBin[j] && std::abs(jetEta[corrJets.idx(lead_jet)]) < etaBin[j+1] ){
	    
	    for(int i=0; i<pt_int; i++){ 
              
              if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){

		
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
		  JetResponseJetHemisphere[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);
		  JetResponseJetHemisphere[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);
		}
		else{
		  JetResponsePhotonHemisphere[i][j][k] -> hResponse -> Fill(response,weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][0] -> hPt       -> Fill(photonPt[0],weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][k] -> hAlpha    -> Fill(alpha,weight*PUWeight);
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
    
    
    // SMALL RESPONSE DIAGNOSTIC::Plots just for small response region
    //if( corrJets.pt(lead_jet)/photonPt[0] <= 0.3 ){
    //hPhotonEta_low_resp_reg -> Fill(photonEta[0],weight);
    
    //if(!filestr.is_open()) cout<<"wrong"<<endl;
    //filestr<<runNum<<":"<<lumiNum<<":"<<eventNum<<endl;	   
    //if(!EventVariables.is_open()) cout<<"EventVariable.txt is not open"<<endl;
    //EventVariables<<photonPt[0]<<":"<<corrJets.pt(lead_jet)<<":"<<photonEta[0]<<":"<<jetEta[corrJets.idx(lead_jet)]<<":"<<photonPhi[0]<<":"<<jetPhi[corrJets.idx(lead_jet)]<<":"<<photonIsoEcal[0]<<":"<<photonIsoHcal[0]<<":"<<photonIsoTrk[0]<<endl;
    //}
    
    //if(( hltPhoton[5] ||  hltPhoton[6] || hltPhoton[7])  && corrJets.pt(lead_jet)/photonPt[0]<=1.0 && photonPt[0]>130.){
    //  hSmallJetPtResponse -> Fill(corrJets.pt(lead_jet),corrJets.pt(lead_jet)/photonPt[0],weight);
    //  hPhotonEta_high_pt_reg_Response -> Fill(photonEta[0],corrJets.pt(lead_jet)/photonPt[0],weight);      
    //}
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
    
    /*// Fill 2d Diagram Response L1 correctio
      hResponseL1 -> Fill(corrJets.pt(lead_jet)/photonPt[0],jetCorrL1[corrJets.idx(lead_jet)],weight);
      hResponseL2L3 -> Fill(corrJets.pt(lead_jet)/photonPt[0],jetCorrL2L3[corrJets.idx(lead_jet)],weight);
      
      //Fill 2d histogram Response agaist Photon Eta
      hRespPhotonEta -> Fill(corrJets.pt(lead_jet)/photonPt[0],photonEta[0],weight);
      
      // Fill 2d Diagram Response against JetEMF JetFHPD JetFRBX
      hResponseJetEMF -> Fill(response,jetEMF[corrJets.idx(lead_jet)],weight);
      hResponseJetFHPD -> Fill(response,jetFHPD[corrJets.idx(lead_jet)],weight);
      hResponseJetFRBX -> Fill(response,jetFRBX[corrJets.idx(lead_jet)],weight);
      
      // Fill 2d Diagram Response against PhotonEMF PhotonFHPD PhotonFRBX
      if(photonidx!=-1){
      hResponsePhotonEMF -> Fill(response,jetEMF[corrJets.idx(photonidx)],weight);
      hResponsePhotonFHPD -> Fill(response,jetFHPD[corrJets.idx(photonidx)],weight);
      hRespoxnsePhotonFRBX -> Fill(response,jetFRBX[corrJets.idx(photonidx)],weight);
      if(response<=1 && corrJets.pt(photonidx)/photonPt[0]<=2)    hResponsePhotonPtRatio -> Fill(response,corrJets.pt(photonidx)/photonPt[0],weight);
      }*/
    
    // Fill Histogram for first Photon pT 
    hPhoton1Pt->Fill(photonPt[0],weight);

    // Fill 2d histogram with first and second Photon of one event
    hPhoton12Pt->Fill(photonPt[0],photonPt[1],weight);
    
    nocut = nocut +1;

  } // End of loop over entries


  TFile *f;
  TString TotFilename;

  for(int i=0;i<pt_int;i++){
    for(int j=0; j<eta_int; j++){

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

      for(int k=0; k<alpha_int; k++){	

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

	TotFilename = RootPath + "response_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";         
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetResponseJetHemisphere[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;

	TotFilename = RootPath + "response_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";           
	f = new TFile(TotFilename,"RECREATE");
	f -> WriteTObject(JetResponsePhotonHemisphere[i][j][k]->hResponse,"histo");
        f->Close();
        delete f;

      }
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
  cout<<"cut11 = "<<cut11<<endl;\
  for(int i=0; i<numTrigger; i++) cout<<"cut12["<<i<<"] = "<<cut12[i]<<endl;
  cout<<"nocut = "<<nocut<<endl;
  cout<<"number of all events = "<<cut1+cut2+cut3+cut4+cut5+cut6+cut7+cut8+cut9+cut10+cut11+cut12[0]+cut12[1]+cut12[2]+cut12[3]+cut12[4]+cut12[5]+cut12[6]+cut12[7]+nocut<<endl<<endl;

  cout<<"testCalc = "<<testCalc<<endl;
  
  delete chain;
  filestr.close();   
  EventVariables.close();

 
  
}


// --------------------------------------------------------------------------------
// Now take the response functions and fit to the core region a gauss function 
// so that one gets the jet energy scale and jet energy resolution


// Calculate sigma and mean for all Response functions (also for different nVtx regions)

void calcScale(){


  std::cout << "Read Response Histograms!" << std::endl<< std::endl;

  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){
      for(int k=0; k<alpha_int; k++){

	JetResponseJetHemisphere[i][j][k]         = new CResponse(1); 		
	JetResponsePhotonHemisphere[i][j][k]      = new CResponse(1); 	
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

      for(int k=0; k<alpha_int; k++){
	
	TotFilename = RootPath + "response_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetResponseJetHemisphere[i][j][k]->hResponse =  (TH1D*) gDirectory->Get("histo");
	JetResponseJetHemisphere[i][j][k]->hResponse -> SetDirectory(0);
	delete file; 
				
	TotFilename = RootPath + "response_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetResponsePhotonHemisphere[i][j][k]->hResponse =  (TH1D*) gDirectory->Get("histo");
	JetResponsePhotonHemisphere[i][j][k]->hResponse -> SetDirectory(0);
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
  

  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){
      
      JetScaleResAlpha[i][j] = new CScaleResAlpha;
    }
    JetScaleRes[j] = new CScaleRes;
  }
  
  
  std::cout << "Calculation of Jet energy scale and resolution" << std::endl;
    
  // Fit a Gaussian to all Response Functions (conducted in function calculate()) 
  for(int i=0; i<pt_int; i++){
    for(int j=0; j<eta_int; j++){

      JetResponseJetHemisphere[i][j][0]          -> calculatePt();
      JetResponsePhotonHemisphere[i][j][0]       -> calculatePt();
      
      for(int k=0; k<alpha_int; k++){

	//cout<<"JetResponseJetHemisphere["<<i<<"]["<<j<<"]["<<k<<"]->hResponse->GetEntries() = "<<JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries()<<endl;
	//cout<<"JetResponsePhotonHemisphere["<<i<<"]["<<j<<"]["<<k<<"]->hResponse->GetEntries() = "<<JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries()<<endl;
	
	if(JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries()>90  && JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries()>90){
	  
	  JetResponseJetHemisphere[i][j][k]       -> calculate(); 
	  JetResponsePhotonHemisphere[i][j][k]    -> calculate();
	  
        }  
	
      }
    }
  }


 
  TString etaRegion[eta_int];
  TString ptRegion[pt_int];
  TString filename;
  TString graphTitle;
  int idx,length;
  
  // some description declarations
  for(int j=0;j<eta_int;j++){
    etaRegion[j].Form("%4.1f < |#eta| < %4.1f",etaBin[j],etaBin[j+1]);  
  }

  for(int i=0;i<pt_int;i++){
    if(i!=pt_int-1)  ptRegion[i].Form("and %4.0f GeV < p_{T}^{#gamma} < %4.0f GeV",bd[i],bd[i+1]);
    else             ptRegion[i].Form("and %4.0f GeV < p_{T}^{#gamma} ",bd[i]);
  }
    
  // Read Residual Imbalances
  TGraphErrors* ResImbScale[eta_int] = {0};
  TGraphErrors* ResImbRes[eta_int]   = {0};
  TString ResImbFilename;
  for(int j=0; j<eta_int; j++){
    if(date==2012){
      
      if(jetType == 1)      ResImbFilename  = "plots_2012/PF_L1FastJet/mc/root_files/Residual_imbalance_Scale_for_";
      else if(jetType == 2) ResImbFilename  = "plots_2012/PF_L1CHS/mc/root_files/Residual_imbalance_Scale_for_";
      
    }
    else if(date==2011) ResImbFilename  = "plots_2011/PF_L1FastJet/mc/root_files/Residual_imbalance_Scale_for_";
    ResImbFilename += j+1;
    if(jetType == 1)      ResImbFilename += "_eta_bin_all_pt_bins_PF_mc.root";
    else if(jetType == 2) ResImbFilename += "_eta_bin_all_pt_bins_PFCHS_mc.root";
    ResImbScale[j]  = readTGraphErrors(ResImbFilename,"Graph;1","Graph;1");  
  }
  for(int j=0; j<eta_int; j++){
    if(date==2012){

      if(jetType == 1)       ResImbFilename  = "plots_2012/PF_L1FastJet/mc/root_files/Residual_imbalance_Resolution_for_";
      else if(jetType == 2)  ResImbFilename  = "plots_2012/PF_L1CHS/mc/root_files/Residual_imbalance_Resolution_for_";
    }
    else if(date==2011) ResImbFilename  = "plots_2011/PF_L1FastJet/mc/root_files/Residual_imbalance_Resolution_for_";
    ResImbFilename += j+1;
    if(jetType == 1)      ResImbFilename += "_eta_bin_all_pt_bins_PF_mc.root";
    else if(jetType == 2) ResImbFilename += "_eta_bin_all_pt_bins_PFCHS_mc.root";
    ResImbRes[j]    = readTGraphErrors(ResImbFilename,"Graph;1","Graph;1");  
  }
  
  cout<<"Residual imbalances are read from following file: "<<ResImbFilename<<endl<<endl;
  
  double_t *YResImbRes,*YResImbScale,*YResImbResError; 
  
  // write in other object all relevant numbers (only if they are different from zero)
  for(int j=0; j<eta_int; j++){

    
    YResImbScale    = ResImbScale[j] -> GetY();        
    YResImbRes      = ResImbRes[j]   -> GetY();
    YResImbResError = ResImbRes[j]   -> GetEY();
    
    
    
    for(int i=0; i<pt_int; i++){
  
      length = alpha_int;     
      idx = 0;
      
      for(int k=0; k<alpha_int;k++){
        if(JetResponseJetHemisphere[i][j][k]->mean_array == 0 ){

	  length = length - 1;       
          continue;

	}

 
	// Calculate weighted mean  
	float weightPhoton = 1/TMath::Power(JetResponsePhotonHemisphere[i][j][k] -> sigma_error,2);
	float weightJet = 1/TMath::Power(JetResponseJetHemisphere[i][j][k] -> sigma_error,2);	
	float sigma = ((JetResponseJetHemisphere[i][j][k]->sigma_array)*weightJet + (JetResponsePhotonHemisphere[i][j][k]->sigma_array)*weightPhoton)/(weightJet+weightPhoton);
	float sigmaError = TMath::Sqrt(1./(weightJet+weightPhoton));
	weightPhoton = 1/TMath::Power(JetResponsePhotonHemisphere[i][j][k] -> alpha_error,2);
	weightJet = 1/TMath::Power(JetResponseJetHemisphere[i][j][k] -> alpha_error,2);	
	float alphaVal = ((JetResponseJetHemisphere[i][j][k]->alpha_array)*weightJet + (JetResponsePhotonHemisphere[i][j][k]->alpha_array)*weightPhoton)/(weightJet+weightPhoton);
	float alphaError = TMath::Sqrt(1./(weightJet+weightPhoton));

	JetScaleResAlpha[i][j]  -> alpha[idx]      = alphaVal;
        JetScaleResAlpha[i][j]  -> alphaError[idx] = alphaError;
	JetScaleResAlpha[i][j]  -> sigma[idx]      = sigma;         
        JetScaleResAlpha[i][j]  -> sigmaError[idx] = sigmaError;
	JetScaleResAlpha[i][j]  -> mean[idx]       = JetResponseJetHemisphere[i][j][k] -> mean_array;         
        JetScaleResAlpha[i][j]  -> meanError[idx]  = JetResponseJetHemisphere[i][j][k] -> mean_error;

	/* OLD
	JetScaleResAlpha[i][j] -> alpha[idx]      = JetResponse[i][j][k] -> alpha_array;
        JetScaleResAlpha[i][j] -> alphaError[idx] = JetResponse[i][j][k] -> alpha_error;
        JetScaleResAlpha[i][j] -> mean[idx]       = JetResponse[i][j][k] -> mean_array;         
        JetScaleResAlpha[i][j] -> meanError[idx]  = JetResponse[i][j][k] -> mean_error;
        JetScaleResAlpha[i][j] -> sigma[idx]      = JetResponse[i][j][k] -> sigma_array;         
        JetScaleResAlpha[i][j] -> sigmaError[idx] = JetResponse[i][j][k] -> sigma_error;
	*/
       
        idx = idx + 1;  

      }
      
      // Calculate weighted mean of Photon Pt

      float weightPhoton = 1/TMath::Power(JetResponsePhotonHemisphere[i][j][0] -> pT_error,2);
      float weightJet = 1/TMath::Power(JetResponseJetHemisphere[i][j][0] -> pT_error,2);	
      float pT = ((JetResponseJetHemisphere[i][j][0]->pT_array)*weightJet + (JetResponsePhotonHemisphere[i][j][0]->pT_array)*weightPhoton)/(weightJet+weightPhoton);
      float pTError = TMath::Sqrt(1./(weightJet+weightPhoton));

      JetScaleResAlpha[i][j] -> pT         = pT;
      JetScaleResAlpha[i][j] -> pTError    = pTError;

      /* OLD
      JetScaleResAlpha[i][j] -> pT         = JetResponse[i][j][0] -> pT_array;
      JetScaleResAlpha[i][j] -> pTError    = JetResponse[i][j][0] -> pT_error;
      */

      // fix q and qprime on residual imbalance value (for systematic Uncertainies take more/less 0.5*q)
      JetScaleResAlpha[i][j] -> q           = YResImbScale[i];   //-0.5*YResImbScale[i];
      JetScaleResAlpha[i][j] -> qprime      = YResImbRes[i]  ;   //-0.5*YResImbRes[i];
      JetScaleResAlpha[i][j] -> qprimeError = YResImbResError[i]  ;
      

      //JetScaleResAlpha[i][j] -> q      = 0;
      //JetScaleResAlpha[i][j] -> qprime = 0; 

      cout<<"ptBoundLow = "<<bd[i]<<endl;
      cout<<"Eta bin = "<<j+1<<endl;
      JetScaleResAlpha[i][j] -> calculate(length,4);

      
      
      if(JetScaleResAlpha[i][j]->scale !=0){
	filename.Form("jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin",j+1,i+1);
	graphTitle = "Jet Energy Scale for " + etaRegion[j] + ptRegion[i];
	plotTGraphErrors(JetScaleResAlpha[i][j] -> gJetScaleAlpha, filename, graphTitle, "alpha","JES",0,20,0.7,1.2,"",0);  

	filename.Form("jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin",j+1,i+1);  
	graphTitle = "Jet Energy Resolution for " + etaRegion[j] + ptRegion[i];
	plotTGraphErrors(JetScaleResAlpha[i][j] -> gJetResolutionAlpha, filename, graphTitle, "alpha","JER",0,20,0.,0.4,"",0); 
      }
    }
  }
   

  // DRAW LAST PLOTS
  char legEntry[100];

  for(int j=0; j<eta_int; j++){

    length = pt_int;
    idx = 0;
    cout<<endl<<"j = "<<j<<endl;
    for(int i=0; i<pt_int; i++){      
      
      if(JetScaleResAlpha[i][j] -> scale == 0 || JetScaleResAlpha[i][j] -> resolution != JetScaleResAlpha[i][j] -> resolution ){
        length = length - 1;  
	//cout<<"JetScaleResAlpha[i][j] -> scale = "<<JetScaleResAlpha[i][j] -> scale<<endl; 
	//cout<<"JetScaleResAlpha[i][j] -> resolution = "<<JetScaleResAlpha[i][j] -> resolution<<endl;    
        continue;
      }
      
      JetScaleRes[j] -> pT[idx]              = JetScaleResAlpha[i][j] -> pT;
      JetScaleRes[j] -> pTError[idx]         = JetScaleResAlpha[i][j] -> pTError;
      JetScaleRes[j] -> scale[idx]           = JetScaleResAlpha[i][j] -> scale;
      JetScaleRes[j] -> scaleError[idx]      = JetScaleResAlpha[i][j] -> scaleError;
      JetScaleRes[j] -> resolution[idx]      = JetScaleResAlpha[i][j] -> resolution;
      JetScaleRes[j] -> resolutionError[idx] = JetScaleResAlpha[i][j] -> resolutionError;
      
      idx = idx +1;

    }

    JetScaleRes[j] ->calculate(length);
    
     
    filename.Form("Scale_for_%i_eta_bin",j+1);
    cout<<"filename = "<<filename<<endl;
    graphTitle = "Jet Energy Scale for " + etaRegion[j];
    sprintf(legEntry,"scale = %4.3f #pm %4.3f ",JetScaleRes[j] -> scalenumber,JetScaleRes[j] ->scalenumberError);
    plotTGraphErrors(JetScaleRes[j] -> gScale, filename, graphTitle, "p_{T}^{#gamma}","JES",0,600,0.7,1.1,legEntry,1); 
     
    filename.Form("Resolution_for_%i_eta_bin",j+1);
    cout<<"filename = "<<filename<<endl;
    graphTitle = "Jet Energy Resolution for " + etaRegion[j];
    //sprintf(legEntry,"#splitline{#splitline{N = %4.3f #pm %4.3f}{S = %4.3f #pm %4.3f}}{#splitline{C = %4.3f #pm %4.3f}{m = %4.3f #pm %4.3f}} ",JetScaleRes[j] -> N,JetScaleRes[j] ->NErr,JetScaleRes[j] -> S,JetScaleRes[j] ->SErr,JetScaleRes[j] -> C,JetScaleRes[j] ->CErr,JetScaleRes[j] -> m,JetScaleRes[j] ->mErr);
    plotTGraphErrors(JetScaleRes[j]->gResolution, filename, graphTitle, "p_{T}^{#gamma}","JER",0,600,-0.1,0.3,"a",0); 
     
    

  }
   

  

  
}

