#include "myDeclarations.h"
#include "myFunctions.h"
#include "myClasses.h"

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>

//---------------------------------------------------------------
// Draw histograms
// --------------------------------------------------------------
void draw() {
  std::cout << "Plotting histograms" << std::endl;
  
  TString texte;
  Int_t oldLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;
  
  //---------- 1d - Histograms----------------------------------
  plotTH1(hDeltaPhi1st2ndJet,"hDeltaPhi1st2ndJet",0);
  plotTH2(hDeltaPhi1st2ndJetDeltaPt,"hDeltaPhi1st2ndJetDeltaPt_2d");

  plotTH1(hAlpha,"hAlpha",0);
  plotTH1(hDeltaPhi,"hDeltaPhi",0);

  //plotTH1(hChiSquareIntrinsic,"hChiSquareIntrinsic",0);
  // plotTH1(hChiSquareImbalance,"hChiSquareImbalance",0);
  //plotTH1(hChiSquareResolution,"hChiSquareResolution",0);
  
  /*// alpha plot (pTJet2/pTPhoton1)
  plotTH1(hRatioPt,"PtRatio_Jet2_Photon1",0);

  // Number of Jets in one Event
  plotTH1(hJetN,"JetN",0); 

  // Eta of first Jet in events with vtxN>1
  // plotTH1(hJetEta_1,"JetEta_with_vtxN_ge_2",0);

  

  // Number of Photons per event
  plotTH1(hPhotonN,"PhotonN",0);

  // Number of tight Photons per event
  plotTH1(hPhotonNtight,"PhotonNtight",0);
 
  // Eta of first Photon
  plotTH1(hPhotonEta,"PhotonEta",0);

  // Eta of first Photon in high Pt region
  plotTH1(hPhotonEta_low_resp_reg,"PhotonEta_low_resp_reg",0);
  */

  // Eta of first Jet
  // plotTH1(hJetEta,"JetEta",0);

  // Phi of Photon
  //plotTH1(hPhotonPhi,"photon_phi_distribution",0);
  
  // First Photon Pt for all Photons
  plotTH1(hPhoton1Pt,"Photon1Pt",1);

  

  // number of pile-ups in a certain eta region
  /*
  plotTH1(hVtxN,"VtxN",0);

  plotTH1(hPhotonIsoEcal,"hPhotonIsoEcal",0);
  plotTH1(hPhotonIsoHcal,"hPhotonIsoHcal",0);
  plotTH1(hPhotonIsoTrk,"hPhotonIsoTrk",0);
  */
  
  for(int i=0; i<numTrigger; i++){
    /*
    texte.Form("VtxN%i",i);
    plotTH1(hVtxPtbinned[i],texte,0);
    
    texte.Form("PhotonPt%i",i);
    plotTH1(hPhotonPt[i],texte,0);

    texte.Form("Rho%i",i);
    plotTH1(hRho[i],texte,0);

    texte.Form("TriggerEffBefore%i",i);
    plotTH1(hTriggerEffBefore[i],texte,0);

    texte.Form("TriggerEffAfter%i",i);
    plotTH1(hTriggerEffAfter[i],texte,0);

    texte.Form("PUgenMC%i",i);
    plotTH1(hPUgenMC[i],texte,0);
    */

  }
  

  /*
    for(int i=0; i<40; i++){
    texte.Form("RhoVtxBinned%i",i);
    plotTH1(hRhoVtxBinned[i],texte,0);
    }
  */

  // number of vertices in different pt bins
  /* 
  plotTH1(hPhotonPtVtxBinned[0],"PhotonPtNVtxBinned0",0);
  plotTH1(hPhotonPtVtxBinned[1],"PhotonPtNVtxBinned1",0);
  plotTH1(hPhotonPtVtxBinned[2],"PhotonPtNVtxBinned2",0);
  plotTH1(hPhotonPtVtxBinned[3],"PhotonPtNVtxBinned3",0);
  plotTH1(hPhotonPtVtxBinned[4],"PhotonPtNVtxBinned4",0);
  if(date == 2012){
    plotTH1(hPhotonPtVtxBinned[5],"PhotonPtNVtxBinned5",0);
    plotTH1(hPhotonPtVtxBinned[6],"PhotonPtNVtxBinned6",0);
  }
  */

  //plotTH1(hWeight,"Weight",1);
  //plotTH1(hWeightWeight,"Weight_weighted",0);
 
  // pT of first tight photon 
  // plotTH1(hPhotonT1Pt[i],"tightPhoton1Pt",1); 

  // Photon pT spectrum in various ptPhoton and Eta Bins
  /*for(int j=0;j<pt_int;j++){
  for(int i=0;i<2;i++){
    texte= "PhotonPtSpectrum_in_";
    texte += j+1;
    texte += "_Pt_bin_and_";
    texte += i+1;
    texte += "_eta_bin";    
    plotTH1(haux[j][i],texte,1);
  }
  }*/
 
    
  /*
  for(int k=0; k<4; k++){
    for(int i=0; i<pt_int; i++){
      for(int j=0; j<eta_int; j++){
	texte.Form("PUsysY%iPU%iPT%iEta",k,i,j);
	plotTH1(hPUsysY[k][i][j],texte,0);
	
      }
    }
  }
  */
  /*
  gStyle->SetOptStat(0);   
  TCanvas* c = 0;  
  TCanvas* cc = 0;  
  for(int k=0; k<alpha_int; k++){
  
    c = new TCanvas("c","c",0,0,500,500);
    TLegend *legend;
    TString legname;
 
    
    c -> cd(); 
    h2ndJetPt1stJetHemisphere[k] -> Draw();
    legname.Form("2nd Jet Pt spectrum (%4.1f GeV < P_{T}^{#gamma} and %4.1f < #alpha < %4.1f and %4.1f < #eta < %4.1f)",bd[pt_int-1],alphaBin[k],alphaBin[k+1],etaBin[0],etaBin[0]);    
    h2ndJetPt1stJetHemisphere[k] -> SetTitle(legname);
    h2ndJetPtPhotonHemisphere[k] -> Scale((h2ndJetPt1stJetHemisphere[k]->Integral())/(h2ndJetPtPhotonHemisphere[k]->Integral()));
    h2ndJetPtPhotonHemisphere[k] -> Draw("same");
    h2ndJetPtPhotonHemisphere[k] -> SetLineColor(3);
    legend  = new TLegend(0.5,0.7,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.030);
    legend -> AddEntry(h2ndJetPt1stJetHemisphere[k],"Jet Hemisphere","l");
    legend -> AddEntry(h2ndJetPtPhotonHemisphere[k],"Photon Hemisphere","l");
    legend -> Draw("same");
    c -> SetLogy();
    texte.Form("h2ndJetPtPhotonJetHemisphere_%i_alphaBin.pdf",k);
    TString TotFilename = PDFPath + texte;
    c -> SaveAs(TotFilename);
    delete legend;
  
    cc = new TCanvas("cc","cc",0,0,500,500);
    cc -> cd(); 
    h2ndJetEta1stJetHemisphere[k] -> Draw();
    h2ndJetEtaPhotonHemisphere[k] -> Scale((h2ndJetEta1stJetHemisphere[k]->Integral())/(h2ndJetEtaPhotonHemisphere[k]->Integral()));
    h2ndJetEtaPhotonHemisphere[k] -> Draw("same");
    h2ndJetEtaPhotonHemisphere[k] -> SetLineColor(3);  
    legend  = new TLegend(0.5,0.1,0.9,0.3);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.030);
    legend -> AddEntry(h2ndJetPt1stJetHemisphere[k],"Jet Hemisphere","l");
    legend -> AddEntry(h2ndJetPtPhotonHemisphere[k],"Photon Hemisphere","l");
    legend -> Draw("same");
    cc -> SetLogy();
    texte.Form("h2ndJetEtaPhotonJetHemisphere_%i_alphaBin.pdf",k);
    TotFilename = PDFPath + texte;
    cc -> SaveAs(TotFilename);

    delete cc;
    delete c;

  } 
  gStyle->SetOptStat("nemr");   
  */
 //---------- 2d - Histograms----------------------------------
  /*
  plotTH2(hJet1PtPhoton1Pt,"Jet1PtPhoton1Pt_2d");
  plotTH2(hImbalanceJetN,"hImbalanceJetN_2d");
  plotTH2(hImbalancePhoton1Pt,"hImbalancePhoton1Pt_2d");
  plotTH2(hImbalanceVtx,"hImbalanceVtx_2d");
  plotTH2(hImbalanceJetPt,"hImbalanceJetPt_2d");
  plotTH2(hImbalanceJetEta,"hImbalanceJetEta_2d");
  plotTH2(hImbalancePhotonEta,"hImbalancePhotonEta_2d");
  plotTH2(hImbalanceWeights,"hImbalanceWeights_2d");
  plotTH2(hImbalanceJet2Pt,"hImbalanceJet2Pt_2d");
  plotTH2(hImbalanceGenJet2Pt,"hImbalanceGenJet2Pt_2d");
  plotTH2(hJet2PtGenJet2Pt,"hJet2PtGenJet2Pt_2d");
  plotTH2(hImbalanceAlpha,"hImbalanceAlpha_2d");
  plotTH2(hImbalanceDeltaRPhoton1stJet,"hImbalanceDeltaRPhoton1stJet_2d");
  plotTH2(hImbalanceDeltaRPhoton2ndJet,"hImbalanceDeltaRPhoton2ndJet_2d");
  plotTH2(hImbalanceDeltaR1stJet2ndJet,"hImbalanceDeltaR1stJet2ndJet_2d");
  plotTH2(hImbalanceEta2ndJet,"hImbalanceEta2ndJet_2d");
  plotTH2(hImbalanceRatioPhotonPt12,"hImbalanceRatioPhotonPt12_2d");
  plotTH2(hImbalanceRatioPhotongenPhoton,"hImbalanceRatioPhotongenPhoton_2d");
  plotTH2(hImbalanceDeltaPhiPhoton2ndJet,"hImbalanceDeltaPhiPhoton2ndJet_2d");
  plotTH2(hImbalancePhoton2Pt,"hImbalancePhoton2Pt_2d");
  plotTH2(hImbalanceDeltaEtaPhoton2ndJet,"hImbalanceDeltaEtaPhoton2ndJet_2d");
  plotTH2(hImbalanceNobjPhoton,"hImbalanceNobjPhoton_2d");
  plotTH2(hImbalance2ndJetCorr,"hImbalance2ndJetCorr_2d");
  plotTH2(hImbalance1stGenJetID,"hImbalance1stGenJetID_2d");
  plotTH2(hPhotonEta2ndJet,"hPhotonEta2ndJet_2d");
  plotTH1(hEnergyBalance,"hEnergyBalance_1d",1);
  plotTH2(hImbalanceRatioPhotonPtJetColPt,"hImbalanceRatioPhotonPtJetColPt_2d");

  plotTH2(hJet1PtPhoton1Pt2,"Jet1PtPhoton1Pt2_2d");
  plotTH2(hImbalanceJetN2,"hImbalanceJetN2_2d");
  plotTH2(hImbalancePhoton1Pt2,"hImbalancePhoton1Pt2_2d");
  plotTH2(hImbalanceVtx2,"hImbalanceVtx2_2d");
  plotTH2(hImbalanceJetPt2,"hImbalanceJetPt2_2d");
  plotTH2(hImbalanceJetEta2,"hImbalanceJetEta2_2d");
  plotTH2(hImbalancePhotonEta2,"hImbalancePhotonEta2_2d");
  plotTH2(hImbalanceWeights2,"hImbalanceWeights2_2d");
  plotTH2(hImbalanceJet2Pt2,"hImbalanceJet2Pt2_2d");
  plotTH2(hImbalanceGenJet2Pt2,"hImbalanceGenJet2Pt2_2d");
  plotTH2(hJet2PtGenJet2Pt2,"hJet2PtGenJet2Pt2_2d");
  plotTH2(hImbalanceAlpha2,"hImbalanceAlpha2_2d");
  plotTH2(hImbalanceDeltaRPhoton1stJet2,"hImbalanceDeltaRPhoton1stJet2_2d");
  plotTH2(hImbalanceDeltaRPhoton2ndJet2,"hImbalanceDeltaRPhoton2ndJet2_2d");
  plotTH2(hImbalanceDeltaR1stJet2ndJet2,"hImbalanceDeltaR1stJet2ndJet2_2d");
  plotTH2(hImbalanceEta2ndJet2,"hImbalanceEta2ndJet2_2d");
  plotTH2(hImbalanceRatioPhotonPt122,"hImbalanceRatioPhotonPt122_2d");
  plotTH2(hImbalanceRatioPhotongenPhoton2,"hImbalanceRatioPhotongenPhoton2_2d");
  plotTH2(hImbalanceDeltaPhiPhoton2ndJet2,"hImbalanceDeltaPhiPhoton2ndJet2_2d");
  plotTH2(hImbalancePhoton2Pt2,"hImbalancePhoton2Pt2_2d");
  plotTH2(hImbalanceDeltaEtaPhoton2ndJet2,"hImbalanceDeltaEtaPhoton2ndJet2_2d");
  plotTH2(hImbalanceNobjPhoton2,"hImbalanceNobjPhoton2_2d");
  plotTH2(hImbalance2ndJetCorr2,"hImbalance2ndJetCorr2_2d");
  plotTH2(hImbalance1stGenJetID2,"hImbalance1stGenJetID2_2d");
  plotTH2(hPhotonEta2ndJet2,"hPhotonEta2ndJet2_2d");
  plotTH1(hEnergyBalance2,"hEnergyBalance2_1d",1);
  plotTH2(hImbalanceRatioPhotonPtJetColPt2,"hImbalanceRatioPhotonPtJetColPt2_2d");
  */

  // Photon Pt against number of vertices
  //plotTH2(hPhotonPtVtx,"PhotonPtNVtx_2d");

 // Photon Pt against number of vertices (before CUTS)
  //plotTH2(hPhotonPtVtxbeforeCUTS,"PhotonPtNVtxbeforeCUTS_2d");
 
 // PhotonPt against Jet Energy Fraction (scalar)
 // plotTH2(hResponseEnergyFraction,"Response_JetEnergyFraction_scalar_2d");
 
 // PhotonPt againstJet Energy Fraction (vectorial)
 //plotTH2(hResponseVecEnergyFraction,"Response_JetEnergyFraction_vectorial_2d");

 // Response against L1 correction
 //plotTH2(hResponseL1,"Response_L1correction_2d");

 // Response against L2L3 correction
 //plotTH2(hResponseL2L3,"Response_L2L3correction_2d");

  

 // Response against PhotonEta
 plotTH2(hResponsePhotonEta,"ResponsePhotonEta_2d");

 // Response against PhotonEta
 //plotTH2(hResponseLeadJet,"hResponseLeadJet");

 // Eta against Pt of first Jet
 //plotTH2(hJetEtaPt1,"JetEtaPt1_2d");

 // Pt of first and second photon
 // plotTH2(hPhoton12Pt,"Photon12Pt_2d");

 // Response against Pt of first Photon
 //plotTH2(hResponsePt,"ResponsePhoton1Pt_2d");

 // Response against Pt of first Photon for delta eta <1.5
 //plotTH2(hResponsePt_1_5_Deta,"ResponsePhoton1Pt_1_5_Deta_2d");

 // Response against Pt of first Photon for delta eta <1.1
 //plotTH2(hResponsePt_1_1_Deta,"ResponsePhoton1Pt_1_1_Deta_2d");

 // Response against Pt of first Photon for delta eta <0.5
 //plotTH2(hResponsePt_0_5_Deta,"ResponsePhoton1Pt_0_5_Deta_2d");

 

 // Response against EMF (Photon)
 //plotTH2(hResponsePhotonEMF,"ResponsePhotonEMF_2d");

 // Response against EMF (leading Jet)
 //plotTH2(hResponseJetEMF,"ResponseJetEMF_2d");

 // Response against Pt of first Photon
 //plotTH2(hResponsePhotonPtRatio,"ResponsePhotonPtRatio_2d");

 /*
 if(file==3){
 // Response against FHPD (Photon)
 plotTH2(hResponsePhotonFHPD,"ResponsePhotonFHPD_2d");

 //  Response against FRBX (Photon)
 plotTH2(hResponsePhotonFRBX,"ResponsePhotonFRBX_2d");

 // Response against FHPD (leading Jet)
 plotTH2(hResponseJetFHPD,"ResponseJetFHPD_2d");

 // Response against FRBX (leading Jet)
 plotTH2(hResponseJetFRBX,"ResponseJetFRBX_2d");
 }
 */

 // Response against Pt of first Photon before dead Ecal cell cut 
 // plotTH2(hResponsePtBeforeEcal,"ResponsePhoton1Pt_before_ecal_cut_2d");

 // Response against Eta of first Photon in high Pt region
 //plotTH2(hPhotonEta_high_pt_reg_Response,"ResponsePhotonEta_high_pt_reg_2d");

 // Response against Jet Pt
 //plotTH2(hSmallJetPtResponse,"ResponseJetPt_small_resp_reg_2d");*/

 

  
  
//----------------------------plots which could not be automized ------------------------


 // some declarations
  char cname[100];
  char tname[100];
  TString legname;
  TString TotFilename;
  
  TFile *f;
  
 
  // Pt Binned Response Function
  TCanvas* cRespBinnedPt[pt_int][2][alpha_int] = {{{0}}};
  TLegend* leg_histogram  = 0; 
  for(int i=0;i<pt_int;i++){
    for(int j=0; j<eta_int; j++){
      for(int k=0; k<alpha_int; k++){

	// if(i==1 && j==1 && k==5 || i==2 && j==0 && k==2 ||i==3 && j==0 && k==3 ||i==5 && j==0 && k==1 ||i==5 && j==1 && k==5 ){
	if(JetResponseJetHemisphere[i][j][k]->hResponse -> GetEntries()<=31 || JetResponsePhotonHemisphere[i][j][k]->hResponse -> GetEntries()<=31) continue;
	JetResponseJetHemisphere[i][j][k]->hResponse -> SetMaximum((JetResponseJetHemisphere[i][j][k]->hResponse->GetMaximum())*10);
	sprintf(cname,"cRespBinnedPtJetHem_%i_%i_%i", i,j,k);

        if(i==pt_int-1) sprintf(tname,"Response function (Jet Hemisphere) with %g GeV < Pt1_gamma (%i. eta bin and %i. alpha bin)",bd[pt_int-1], j+1, k+1);
        else sprintf(tname,"Response function (Jet Hemisphere) with %g GeV <Pt1Gamma< %g GeV  (%i. eta bin and %i. alpha bin)",bd[i],bd[i+1], j+1,k+1);
	
        cRespBinnedPt[i][j][k] = new TCanvas(cname,tname,0,0,500,500);
        cRespBinnedPt[i][j][k] -> cd(); 
        //JetResponseJetHemisphere[i][j][k]->hResponse -> Draw();
	
	
	
        if(type == 2){
	  JetIntrinsic[i][j][k]->hResponse -> GetXaxis()->SetLabelSize(0.04);
	  JetIntrinsic[i][j][k]->hResponse -> GetXaxis()->SetTitleOffset(1.2);
	  JetIntrinsic[i][j][k]->hResponse -> GetYaxis()->SetLabelSize(0.04);
	  JetIntrinsic[i][j][k]->hResponse -> GetYaxis()->SetTitleOffset(1.3);
	  JetIntrinsic[i][j][k]->hResponse -> GetYaxis()->SetTitle("# Ereignisse");
	  //JetImbalance[i][j][k]->hResponse -> Draw();
	  //JetImbalance[i][j][k]->hResponse -> Draw("same");
	  JetIntrinsic[i][j][k]->hResponse -> Draw(); 
	  //JetIntrinsic[i][j][k]->hResponse -> Draw("same");
	  //JetImbalance[i][j][k]->hResponse -> SetLineColor(3);  
	  //JetIntrinsic[i][j][k]->hResponse -> SetLineColor(6); 
	} 
      



        leg_histogram = new TLegend(0.10,0.75,0.78,0.9);
        leg_histogram -> SetFillColor(0);
        leg_histogram -> SetFillStyle(0);
        leg_histogram->SetTextSize(0.033);
        if(i==pt_int-1) legname.Form("#splitline{%4.1f GeV < P_{T}^{#gamma} and %4.1f < #alpha < %4.1f}{%4.1f < #eta < %4.1f}",bd[i],alphaBin[k],alphaBin[k+1],etaBin[j],etaBin[j+1]);
        else legname.Form("#splitline{%4.1f GeV < P_{T}^{#gamma} < %4.1f GeV and %4.1f < #alpha < %4.1f}{%4.1f < #eta < %4.1f}",bd[i],bd[i+1],alphaBin[k],alphaBin[k+1],etaBin[j],etaBin[j+1]);       
        //leg_histogram->SetHeader(legname);
        //leg_histogram->Draw("same");     
        cRespBinnedPt[i][j][k] -> SetLogy();
	
        // Saving plot as PDf
        TotFilename = PDFPath + "response_jet_in_";
        TotFilename+= (i+1);
        TotFilename+= "_Pt_bin_";
        TotFilename+= (j+1);
        TotFilename+= "_eta_bin_";
        TotFilename+= (k+1);
        TotFilename+= "_alpha_bin_" + DataType + ".pdf";      
        cRespBinnedPt[i][j][k] -> SaveAs(TotFilename);

        // Writing histogram in root file
        TotFilename = RootPath + "response_jet_in_";
        TotFilename+= (i+1);
        TotFilename+= "_Pt_bin_";
        TotFilename+= (j+1);
        TotFilename+= "_eta_bin_";
        TotFilename+= (k+1);
        TotFilename+= "_alpha_bin_" + DataType + ".root";         
        f = new TFile(TotFilename,"RECREATE");
        f -> WriteTObject(JetResponseJetHemisphere[i][j][k]->hResponse);
        f->Close();
        delete f;
	

	TotFilename = RootPath + "response_intrinsic_jet_in_";
        TotFilename+= (i+1);
        TotFilename+= "_Pt_bin_";
        TotFilename+= (j+1);
        TotFilename+= "_eta_bin_";
        TotFilename+= (k+1);
        TotFilename+= "_alpha_bin_" + DataType + ".root";    
	f = new TFile(TotFilename,"RECREATE");
        f -> WriteTObject(JetIntrinsic[i][j][k]->hResponse);
        f->Close();
        delete f;

	TotFilename = RootPath + "response_imbalance_jet_in_";
        TotFilename+= (i+1);
        TotFilename+= "_Pt_bin_";
        TotFilename+= (j+1);
        TotFilename+= "_eta_bin_";
        TotFilename+= (k+1);
        TotFilename+= "_alpha_bin_" + DataType + ".root";
	f = new TFile(TotFilename,"RECREATE");
        f -> WriteTObject(JetImbalanceJetHemisphere[i][j][k]->hResponse);
        f->Close();
        delete f;
	delete  cRespBinnedPt[i][j][k];
      
      }
    }
  }

  
  for(int i=0;i<pt_int;i++){
    for(int j=0; j<eta_int; j++){
      for(int k=0; k<alpha_int; k++){
	
	
	if(JetResponsePhotonHemisphere[i][j][k]->hResponse -> GetEntries()<=31 || JetResponseJetHemisphere[i][j][k]->hResponse -> GetEntries()<=31 ) continue;

	JetResponsePhotonHemisphere[i][j][k]->hResponse -> SetMaximum((JetResponsePhotonHemisphere[i][j][k]->hResponse->GetMaximum())*10);
	sprintf(cname,"cRespBinnedPtPhotonHem_%i_%i_%i", i,j,k);

        if(i==pt_int-1) sprintf(tname,"Response function (Photon Hemisphere) with %g GeV < Pt1_gamma (%i. eta bin and %i. alpha bin)",bd[pt_int-1], j+1, k+1);
        else sprintf(tname,"Response function (Photon Hemisphere) with %g GeV <Pt1Gamma< %g GeV  (%i. eta bin and %i. alpha bin)",bd[i],bd[i+1], j+1,k+1);
  
        cRespBinnedPt[i][j][k] = new TCanvas(cname,tname,0,0,500,500);
        cRespBinnedPt[i][j][k] -> cd(); 
        JetResponsePhotonHemisphere[i][j][k]->hResponse -> Draw();
        if(type == 2){
	  //JetImbalance[i][j][k]->hResponse -> Draw();
	  //JetImbalance[i][j][k]->hResponse -> Draw("same");
	  //JetIntrinsic[i][j][k]->hResponse -> Draw(); 
	  //JetIntrinsic[i][j][k]->hResponse -> Draw("same");
	  //JetImbalance[i][j][k]->hResponse -> SetLineColor(3);  
	  //JetIntrinsic[i][j][k]->hResponse -> SetLineColor(6); 
	} 
      
        leg_histogram = new TLegend(0.10,0.75,0.78,0.9);
        leg_histogram -> SetFillColor(0);
        leg_histogram -> SetFillStyle(0);
        leg_histogram->SetTextSize(0.033);
        if(i==pt_int-1) legname.Form("#splitline{%4.1f GeV < P_{T}^{#gamma} and %4.1f < #alpha < %4.1f}{%4.1f < #eta < %4.1f}",bd[i],alphaBin[k],alphaBin[k+1],etaBin[j],etaBin[j+1]);
        else legname.Form("#splitline{%4.1f GeV < P_{T}^{#gamma} < %4.1f GeV and %4.1f < #alpha < %4.1f}{%4.1f < #eta < %4.1f}",bd[i],bd[i+1],alphaBin[k],alphaBin[k+1],etaBin[j],etaBin[j+1]);       
        leg_histogram->SetHeader(legname);
        leg_histogram->Draw("same");     
        cRespBinnedPt[i][j][k] -> SetLogy();

        // Saving plot as PDf
        TotFilename = PDFPath + "response_photon_in_";
        TotFilename+= (i+1);
        TotFilename+= "_Pt_bin_";
        TotFilename+= (j+1);
        TotFilename+= "_eta_bin_";
        TotFilename+= (k+1);
        TotFilename+= "_alpha_bin_" + DataType + ".pdf";      
        cRespBinnedPt[i][j][k] -> SaveAs(TotFilename);

        // Writing histogram in root file
        TotFilename = RootPath + "response_photon_in_";
        TotFilename+= (i+1);
        TotFilename+= "_Pt_bin_";
        TotFilename+= (j+1);
        TotFilename+= "_eta_bin_";
        TotFilename+= (k+1);
        TotFilename+= "_alpha_bin_" + DataType + ".root";         
        f = new TFile(TotFilename,"RECREATE");
        f -> WriteTObject(JetResponsePhotonHemisphere[i][j][k]->hResponse);
        f->Close();
        delete f;

	

	TotFilename = RootPath + "response_imbalance_photon_in_";
        TotFilename+= (i+1);
        TotFilename+= "_Pt_bin_";
        TotFilename+= (j+1);
        TotFilename+= "_eta_bin_";
        TotFilename+= (k+1);
        TotFilename+= "_alpha_bin_" + DataType + ".root";
	f = new TFile(TotFilename,"RECREATE");
        f -> WriteTObject(JetImbalancePhotonHemisphere[i][j][k]->hResponse);
        f->Close();
        delete f;
	delete  cRespBinnedPt[i][j][k];
	
      
      }
    }
  }  
  
  
  /*
  // Delta Phi
  TCanvas* cDphi = new TCanvas("cDphi","Delta Phi",0,0,500,500);
  cDphi->SetLeftMargin(0.17);
  hDphi->GetYaxis()->SetTitleOffset(1.5); 
  cDphi->cd();
  hDphi->Draw();
  TLegend* leg_hist       =  0 ;
  double_t inte = 0;
  int bmin = hDphi->FindBin(2.95); 
  //cout<<"bmin="<<bmin<<endl;
  int bmax = hDphi->FindBin(3.14);
  //cout<<"bmax="<<bmax<<endl;
  inte = (hDphi->Integral(bmin,bmax)/hDphi->Integral())*100;
  sprintf(legname,"in %g %% of events is #Delta #Phi >= 2.95",inte );
  leg_hist = new TLegend(0.2,0.6,0.81,0.77);
  leg_hist->SetTextSize(0.033);
  leg_hist->SetHeader(legname);
  leg_hist->Draw("same");
  sprintf(tot_filename,"%sdelta_phi.pdf",savename);
  cDphi->SaveAs(tot_filename);  
  sprintf(tot_filename,"%sdelta_phi.root",save_root);
  f = new TFile(tot_filename,"RECREATE");
  hDphi->WriteTObject();
  f->Close();
  delete f;*/

  /*
  
  // Pt of first (gen and reco) and second (reco) Jet
  TLegend* leg_histo  = 0;  
  TCanvas* cJet1Pt = new TCanvas("cJet1Pt","Pt of first Jet",0,0,500,500);
  cJet1Pt -> cd();
  hJet2Pt -> Draw();
  hJet1Pt -> Draw("same");
  hJet2Pt -> SetLineColor(2);  
  leg_histo = new TLegend(0.21,0.3,0.5,0.49);
  leg_histo -> AddEntry(hJet1Pt,"first Jet","l");
  leg_histo -> AddEntry(hJet2Pt,"second Jet","l");
  leg_histo -> Draw("same");
  sprintf(sname_1,"Jet12Pt");
  sprintf(tot_filename,"%s%s.pdf",savename,sname_1);
  cJet1Pt ->SaveAs(tot_filename);
  sprintf(tot_filename,"%s%s.root",save_root,sname_1);
  f = new TFile(tot_filename,"RECREATE");
  hJet1Pt->WriteTObject();
  f->Close();
  delete f;*/
  gErrorIgnoreLevel = oldLevel;

}
