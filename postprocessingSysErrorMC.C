// ------------------------------------------------------------------------------------------------
// -------  Script to compare Jet Energy Resolutions and Scales between Data and MC (T.L.)  ------- 
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TChain.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TObject.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"

//#include "utils.h"
//#include "HistOps.h"

TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName);
TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName);
TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName);


const int numEntries = 9;
int postprocessing(int etaBin){

  // For looking at different systematic uncertainties independently
  bool flav = false;
  bool PU = false;
  bool PES = false;
  bool MC = false;
  bool JEC = true;

  TString etaString;   
  char pdfFile[100];

  TString path = "plots_2012/PF_L1CHS/";
  TString rootFiles;
  
  
  // Read the MC files and all systematic error files
  rootFiles.Form("mc/root_files_3sigma/Resolution_for_%i_eta_bin_PFCHS_mc.root",etaBin);
  TGraphErrors* JERMC = readTGraphErrors(path+rootFiles,"Graph;1","Graph;1");
  
  path = "plots_2012/PF_L1CHS/systematicUncertaintiesMC/";
  /*
  rootFiles.Form("JetPartonFlavorUncer/root_files/Resolution_for_%i_eta_bin_quark_PF_mc.root",etaBin);  
  TGraph* sysJetFlavorQuark = readTGraph(path+rootFiles,"Graph;1","Graph;1");   
  rootFiles.Form("JetPartonFlavorUncer/root_files/Resolution_for_%i_eta_bin_gluon_PF_mc.root",etaBin);  
  TGraph* sysJetFlavorGluon = readTGraph(path+rootFiles,"Graph;1","Graph;1");
 
  
  rootFiles.Form("PUuncertanties/root_files_data_le_10/Resolution_for_%i_eta_bin_PF_data.root",etaBin);  
  TGraph* sysJetPUle10 = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("PUuncertanties/root_files_data_gt_10_le_15/Resolution_for_%i_eta_bin_PF_data.root",etaBin);  
  TGraph* sysJetPUgt10 = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("PUuncertanties/root_files_data_gt_15_le_20/Resolution_for_%i_eta_bin_PF_data.root",etaBin);  
  TGraph* sysJetPUgt15 = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("PUuncertanties/root_files_data_gt_20/Resolution_for_%i_eta_bin_PF_data.root",etaBin);  
  TGraph* sysJetPUgt20 = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  */
  /*
  rootFiles.Form("plots_2012/PF_L1FastJet/mc/root_files/Resolution_for_%i_eta_bin_PUle10_PF_mc.root",etaBin);  
  TGraph* sysJetPUle10 = readTGraph(rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("plots_2012/PF_L1FastJet/mc/root_files/Resolution_for_%i_eta_bin_PUgt10le15_PF_mc.root",etaBin);  
  TGraph* sysJetPUgt10 = readTGraph(rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("plots_2012/PF_L1FastJet/mc/root_files/Resolution_for_%i_eta_bin_PUgt15le20_PF_mc.root",etaBin);  
  TGraph* sysJetPUgt15 = readTGraph(rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("plots_2012/PF_L1FastJet/mc/root_files/Resolution_for_%i_eta_bin_PUgt20_PF_mc.root",etaBin);  
  TGraph* sysJetPUgt20 = readTGraph(rootFiles,"Graph;1","Graph;1");
  */
  /*
  
  rootFiles.Form("ResImbUncertainties/root_files_lower_bound/Resolution_for_%i_eta_bin_PF_data.root",etaBin);  
  TGraph* sysMClow = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("ResImbUncertainties/root_files_upper_bound/Resolution_for_%i_eta_bin_PF_data.root",etaBin);    
  TGraph* sysMCup = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  */  
  rootFiles.Form("JECUncertainties/root_files_mc_lower_bound/Resolution_for_%i_eta_bin_PFCHS_mc.root",etaBin);    
  TGraphErrors* sysJEClow = readTGraphErrors(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("JECUncertainties/root_files_mc_upper_bound/Resolution_for_%i_eta_bin_PFCHS_mc.root",etaBin);    
  TGraphErrors* sysJECup = readTGraphErrors(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("PESUncertainties/root_files_mc_lower_bound/Resolution_for_%i_eta_bin_PFCHS_mc.root",etaBin);  
  TGraph* sysPESlow = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("PESUncertainties/root_files_mc_upper_bound/Resolution_for_%i_eta_bin_PFCHS_mc.root",etaBin);  
  TGraph* sysPESup = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  
  sprintf(pdfFile,"Resolution_for_%d_eta_bin_PF_data_mc_comparison.pdf",etaBin);
  if(etaBin == 1)      etaString = "Jet Energy Resolution for |#eta|<0.5";
  else if(etaBin == 2) etaString = "Jet Energy Resolution for 0.5<|#eta|<1.1";
  else if(etaBin == 3) etaString = "Jet Energy Resolution for 1.1<|#eta|<1.7";
  else if(etaBin == 4) etaString = "Jet Energy Resolution for 1.7<|#eta|<2.3";
  else if(etaBin == 5) etaString = "Jet Energy Resolution for 2.3<|#eta|<5.0";
  
  
 
  //const int numEntries = JERMC->GetN();
  cout<<endl<<"Number of entries in JERMC = "<<numEntries<<endl<<endl;
  
  
  double xMC[numEntries] = {0};
  double yMC[numEntries] = {0};
  double xMCError[numEntries] = {0};
  double yMCError[numEntries]  = {0};
  
   
  // loop over pt bins for every systematic uncertainty
  int idxMC = 0;
  for(int i=0; i<numEntries; i++){

    JERMC     -> GetPoint(i,xMC[idxMC],yMC[idxMC]);
    xMCError[idxMC] = JERMC -> GetErrorX(i);
    yMCError[idxMC] = JERMC -> GetErrorY(i);

    cout<<"xMC["<<idxMC<<"] = "<<xMC[idxMC]<<endl;
    cout<<"yMC["<<idxMC<<"] = "<<yMC[idxMC]<<endl;
    cout<<"xMCError["<<idxMC<<"] = "<<xMCError[idxMC]<<endl;
    cout<<"yMCError["<<idxMC<<"] = "<<yMCError[idxMC]<<endl;
    idxMC += 1;
  }
  cout<<endl;

  double xJECUp[numEntries] = {0};
  double yJECUp[numEntries] = {0};
  double xJECErrorUp[numEntries] = {0};
  double yJECErrorUp[numEntries] = {0};
  double xJECLow[numEntries] = {0};
  double yJECLow[numEntries] = {0};
  double xJECErrorLow[numEntries] = {0};
  double yJECErrorLow[numEntries] = {0};
  double yJECsysUp[numEntries] = {0};
  double yJECErrorsysUp[numEntries] = {0};
  double yJECsysLow[numEntries] = {0};
  double yJECErrorsysLow[numEntries] = {0};
  double sysUncertaintyJEC[numEntries] = {0};
  double sysUncertaintyJECError[numEntries] = {0};
  
    
  if(JEC){


    // Plot upper and lower uncertainty and compare if they are compatible
    for(int i=0; i<numEntries; i++){
      sysJECup  -> GetPoint(i,xJECUp[i],yJECUp[i]);
      sysJEClow -> GetPoint(i,xJECLow[i],yJECLow[i]);
   
      xJECErrorUp[i]  = sysJECup -> GetErrorX(i);
      yJECErrorUp[i]  = sysJECup -> GetErrorY(i);
      xJECErrorLow[i] = sysJEClow -> GetErrorX(i);
      yJECErrorLow[i]  = sysJEClow -> GetErrorY(i);

      cout<<"xJECUp["<<i<<"] = "<<xJECUp[i]<<endl;
      cout<<"xJECLow["<<i<<"] = "<<xJECLow[i]<<endl;  
      cout<<"xMC["<<i<<"] = "<<xMC[i]<<endl;     
      cout<<"yMC["<<i<<"] = "<<yMC[i]<<endl;
      cout<<"yJECUp["<<i<<"] = "<<yJECUp[i]<<endl;
      cout<<"yJECLow["<<i<<"] = "<<yJECLow[i]<<endl<<endl;

      // Calculate delta with JEC Res*(1+delta)
      yJECsysUp[i]  = std::abs(yJECUp[i]/yMC[i] - 1.);
      yJECsysLow[i] = std::abs(yJECLow[i]/yMC[i] - 1.);
      yJECErrorsysUp[i] = TMath::Sqrt(TMath::Power(1./yMC[i],2)*TMath::Power(yJECErrorUp[i],2));
      yJECErrorsysLow[i] = TMath::Sqrt(TMath::Power(1./yMC[i],2)*TMath::Power(yJECErrorLow[i],2));
      

      if(yJECsysUp[i] != yJECsysUp[i]) cout<<"yJECsysUp:  i = "<<i<<endl;
      if(yJECErrorsysUp[i] != yJECErrorsysUp[i]) cout<<"yJECErrorsysUp:  i = "<<i<<endl;
      if(yJECsysLow[i] != yJECsysLow[i]) cout<<"yJECsysLow:  i = "<<i<<endl;
      if(yJECErrorsysLow[i] != yJECErrorsysLow[i]) cout<<"yJECErrorsysLow:  i = "<<i<<endl;
      if(xMC[i] != xMC[i]) cout<<"xMC: i = "<<i<<endl;
      if(xMCError[i] != xMCError[i]) cout<<"xMCError: i = "<<i<<endl;


      //Calculate for every pt bin a symmetric error for JEC
      if((yJECUp[i]<yMC[i] && yJECLow[i]>yMC[i]) || (yJECUp[i]>yMC[i] && yJECLow[i]<yMC[i])){
	
	sysUncertaintyJEC[i] = (std::abs(yJECUp[i]-yJECLow[i])/2.)/yMC[i];
	sysUncertaintyJECError[i] = TMath::Sqrt(TMath::Power(1./(2.*yMC[i]),2)*(TMath::Power(yJECErrorUp[i],2)+TMath::Power(yJECErrorLow[i],2))+TMath::Power(std::abs(yJECUp[i]-yJECLow[i])/(2.*TMath::Power(yMC[i],2)),2)*TMath::Power(yMCError[i],2));
	cout<<"enclosing JER!"<<endl<<endl<<endl;
      }
      else{

	if(std::abs(1-yJECUp[i]/yMC[i]) > std::abs(1-yJECLow[i]/yMC[i])) {
	  sysUncertaintyJEC[i] = std::abs((yJECUp[i]-yMC[i])/2.)/yMC[i];
	  
	  sysUncertaintyJECError[i] = TMath::Sqrt(TMath::Power(1./(2.*yMC[i]),2)*TMath::Power(yJECErrorUp[i],2)+TMath::Power(yJECUp[i]/(2*TMath::Power(yMC[i],2)),2)*TMath::Power(yMCError[i],2));

	}
	else{
	  sysUncertaintyJEC[i] = std::abs((yJECLow[i]-yMC[i]))/2./yMC[i];
	  sysUncertaintyJECError[i] = TMath::Sqrt(TMath::Power(1./(2.*yMC[i]),2)*TMath::Power(yJECErrorLow[i],2)+TMath::Power(yJECLow[i]/(2*TMath::Power(yMC[i],2)),2)*TMath::Power(yMCError[i],2));

	}
	cout<<"NOT enclosing JER!"<<endl<<endl;

	cout<<"yMCError["<<i<<"] = "<<yMCError[i]<<endl;
	cout<<"yJECErrorUp["<<i<<"] = "<<yJECErrorUp[i]<<endl;
	cout<<"yJECErrorLow["<<i<<"] = "<<yJECErrorLow[i]<<endl<<endl<<endl<<endl;
      }
      
      cout<<"sysUncertaintyJECError[i] = "<<sysUncertaintyJECError[i]<<endl;
    }

    TCanvas *cJEC3 = new TCanvas("cJEC3",etaString,200,10,500,500);
    cJEC3 -> SetLeftMargin(0.17);
    cJEC3 -> cd();
    TGraphErrors* sysUncertaintyFinal = new TGraphErrors(numEntries,xJECUp,sysUncertaintyJEC,xJECErrorUp,sysUncertaintyJECError);
    TF1* f1 = new TF1("name","pol0",0,600);
    sysUncertaintyFinal -> Fit("name","");
    cout<<endl<<"chi2 = "<<f1 -> GetChisquare()<<endl<<endl; 
    sysUncertaintyFinal -> SetTitle("Final Relative JEC Uncertainty");
    sysUncertaintyFinal -> GetXaxis()->SetTitle("p_{T}^{#gamma}");
    sysUncertaintyFinal -> GetYaxis()->SetTitleOffset(1.5);
    sysUncertaintyFinal -> GetYaxis()->SetTitle("rel. Uncertainty");
    sysUncertaintyFinal -> SetMarkerStyle(22);
    sysUncertaintyFinal -> SetMaximum(0.1);
    sysUncertaintyFinal -> SetMarkerSize(1);
    sysUncertaintyFinal -> SetMarkerColor(1);
    sysUncertaintyFinal -> SetLineColor(1);
    sysUncertaintyFinal -> Draw("AP");

    TLatex*  info   = new TLatex();
    char legname[100] = {0};
    info->SetTextFont(132);
    info-> SetNDC();
    info->SetTextSize(0.040);
    sprintf(legname,"#splitline{#chi^{2} = %4.2f}{dof = %i}",f1 -> GetChisquare(),f1 -> GetNDF());
    info->DrawLatex(0.60,0.80,legname);
    
    delete info;
    info   = new TLatex();
    info->SetTextFont(132);
    info-> SetNDC();
    info->SetTextSize(0.040);
    sprintf(legname,"#splitline{sys Uncer. = %4.3f}{stat Error = %4.3f}",f1 -> GetParameter(0),f1 -> GetParError(0));
    info->DrawLatex(0.60,0.6,legname);
    
    
    TGraphErrors* JECsysUp = new TGraphErrors(numEntries,xJECUp,yJECsysUp,xJECErrorUp,yJECErrorsysUp);
    TGraphErrors* JECsysLow = new TGraphErrors(numEntries,xJECLow,yJECsysLow,xJECErrorLow,yJECErrorsysLow);

    TCanvas *cJEC = new TCanvas("cJEC",etaString,200,10,500,500);
    cJEC -> SetLeftMargin(0.17);
    cJEC -> cd();

    JECsysUp->SetTitle("JEC Uncertainty");
    JECsysUp->GetYaxis()->SetTitle("|JER(#pm #Delta)/JER - 1|");
    JECsysUp->GetYaxis()->SetTitleOffset(1.5);
    JECsysUp->GetXaxis()->SetTitle("photonPt");

    JECsysUp -> SetMarkerStyle(22);
    JECsysUp -> SetMarkerSize(1);
    JECsysUp -> SetMarkerColor(3);
    JECsysUp -> SetLineColor(3);
    JECsysUp -> Draw("AP");
    JECsysLow -> SetMarkerStyle(20);
    JECsysLow -> SetMarkerColor(2);
    JECsysLow -> SetLineColor(2);
    JECsysLow -> Draw("sameP");

    TLegend *legend  = new TLegend(0.6,0.8,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.033);
    legend -> AddEntry(JECsysUp,"JEC(1+#Delta)","l");
    legend -> AddEntry(JECsysLow,"JEC(1-#Delta)","l");
    legend -> Draw("same");
    sprintf(pdfFile,"Delta_JEC_Uncertainties_for_%d_eta_bin_PFCHS_mc.pdf",etaBin);
    cJEC->SaveAs(path + "JECUncertainties/plots_mc/" +pdfFile);

    TCanvas *cJEC2 = new TCanvas("cJEC2",etaString,200,10,500,500);
    cJEC2 -> SetLeftMargin(0.17);
    cJEC2 -> cd();

    JERMC -> GetFunction("fResolution") -> SetLineColor(1);
    JERMC->SetMarkerStyle(21);
    JERMC->SetMarkerSize(1);
    JERMC ->SetMinimum(0.04);
    JERMC->SetMaximum(0.11);
    JERMC->Draw("AP");
    sysJECup -> SetMarkerStyle(22);
    sysJECup -> SetMarkerColor(3);
    sysJECup -> SetLineColor(3);
    sysJECup -> GetFunction("fResolution") -> SetLineColor(3);
    sysJECup -> Draw("sameP");
    sysJEClow -> GetFunction("fResolution") -> SetLineColor(2);
    sysJEClow -> SetMarkerStyle(20);
    sysJEClow -> SetMarkerColor(2);
    sysJEClow -> SetLineColor(2);
    sysJEClow -> Draw("sameP");
    delete legend;
    legend  = new TLegend(0.6,0.8,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.033);
    legend -> AddEntry(JERMC,"JEC","l");
    legend -> AddEntry(sysJECup,"JEC(1+#Delta)","l");
    legend -> AddEntry(sysJEClow,"JEC(1-#Delta)","l");
    legend -> Draw("same");

    sprintf(pdfFile,"Resolutions_JEC_Uncertainties_for_%d_eta_bin_PFCHS_mc.pdf",etaBin);
    cJEC2->SaveAs(path + "JECUncertainties/plots_mc/" +pdfFile);


 
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // PES Uncertainty

  double xPESUp[numEntries] = {0};
  double yPESUp[numEntries] = {0};
  double xPESErrorUp[numEntries] = {0};
  double yPESErrorUp[numEntries] = {0};
  double xPESLow[numEntries] = {0};
  double yPESLow[numEntries] = {0};
  double xPESErrorLow[numEntries] = {0};
  double yPESErrorLow[numEntries] = {0};
  double yPESsysUp[numEntries] = {0};
  double yPESErrorsysUp[numEntries] = {0};
  double yPESsysLow[numEntries] = {0};
  double yPESErrorsysLow[numEntries] = {0};
  
    
  if(PES){


    // Plot upper and lower uncertainty and compare if they are compatible

    for(int i=0; i<numEntries; i++){
      sysPESup  -> GetPoint(i,xPESUp[i],yPESUp[i]);
      sysPESlow -> GetPoint(i,xPESLow[i],yPESLow[i]);
   
      xPESErrorUp[i]  = sysPESup -> GetErrorX(i);
      yPESErrorUp[i]  = sysPESup -> GetErrorY(i);
      xPESErrorLow[i] = sysPESlow -> GetErrorX(i);
      yPESErrorLow[i]  = sysPESlow -> GetErrorY(i);

      cout<<"xPESUp["<<i<<"] = "<<xPESUp[i]<<endl;
      cout<<"xPESLow["<<i<<"] = "<<xPESLow[i]<<endl;     
      cout<<"yMC["<<i<<"] = "<<yMC[i]<<endl;
      cout<<"yPESUp["<<i<<"] = "<<yPESUp[i]<<endl;
      cout<<"yPESLow["<<i<<"] = "<<yPESLow[i]<<endl<<endl;

      // Calculate delta with PES Res*(1+delta)
      yPESsysUp[i]  = std::abs(yPESUp[i]/yMC[i] - 1.);
      yPESsysLow[i] = std::abs(yPESLow[i]/yMC[i] - 1.);
      yPESErrorsysUp[i] = TMath::Sqrt(TMath::Power(1./yMC[i],2)*TMath::Power(yPESErrorUp[i],2));
      yPESErrorsysLow[i] = TMath::Sqrt(TMath::Power(1./yMC[i],2)*TMath::Power(yPESErrorLow[i],2));
      

      if(yPESsysUp[i] != yPESsysUp[i]) cout<<"yPESsysUp:  i = "<<i<<endl;
      if(yPESErrorsysUp[i] != yPESErrorsysUp[i]) cout<<"yPESErrorsysUp:  i = "<<i<<endl;
      if(yPESsysLow[i] != yPESsysLow[i]) cout<<"yPESsysLow:  i = "<<i<<endl;
      if(yPESErrorsysLow[i] != yPESErrorsysLow[i]) cout<<"yPESErrorsysLow:  i = "<<i<<endl;
      if(xMC[i] != xMC[i]) cout<<"xMC: i = "<<i<<endl;
      if(xMCError[i] != xMCError[i]) cout<<"xMCError: i = "<<i<<endl;
      

    }
    
    TGraphErrors* PESsysUp = new TGraphErrors(numEntries,xPESUp,yPESsysUp,xPESErrorUp,yPESErrorsysUp);
    TGraphErrors* PESsysLow = new TGraphErrors(numEntries,xPESLow,yPESsysLow,xPESErrorLow,yPESErrorsysLow);

    TCanvas *cPES = new TCanvas("cPES",etaString,200,10,500,500);
    cPES -> SetLeftMargin(0.17);
    cPES -> cd();

    PESsysUp->SetTitle("PES Uncertainty");
    PESsysUp->GetYaxis()->SetTitle("|JER(#pm #Delta)/JER - 1|");
    PESsysUp->GetYaxis()->SetTitleOffset(1.1);
    PESsysUp->GetXaxis()->SetTitle("photonPt");

    PESsysUp -> SetMarkerStyle(22);
    PESsysUp -> SetMarkerSize(1);
    PESsysUp -> SetMarkerColor(3);
    PESsysUp -> SetLineColor(3);
    PESsysUp -> Draw("AP");
    PESsysLow -> SetMarkerStyle(20);
    PESsysLow -> SetMarkerColor(2);
    PESsysLow -> SetLineColor(2);
    PESsysLow -> Draw("sameP");

    TLegend *legend  = new TLegend(0.6,0.8,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.033);
    legend -> AddEntry(PESsysUp,"PES(1+#Delta)","l");
    legend -> AddEntry(PESsysLow,"PES(1-#Delta)","l");
    legend -> Draw("same");
    sprintf(pdfFile,"Delta_PES_Uncertainties_for_%d_eta_bin_PFCHS_mc.pdf",etaBin);
    cPES->SaveAs(path + "PESUncertainties/plots_mc/" +pdfFile);

    TCanvas *cPES2 = new TCanvas("cPES2",etaString,200,10,500,500);
    cPES2 -> SetLeftMargin(0.17);
    cPES2 -> cd();

    JERMC -> GetFunction("fResolution") -> SetLineColor(1);
    JERMC->SetMarkerStyle(21);
    JERMC->SetMarkerSize(1);
    JERMC ->SetMinimum(0.05);
    JERMC->SetMaximum(0.13);
    JERMC->Draw("AP");
    sysPESup -> SetMarkerStyle(22);
    sysPESup -> SetMarkerColor(3);
    sysPESup -> SetLineColor(3);
    sysPESup -> GetFunction("fResolution") -> SetLineColor(3);
    sysPESup -> Draw("sameP");
    sysPESlow -> GetFunction("fResolution") -> SetLineColor(2);
    sysPESlow -> SetMarkerStyle(20);
    sysPESlow -> SetMarkerColor(2);
    sysPESlow -> SetLineColor(2);
    sysPESlow -> Draw("sameP");
    delete legend;
    legend  = new TLegend(0.6,0.8,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.033);
    legend -> AddEntry(JERMC,"PES","l");
    legend -> AddEntry(sysPESup,"PES(1+#Delta)","l");
    legend -> AddEntry(sysPESlow,"PES(1-#Delta)","l");
    legend -> Draw("same");

    sprintf(pdfFile,"Resolutions_PES_Uncertainties_for_%d_eta_bin_PFCHS_mc.pdf",etaBin);
    cPES2->SaveAs(path + "PESUncertainties/plots_mc/" +pdfFile);


 
  }
 
  return 0;
    /*
    //Plot

  // 5.) Symmetrize JEC uncertanty 
  
  idxData=0;
  idx1 = 0;
  idx2 = 0;
  while(idxData<numEntries){
    sysJECup -> GetPoint(idx1,xAux1,yAux1);
    sysJEClow -> GetPoint(idx2,xAux2,yAux2);
   
    cout<<"idx1 = "<<idx1<<endl;
    cout<<"idx2 = "<<idx2<<endl;
    cout<<"xJECup["<<idx1<<"] = "<<xAux1<<endl;
    cout<<"xJEClow["<<idx2<<"] = "<<xAux2<<endl;
    if(abs(xAux1-xAux2)/xAux1>0.1){
      if(xAux1>xAux2) idx2 += 1;
      else idx1 += 1;
      cout<<"Error JEC sys!!"<<endl;
      continue;
    }
    if(abs(xAux1-xData[idxData])/xData[idxData]<0.1){
      yJEC[idxData] = abs(yAux1-yAux2)/2.;
      cout<<"yJEC["<<idxData<<"] = "<<yJEC[idxData]<<endl;
      idxData += 1;
      idx1 += 1;
      idx2 += 1;
    }
    else if(abs(xAux1-xData[idxData])/xData[idxData]>0.1 && xAux1<xData[idxData]){
      idx1 += 1;
      idx2 += 1;
      continue;
    }
    else{
      cout<<"Error in JEC!!!!"<<endl;
      return 0;
    }
  }
  }
    

  //Calculate all systematic uncertainties (add in quadrature)
  for(int i=0; i<numEntries; i++){
    xSys[i] = 0;
    ySys[i] = sqrt(yPU[i]*yPU[i] + yErrorJetFlavor[i]*yErrorJetFlavor[i] +  yPES[i]*yPES[i] + yResImb[i]*yResImb[i] + yJEC[i]*yJEC[i]);
    cout<<endl<<"xData["<<i<<"] = "<<xData[i]<<endl;
    cout<<"yData["<<i<<"] = "<<yData[i]<<endl;
    cout<<"xSys["<<i<<"] = "<<xSys[i]<<endl;
    cout<<"ySys["<<i<<"] = "<<ySys[i]<<endl;
  }

   

  TGraphErrors* JERDataSys = new TGraphErrors(numEntries,xData,yData,xSys,ySys);
  TGraphErrors* JERMCSys   = new TGraphErrors(numEntries,xMC,yMC,xMCSys,yMCSys);

  

  // Plot resolution graph for one eta bin with sys error
  TCanvas *c = new TCanvas("c",etaString,200,10,500,500);
  c -> SetLeftMargin(0.17);
  c -> cd();
  
  //TMultiGraph* mg = new TMultiGraph();
  //mg -> SetTitle(etaString);


   //  Res_2012ABC->SetMarkerColor(kBlue-3);
   
  //Res_2012AB->SetXTitle("|#eta|");
  //Res_2012AB->GetXaxis()->SetRangeUser(0., 5.);
  //Res_2012AB->SetYTitle("Data/MC ratio (const fit)");
  //Res_2012AB->GetYaxis()->SetRangeUser(0.8, 1.5);
  //Res_2012AB->SetMarkerStyle(21);
  //Res_2012AB->SetMarkerSize(1.4);
  //Res_2012AB->SetLineColor(kPink-3);
  //Res_2012AB->SetMarkerColor(kPink-3);
  //Res_2012AB->Draw("e1p");

  if(PU && flav && PES && JEC && MC) etaString += " - all sys Uncertainties";
  if(PU && !flav && !PES && !JEC && !MC) etaString += " - only PU uncert.";
  if(!PU && flav && !PES && !JEC && !MC) etaString += " - only flavor uncert.";
  if(!PU && !flav && PES && !JEC && !MC) etaString += " - only PES uncert.";
  if(!PU && !flav && !PES && JEC && !MC) etaString += " - only JEC uncert.";
  if(!PU && !flav && !PES && !JEC && MC) etaString += " - only MC uncert.";

  JERData -> SetTitle(etaString);

  JERData -> GetXaxis() -> SetTitle("p_{T}^{#gamma}");
  JERData -> GetYaxis() -> SetTitleOffset(2.0); 
  JERData -> GetYaxis() -> SetTitle("JER");
  JERDataSys -> SetMinimum(-0.1);
  JERDataSys -> SetMaximum(0.3);   
  JERData -> SetMarkerStyle(20); 
  JERData->SetMarkerSize(1.4);
  JERData -> SetMinimum(-0.05);
  JERData -> SetMaximum(0.25);   
  //JERData -> SetMarkerColor(2); 
  //JERData -> SetLineColor(2); 
  JERDataSys->SetMarkerStyle(20);
  JERDataSys->SetMarkerSize(1.4);
  JERDataSys->SetFillColor(kGray);
  JERDataSys->SetFillStyle(3001);
  JERDataSys->SetLineColor(kGray);
  
  //JERDataSys -> Draw("a3");
  JERData -> Draw("AP"); 
  JERDataSys->DrawClone("e3p");
  JERData -> Draw("e1psame");  
  c -> SaveAs(pdfFile);

  //TLegend *legend  = new TLegend(0.6,0.8,0.9,0.9);
  //legend -> SetFillColor(0);
  //legend -> SetTextSize(0.033);
  //legend -> AddEntry(JERData,"Data","l");
  //legend -> AddEntry(JERMC,"MC","l");
  
  //mg -> Add(JERData);
  //mg -> Add(JERDataSys);
  //mg -> Draw("AP");
  
  //mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma}");
  //mg -> GetYaxis() -> SetTitleOffset(2.0); 
  //mg -> GetYaxis()->SetLimits(-0.4,0.3); 
  //mg -> GetXaxis()->SetLimits(0.,600.);
  
  
  //mg -> GetYaxis() -> SetTitle("JER");
  //mg -> SetMinimum(-0.1);
  //mg -> SetMaximum(0.3);   
  
  //legend -> Draw("same");
  //c -> SaveAs(pdfFile);
  
  //return 0;
  // --------------------------------  Ratio Plot for Resolution  ----------------------------------------
  
  double *xRatio = new double[numEntries];
  double *yRatio = new double[numEntries];
  double *xRatioError = new double[numEntries];
  double *yRatioError = new double[numEntries];
  

  // Get points in MC
  idxData = 0;
  idx = 0;
  while(idxData<numEntries){
    
    JERMC   -> GetPoint(idx,xAux1,yAux1);
    
    if(abs(xAux1-xData[idxData])/xData[idxData]<0.1){
      xMC[idxData] = xAux1;
      yMC[idxData] = yAux1;
      xMCError[idxData] = JERMC -> GetErrorX(idx);
      yMCError[idxData] = JERMC -> GetErrorY(idx);
      cout<<"xMC["<<idxData<<"] = "<<xMC[idxData]<<endl;
      cout<<"yMC["<<idxData<<"] = "<<yMC[idxData]<<endl;
      idxData += 1;
    }

    idx += 1;
    
  }

  cout<<endl;
   
 
  for(int i=0; i<numEntries; i++){

    cout<<endl<<"xMC["<<i<<"] = "<<xMC[i]<<endl;
    cout<<"yMC["<<i<<"] = "<<yMC[i]<<endl;
    cout<<"yMCError["<<i<<"] = "<<yMCError[i]<<endl;
    cout<<"xData["<<i<<"] = "<<xData[i]<<endl;
    cout<<"yData["<<i<<"] = "<<yData[i]<<endl;
    cout<<"yDataError["<<i<<"] = "<<yDataError[i]<<endl;
  
    xRatio[i] = 1./2.*(xData[i] + xMC[i]);
    yRatio[i] = yData[i]/yMC[i];
    

    // Statistical Error
    xRatioError[i] = 1./2.*TMath::Sqrt(TMath::Power(xDataError[i],2)+TMath::Power(xMCError[i],2));
    yRatioError[i] = TMath::Sqrt(TMath::Power((1./yMC[i]),2)*TMath::Power(yDataError[i],2)+TMath::Power((yData[i]/(TMath::Power(yMC[i],2))),2)*TMath::Power(yMCError[i],2));
    
  }

  TString ratioName;

  //TGraphErrors *Ratio = new TGraphErrors(10,axcomp,aycomp,axcompError,aycompError);
  TGraphErrors *Ratio = new TGraphErrors(numEntries,xRatio,yRatio,xRatioError,yRatioError);
  Ratio -> GetXaxis()->SetLimits(0,600);
  if(etaBin == 1 )     ratioName = "Ratio between Data and MC for |#eta|<0.5";
  else if(etaBin == 2) ratioName = "Ratio between Data and MC for 0.5<|#eta|<1.1";
  else if(etaBin == 3) ratioName = "Ratio between Data and MC for 1.1<|#eta|<1.7";
  else if(etaBin == 4) ratioName = "Ratio between Data and MC for 1.7<|#eta|<2.3";
  else if(etaBin == 5) ratioName = "Ratio between Data and MC for 2.3<|#eta|<5.0";
 
  Ratio -> SetTitle(ratioName); 
  Ratio -> GetXaxis() -> SetTitle("Photon pT");
  Ratio -> GetXaxis() -> SetTitleOffset(1.1); 
  Ratio -> GetYaxis() -> SetTitle("Ratio of JER (DATA/MC)");
  Ratio -> GetYaxis() -> SetTitleOffset(1.2);   

  TF1* f1 = new TF1("name","pol0",0,600);   
  Ratio -> Fit("name","");
  cout<<"ChiSquare = "<<f1 -> GetChisquare()<<endl;
  TLegend *legend  = 0;
  legend = new TLegend(0.55,0.8,0.9,0.9);
  legend -> SetFillColor(0);
  char legname[100];
  double fitPar = f1 -> GetParameter(0);
  sprintf(legname," %4.3f #pm %4.3f", f1 -> GetParameter(0), f1->GetParError(0));
  legend -> SetHeader(legname);
  TCanvas *c11 = new TCanvas("c11",ratioName,200,10,500,500);
  c11 -> cd();
  Ratio -> SetMinimum(0.5);
  Ratio -> SetMaximum(2.0);
  

  // Draw info boxes
  Ratio  -> Draw("AP"); 
  legend -> Draw("same");
  TLatex*  info   = new TLatex();
  info->SetTextFont(132);
  info-> SetNDC();
  info->SetTextSize(0.035);
  sprintf(legname,"#splitline{#chi^{2} = %4.2f}{dof = %i}",f1 -> GetChisquare(),f1 -> GetNDF());
  info->DrawLatex(0.22,0.84,legname);
  
  sprintf(pdfFile,"Ratio_Resolution_for_%d_eta_bin_PF_data_comparison.pdf",etaBin);
  c11 -> SaveAs(pdfFile);

  //delete JERMC;
  //delete JERData;
  //delete f1;
  //delete Ratio;
  //-------------------------------------------------- Ratio Plot with systematic Uncertainties -------------------------------------------------------------------
  
  // Make the fit with all points moved by systematic uncertainties up / down
  double *yRatioUp = new double[numEntries];
  double *yRatioLow = new double[numEntries];
  double *yRatioErrorUp = new double[numEntries];
  double *yRatioErrorLow = new double[numEntries];

  for(int i=0; i<numEntries; i++){

    yRatioUp[i] = (yData[i] + ySys[i])/yMC[i];
    yRatioLow[i] = (yData[i] - ySys[i])/yMC[i];

    // Statistical Error
    xRatioError[i] = 1./2.*TMath::Sqrt(TMath::Power(xDataError[i],2)+TMath::Power(xMCError[i],2));
    yRatioError[i] = TMath::Sqrt(TMath::Power((1./yMC[i]),2)*TMath::Power(yDataError[i],2)+TMath::Power((yData[i]/(TMath::Power(yMC[i],2))),2)*TMath::Power(yMCError[i],2));
    yRatioErrorUp[i] = TMath::Sqrt(TMath::Power((1./yMC[i]),2)*TMath::Power(yDataError[i],2)+TMath::Power(((yData[i]-ySys[i])/(TMath::Power(yMC[i],2))),2)*TMath::Power(yMCError[i],2));
    yRatioErrorLow[i] = TMath::Sqrt(TMath::Power((1./yMC[i]),2)*TMath::Power(yDataError[i],2)+TMath::Power(((yData[i]-ySys[i])/(TMath::Power(yMC[i],2))),2)*TMath::Power(yMCError[i],2));

  }

  TGraphErrors *RatioUp = new TGraphErrors(numEntries,xRatio,yRatioUp,xRatioError,yRatioErrorUp);
  TGraphErrors *RatioLow = new TGraphErrors(numEntries,xRatio,yRatioLow,xRatioError,yRatioErrorLow);
      
  RatioUp -> Fit("name","");
  double fitParUp = f1->GetParameter(0);
  cout<<"fitParUp = "<<fitParUp<<endl;
  RatioLow -> Fit("name","");
  double fitParLow = f1->GetParameter(0);
  cout<<"fitParLow = "<<fitParLow<<endl;

  double *x = new double[numEntries];
  double *y = new double[numEntries];
  double *exl = new double[numEntries];
  double *exh = new double[numEntries];
  double *eyl = new double[numEntries];
  double *eyh = new double[numEntries];

  for(int i=0; i<numEntries; i++){
    x[i] = xRatio[i];
    y[i] = fitPar;
    exl[i] = 0;
    exh[i] = 0;
    eyl[i] = fitPar-fitParLow;
    eyh[i] = fitParUp-fitPar;

  }
  cout<<endl<<"sys ERROR = "<<(fitParUp-fitParLow)/fitPar<<endl;

  TGraphAsymmErrors *sysBorders = new TGraphAsymmErrors(numEntries,x,y,exl,exh,eyl,eyh);
  TCanvas *c12 = new TCanvas("c12",ratioName,200,10,500,500);
  c12 -> cd();
  Ratio -> SetMinimum(0.6);
  Ratio -> SetMaximum(1.3);
  sysBorders ->SetMarkerStyle(20);
  sysBorders ->SetMarkerSize(1.4);
  sysBorders ->SetFillColor(kGray);
  sysBorders ->SetFillStyle(3001);
  sysBorders ->SetLineColor(kGray);

  Ratio  -> Draw("AP"); 
  sysBorders ->DrawClone("e3");  
 
  sprintf(pdfFile,"Ratio_Resolution_for_%d_eta_bin_PF_data_comparison_With_Sys_Uncert.pdf",etaBin);
  c12 -> SaveAs(pdfFile);

  return 0;
    */
}

//! Read TGraphErrors from file
// -------------------------------------------------------------------------------------
TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName) {
  TFile file(fileName,"READ");
  TGraphErrors *g = 0;
  file.GetObject(gName,g);
  if( g ) {
    if( newGName.Length() ) g->SetName(newGName);
  } else {
    std::cerr << "ERROR in FileOps::readTGraph: TGraphAsymmErrors with name '" << gName << "' does not exist in file '" << fileName << "'\n.";
    file.Close();
    exit(-1);
  }
  file.Close();
    
  return g;
} 

//! Read TGraph from file
// -------------------------------------------------------------------------------------
TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName) {
  TFile file(fileName,"READ");
  TGraph *g = 0;
  file.GetObject(gName,g);
  if( g ) {
    if( newGName.Length() ) {
    }  
    g->SetName(newGName);
  } else {
    std::cerr << "ERROR in FileOps::readTGraph: TGraph with name '" << gName << "' does not exist in file '" << fileName << "'\n.";
    file.Close();
    exit(-1);
  }
  file.Close();
    
  return g;
}

//! Read TGraphAsymmErrors from file
// -------------------------------------------------------------------------------------
TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName) {
  TFile file(fileName,"READ");
  TGraphAsymmErrors *g = 0;
  file.GetObject(gName,g);
  if( g ) {
    if( newGName.Length() ) g->SetName(newGName);
  } else {
    std::cerr << "ERROR in FileOps::readTGraph: TGraphAsymmErrors with name '" << gName << "' does not exist in file '" << fileName << "'\n.";
    file.Close();
    exit(-1);
  }
  file.Close();
    
  return g;
}



/*
  namespace util
  {
  class FileOps
  {
  public:
  static TF1* readTF1(const TString &fileName, const TString &fName, const TString &newFName = "");
  static TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName = "");
  static TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName = "");
  static TH2* readTH2(const TString &fileName, const TString &hName, const TString &newHName = "", bool useCurrentStyle = true);
  static TH1* readTH1(const TString &fileName, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  static util::HistVec readTH1(const std::vector<TString> &fileNames, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  static util::HistVec readHistVec(const TString &fileName, const std::vector<TString> &histNames, const TString &newHistNameSuffix = "", bool useCurrentStyle = true);
  static util::HistVec readHistVec(const TString &fileName, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  static std::vector<TF1*> readTF1Vec(const TString &fileName, const TString &fName, const TString &newFName = "");
  static THStack* readTHStack(const TString &fileName, const TString &stackName, const TString &newStackName = "");
  static TChain* createTChain(const TString &fileName, const TString &treeName = "", unsigned int verbosity = 1);
  };


  // -------------------------------------------------------------------------------------
  TH2* FileOps::readTH2(const TString &fileName, const TString &hName, const TString &newHName, bool useCurrentStyle) {
  TFile file(fileName,"READ");
  TH2* h = 0;
  file.GetObject(hName,h);
  if( h ) {
  h->SetDirectory(0);
  if( useCurrentStyle ) h->UseCurrentStyle();
  if( newHName.Length() ) h->SetName(newHName);
  } else {
  std::cerr << "ERROR in FileOps::readTH2: No TH2 with name '" << hName << "' in file '" << fileName << "'\n.";
  file.Close();
  exit(-1);
  }
  file.Close();
    
  return h;
  }


  //! Read TH1 histogram from file
  // -------------------------------------------------------------------------------------
  TH1* FileOps::readTH1(const TString &fileName, const TString &histName, const TString &newHistName, bool useCurrentStyle) {
  TFile file(fileName,"READ");
  TH1 *h = 0;
  file.GetObject(histName,h);
  if( h ) {
  h->SetDirectory(0);
  if( useCurrentStyle ) h->UseCurrentStyle();
  if( newHistName.Length() ) h->SetName(newHistName);
  } else {
  std::cerr << "ERROR in FileOps::readTH1: Histogram with name '" << histName << "' does not exist in file '" << fileName << "'\n.";
  file.Close();
  exit(-1);
  }
  file.Close();
    
  return h;
  }


  //! Read THStack from file
  // -------------------------------------------------------------------------------------
  THStack* FileOps::readTHStack(const TString &fileName, const TString &stackName, const TString &newStackName) {
  TFile file(fileName,"READ");
  THStack *s = 0;
  file.GetObject(stackName,s);
  if( s ) {
  //      s->SetDirectory(0);
  if( newStackName.Length() ) s->SetName(newStackName);
  } else {
  std::cerr << "ERROR in FileOps::readTHStack: THStack with name '" << stackName << "' does not exist in file '" << fileName << "'\n.";
  file.Close();
  exit(-1);
  }
  file.Close();
    
  return s;
  }



  //! Read TGraph from file
  // -------------------------------------------------------------------------------------
  TGraph* FileOps::readTGraph(const TString &fileName, const TString &gName, const TString &newGName) {
  TFile file(fileName,"READ");
  TGraph *g = 0;
  file.GetObject(gName,g);
  if( g ) {
  if( newGName.Length() ) g->SetName(newGName);
  } else {
  std::cerr << "ERROR in FileOps::readTGraph: TGraph with name '" << gName << "' does not exist in file '" << fileName << "'\n.";
  file.Close();
  exit(-1);
  }
  file.Close();
    
  return g;
  }
  
  
  
  
  
  //! Read TH1 histograms from different files
  // -------------------------------------------------------------------------------------
  util::HistVec FileOps::readTH1(const std::vector<TString> &fileNames, const TString &histName, const TString &newHistName, bool useCurrentStyle) {
  util::HistVec v(fileNames.size());
  for(unsigned int i = 0; i < fileNames.size(); ++i) {
  v[i] = readTH1(fileNames[i],histName,newHistName+util::toTString(i),useCurrentStyle);
  }
  std::cout << "Done\n";
    
  return v;
  }
  

  //! Read TH1 histograms from one file
  // -------------------------------------------------------------------------------------
  util::HistVec FileOps::readHistVec(const TString &fileName, const std::vector<TString> &histNames, const TString &newHistNameSuffix, bool useCurrentStyle) {
  util::HistVec v;
  TFile file(fileName,"READ");
  for(std::vector<TString>::const_iterator it = histNames.begin();
  it != histNames.end(); ++it) {
  TH1 *h = 0;
  file.GetObject(*it,h);
  if( h ) {
  h->SetDirectory(0);
  if( useCurrentStyle ) h->UseCurrentStyle();
  h->SetName((*it)+newHistNameSuffix);
  v.push_back(h);
  } else {
  std::cerr << "ERROR in FileOps::readHistVec: Histogram with name '" << *it << "' does not exist in file '" << fileName << "'\n.";
  file.Close();
  exit(-1);
  }
  }
  file.Close();
    
  return v;
  }


  
  //! Read TH1 histograms from one file
  // -------------------------------------------------------------------------------------
  util::HistVec FileOps::readHistVec(const TString &fileName, const TString &histName, const TString &newHistName, bool useCurrentStyle) {
  util::HistVec v;
  TFile file(fileName,"READ");
  bool binExists = true;
  unsigned int bin = 0;
  while( binExists ) {
  TH1 *h = 0;
  file.GetObject(histName+util::toTString(bin),h);
  if( h ) {
  h->SetDirectory(0);
  if( useCurrentStyle ) h->UseCurrentStyle();
  if( newHistName.Length() ) h->SetName(newHistName+util::toTString(bin));
  v.push_back(h);
  ++bin;
  } else {
  binExists = false;
  }
  }
  file.Close();
    
  if( v.size() == 0 ) std::cerr << "WARNING in util::FileOps::readHistVec(): No histogram read!\n";
    
  return v;
  }


  // -------------------------------------------------------------------------------------
  TF1* FileOps::readTF1(const TString &fileName, const TString &fName, const TString &newFName) {
  TFile file(fileName,"READ");
  TF1 *f = 0;
  file.GetObject(fName,f);
  if( f ) {
  //f->SetDirectory(0);
  if( newFName.Length() ) f->SetName(newFName);
  } else {
  std::cerr << "ERROR in FileOps::readTF1: TF1 with name '" << fName << "' does not exist in file '" << fileName << "'\n.";
  }
  file.Close();
    
  return f;
  }


  // -------------------------------------------------------------------------------------
  std::vector<TF1*> FileOps::readTF1Vec(const TString &fileName, const TString &fName, const TString &newFName) {
  std::vector<TF1*> v;
  TFile file(fileName,"READ");
  bool binExists = true;
  unsigned int bin = 0;
  while( binExists ) {
  TF1* f = 0;
  file.GetObject(fName+util::toTString(bin),f);
  if( f ) {
  //f->SetDirectory(0);
  if( newFName.Length() ) f->SetName(newFName+util::toTString(bin));
  v.push_back(f);
  ++bin;
  } else {
  binExists = false;
  }
  }
  file.Close();
    
  if( v.size() == 0 ) std::cerr << "WARNING in util::FileOps::readTF1Vec(): No TF1 read!\n";
    
  return v;
  }



  //! Create TChain from input root files. The root
  //! files are expected to contain a TTree "DiJetTree".
  //! There are two possible input options:
  //!
  //! 1) 'fileName' specifies a single root file; it ends
  //!    with '.root';
  //! 2) 'fileName' contains a list of root file names.
  // --------------------------------------------------
  TChain* FileOps::createTChain(const TString &fileName, const TString &treeName, unsigned int verbosity) {
  if( verbosity >= 1 ) std::cout << "Creating TChain" << std::endl;

  TString tree = (treeName=="" ? "DiJetTree" : treeName);
  TChain* chain = new TChain(tree); 
    
  // Option 1: single root file
  if( fileName.EndsWith(".root") ) {
  if( verbosity >= 1 ) std::cout << "  Adding '" << tree << "' from file '" << fileName << "'" << std::endl;
  chain->Add(fileName);
  }
  // Option 2: list of root files
  else {
  if( verbosity >= 1 ) std::cout << "  Opening files from list '" << fileName << "'" << std::endl;
  std::ifstream filelist;
  filelist.open(fileName);
  int nOpenedFiles = 0;
  if( filelist.is_open() ) {
  TString name = "";
  while( !filelist.eof() ) {
  filelist >> name;
  if( filelist.eof() ) break;
  if( verbosity >= 1 ) std::cout << "  Adding '" << tree << "' from file '" << name << "'" << std::endl;
  chain->Add(name);
  nOpenedFiles++;
  }
  } else {
  std::cerr << "ERROR opening file '" << fileName << "'\n";
  exit(1);
  }
  filelist.close();
  }
    
  return chain;
  }
  }
*/







