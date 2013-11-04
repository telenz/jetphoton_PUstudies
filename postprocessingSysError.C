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

int postprocessing(int etaBin){

  // For looking at different systematic uncertainties independently
  bool flav = false;
  bool PU = false;
  bool PES = false;
  bool MC = false;
  bool JEC = true;

  TString etaString;   
  char pdfFile[100];

  TString path = "plots_2012/PF_L1FastJet/";
  TString rootFiles;
  
  
  // Read the data files and all systematic error files
  rootFiles.Form("data/root_files/Resolution_for_%i_eta_bin_PF_data.root",etaBin);
  TGraphErrors* JERData = readTGraphErrors(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("mc/root_files/Resolution_for_%i_eta_bin_PF_mc.root",etaBin);
  TGraphErrors* JERMC = readTGraphErrors(path+rootFiles,"Graph;1","Graph;1");

  path = "plots_2012/PF_L1FastJet/systematicUncertainties/";
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
  rootFiles.Form("PESUncertainties/root_files_lower_bound/Resolution_for_%i_eta_bin_PF_data.root",etaBin);  
  TGraph* sysPESlow = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("PESUncertainties/root_files_upper_bound/Resolution_for_%i_eta_bin_PF_data.root",etaBin);  
  TGraph* sysPESup = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("ResImbUncertainties/root_files_lower_bound/Resolution_for_%i_eta_bin_PF_data.root",etaBin);  
  TGraph* sysMClow = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("ResImbUncertainties/root_files_upper_bound/Resolution_for_%i_eta_bin_PF_data.root",etaBin);    
  TGraph* sysMCup = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("JECUncertainties/root_files_mc_lower_bound/Resolution_for_%i_eta_bin_PF_mc.root",etaBin);    
  TGraph* sysJEClow = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  rootFiles.Form("JECUncertainties/root_files_mc_upper_bound/Resolution_for_%i_eta_bin_PF_mc.root",etaBin);    
  TGraph* sysJECup = readTGraph(path+rootFiles,"Graph;1","Graph;1");
  
  
  sprintf(pdfFile,"Resolution_for_%d_eta_bin_PF_data_mc_comparison.pdf",etaBin);
  if(etaBin == 1)      etaString = "Jet Energy Resolution for |#eta|<0.5";
  else if(etaBin == 2) etaString = "Jet Energy Resolution for 0.5<|#eta|<1.1";
  else if(etaBin == 3) etaString = "Jet Energy Resolution for 1.1<|#eta|<1.7";
  else if(etaBin == 4) etaString = "Jet Energy Resolution for 1.7<|#eta|<2.3";
  else if(etaBin == 5) etaString = "Jet Energy Resolution for 2.3<|#eta|<5.0";
  
  
  //TGraphErrors* JERMC = readTGraphErrors(rootFileMC,"Graph;1","Graph;1");  
  //JERMC -> SetMarkerColor(3);
  //JERMC -> SetLineColor(3);
  //TCanvas* c1 = new TCanvas("c1","c1",100,100,500,500);
  //c1->cd();
  //JERMC -> Draw("AP");
  
  //JERData -> SetMarkerColor(4);
  //JERData -> SetLineColor(4);
  
  // Evaluate systematic jetFlavor uncertainty for every pt point  
  int numEntries = JERData->GetN();
  cout<<endl<<"Number of entries in JERData = "<<numEntries<<endl<<endl;
  
  double auxL=1000;
  double auxG=0;
  double *xData = new double[numEntries];
  double *yData = new double[numEntries];
  double *xDataError = new double[numEntries];
  double *yDataError = new double[numEntries];
  
  double xAux1 = 0;
  double xAux2 = 0;
  double yAux1 = 0;
  double yAux2 = 0;
   
  // loop over pt bins for every systematic uncertainty
  int idxData = 0;
  for(int i=0; i<numEntries; i++){

    JERData             -> GetPoint(i,xData[idxData],yData[idxData]);
    if(yData[idxData]<0.001){
      cout<<"For "<<i<<". pT Bin Resolution was too small to take into account!!"<<endl;
      continue;
    }
    xDataError[idxData] = JERData -> GetErrorX(i);
    yDataError[idxData] = JERData -> GetErrorY(i);
    cout<<"xData["<<idxData<<"] = "<<xData[idxData]<<endl;
    cout<<"yData["<<idxData<<"] = "<<yData[idxData]<<endl;
    idxData += 1;
  }
  cout<<endl;

  numEntries = idxData;
  cout<<"New number of Entries is now "<<numEntries<<endl;

  delete JERData;
  JERData = new TGraphErrors(numEntries,xData,yData,xDataError,yDataError);

  double *xMC = new double[numEntries];
  double *yMC = new double[numEntries];
  double *xMCError = new double[numEntries];
  double *yMCError = new double[numEntries];
  
  double *yPU = new double[numEntries];
  double *yResImb = new double[numEntries];
  double *yJEC = new double[numEntries];
  double *yPES = new double[numEntries];
  double *yErrorJetFlavor = new double[numEntries];

  double *xSys = new double[numEntries];
  double *ySys = new double[numEntries];
  double *xMCSys = new double[numEntries];
  double *yMCSys = new double[numEntries];
  double xAux[4] = {0};
  double yAux[4] = {0};  
  
  for(int i=0;i<numEntries;i++){
    xMC[i] = 0;
    yMC[i] = 0;
    xMCError[i] = 0;
    yMCError[i] = 0;

    yPU[i] = 0;
    yResImb[i] = 0;
    yJEC[i] = 0;
    yPES[i] = 0;
    yErrorJetFlavor[i] = 0;
    xSys[i] = 0;
    ySys[i] = 0;
  }
  //return 0;
  // 1.) Pileup: Find the largest and smallest value of PU uncertainties for every pt bin

  int idx0 = 0;
  int idx1 = 0;
  int idx2 = 0;
  int idx3 = 0;
  idxData = 0;  
  int itr = 0;
  bool test;
  if(PU){
  while(idxData<numEntries && itr<30){
    itr += 1;
    cout<<"itr = "<<itr<<endl;
    yPU[idxData] = 0;
    cout<<"idxData = "<<idxData<<endl;
    cout<<"idx0 = "<<idx0<<endl;
    cout<<"idx1 = "<<idx1<<endl;
    cout<<"idx2 = "<<idx2<<endl;
    cout<<"idx3 = "<<idx3<<endl;
    cout<<"N = "<<sysJetPUle10   -> GetN()<<endl;
    cout<<"N = "<<sysJetPUgt10   -> GetN()<<endl;
    cout<<"N = "<<sysJetPUgt15   -> GetN()<<endl;
    cout<<"N = "<<sysJetPUgt20   -> GetN()<<endl;
    sysJetPUle10   -> GetPoint(idx0,xAux[0],yAux[0]);
    sysJetPUgt10   -> GetPoint(idx1,xAux[1],yAux[1]);
    sysJetPUgt15   -> GetPoint(idx2,xAux[2],yAux[2]);
    sysJetPUgt20   -> GetPoint(idx3,xAux[3],yAux[3]);

    cout<<"xAux[0] = "<<xAux[0]<<endl;
    cout<<"xAux[1] = "<<xAux[1]<<endl;
    cout<<"xAux[2] = "<<xAux[2]<<endl;
    cout<<"xAux[3] = "<<xAux[3]<<endl;
    cout<<"xData["<<idxData<<"] = "<<xData[idxData]<<endl;

    auxL=1000;
    auxG=0;
    test = false;

    double Pileup[4] = {5,12.5,17.5,25};
    double JERpu[4] = {0};
    TCanvas *c = new TCanvas("c",etaString,200,10,500,500);
 
    for(int j=0; j<4; j++){
      cout<<"j = "<<j<<endl;
      if(abs(xAux[0]-xData[idxData])/xData[idxData]>0.1 && xAux[0]<xData[idxData] && idxData<5){
	idx0 += 1;
	break;
      }
      if(abs(xAux[1]-xData[idxData])/xData[idxData]>0.1 && xAux[1]<xData[idxData]  && idxData<5){
	idx1 += 1;
	break;
      }
      
      if(abs(xAux[2]-xData[idxData])/xData[idxData]>0.1 && xAux[2]<xData[idxData]  && idxData<5){
	idx2 += 1;
	break;
      }
      if(abs(xAux[3]-xData[idxData])/xData[idxData]>0.1 && xAux[3]<xData[idxData] && idxData<5){
      	cout<<"IN"<<endl;
      idx3 += 1;
      break;
      }
      
      
      if(abs(xAux[j]-xData[idxData])/xData[idxData]<0.1){
	

	if(yAux[j]>auxG) auxG = yAux[j]; 
	if(yAux[j]<auxL) auxL = yAux[j];
	if(j==0) idx0 += 1;
	else if(j==1) idx1 += 1;
	else if(j==2) idx2 += 1;
	else if(j==3) idx3 += 1;
	//if(idxData==0){
	  cout<<"j = "<<j<<endl;
	  cout<<"yAux[j] = "<<yAux[j]<<endl;
	  cout<<"auxG = "<<auxG<<endl;
	  cout<<"auxL = "<<auxL<<endl;
	  //}  
      }
      if(j==3) test = true;
      
    }
    
    if(auxG != 0 && auxL != 1000 && test){
      yPU[idxData] = abs(auxG - auxL)/2;
      cout<<"yPU["<<idxData<<"] = "<<yPU[idxData]<<endl;
      //cout<<"yPU["<<idxData<<"] = "<<yPU[idxData]<<endl;
      idxData += 1;
    }

    

  }

  if(itr>=30){
    cout<<"PU not correct !!"<<endl;
    return 0;

  }
  cout<<endl;
  }

  //yPU[0] = 0;
  //yPU[1] = 0;
  //yPU[2] = 0;
  //yPU[3] = 0;
  //yPU[4] = 0;
  //yPU[5] = 0;
  //yPU[6] = 0;
  //yPU[7] = 0;
  //yPU[8] = 0;
  //yPU[9] = 0;

  // 2.) Flavor uncertainty
  int idx=0;
  idxData = 0;
  itr=0;
  if(flav){
  while(idxData<numEntries && itr<20){
    itr += 1;
    sysJetFlavorQuark   -> GetPoint(idx,xAux1,yAux1);
    sysJetFlavorGluon   -> GetPoint(idx,xAux2,yAux2);
    cout<<"xQuark["<<idx<<"] = "<<xAux1<<endl;
    cout<<"xGluon["<<idx<<"] = "<<xAux2<<endl;
    cout<<"xData["<<idxData<<"] = "<<xData[idxData]<<endl;
    idx += 1;
    if(abs(xAux1-xAux2)/xAux1>0.1){
      cout<<"Error flav sys!!"<<endl;
      return 0;
    }
    if(abs(xAux1-xData[idxData])/xData[idxData]<0.1){
      yErrorJetFlavor[idxData] = abs(yAux1-yAux2)/2;
      cout<<"yErrorJetFlavor["<<idxData<<"] = "<<yErrorJetFlavor[idxData]<<endl;
      idxData += 1;
    }
    else{
      continue;
    }

    //cout<<"yData["<<i<<"] = "<<yData[i]<<endl;
    //cout<<"yQuark["<<i<<"] = "<<yQuark[i]<<endl;
    //cout<<"yGluon["<<i<<"] = "<<yGluon[i]<<endl;
  }
  //return 0;
  
  cout<<endl;
  }
    
  // 3.) Symmetrize PES uncertanty
  idx = 0;
  idxData = 0;
  if(PES){
  while(idxData<numEntries){
    sysPESup  -> GetPoint(idx,xAux1,yAux1);
    sysPESlow -> GetPoint(idx,xAux2,yAux2);
    cout<<"xPESup["<<idx<<"] = "<<xAux1<<endl;
    cout<<"xPESlow["<<idx<<"] = "<<xAux2<<endl;
    idx += 1;
    if(abs(xAux1-xAux2)/xAux1>0.1){
      cout<<"Error PES sys!!"<<endl;
      return 0;
    }
    if(abs(xAux1-xData[idxData])/xData[idxData]<0.1){
      yPES[idxData] = abs(yAux1-yAux2)/2.;
      cout<<"yPES["<<idxData<<"] = "<<yPES[idxData]<<endl;
      idxData += 1;
    }
    else{
      continue;
    }
    
  }
  cout<<endl;
  }

  // 4.) Symmetrize MC uncertainty
  if(MC){
  idx=0;
  idxData = 0;
  while(idxData<numEntries){
    sysMCup -> GetPoint(idx,xAux1,yAux1);
    sysMClow -> GetPoint(idx,xAux2,yAux2);
    cout<<"xMCup["<<idx<<"] = "<<xAux1<<endl;
    cout<<"xMClow["<<idx<<"] = "<<xAux2<<endl;
    cout<<"xData["<<idxData<<"] = "<<xData[idxData]<<endl;
    idx += 1;
    if(abs(xAux1-xAux2)/xAux1>0.1){
      cout<<"Error MC sys!!"<<endl;
      return 0;
    }
    if(abs(xAux1-xData[idxData])/xData[idxData]<0.1){
      yResImb[idxData] = abs(yAux1-yAux2)/2.;
      cout<<"yResImb["<<idxData<<"] = "<<yResImb[idxData]<<endl;
      idxData += 1;
    }
    else{
      
      cout<<"Take care in MC!!!"<<endl;
      continue;
    }
  }
  cout<<endl;
  }
  

  // 5.) Symmetrize JEC uncertanty
  
  idxData=0;
  idx1 = 0;
  idx2 = 0;
  if(JEC){
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







