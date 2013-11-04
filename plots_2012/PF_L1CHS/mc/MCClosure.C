
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
#include "TH2.h"
#include "TObject.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TColorWheel.h"
#include "TLatex.h"

TCanvas* DrawComparison(TH1D* prediction, TH1D* selection, TString Title, TString LumiTitle, TString xTitle, bool isData);

int makeGraphData(){

 
  
  TGraphErrors* histoJet, *histoPhoton;
  
  
  TCanvas* canvas[13];
  TFile *file;
  TString fileName;

  for(int i=4; i<13; i++){

    fileName.Form("c%i",i);
    canvas[i] = new TCanvas(fileName,fileName,0,0,500,500);

    canvas[i] ->cd();
    fileName.Form("root_files_JetHemisphere/jet_energy_resolution_for_1_eta_bin_%i_pTGamma_bin_PFCHS_mc.root",i);
    file = TFile::Open(fileName);
    file->GetObject("Graph",histoJet);
    histoJet -> SetMarkerColor(2);
    histoJet -> SetLineColor(2);
    histoJet -> GetFunction("fResolutionAlpha")->SetLineColor(2);
    histoJet ->Draw("AP"); 
    delete file;
    
    fileName.Form("root_files/jet_energy_resolution_for_1_eta_bin_%i_pTGamma_bin_PFCHS_mc.root",i);
    file = TFile::Open(fileName);
    file->GetObject("Graph",histoPhoton);
    histoPhoton -> SetMarkerColor(3);
    histoPhoton -> SetLineColor(3);
    histoPhoton -> GetFunction("fResolutionAlpha")->SetLineColor(3);
    histoPhoton ->Draw("Psame");
    delete file;
  }


  return 0;

}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int MCClosure(){

 
  
  TGraphErrors* intrinsic, *full;
  TGraphErrors* Ratio;
  TGraphErrors* RatioAbs;
  
  
  TCanvas* canvas[4];
  TCanvas* canvas2[4];
  TCanvas* canvas3[4];
  TFile *file;
  TString fileName;

  const int eta_int = 4;
  double etaBorders[eta_int+1] = {0};
  TString etaRegion; 
  etaBorders[0] = 0.0;
  etaBorders[1] = 0.5;
  etaBorders[2] = 1.1; 
  etaBorders[3] = 1.7; 
  etaBorders[4] = 2.3; 

  for(int j=0; j<4; j++){

    if(j == 0) etaRegion.Form("|#eta| < %4.1f", etaBorders[j+1]);
    else       etaRegion.Form("%4.1f < #eta < %4.1f",etaBorders[j], etaBorders[j+1]);

    fileName.Form("c%i",j);
    canvas[j] = new TCanvas(fileName,fileName,0,0,500,500);

    canvas[j] ->cd();

    //char file1[100] = "_JetHemisphere_widthDivMean";
    char file1[100] = "";
    char file2[100] = "";

    fileName.Form("root_files%s/Resolution_for_%i_eta_bin_PFCHS_mc.root",file1,j+1);
    //fileName.Form("root_files/Scale_for_%i_eta_bin_PFCHS_mc.root",j+1);
    file = TFile::Open(fileName);
    file->GetObject("Graph",full);
    full -> SetMarkerColor(3);
    full -> SetLineColor(3);
    full -> GetFunction("fResolution")->SetLineColor(3);
    full ->SetMaximum(0.12);
    full ->SetMinimum(0.04);
    full ->Draw("AP");
    delete file;


    fileName.Form("root_files%s/Resolution_for_%i_eta_bin_intrinsic_PFCHS_mc.root",file1,j+1);
    //fileName.Form("root_files/Scale_for_%i_eta_bin_intrinsic_PFCHS_mc.root",j+1);
    file = TFile::Open(fileName);
    file->GetObject("Graph",intrinsic);
    intrinsic -> SetMarkerColor(2);
    intrinsic -> SetLineColor(2);
    
    intrinsic -> GetFunction("fResolution")->SetLineColor(2);
    intrinsic ->Draw("Psame"); 
    
   
    TLegend *legend  = new TLegend(0.6,0.7,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.033);
    legend -> AddEntry(intrinsic,"intrinsic","l");
    legend -> AddEntry(full,"full resolution","l");

    legend->Draw("same");

    fileName.Form("Resolution_intrinsic_full_for_%i_eta_bin%s.pdf",j+1,file1);
    canvas[j]->SaveAs(fileName);



    delete file;

    const int numEntries = 9;
    double xIntrinsic[numEntries] = {0};
    double yIntrinsic[numEntries] = {0};
    double xIntrinsicError[numEntries] = {0};
    double yIntrinsicError[numEntries] = {0};
    double xFull[numEntries] = {0};
    double yFull[numEntries] = {0};
    double xFullError[numEntries] = {0};
    double yFullError[numEntries] = {0};

    double x[numEntries] = {0};
    double y[numEntries] = {0};
    double xError[numEntries] = {0};
    double yError[numEntries] = {0};

    for(int i=0; i<numEntries; i++){
      intrinsic   -> GetPoint(i+0,xIntrinsic[i],yIntrinsic[i]);
      full        -> GetPoint(i+0,xFull[i],yFull[i]);

      cout<<"xFull["<<i<<"] = "<<xFull[i]<<endl;
      cout<<"xIntrinsic["<<i<<"] = "<<xIntrinsic[i]<<endl;

      xIntrinsicError[i] = intrinsic -> GetErrorX(i+0);
      yIntrinsicError[i] = intrinsic -> GetErrorY(i+0);
      xFullError[i]      = full      -> GetErrorX(i+0);
      yFullError[i]      = full      -> GetErrorY(i+0);
      
      y[i] = yFull[i]/yIntrinsic[i] -1 ;
      x[i] = 1./2.*(xFull[i]+xIntrinsic[i]);
      
      yError[i] =  TMath::Sqrt(TMath::Power((1./yIntrinsic[i]),2)*TMath::Power(yFullError[i],2)+TMath::Power((yFull[i]/(TMath::Power(yIntrinsic[i],2))),2)*TMath::Power(yIntrinsicError[i],2));
      //cout<<"y["<<i<<"] = "<<y[i]*100<<endl<<endl;
      //cout<<"yError["<<i<<"] = "<<yError[i]*100<<endl<<endl;

      xError[i] = 1./2.*TMath::Sqrt(TMath::Power(xFullError[i],2)+TMath::Power(xIntrinsicError[i],2));
    }

    fileName.Form("c2%i",j);
    canvas2[j] = new TCanvas(fileName,fileName,0,0,500,500);
    canvas2[j] ->SetBottomMargin(0.15);
    canvas2[j] ->SetLeftMargin(0.17);
    canvas2[j] ->cd();

    Ratio = new TGraphErrors(numEntries,x,y,xError,yError);

    
    char titleName[100] = {0};
    sprintf(titleName,"MC Closure Test  %s",file2);
    Ratio->SetTitle(titleName);
    Ratio->GetXaxis()->SetTitleSize(0.06);
    Ratio->GetYaxis()->SetTitleSize(0.04);
    Ratio-> GetXaxis()->SetTitle("p_{T}^{#gamma}");
    Ratio -> GetXaxis()->SetLimits(0,600);
    Ratio-> GetYaxis()->SetTitle("#frac{predic. intrinsic}{intrinsic} -1");
    Ratio->GetYaxis()->SetTitleOffset(1.5); 
    Ratio->SetMarkerStyle(20);
    Ratio->SetMarkerSize(1.1);
    Ratio ->SetMaximum(0.3);
    Ratio ->SetMinimum(-0.3);
    Ratio->Draw("AP");

    TF1* line = new TF1("line","0.",0,600);
    
    line->SetLineColor(2);
    line->Draw("same");

    TLatex*  info   = new TLatex();
    char legname[100] = {0};
    info->SetTextFont(132);
    info-> SetNDC();
    info->SetTextSize(0.040);
    info->DrawLatex(0.60,0.80,etaRegion);

    fileName.Form("MCClosure_for_%i_eta_bin%s.pdf",j+1,file1);
    canvas2[j]->SaveAs(fileName);


    canvas3[j] = new TCanvas(fileName,fileName,0,0,500,500);
    canvas3[j] ->cd();

    for(int i=0; i<numEntries; i++) y[i] = std::abs(y[i]);
    RatioAbs = new TGraphErrors(numEntries,x,y,xError,yError);
  
    TF1* fit = new TF1("fit","pol0",0,600);
    //TF1* fit = new TF1("fit","expo(0)",0,600);
    sprintf(titleName,"MCClosure - %s absolute Values",file2);
    RatioAbs -> SetTitle(titleName);
    RatioAbs -> GetXaxis()->SetTitle("p_{T}^{#gamma}");
    RatioAbs -> GetXaxis()->SetLimits(0,600);
    RatioAbs -> GetYaxis()->SetTitle("JER/intrinsic -1");
    RatioAbs -> GetYaxis()->SetTitleOffset(1.5); 
    RatioAbs -> SetMarkerStyle(20);
    RatioAbs -> SetMarkerSize(1.1);
    RatioAbs -> SetMaximum(0.25);
    RatioAbs -> SetMinimum(-0.05);
    RatioAbs -> Fit("fit","Q");
    RatioAbs -> Draw("AP");

    line ->Draw("same");
    fit -> SetLineColor(3);
    fit -> Draw("esame");

    delete info;
    info   = new TLatex();   
    info->SetTextFont(132);
    info-> SetNDC();
    info->SetTextSize(0.040);
    sprintf(legname,"%i. eta Bin",j+1);
    info->DrawLatex(0.60,0.80,legname);

    
    sprintf(legname,"fit = %4.3f #pm %4.3f",fit->GetParameter(0),fit->GetParError(0));
    info -> DrawLatex(0.60,0.70,legname);
    info -> Draw("same");

    sprintf(legname,"Chi^2 = %4.3f",fit->GetChisquare());
    info -> DrawLatex(0.60,0.63,legname);
    info -> Draw("same");
    sprintf(legname,"ndof = %i",fit->GetNDF());
    info -> DrawLatex(0.60,0.56,legname);
    info -> Draw("same");

    fileName.Form("MCClosure_for_%i_eta_bin%s_absoluteValues.pdf",j+1,file1);
    canvas2[j]->SaveAs(fileName);

  }

  return 0;

}

// -------------------------------------------------------
int plot(){

 
  
  TGraphErrors* intrinsic, *full;
  TGraphErrors* Ratio;
  
  
  TCanvas* canvas[20];
  TCanvas* canvas2[4];
  TFile *file;
  TString fileName;

  for(int j=10; j<12; j++){

    fileName.Form("c%i",j);
    canvas[j] = new TCanvas(fileName,fileName,0,0,500,500);

    canvas[j] ->cd();

    char file1[100] = "";
    char file2[100] = "";

    fileName.Form("root_files/jet_energy_resolution_for_1_eta_bin_%i_pTGamma_bin_imbalance_PFCHS_mc.root",j+1);
    //fileName.Form("root_files/Scale_for_%i_eta_bin_PFCHS_mc.root",j+1);
    file = TFile::Open(fileName);
    file->GetObject("Graph",full);
    full -> SetMarkerColor(1);
    full -> SetLineColor(1);
    full -> GetFunction("fResolutionAlpha")->SetLineColor(3);
    full ->SetMaximum(0.12);
    full ->SetMinimum(0.00);
    full ->Draw("AP");
    delete file;

    TLatex*  info   = new TLatex();
    char legname[100] = {0};
    info->SetTextFont(132);
    info-> SetNDC();
    info->SetTextSize(0.040);
    sprintf(legname,"#Chi^{2} = %f",full->GetFunction("fResolutionAlpha")->GetChisquare());
    info->DrawLatex(0.60,0.80,legname);

    delete info; 

    info   = new TLatex();
    
    info->SetTextFont(132);
    info-> SetNDC();
    info->SetTextSize(0.040);
    sprintf(legname,"ndof = %i",full->GetFunction("fResolutionAlpha")->GetNDF());
    info->DrawLatex(0.60,0.70,legname);
   
   
    

    fileName.Form("jet_energy_resolution_imbalance_for_%i_pt_bin.pdf",j+1);
    canvas[j]->SaveAs(fileName);



    

    /*

    const int numEntries = 9;
    double xIntrinsic[numEntries] = {0};
    double yIntrinsic[numEntries] = {0};
    double xIntrinsicError[numEntries] = {0};
    double yIntrinsicError[numEntries] = {0};
    double xFull[numEntries] = {0};
    double yFull[numEntries] = {0};
    double xFullError[numEntries] = {0};
    double yFullError[numEntries] = {0};

    double x[numEntries] = {0};
    double y[numEntries] = {0};
    double xError[numEntries] = {0};
    double yError[numEntries] = {0};

    for(int i=0; i<numEntries; i++){
      intrinsic   -> GetPoint(i+0,xIntrinsic[i],yIntrinsic[i]);
      full        -> GetPoint(i+0,xFull[i],yFull[i]);

      cout<<"xFull["<<i<<"] = "<<xFull[i]<<endl;
      cout<<"xIntrinsic["<<i<<"] = "<<xIntrinsic[i]<<endl;

      xIntrinsicError[i] = intrinsic -> GetErrorX(i+0);
      yIntrinsicError[i] = intrinsic -> GetErrorY(i+0);
      xFullError[i]      = full      -> GetErrorX(i+0);
      yFullError[i]      = full      -> GetErrorY(i+0);
      
      y[i] = yFull[i]/yIntrinsic[i] -1 ;
      x[i] = 1./2.*(xFull[i]+xIntrinsic[i]);
      
      yError[i] =  TMath::Sqrt(TMath::Power((1./yIntrinsic[i]),2)*TMath::Power(yFullError[i],2)+TMath::Power((yFull[i]/(TMath::Power(yIntrinsic[i],2))),2)*TMath::Power(yIntrinsicError[i],2));
      cout<<"y["<<i<<"] = "<<y[i]*100<<endl<<endl;
      cout<<"yError["<<i<<"] = "<<yError[i]*100<<endl<<endl;

      xError[i] = 1./2.*TMath::Sqrt(TMath::Power(xFullError[i],2)+TMath::Power(xIntrinsicError[i],2));
    }

    fileName.Form("c2%i",j);
    canvas2[j] = new TCanvas(fileName,fileName,0,0,500,500);
    canvas2[j] ->cd();

    Ratio = new TGraphErrors(numEntries,x,y,xError,yError);
    char titleName[100] = {0};
    sprintf(titleName,"MCClosure - %s",file2);
    Ratio->SetTitle(titleName);
    Ratio-> GetYaxis()->SetTitle("JER/intrinsic -1");
    Ratio->GetYaxis()->SetTitleOffset(1.5); 
    Ratio->SetMarkerStyle(20);
    Ratio->SetMarkerSize(1.1);
    Ratio ->SetMaximum(0.1);
    Ratio ->SetMinimum(-0.1);
    Ratio->Draw("AP");

    TF1* line = new TF1("line","0.",0,600);
    line->SetLineColor(2);
    line->Draw("same");

    TLatex*  info   = new TLatex();
    char legname[100] = {0};
    info->SetTextFont(132);
    info-> SetNDC();
    info->SetTextSize(0.040);
    sprintf(legname,"%i. eta Bin",j+1);
    info->DrawLatex(0.60,0.80,legname);

    fileName.Form("MCClosure_for_%i_eta_bin%s.pdf",j+1,file1);
    canvas2[j]->SaveAs(fileName);

    */
  }


  return 0;

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int makeGraphMC(){

 
  TGraphErrors * graph;
  double x[7],y[7],xErr[7],yErr[7];
  TH2D* histo2D;
  

  TFile *file = TFile::Open("PhotonPtNVtx_2d_PF_mc.root");
  file->GetObject("Photon Pt against Number of Vertices",histo2D);

  TH1D * histoVtx = (TH1D*)histo2D->ProjectionY("histoVtx",0,1200);

  
  TH1D* histo1D = (TH1D*)histo2D->ProjectionX("histo1D",1,5);
  TH1D* histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histoVtxPart->SetBinContent(1,histoVtx->GetBinContent(1));
  histoVtxPart->SetBinContent(2,histoVtx->GetBinContent(2));
  histoVtxPart->SetBinContent(3,histoVtx->GetBinContent(3));
  histoVtxPart->SetBinContent(4,histoVtx->GetBinContent(4));
  histoVtxPart->SetBinContent(5,histoVtx->GetBinContent(5));
  y[0] = histo1D->GetMean();
  yErr[0] = histo1D->GetMeanError();  
  x[0] = histoVtxPart->GetMean();
  xErr[0] = histoVtxPart->GetMeanError();

  cout<<"x0 = "<<x[0]<<endl;
  cout<<"x0Err = "<<xErr[0]<<endl;

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",6,10);
  histoVtxPart->SetBinContent(6,histoVtx->GetBinContent(6));
  histoVtxPart->SetBinContent(7,histoVtx->GetBinContent(7));
  histoVtxPart->SetBinContent(8,histoVtx->GetBinContent(8));
  histoVtxPart->SetBinContent(9,histoVtx->GetBinContent(9));
  histoVtxPart->SetBinContent(10,histoVtx->GetBinContent(10));
  y[1] = histo1D->GetMean();
  yErr[1] = histo1D->GetMeanError();
  x[1] = histoVtxPart->GetMean();
  xErr[1] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",11,15);
  histoVtxPart->SetBinContent(11,histoVtx->GetBinContent(11));
  histoVtxPart->SetBinContent(12,histoVtx->GetBinContent(12));
  histoVtxPart->SetBinContent(13,histoVtx->GetBinContent(13));
  histoVtxPart->SetBinContent(14,histoVtx->GetBinContent(14));
  histoVtxPart->SetBinContent(15,histoVtx->GetBinContent(15));
  y[2] = histo1D->GetMean();
  yErr[2] = histo1D->GetMeanError();
  x[2] = histoVtxPart->GetMean();
  xErr[2] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",16,20);
  histoVtxPart->SetBinContent(16,histoVtx->GetBinContent(16));
  histoVtxPart->SetBinContent(17,histoVtx->GetBinContent(17));
  histoVtxPart->SetBinContent(18,histoVtx->GetBinContent(18));
  histoVtxPart->SetBinContent(19,histoVtx->GetBinContent(19));
  histoVtxPart->SetBinContent(20,histoVtx->GetBinContent(20));
  y[3] = histo1D->GetMean();
  yErr[3] = histo1D->GetMeanError();
  x[3] = histoVtxPart->GetMean();
  xErr[3] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",21,25);
  histoVtxPart->SetBinContent(21,histoVtx->GetBinContent(21));
  histoVtxPart->SetBinContent(22,histoVtx->GetBinContent(22));
  histoVtxPart->SetBinContent(23,histoVtx->GetBinContent(23));
  histoVtxPart->SetBinContent(24,histoVtx->GetBinContent(24));
  histoVtxPart->SetBinContent(25,histoVtx->GetBinContent(25));
  y[4] = histo1D->GetMean();
  yErr[4] = histo1D->GetMeanError();
  x[4] = histoVtxPart->GetMean();
  xErr[4] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",26,30);
  histoVtxPart->SetBinContent(26,histoVtx->GetBinContent(26));
  histoVtxPart->SetBinContent(27,histoVtx->GetBinContent(27));
  histoVtxPart->SetBinContent(28,histoVtx->GetBinContent(28));
  histoVtxPart->SetBinContent(29,histoVtx->GetBinContent(29));
  histoVtxPart->SetBinContent(30,histoVtx->GetBinContent(30));
  y[5] = histo1D->GetMean();
  yErr[5] = histo1D->GetMeanError();
  x[5] = histoVtxPart->GetMean();
  xErr[5] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",31,35);
  histoVtxPart->SetBinContent(31,histoVtx->GetBinContent(31));
  histoVtxPart->SetBinContent(32,histoVtx->GetBinContent(32));
  histoVtxPart->SetBinContent(33,histoVtx->GetBinContent(33));
  histoVtxPart->SetBinContent(34,histoVtx->GetBinContent(34));
  histoVtxPart->SetBinContent(35,histoVtx->GetBinContent(35));
  y[6] = histo1D->GetMean();
  yErr[6] = histo1D->GetMeanError();
  x[6] = histoVtxPart->GetMean();
  xErr[6] = histoVtxPart->GetMeanError();
 

  graph = new TGraphErrors(7,x,y,xErr,yErr);
  TF1* f1 = new TF1("pol1","pol1"); 
  graph -> Fit(f1,"Q","",5.0,35.0);

  char legEntry[100];


  cout<<"par0 = "<<f1->GetParameter(0)<<endl;
  cout<<"par1 = "<<f1->GetParameter(1)<<endl;
  sprintf(legEntry,"%4.2f + %4.2f (#pm %4.2f) * x",f1->GetParameter(0),f1->GetParameter(1),f1->GetParError(1));

  TLegend* leg_hist       =  0 ;
  leg_hist = new TLegend(0.5,0.75,0.9,0.9);
  leg_hist->SetTextSize(0.033);
  leg_hist->SetHeader(legEntry);
  leg_hist -> SetFillColor(0);
  

  graph -> SetTitle("Pile-up dependence of Photon pT (MC)");
  graph -> GetXaxis() -> SetTitle("#Vtx"); 
  graph -> GetYaxis() -> SetTitle("Photon pT");
  graph -> SetMaximum(230);
  graph  -> GetYaxis() -> SetTitleOffset(1.3); 
  

  TCanvas* canvas = new TCanvas("c","c",0,0,500,500);
  canvas ->cd();
  graph->Draw("AP");
  leg_hist->Draw("same");

  canvas->SaveAs("PUDependence_MC.pdf");

  return 0;

}


int makeOtherGraph(){


  TH1D *hPhotonPtData;
  TH1D *hPhotonPtMC;
  TFile *file, *file1;  
  /*
  file = TFile::Open("Photon1Pt_PF_data.root");
  file->GetObject("4",hPhotonPtData);
    
  file1 = TFile::Open("Photon1Pt_PF_mc.root");
  file1->GetObject("4",hPhotonPtMC);

  hPhotonPtMC->Scale(hPhotonPtData->Integral()/hPhotonPtMC->Integral()); 

  //TCanvas* canvas = DrawComparison(hPhotonPtData,hPhotonPtMC,"Photon Pt spectrum (data and MC) New","Lumititle","Photon Pt", 1);
  
  //canvas ->Print("PhotonPtComparison.pdf");

  
  delete file;
  delete file1; 
  */ 

  TH1D* hNVtxData;
  TH1D* hNVtxMC;

  file = TFile::Open("VtxN_PF_data.root");
  file->GetObject("hVtxN",hNVtxData);

  
  file1 = TFile::Open("VtxN_PF_mc.root");
  file1->GetObject("hVtxN",hNVtxMC);

  hNVtxMC->Scale(hNVtxData->Integral()/hNVtxMC->Integral()); 

  //canvas = DrawComparison(hNVtxData, hNVtxMC,"Number of Vertices (data and MC) New PU Dist","Lumititle","#Vtx", 1);
  //canvas ->Print("NVtxComparison.pdf");

  

  char title[100];
  TH1D* hNVtxD[8];
  TH1D* hNVtxM[8];
  TCanvas* canv[8];
  TCanvas* canvas;
  for(int i=0; i<8;i++){

  delete file;
  delete file1;
 
  sprintf(title,"VtxN%i_PF_data.root",i);
  file = TFile::Open(title);
  
  sprintf(title,"hVtxPtBinned%i",i);
  file->GetObject(title,hNVtxD[i]);
  
  sprintf(title,"VtxN%i_PF_mc.root",i);
  file1 = TFile::Open(title);

  sprintf(title,"hVtxPtBinned%i",i);
  file1->GetObject(title,hNVtxM[i]);
  
  
  hNVtxM[i]->Scale(hNVtxD[i]->Integral()/hNVtxM[i]->Integral()); 
  
  sprintf(title,"Number of Vertices (data and MC) pt %i  NEW PU Dist",i);
  canv[i] = DrawComparison(hNVtxD[i], hNVtxM[i],title,"","#Vtx", 1);
  
  sprintf(title,"NVtxComparison%i.pdf",i);
  canv[i] ->Print(title);

  }

  delete file;
  delete file1;
  TH1D* hIsoData;
  TH1D* hIsoMC;
 
  sprintf(title,"hPhotonIsoEcal_PF_data.root");
  file = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoEcal");
  file->GetObject(title,hIsoData);
  
  sprintf(title,"hPhotonIsoEcal_PF_mc.root");
  file1 = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoEcal");
  file1->GetObject(title,hIsoMC);
  
  
  hIsoMC->Scale(hIsoData->Integral()/hIsoMC->Integral()); 
  
  sprintf(title,"Ecal Isolation");
  canvas = DrawComparison(hIsoData,hIsoMC,title,"","Iso", 1);
  
  sprintf(title,"PhotonIsoEcal.pdf");
  canvas ->Print(title);

  delete file;
  delete file1;
  
 
  sprintf(title,"hPhotonIsoHcal_PF_data.root");
  file = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoHcal");
  file->GetObject(title,hIsoData);
  
  sprintf(title,"hPhotonIsoHcal_PF_mc.root");
  file1 = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoHcal");
  file1->GetObject(title,hIsoMC);
  
  
  hIsoMC->Scale(hIsoData->Integral()/hIsoMC->Integral()); 
  
  sprintf(title,"Hcal Isolation");
  canvas = DrawComparison(hIsoData, hIsoMC,title,"","Iso", 1);
  
  sprintf(title,"PhotonIsoHcal.pdf");
  canvas ->Print(title);

  delete file;
  delete file1;
  
  sprintf(title,"hPhotonIsoTrk_PF_data.root");
  file = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoTrk");
  file->GetObject(title,hIsoData);
  
  sprintf(title,"hPhotonIsoTrk_PF_mc.root");
  file1 = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoTrk");
  file1->GetObject(title,hIsoMC);
  
  
  hIsoMC->Scale(hIsoData->Integral()/hIsoMC->Integral()); 
  
  sprintf(title,"Trk Isolation");
  canvas = DrawComparison(hIsoData, hIsoMC,title,"","Iso", 1);
  
  sprintf(title,"PhotonIsoTrk.pdf");
  canvas ->Print(title);
  

  


  return 0;
}

int makeRhoGraph(){

  char title[100];
  TH1D* hRhoD[8];
  TH1D* hRhoM[8];
  TCanvas* canv[8];
  TFile *file;
  TFile *file1;
  for(int i=0; i<8;i++){

    sprintf(title,"Rho%i_PF_data.root",i);
    file = TFile::Open(title);
  
    sprintf(title,"hRhoPtBinned%i",i);
    file->GetObject(title,hRhoD[i]);
    
    sprintf(title,"Rho%i_PF_mc.root",i);
    file1 = TFile::Open(title);
    
    sprintf(title,"hRhoPtBinned%i",i);
    file1->GetObject(title,hRhoM[i]);
    
 
    hRhoM[i]->Scale(hRhoD[i]->Integral()/hRhoM[i]->Integral()); 
    
    sprintf(title,"Rho (data and MC) pt %i ",i);
    canv[i] = DrawComparison(hRhoD[i], hRhoM[i],title,"","Rho", 1);
    
    sprintf(title,"RhoComparison%i.pdf",i);
    
    canv[i]->Print(title);
   
    delete file;
    delete file1;
    
  }

  return 0;
}

int makeOtherRhoGraph(){

  char title[100];
  TH1D* hRhoD[40];
  TH1D* hRhoM[40];
  TCanvas* canv[40];
  TFile *file;
  TFile *file1;
  for(int i=1; i<40;i++){
    
    sprintf(title,"RhoVtxBinned/RhoVtxBinned%i_PF_data.root",i);
    file = TFile::Open(title);
   
    sprintf(title,"hRhoVtxBinned%i",i);
    file->GetObject(title,hRhoD[i]);
     
    sprintf(title,"RhoVtxBinned/RhoVtxBinned%i_PF_mc.root",i);
    file1 = TFile::Open(title);
     
    sprintf(title,"hRhoVtxBinned%i",i);
    file1->GetObject(title,hRhoM[i]);
  
 
    hRhoM[i]->Scale(hRhoD[i]->Integral()/hRhoM[i]->Integral()); 
     
    sprintf(title,"Rho (data and MC) for #Vtx = %i ",i);
    canv[i] = DrawComparison(hRhoD[i], hRhoM[i],title,"","Rho", 1);
     
    sprintf(title,"RhoVtxBinned/RhoComparisonVtxBinned%i.pdf",i);
      
    canv[i]->Print(title);
     
    delete file;
    delete file1;
    
  }

  return 0;
}

int makeTriggerEffGraph(){

  char title[100];
  TH1D* hTriggerBefore[8];
  TH1D* hTriggerAfter[8];
  TH1D* hTrigger[8];
  TCanvas* canv[8];
  TFile *file;
  TFile *file1;
  for(int i=0; i<8;i++){


    sprintf(title,"TriggerEffBefore%i_PF_mc.root",i);
    file = TFile::Open(title);
    
    sprintf(title,"hTriggerEffPtBinnedBefore%i",i);
    file->GetObject(title,hTriggerBefore[i]);
     
    sprintf(title,"TriggerEffAfter%i_PF_mc.root",i);
    file1 = TFile::Open(title);
    
    sprintf(title,"hTriggerEffPtBinnedAfter%i",i);
    file1->GetObject(title,hTriggerAfter[i]);
    
    hTrigger[i] = new TH1D("ratio","ratio",60,0,60);

    hTrigger[i]->Divide(hTriggerAfter[i],hTriggerBefore[i]);
    canv[i] = new TCanvas("ratio","ratio",500,500);
    canv[i] -> cd();

    hTrigger[i]->SetTitle("Trigger efficiency in MC");

    hTrigger[i]->SetXTitle("#Vtx");

    hTrigger[i]->SetYTitle("#events pass Trigger / # all events");

    hTrigger[i]->Draw();
   
    sprintf(title,"TriggerEff%i.pdf",i);
    canv[i]->Print(title);
    
    delete file;
    delete file1;
    
  }

  return 0;
}


TCanvas* DrawComparison(TH1D* prediction, TH1D* selection, TString Title, TString LumiTitle, TString xTitle, bool isData)
{
   double MinX = selection->GetXaxis()->GetXmin();
   double MaxX = selection->GetXaxis()->GetXmax();
   double MaxY = selection->GetMaximum();
   double YRangeMax = MaxY;
   TString titlePrediction;
   TString titleSelection;
   TString RatioTitle;
   
   if( isData ){
      titlePrediction = "Data";
      titleSelection = "MC";
      RatioTitle = "(Data-MC)/MC";
   }
   else {
      titlePrediction = "Data-driven Pred. from MC";
      titleSelection = "MC Expectation";
      RatioTitle = "(Pred-MC)/MC";
   }

   //static Int_t c_LightBrown   = TColor::GetColor( "#D9D9CC" );
   static Int_t c_LightGray    = TColor::GetColor( "#DDDDDD" );

   prediction->SetAxisRange(MinX, MaxX, "X");
   prediction->GetYaxis()->SetRangeUser(0.005, YRangeMax);
   prediction->SetMarkerStyle(20);
   prediction->SetMarkerSize(0.9);
   prediction->SetMarkerColor(kBlack);
   prediction->SetXTitle(xTitle);
   prediction->SetYTitle("Events");
   selection->SetAxisRange(MinX, MaxX, "X");
   selection->GetYaxis()->SetRangeUser(0.05, YRangeMax);
   // selection->SetFillColor(c_LightBrown);
   selection->SetFillColor(c_LightGray);
   selection->SetTitle("");
   selection->SetXTitle(xTitle);
   selection->SetYTitle("Events");
   TCanvas *c = new TCanvas("ca", "Comparison and ratio of two histos", 700, 700);
   TPad *pad1 = new TPad("pad1a", "pad1a", 0, 0.35, 1, 1);
   //pad1->SetLogy();
   pad1->SetBottomMargin(0);
   pad1->Draw();
   pad1->cd();
  
   selection->DrawCopy("hist");
   prediction->Draw("same");
   selection->SetFillColor(kAzure-3);
   selection->SetFillStyle(3354);
   selection->DrawCopy("e2same");
   selection->SetFillStyle(1001);
   //  selection->SetFillColor(c_LightBrown);
   selection->SetFillColor(c_LightGray);

   TLegend* leg1 = new TLegend(0.48, 0.63, 0.95, 0.83);
   leg1->SetFillStyle(0);
   leg1->SetLineStyle(1);
   leg1->SetTextFont(42);
   leg1->SetTextSize(0.04);
   leg1->AddEntry(selection, titleSelection, "lf");
   leg1->AddEntry(prediction, titlePrediction, "lep");
   leg1->Draw("same");
 
   TPaveText* pt = new TPaveText(0.11, 0.98, 0.95, 0.86, "NDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextSize(0.045);
   pt->AddText(Title);
   pt->AddText(LumiTitle);
   pt->Draw();
   c->cd();
   TPad *pad2 = new TPad("pad2a", "pad2a", 0, 0, 1, 0.35);
   pad2->SetTopMargin(0);
   pad2->Draw();
   pad2->cd();
   TH1D* r = new TH1D(*selection);
   r->SetTitle("");
   r->SetLabelSize(0.08, "XYZ");
   r->SetLabelOffset(0.01, "XYZ");
   r->SetTitleSize(0.09, "XYZ");
   r->SetTitleOffset(0.65, "Y");
   r->SetTickLength(0.05);
   r->SetYTitle(RatioTitle);
   r->SetStats(0);
   r->SetMarkerStyle(20);
   r->SetMarkerSize(0.9);
   r->SetMarkerColor(kBlack);
   r->Reset();
   r->Add(prediction, 1);
   r->Add(selection, -1);
   r->Divide(selection);
   r->SetMaximum(1.2);
   r->SetMinimum(-1.2);
   r->Draw("ep");
   TLine l;
   l.DrawLine(MinX, 0., MaxX, 0.);
   c->cd();
   return c;
}
