// This file is supposed to plot TGraphs also more than 1 in one plot

#ifndef FileOps_h
#define FileOps_h

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
#include "TText.h"
#include "TROOT.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraph.h"
#include "TRandom.h"
//#include "FitResult.h"

//#include "utils.h"
//#include "HistOps.h"

TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName);
//TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName);
TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName);
TF1* readTF1(const TString &fileName, const TString &fName, const TString &newFName);

const int pt_int = 10;
const int eta_int = 3;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Double_t fitfunc(Double_t *x, Double_t *p){

      double Npu = int((*x)/1000.);
      double pTgen = (*x) - Npu*1000.;
      double N = p[0] + p[1]*sqrt(Npu); // or p[0] + p[1]*Npu
      double S = p[2];
      double m = p[3];
      double C = 0;

      return sqrt(N*N/(pTgen*pTgen) + S*S*pow(pTgen,m-1) + C*C);
}

int plotPUdependence_simultaneousFit(){


  const int NPUbins     = 5;
  const bool linear     = 0;

  gDirectory->Delete(); 
  gROOT->GetListOfCanvases()->Delete();

  TString pdfFile, etaRegion, ptRegion, title, EtaPtRegion; 
  TString sourcePUle5, sourcePUle10gt5, sourcePUle15gt10, sourcePUle20gt15, sourcePUgt20; 
  double ptBorders[pt_int+1]   = {0}; 
  double etaBorders[eta_int+1] = {0};
  TString sourcePU[NPUbins];

  etaBorders[0] = 0.0;
  etaBorders[1] = 1.3;
  etaBorders[2] = 2.0; 
  etaBorders[3] = 3.0; \

  const double bd[pt_int+1]          = {30,40,60,80,120,160,200,250,300,400,1000000}; 
  

  ptBorders[0]  = 30;
  ptBorders[1]  = 40;
  ptBorders[2]  = 60;
  ptBorders[3]  = 80;
  ptBorders[4]  = 120;
  ptBorders[5]  = 160;
  ptBorders[6]  = 200;
  ptBorders[7]  = 250;
  ptBorders[8]  = 300;
  ptBorders[9]  = 400;

  TLatex*  info;
  TLegend *legend, *legend1;
  TGraphErrors* Add;
  TF1 *fScaleAlpha;
  TF1 *total;

  TString text;

  TF1* deltaJER[eta_int];
  TF1* test[eta_int];

  float NerrorNPU0[eta_int] = {0};
  float NerrorNPU140[eta_int] = {0};
  float Nresult[eta_int] = {0};
  float NresultError[eta_int] = {0};

  TGraphErrors *PU[NPUbins];
  TH1 *vtx[NPUbins];
  TH1 *vtxEta[NPUbins][pt_int];

  TGraph* chi2_ndf[eta_int];
  double NPU[NPUbins];
  double chi2NDF[NPUbins];

  for(int etaBin=0; etaBin<eta_int; etaBin++){
    //for(int etaBin=0; etaBin<1; etaBin++){

    if(etaBin == 0) etaRegion.Form("|#eta| < %4.1f", etaBorders[etaBin+1]);
    else       etaRegion.Form("%4.1f < |#eta| < %4.1f",etaBorders[etaBin], etaBorders[etaBin+1]);
  
    title = "Jet Energy Resolution";
  
    sourcePU[0].Form("plots_2012/PF_L1CHS/mc/Resolution_for_%i_eta_bin_intrinsic_PFCHS_mc_le5.root",etaBin+1);
    sourcePU[1].Form("plots_2012/PF_L1CHS/mc/Resolution_for_%i_eta_bin_intrinsic_PFCHS_mc_le10gt5.root",etaBin+1);
    sourcePU[2].Form("plots_2012/PF_L1CHS/mc/Resolution_for_%i_eta_bin_intrinsic_PFCHS_mc_le15gt10.root",etaBin+1);
    sourcePU[3].Form("plots_2012/PF_L1CHS/mc/Resolution_for_%i_eta_bin_intrinsic_PFCHS_mc_le20gt15.root",etaBin+1);
    sourcePU[4].Form("plots_2012/PF_L1CHS/mc/Resolution_for_%i_eta_bin_intrinsic_PFCHS_mc_gt20.root",etaBin+1);

    for(int puBin=0; puBin<NPUbins; puBin++)  PU[puBin] = readTGraphErrors(sourcePU[puBin],"Graph","Graph");
  
    TCanvas *c = new TCanvas("c",title,200,10,450,450);
    c -> SetLeftMargin(0.12);
    c -> cd();
    c->SetBottomMargin(0.12);
    c->SetLeftMargin(0.12);
  
    legend1  = new TLegend(0.4,0.65,0.9,0.9);
    legend1 -> SetFillColor(0);
    legend1 -> SetTextSize(0.030);
    PU[0]->GetFunction("fResolution")->SetLineColor(3);  
    chi2NDF[0] = PU[0]->GetFunction("fResolution")->GetChisquare()/PU[0]->GetFunction("fResolution")->GetNDF();
    //cout<<endl<<"Chi2_1/NDF = "<<chi2NDF[0]<<endl<<endl;
  
    PU[0]->GetXaxis()->SetTitle("p_{T}^{gen. Jet}");
    PU[0]->SetMarkerColor(3);
    PU[0]->SetLineColor(3);
    PU[0]->SetMinimum(0.0);
    PU[0]->SetMaximum(0.70);
    PU[0]->GetXaxis()->SetLimits(0,700);
    PU[0]->Draw("AP");
    legend1 -> AddEntry(PU[0],"NPU #leq 8","l");

    PU[1] -> GetFunction("fResolution")->SetLineColor(4);
    PU[1] -> SetMarkerColor(4);
    PU[1] -> SetLineColor(4);
    legend1 -> AddEntry(PU[1],"8 < NPU #leq 16","l");

    PU[2]->GetFunction("fResolution")->SetLineColor(5);
    PU[2]->SetMarkerColor(5);
    PU[2]->SetLineColor(5);
    legend1 -> AddEntry(PU[2],"16 < NPU #leq 24","l");

    PU[3] -> GetFunction("fResolution")->SetLineColor(6);
    PU[3] -> SetMarkerColor(6);
    PU[3] -> SetLineColor(6);
    legend1 -> AddEntry(PU[3],"24 < NPU #leq 32","l");

    PU[4]->GetFunction("fResolution")->SetLineColor(7);
    PU[4]->SetMarkerColor(7);
    PU[4]->SetLineColor(7);
    legend1 -> AddEntry(PU[4],"32 < NPU","l");

    for(int puBin=1; puBin<NPUbins; puBin++) PU[puBin]->Draw("Psame");

    
    double xN[NPUbins], xS[NPUbins], xm[NPUbins],xC[NPUbins], xEta[NPUbins], xEtaError[NPUbins];
    //int xEta[NPUbins];

    double xNError[NPUbins] = {0};
    double xSError[NPUbins] = {0};
    double xmError[NPUbins] = {0};
    double xCError[NPUbins] = {0};

 
    for(int ptBin=0; ptBin<10; ptBin++){

      sourcePU[0].Form("plots_2012/PF_L1CHS/mc/root_files/hVtxN1_%i_Pt_bin_%i_eta_bin_PFCHS_mc.root",ptBin+1,etaBin+1);
      sourcePU[1].Form("plots_2012/PF_L1CHS/mc/root_files/hVtxN2_%i_Pt_bin_%i_eta_bin_PFCHS_mc.root",ptBin+1,etaBin+1);
      sourcePU[2].Form("plots_2012/PF_L1CHS/mc/root_files/hVtxN3_%i_Pt_bin_%i_eta_bin_PFCHS_mc.root",ptBin+1,etaBin+1);
      sourcePU[3].Form("plots_2012/PF_L1CHS/mc/root_files/hVtxN4_%i_Pt_bin_%i_eta_bin_PFCHS_mc.root",ptBin+1,etaBin+1);
      sourcePU[4].Form("plots_2012/PF_L1CHS/mc/root_files/hVtxN5_%i_Pt_bin_%i_eta_bin_PFCHS_mc.root",ptBin+1,etaBin+1);
    
      for(int puBin = 0; puBin<NPUbins; puBin++){

	TFile * file = TFile::Open(sourcePU[puBin]);
	file->GetObject("histo",vtxEta[puBin][ptBin]);
	vtxEta[puBin][ptBin] -> SetDirectory(0);
	delete file;
   
      }
    }


    for(int puBin = 0; puBin<NPUbins; puBin++){

      for(int ptBin=1; ptBin<10; ptBin++) vtxEta[puBin][0]->Add(vtxEta[puBin][ptBin]);
      xEta[puBin] = vtxEta[puBin][0] -> GetMean(1);
      cout<<"xEta["<<puBin<<"] ="<<xEta[puBin]<<endl;
      xEtaError[puBin] = vtxEta[puBin][0] -> GetMeanError(1);

    }
  

    // Try Multigraph Fit
    double *x    = new double[pt_int*NPUbins];
    double *y    = new double[pt_int*NPUbins];
    double *yErr = new double[pt_int*NPUbins];
    double *xErr = new double[pt_int*NPUbins];
    int l = 0;
    int i = 0;
   
  
    for(int j=0; j< NPUbins; j++){
      //cout<<"xEta["<<j<<"] = "<<xEta[j]<<endl;
      i = 0;
      while(i < PU[j]->GetN()){
     
	PU[j]->GetPoint(i,x[l],y[l]);
      	
	x[l] = x[l] + 1000.*(int)(xEta[j]);
	//x[l] = x[l] + 1000.*(xEta[j]);
	
	
	yErr[l] = PU[j] -> GetErrorY(i);
	
	xErr[l] = PU[j] -> GetErrorX(i);
	
	if(x[l] < ptBorders[i] &&  x[l] > ptBorders[i+1]) cout<< "not every pt bin filled!!!!!!!!!!!!!!!!!!"<<endl;
      
      i += 1;
      l += 1;
      }
    }

   
    TGraphErrors *gr = new TGraphErrors(l,x,y,xErr,yErr);
   
    double *p = new double[4];
       
    TCanvas *cgr = new TCanvas("cgr","gr",200,10,450,450);
    cgr->cd();
    gr->Draw("AP");
    cgr->SaveAs("PUdependence_simultaneousFit/test.pdf");


    
    TF1 *fit = new TF1("fit",fitfunc,0.,40*1000,4);
    fit->SetParName(0,"a");
    fit->SetParName(1,"b");
    fit->SetParName(2,"S");
    fit->SetParName(3,"m");
 
    fit->SetParameter(0,2.0);
    fit->SetParameter(1,1.0);
    fit->SetParameter(2,0.369);
    fit->SetParameter(3,0.386);

    gr->Fit(fit, "QRN"); 
   
    cout<<"Etabin Nr:  "<<etaBin<<endl;
    cout<<"fit->GetChisquare() = "<<fit->GetChisquare()<<endl;
    cout<<"fit->GetNDF() = "<<fit->GetNDF()<<endl;
    cout<<"S = "<<fit->GetParameter("S")<<endl;
    cout<<"SErr = "<<fit->GetParError(2)<<endl;
    cout<<"m = "<<fit->GetParameter("m")<<endl;
    cout<<"mErr = "<<fit->GetParError(3)<<endl;
    cout<<"a = "<<fit->GetParameter("a")<<endl;
    cout<<"aErr = "<<fit->GetParError(0)<<endl;
    cout<<"b = "<<fit->GetParameter("b")<<endl;
    cout<<"bErr = "<<fit->GetParError(1)<<endl<<endl;

 
    TCanvas *c2 = new TCanvas("c2",title,200,10,450,450);
    c2 -> SetLeftMargin(0.12);
    c2 -> cd();
    c2 -> SetBottomMargin(0.12);
    c2 -> SetLeftMargin(0.12);

    legend  = new TLegend(0.5,0.75,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.030);
  
    TF1* f1;

    TF1* fResolution = new TF1("fResolution","TMath::Sqrt(TMath::Sign(1,[0])*TMath::Power([0]/x,2)+TMath::Power([1],2)*TMath::Power(x,[3]-1)+TMath::Power([2],2))", ptBorders[0], 1000);
    fResolution->SetParName(0,"N");
    fResolution->SetParName(1,"S");
    fResolution->SetParName(2,"C");
    fResolution->SetParName(3,"m");
 
    fResolution->FixParameter(1,fit->GetParameter("S"));
    fResolution->FixParameter(3,fit->GetParameter("m"));
    fResolution->FixParameter(2,0);

    // Fit again the NmSC Fit to JER(ptGen)
    
    for(int puBin = 0; puBin<NPUbins; puBin++){

      PU[puBin] -> Fit("fResolution","QR");
      chi2NDF[puBin] = PU[puBin] -> GetFunction("fResolution") -> GetChisquare()/PU[puBin] -> GetFunction("fResolution") -> GetNDF();
      cout<<endl<<"Chi2_1/NDF = "<<chi2NDF[puBin]<<endl<<endl;
      xN[puBin] = fResolution->GetParameter("N");
      cout<<"xN["<<puBin<<"] = "<<xN[puBin]<<endl;
      xNError[puBin] = fResolution->GetParError(0);
      xS[puBin] = fResolution->GetParameter("S");
      xSError[puBin] = fResolution->GetParError(1);
      cout<<"S = "<<xS[puBin]<<endl;
      xC[puBin] = fResolution->GetParameter("C");
      xCError[puBin] = fResolution->GetParError(2);
      xm[puBin] = fResolution->GetParameter("m");
      cout<<"m = "<<xm[puBin]<<endl;
      xmError[puBin] = fResolution->GetParError(3);

    }
        
    PU[0]->GetFunction("fResolution")->SetLineColor(3);
    PU[1]->GetFunction("fResolution")->SetLineColor(4);
    PU[2]->GetFunction("fResolution")->SetLineColor(5);
    PU[3]->GetFunction("fResolution")->SetLineColor(6);
    PU[4]->GetFunction("fResolution")->SetLineColor(7);
    

    //TGraphErrors *N = new TGraphErrors(NPUbins,xEta,xN,xEtaError,xNError);
    TGraphErrors *S = new TGraphErrors(NPUbins,xEta,xS,xEtaError,xSError);
    TGraphErrors *m = new TGraphErrors(NPUbins,xEta,xm,xEtaError,xmError);
    TGraphErrors *C = new TGraphErrors(NPUbins,xEta,xC,xEtaError,xCError);

    TGraphErrors *N = new TGraphErrors(NPUbins,xEta,xN,xEtaError,xNError);
    chi2_ndf[etaBin] = new TGraph(NPUbins,xEta,chi2NDF);

    if(linear == 1) f1 = new TF1("name","pol1",0,30);
    else            f1 = new TF1("name","[0] + [1]*TMath::Sqrt(x)",0,30); 
  
    N->Fit("name","SQ");

    cout<<"Chisquare = "<<f1->GetChisquare()<<endl;
    cout<<"NDF = "<<f1->GetNDF()<<endl;
  
    N -> SetTitle("NPU Dependency of JER Fit Parameter");
    N -> GetXaxis()->SetTitle("NPU");
    N -> GetYaxis()->SetTitle("N/S/m/C");
    N -> SetMinimum(-5.0);
    N -> SetMaximum(25);
    N -> GetYaxis() -> SetTitleOffset(1.5);
    N -> SetMarkerStyle(20); 
    N -> SetMarkerColor(2);
    N -> SetLineColor(2);
    legend -> AddEntry(N,"N","l");

    S -> SetMarkerColor(3);
    S -> SetLineColor(3);
    S -> SetMarkerStyle(20); 
    //S -> GetFunction("horizontalS")->SetLineColor(3);
    legend -> AddEntry(S,"S","l");
    
    m -> SetMarkerColor(1);
    m -> SetLineColor(1);
    m -> SetMarkerStyle(20); 
    //m -> GetFunction("horizontalm")->SetLineColor(1);
    legend -> AddEntry(m,"m","l");
  
    C -> SetMarkerColor(5);
    C -> SetLineColor(5);
    C -> SetMarkerStyle(20); 
    legend -> AddEntry(C,"C","l");
  
    N -> Draw("AP");
    S -> Draw("Psame");
    m -> Draw("Psame");
    C -> Draw("Psame");
    legend->Draw("same");

    info   = new TLatex();
    info->SetTextFont(132);
    info-> SetNDC();
    info->SetTextSize(0.042);
  
    info->DrawLatex(0.14,0.83,  etaRegion.Data());
    info->DrawLatex(0.14,0.75, "Anti-k_{T} 0.5 PFCHSJets");
 
    if(linear == 1) text.Form("f(NPU) = %4.2f(#pm %4.2f) #upoint NPU + %4.2f(#pm %4.2f)",f1->GetParameter(1),f1->GetParError(1),f1->GetParameter(0),f1->GetParError(0));
    else text.Form("f(NPU) = %4.2f(#pm %4.2f) #upoint #sqrt{NPU} + %4.2f(#pm %4.2f)",f1->GetParameter(1),f1->GetParError(1),f1->GetParameter(0),f1->GetParError(0));
    info->DrawLatex(0.14,0.65, text);

    text.Form("Chi2 = %5.2f, ndof = %i",f1->GetChisquare(),f1->GetNDF());
    info->DrawLatex(0.50,0.59, text);

    if(linear == 1) pdfFile.Form("PUdependence_simultaneousFit/PUDependecy_%i_EtaBin_FitParameter_linear.pdf",etaBin+1);
    else pdfFile.Form("PUdependence_simultaneousFit/PUDependecy_%i_EtaBin_FitParameter_sqrt.pdf",etaBin+1);
    c2->SaveAs(pdfFile);
    delete c2;
    delete info;
    //====================================================================================================================================================================================

    c->cd();
    // Now Fill again the JER(ptGen) plot with the results of Fit parameter N for NPU=0 and NPU =140
    TF1* extra0, *extra1;
    fResolution->FixParameter(0,f1->GetParameter(0));
    extra0 = (TF1*)fResolution->Clone();  
 
    extra0->SetLineColor(9);
    extra0->Draw("same");
    legend1 -> AddEntry(extra0,"NPU = 0 (extrapolation)","l");
    legend1->Draw("same");

    info   = new TLatex();
    info -> SetTextFont(132);
    info -> SetNDC();
    info -> SetTextSize(0.042);
  
    NerrorNPU0[etaBin] = TMath::Sqrt(TMath::Power(0.,2)*TMath::Power(f1->GetParError(1),2)+TMath::Power(f1->GetParError(0),2)); 

    text.Form("NPU = 0    : N = %4.2f (#pm %4.2f) GeV",fResolution->GetParameter("N"),NerrorNPU0[etaBin]);
    info->DrawLatex(0.25,0.60, text);
 
    if(linear == 1) fResolution->FixParameter(0,f1->GetParameter(0)+ 140.* (f1->GetParameter(1)));
    else fResolution->FixParameter(0,f1->GetParameter(0)+ TMath::Sqrt(140.)* (f1->GetParameter(1)));
    extra1 = fResolution;
    extra1 -> SetLineColor(44);
    extra1 -> Draw("same");
    legend1 -> AddEntry(extra1,"NPU = 140 (extrapolation)","l");
    legend1->Draw("same");
    

    if(linear == 1)  NerrorNPU140[etaBin] = TMath::Sqrt(TMath::Power(140.,2)*TMath::Power(f1->GetParError(1),2)+TMath::Power(f1->GetParError(0),2)); 
    else             NerrorNPU140[etaBin] = TMath::Sqrt(140.*TMath::Power(f1->GetParError(1),2)+TMath::Power(f1->GetParError(0),2)); 
    text.Form("NPU = 140: N = %4.2f (#pm %4.2f) GeV",fResolution->GetParameter("N"),NerrorNPU140[etaBin]);
    info->DrawLatex(0.25,0.55, text);
      text.Form("S = %4.3f, m = %4.3f",fResolution->GetParameter("S"),fResolution->GetParameter("m"));
      info->DrawLatex(0.25,0.50, text);
      info->DrawLatex(0.50,0.3, "Anti-k_{T} 0.5 PFCHSJets");
  
      if(linear == 1)      pdfFile.Form("PUdependence_simultaneousFit/PUDependecy_%i_EtaBin_linear.pdf",etaBin+1);
      else if(linear == 0) pdfFile.Form("PUdependence_simultaneousFit/PUDependecy_%i_EtaBin_sqrt.pdf",etaBin+1);
      c -> SaveAs(pdfFile);
      delete info;
      delete legend1;
      delete c;
      //==================================================================================================================================================================================
  
      // deltaJER = Sqrt(JER(NPU=140)^2 - JER(NPU=0)^2)
      deltaJER[etaBin] = new TF1("deltaJER","(TMath::Sqrt(TMath::Sign(1,[0])*TMath::Power([0]/x,2)+TMath::Power([1],2)*TMath::Power(x,[3]-1)+TMath::Power([2],2) - (TMath::Sign(1,[4])*TMath::Power([4]/x,2)+TMath::Power([5],2)*TMath::Power(x,[7]-1)+TMath::Power([6],2))))*x",30,1000);
   
      deltaJER[etaBin]->FixParameter(0,f1->Eval(140));
      deltaJER[etaBin]->FixParameter(1,xS[0]);
      deltaJER[etaBin]->FixParameter(2,0.0);
      deltaJER[etaBin]->FixParameter(3,xm[0]);
      deltaJER[etaBin]->FixParameter(4,f1->Eval(0));
      deltaJER[etaBin]->FixParameter(5,xS[0]);
      deltaJER[etaBin]->FixParameter(6,0.0);
      deltaJER[etaBin]->FixParameter(7,xm[0]);

      test[etaBin] = new TF1("deltaJER","TMath::Sqrt((TMath::Sign(1,[4])*TMath::Power([4]/x,2)+TMath::Power([5],2)*TMath::Power(x,[7]-1)+TMath::Power([6],2)))*x",30,1000);

      test[etaBin]->FixParameter(4,23.9);
      test[etaBin]->FixParameter(5,0.00045);
      test[etaBin]->FixParameter(6,0.0);
      test[etaBin]->FixParameter(7,2.6); 

      Nresult[etaBin] = TMath::Sqrt(TMath::Power(f1->Eval(140),2) - TMath::Power(f1->Eval(0),2));
      NresultError[etaBin] = 2./TMath::Sqrt(TMath::Power(f1->Eval(140),2) - TMath::Power(f1->Eval(0),2)) * TMath::Sqrt(TMath::Power(f1->Eval(140),2)*TMath::Power(NerrorNPU140[etaBin],2) + TMath::Power(f1->Eval(0),2)* TMath::Power(NerrorNPU0[etaBin],2));
      cout<<"NerrorNPU0 ="<<NerrorNPU0[etaBin]<<endl;
      cout<<"Nresult["<<etaBin<<"] = "<<Nresult[etaBin]<<endl;
      cout<<"NresultError["<<etaBin<<"] = "<<NresultError[etaBin]<<endl;

      delete fResolution;
      delete N;
      delete m;
      delete C;
      delete S;
      delete f1;

  }

  TCanvas *c5 = new TCanvas("c5",title,200,10,450,450);
  c5->cd();

  legend  = new TLegend(0.5,0.75,0.9,0.9);
  legend -> SetFillColor(0);
  legend -> SetTextSize(0.030);

  deltaJER[0]->SetLineColor(2);
  deltaJER[1]->SetLineColor(3);
  deltaJER[2]->SetLineColor(4);

  deltaJER[0] ->SetMinimum(00.0);
  deltaJER[0] ->SetMaximum(30.0);

  for(int etaBin=0;etaBin<3;etaBin++){
    if(etaBin == 0) etaRegion.Form("|#eta| < %4.1f", etaBorders[etaBin+1]);
    else            etaRegion.Form("%4.1f < |#eta| < %4.1f",etaBorders[etaBin], etaBorders[etaBin+1]);
    
    deltaJER[etaBin]->SetTitle("#Delta JER [GeV]");
    
    legend->AddEntry(deltaJER[etaBin],etaRegion.Data(),"l");
    
    if(etaBin==0)  deltaJER[etaBin]->Draw();  
    else deltaJER[etaBin]->Draw("same");  

    //test[0]->Draw("same");  

    
  }

  TCanvas *c6 = new TCanvas("c6",title,200,10,450,450);
  c6->cd();

  
  legend1  = new TLegend(0.5,0.75,0.9,0.9);
  legend1 -> SetFillColor(0);
  legend1 -> SetTextSize(0.030);

  chi2_ndf[0]->SetMarkerStyle(20);
  chi2_ndf[1]->SetMarkerStyle(20);
  chi2_ndf[2]->SetMarkerStyle(20);
  chi2_ndf[0]->SetMarkerColor(1);
  chi2_ndf[1]->SetMarkerColor(2);
  chi2_ndf[2]->SetMarkerColor(3);
  chi2_ndf[0]->SetMaximum(60);
  chi2_ndf[0]->SetTitle("#chi^{2}/NDF for the JER over ptGen fits (C fixed)");
  chi2_ndf[0]->Draw("AP");
  legend1->AddEntry(chi2_ndf[0],"1. eta Bin","P");
  chi2_ndf[1]->Draw("Psame");
  legend1->AddEntry(chi2_ndf[1],"2. eta Bin","P");
  chi2_ndf[2]->Draw("Psame");
  legend1->AddEntry(chi2_ndf[2],"3. eta Bin","P");
  legend1->Draw("same");
  if(linear == 1) pdfFile.Form("PUdependence_simultaneousFit/Chisquare_linear.pdf");
  else pdfFile.Form("PUdependence_simultaneousFit/Chisquare_sqrt.pdf");
  c6->SaveAs(pdfFile);
  legend1->Draw("same");
  delete legend1;

  c5->cd();
  info   = new TLatex();
  info->SetTextFont(132);
  info-> SetNDC();
  info->SetTextSize(0.042);
  info->DrawLatex(0.13,0.7, "Anti-k_{T} 0.5 PFCHSJets");
  info->DrawLatex(0.13,0.15, "f(p_{T}) = #sqrt{JER(p_{T},N_{PU}=140)^{2} - JER(p_{T},N_{PU}=0)^{2}} #upoint p_{T}");
  legend->Draw("same");

  text.Form("|#eta| < %4.1f: N = %4.2f (#pm %4.2f) GeV", etaBorders[1],Nresult[0] , NresultError[0]);
  info->DrawLatex(0.25,0.42,text);
  text.Form("%4.1f < |#eta| < %4.1f: N = %4.2f (#pm %4.2f) GeV",etaBorders[1], etaBorders[2],Nresult[1], NresultError[1]);
  info->DrawLatex(0.25,0.35,text);
  text.Form("%4.1f < |#eta| < %4.1f: N = %4.2f (#pm %4.2f) GeV",etaBorders[2], etaBorders[3],Nresult[2] , NresultError[2]);
  info->DrawLatex(0.25,0.28,text);

  deltaJER[0]->GetHistogram()->GetXaxis()->SetTitle("p_{T}^{gen. Jet}");
  deltaJER[0]->GetHistogram()->GetYaxis()->SetTitle("JER #upoint p_{T}^{gen. Jet}");
  c5->Modified();

  if(linear == 1) pdfFile.Form("PUdependence_simultaneousFit/DeltaJER_linear.pdf");
  else            pdfFile.Form("PUdependence_simultaneousFit/DeltaJER_sqrt.pdf");
  c5 -> SaveAs(pdfFile);

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
      cout<<"in"<<endl;
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


// -------------------------------------------------------------------------------------
TF1* readTF1(const TString &fileName, const TString &fName, const TString &newFName) {
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

//! Read TGraphAsymmErrors from file
// -------------------------------------------------------------------------------------
/*  TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName) {
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
    }*/



/*
  namespace util
  {
  class FileOps
  {
  public:
  static TF1* readTF1(const TString &fileName, const TString &fName, const TString &newFName = "");
  static TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName = "");
  static TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName = "");
  //static TH2* readTH2(const TString &fileName, const TString &hName, const TString &newHName = "", bool useCurrentStyle = true);
  //static TH1* readTH1(const TString &fileName, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  //static util::HistVec readTH1(const std::vector<TString> &fileNames, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  //static util::HistVec readHistVec(const TString &fileName, const std::vector<TString> &histNames, const TString &newHistNameSuffix = "", bool useCurrentStyle = true);
  //static util::HistVec readHistVec(const TString &fileName, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  //static std::vector<TF1*> readTF1Vec(const TString &fileName, const TString &fName, const TString &newFName = "");
  //static THStack* readTHStack(const TString &fileName, const TString &stackName, const TString &newStackName = "");
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
#endif





