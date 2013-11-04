// $Id: FileOps.h,v 1.15 2012/04/12 12:55:18 mschrode Exp $

#ifndef FileOps_h
#define FileOps_h

#include <fstream>
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TObject.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"

//#include "utils.h"
//#include "HistOps.h"

TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName);
TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName);
TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName);

int readHistograms(){
  

  TMultiGraph *mg = 0;
  TGraph *JERSigma = 0;
  TGraph *JERRMS = 0;
  TLegend *legend = 0;
  TMultiGraph *mg1 = 0;
  TGraph *JERSigma1 = 0;
  TGraph *JERRMS1 = 0;
  TLegend *legend1 = 0;

  TCanvas *c1 = new TCanvas("c","Jet Energy Resolution for |eta|<0.5",200,10,800,800);
  c1 -> SetLeftMargin(0.17);
  c1 -> cd();

  mg = new TMultiGraph();
  mg -> SetTitle("Jet Energy Resolution (|#eta| < 0.5)");
  

  // mg -> GetTitle() -> SetTitleFontSize(0.1);
  
  JERSigma = readTGraph("plots/data/PF_L1FastJet/root_files/jet_energy_resolution_for_1_eta_bin.root","Graph;1","JERSigma1");  
  JERSigma -> SetMarkerColor(3);
  JERSigma -> SetLineColor(3);

  JERRMS = readTGraph("plots/data/PF_L1FastJet/root_files/jet_energy_resolution_RMS_for_1_eta_bin.root","Graph;1","JERSigma1");
  JERRMS -> SetMarkerColor(4);
  JERRMS -> SetLineColor(4);

  legend  = 0;
  legend = new TLegend(0.55,0.8,0.9,0.9);
  legend -> SetFillColor(0);
  legend -> SetTextSize(0.033);
  legend -> AddEntry(JERSigma,"#sigma_{Gauss}","l");
  legend -> AddEntry(JERRMS,"RMS","l");
  
  mg -> Add(JERSigma);
  mg -> Add(JERRMS);
  mg -> Draw("AP");
  mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma}");
  mg -> GetYaxis() -> SetTitle("JER");
  mg -> GetYaxis() -> SetTitleOffset(2.0); 
  mg -> GetYaxis()->SetRangeUser(0.11,0.3); 
  mg -> GetXaxis()->SetRangeUser(0,300); 
  legend -> Draw("same");

  c1 -> SaveAs("plots/data/PF_L1FastJet/pdfs/jet_energy_resolution_for_1_eta_bin_MG.pdf");

  

  //-------------------------------------------------------

  TCanvas *c2 = new TCanvas("c2","Jet Energy Scale for 0.5<|eta|<1.3",200,10,800,800);
  c2 -> SetLeftMargin(0.17);
  c2 -> cd();
  

  mg1 = new TMultiGraph();
  mg1 -> SetTitle("Jet Energy Resolution (0.5 < |#eta| < 1.3)");
  

  JERSigma1 = readTGraph("plots/data/PF_L1FastJet/root_files/jet_energy_resolution_for_2_eta_bin.root","Graph;1","JERSigma2");  
  JERSigma1 -> SetMarkerColor(3);
  JERSigma1 -> SetLineColor(3);

  JERRMS1 = readTGraph("plots/data/PF_L1FastJet/root_files/jet_energy_resolution_RMS_for_2_eta_bin.root","Graph;1","JERSigma2");
  JERRMS1 -> SetMarkerColor(4);
  JERRMS1 -> SetLineColor(4);

  legend1  = 0;
  legend1 = new TLegend(0.55,0.8,0.9,0.9);
  legend1 -> SetFillColor(0);
  legend1 -> SetTextSize(0.033);
  legend1 -> AddEntry(JERSigma1,"#sigma_{Gauss}","l");
  legend1 -> AddEntry(JERRMS1,"RMS","l");
 
  mg1 -> Add(JERSigma1);
  mg1 -> Add(JERRMS1);
  mg1 -> Draw("AP");
  mg1 -> GetXaxis() -> SetTitle("p_{T}^{#gamma}");
  mg1 -> GetYaxis() -> SetTitle("JER");
  mg1 -> GetYaxis() -> SetTitleOffset(2.0); 
  mg1 -> GetYaxis()->SetRangeUser(0.11,0.3); 
  mg1 -> GetXaxis()->SetRangeUser(0,300); 
  legend1 -> Draw("same");

  c2 -> SaveAs("plots/data/PF_L1FastJet/pdfs/jet_energy_resolution_for_2_eta_bin_MG.pdf");
  
  
    

  return 0;
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
}*/
#endif





