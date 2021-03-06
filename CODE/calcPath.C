#include "myDeclarations.h"
#include "TROOT.h"
#include "TString.h"
#include <iostream>    

void calcPath(int step){

  TString pathYear, pathType, pathSample; 
  TString command; 
    
   // Determine the path where to save all files
   if(date==2012)      pathYear = "plots_2012/";
   else if(date==2011) pathYear = "plots_2011/";
   else cout<<"Please choose a valid number for the integer \"date\"! (or change the allowed range in the file \"calcPath.C\")"<<endl;
   
   if(type==1)      pathType = "data";
   else if(type==2) pathType = "mc";
   else cout<<"Please choose a valid number for the integer \"type\"! (or change the allowed range in the file \"calcPath.C\")"<<endl;
   
   if(jetType==1)      pathSample  = "PF_L1FastJet/";
   else if(jetType==2) pathSample  = "PF_L1CHS/";
   else if(jetType==3) pathSample  = "Calo_L1FastJet/";
   else cout<<"Please choose a valid number for the integer \"file\"! (or change the allowed range in the file \"calcPath.C\")"<<endl;
   
   if(jetType == 1)      DataType = "_PF_";
   else if(jetType == 2) DataType = "_PFCHS_";
   else if(jetType == 3) DataType = "_Calo_";
   else cout<<"Please choose a valid number for the integer \"file\"! (or change the allowed range in the file \"calcPath.C\")"<<endl;
 
   // The whole path to save in the End the .root and .pdf files 
   PDFPath  = pathYear + pathSample + pathType + "/pdfs/"; 
   RootPath = pathYear + pathSample + pathType + "/root_files/"; 
   
   // Suffix to be determined for the filenames in the end
   DataType += pathType;
  
   // Special files to write strange events in extra .txt files
   filestr.open(pathYear + pathSample + pathType + "RunLumiEventNum.txt");
   EventVariables.open(pathYear + pathSample + pathType + "EventVariables.txt");
   EventVariables<<"GammaPt:JetPt:GammaEta:JetEta:GammaPhi:JetPhi:GammaIsoEcal:GammaIsoHcal:GammaIsoTrk"<<endl;

   // Delete all histogram 
   gDirectory->Delete(); 
   // Delete all canvas
   gROOT->GetListOfCanvases()->Delete();

   
   // Delete former files in the root_files-and pdfs-folders   
   if(step == 1){
     command = ".! rm -r " + pathYear + pathSample + pathType + "/root_files/*";
     gROOT->ProcessLine(command);
     command = ".! rm -r " + pathYear + pathSample + pathType + "/pdfs/*";
     gROOT->ProcessLine(command); 
   }
}
