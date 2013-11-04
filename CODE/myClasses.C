#include "myDeclarations.h"
#include "myClasses.h"
#include "myFunctions.h"
#include "TF1.h"
#include "TGraphErrors.h"


    
CResponse::CResponse(int typeResponse){
     
  sigma_array = 0.5;
  sigmaAUX = 0.5;
  sigma_error = 0;
  mean_array = 0.0;
  mean_error = 0;
  alpha_array = 0;
  alpha_error = 0;
  pT_array = 0;
  pT_error = 0;
  RMS_array = 0;
  RMS_error = 0;
  chi2 = 0;
  ndof = 0;
  chi2_ndof = 0;
  meanAUX = 1.;
  meanAUX2 = 1.;
  sigma = 0;
  
  if(typeResponse==1)      hResponse  = createTH1(hResponse,"","Response",200,0,2,"p_{T}^{JetReco}/p_{T}^{#gamma}");
  else if(typeResponse==2) hResponse  = createTH1(hResponse,"","Intrinsic Response",200,0,2,"p_{T}^{rekon. Jet}/p_{T}^{wahrer Jet}");
  else if(typeResponse==3) hResponse  = createTH1(hResponse,"","Imbalance Response",200,0,2,"p_{T}^{JetGen}/p_{T}^{#gamma}");
  hAlpha = createTH1(hAlpha,"","p_{T}^(Jet_{2})/p_{T}^{#gamma}",1000,0,alphaBin[alpha_int],"alpha");
  hPt = createTH1(hPt,"","p_{T}^{#gamma} - Distribution",10000,0,1500,"p^{#gamma}_{T}");
  hPtLeadingJet = createTH1(hPt,"","p_{T}^{#gamma} - Distribution",10000,0,1500,"p^{Jet1}_{T}");
  
};
 

CResponse::~CResponse(){
  delete hResponse;
  delete hAlpha;
  delete hPt;   
};

void CResponse::calculatePt(){
  
  pT_array    = hPt    -> GetMean(1);
  pT_error    = hPt    -> GetMeanError(1);
}

void CResponse::calculate(){
   
  alpha_array = hAlpha -> GetMean(1);
  alpha_error = hAlpha -> GetMeanError(1);
  
  RMS_array = hResponse -> GetRMS(1);
  RMS_error = hResponse -> GetRMSError(1);
  
  TF1* f1 = new TF1("gauss","gaus",0,2); 
  
  // meanAUX  = hResponse -> GetBinCenter(hResponse->GetMaximumBin());

  meanAUX  = hResponse -> GetMean(1);
  sigmaAUX = hResponse -> GetRMS(1);

  float eps = 1.0;
  int l = 0;
  while(eps > 0.0001 && l < 3){
    
    sigma_array = sigmaAUX;
       
    hResponse -> Fit("gauss","QIL","",meanAUX-1.5*sigma_array,meanAUX+1.5*sigma_array);  
    //hResponse -> Fit("gauss","QIL","",meanAUX-2.0*sigma_array,meanAUX+2.0*sigma_array);  
      
    sigmaAUX = f1 -> GetParameter(2);
    meanAUX  = f1 -> GetParameter(1);
    
    if (sigma_array>=sigmaAUX) eps = sigma_array/sigmaAUX-1.;
    else                       eps = sigmaAUX/sigma_array -1.;
   
    l++;    
  }

  int lowBin  = 0;
  int highBin = 0;

  
  // 2.) Calculate the truncated standard deviation 

  for(int i = 0; i<hResponse->GetNbinsX(); i++){
    
    if(hResponse->Integral(hResponse->FindBin(meanAUX)-i, hResponse->FindBin(meanAUX)+i)/hResponse->Integral() >= 0.99){

      lowBin  =  hResponse->FindBin(meanAUX)-i;
      highBin =  hResponse->FindBin(meanAUX)+i;
      break;
    }
  }
  
  hResponse -> GetXaxis() -> SetRange(lowBin,highBin);

  sigma_array = (hResponse->GetRMS(1))/(hResponse->GetMean(1));
  sigma_error = TMath::Sqrt(TMath::Power(1./hResponse->GetMean(1),2)*TMath::Power(hResponse->GetRMSError(1),2) + TMath::Power(hResponse->GetRMS(1)/TMath::Power(hResponse->GetMean(1),2),2)*TMath::Power(hResponse->GetMeanError(1),2));
    
  mean_array = hResponse -> GetMean(1);
  mean_error = hResponse -> GetMeanError(1);

  //sigma_array = (f1 -> GetParameter(2))/(f1 -> GetParameter(1));
  //sigma_error = TMath::Sqrt(TMath::Power(1./f1->GetParameter(1),2)*TMath::Power(f1->GetParError(2),2) + TMath::Power(f1->GetParameter(2)/TMath::Power(f1->GetParameter(1),2),2)*TMath::Power(f1->GetParError(1),2));
    
  //mean_array  = f1 -> GetParameter(1);
  //mean_error  = f1 -> GetParError(1);
  
  //chi2 = f1 -> GetChisquare();
  //ndof = f1 -> GetNDF();  
  //chi2_ndof = chi2/(1.0*ndof);
  
  //delete f1; 
};

//----
CScaleResAlpha::CScaleResAlpha(){
  
  pT = 0;
  pTError = 0;  
  scale = 0;
  scaleError = 0;
  resolution = 0;
  resolutionError = 0;
  q = 0;
  c = 0;
  m = 0;
  qprime = 0;
  qprimeError = 0;
  cprime = 0;
  cprimeError = 0;
  mprime = 0;
  mprimeError = 0;
  
  for(int i=0; i<2*alpha_int; i++){
    alpha[i] = 0;
    alphaError[i] = 0;
    mean[i] = 0;
    meanError[i] = 0;
    sigma[i] = 0;
    sigmaError[i] = 0; 
  }
}


void CScaleResAlpha::calculate(int length, int fit){
  
 
  if(length>=1){
   
    gJetScaleAlpha = new TGraphErrors(length, alpha , mean, alphaError, meanError);
    gJetResolutionAlpha = new TGraphErrors(length, alpha , sigma, alphaError, sigmaError);
    
    if(fit == 2){  //mc intrinsic
      fScaleAlpha = new TF1("fScaleAlpha","[0]",0,alphaBin[alpha_int]);
      gJetScaleAlpha -> Fit("fScaleAlpha","QR");
      c = fScaleAlpha -> GetParameter(0);
      scale      = fScaleAlpha -> GetParameter(0);
      scaleError = fScaleAlpha -> GetParError(0);
      
      fResolutionAlpha = new TF1("fResolutionAlpha","[0]",0,alphaBin[alpha_int]); 
      gJetResolutionAlpha -> Fit("fResolutionAlpha","QR");
      cprime      = fResolutionAlpha -> GetParameter(0);
      cprimeError = fResolutionAlpha -> GetParError(0);
      resolution      = fResolutionAlpha -> GetParameter(0);
      resolutionError = fResolutionAlpha -> GetParError(0);

      cout<<"Intrinsic = "<<endl;
      cout<<"cprime = "<<cprime<<endl;
      cout<<"cprimeError = "<<resolutionError<<endl<<endl;      
    }
    else if(fit == 3){   //mc imbalance
      
      fScaleAlpha = new TF1("fScaleAlpha","1. - [0] - [1]*TMath::Power(x,2)",0,alphaBin[alpha_int]); 
      gJetScaleAlpha -> Fit("fScaleAlpha","QR");
      q = fScaleAlpha -> GetParameter(0);
      qError = fScaleAlpha -> GetParError(0);
      m = fScaleAlpha -> GetParameter(1);
      scale      = 1 - fScaleAlpha -> GetParameter(0);
      scaleError = fScaleAlpha -> GetParError(0);
     
      fResolutionAlpha = new TF1("fResolutionAlpha","[0] + [1]*x",0,alphaBin[alpha_int]); 
      gJetResolutionAlpha -> Fit("fResolutionAlpha","QR");
      qprime = fResolutionAlpha -> GetParameter(0);
      qprimeError = fResolutionAlpha -> GetParError(0);
      mprime = fResolutionAlpha -> GetParameter(1);
      resolution      = fResolutionAlpha -> GetParameter(0);
      resolutionError = fResolutionAlpha -> GetParError(0);
      
      cout<<"qprime = "<<qprime<<endl;
      cout<<"qprimeError = "<<qprimeError<<endl<<endl;;

    }
    else if(fit == 4){ //data and MC
      
      fScaleAlpha      = new TF1("fScaleAlpha","[0]*(1. - [1] - [2]*TMath::Power(x,2))",0,alphaBin[alpha_int]); 
      fScaleAlpha -> FixParameter(1,q);
      gJetScaleAlpha -> Fit("fScaleAlpha","QR");    
      scale      = fScaleAlpha -> GetParameter(0);
      scaleError = fScaleAlpha -> GetParError(0);
      
      fResolutionAlpha = new TF1("fResolutionAlpha","TMath::Sqrt(TMath::Power([0],2) +TMath::Power([1],2) +2*[1]*[2]*x + TMath::Power(([2]*x),2) )",0,alphaBin[alpha_int]); 
      
      if(qprime<0){

	cout<<"----------------------------------------------------- qprime negative!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	cout<<"qprime = "<<qprime<<endl;
	cout<<"qprimeError = "<<qprimeError<<endl;
	qprime = 0;
      }

      fResolutionAlpha ->FixParameter(1,qprime);        
      gJetResolutionAlpha -> Fit("fResolutionAlpha","QRB");

      
      if(fResolutionAlpha -> GetParameter(2) < 0){	
	if(fResolutionAlpha -> GetParameter(2) + fResolutionAlpha -> GetParError(2) > 0){

	  fResolutionAlpha ->FixParameter(2,0.0);   

	}
	else{
	  cout<<"----------------------------------------------------Fit Paramter mprime < 0 -> Set bounds on it!!!!"<<endl;
	  cout<<"mprime = "<<fResolutionAlpha -> GetParameter(2)<<endl;
	  cout<<"mprimeError = "<<fResolutionAlpha -> GetParError(2)<<endl;
	  fResolutionAlpha -> SetParLimits(2,0,10000000);
	  fResolutionAlpha -> SetParameter(2,0.001);
	  if(type == 2){
	    fResolutionAlpha -> SetParameter(0,cprime);
	  }
	  else if(type == 1){
	    fResolutionAlpha -> SetParameter(0,0.05);
	  }
	  gJetResolutionAlpha -> Fit("fResolutionAlpha","QRB");	
	}
      }

      if(fResolutionAlpha -> GetParameter(0) < 0){

	cout<<"----------------------------------------------------Fit Paramter cprime < 0 -> Set bounds on it!!!!"<<endl;

	fResolutionAlpha -> SetParLimits(2,0,10000000);
	fResolutionAlpha -> SetParLimits(0,0,10000000);
	
	if(type == 2){
	  fResolutionAlpha -> SetParameter(0,cprime);
	  fResolutionAlpha -> SetParameter(2,mprime);
	}
	else if(type == 1){
	  fResolutionAlpha -> SetParameter(0,0.05);
	  fResolutionAlpha -> SetParameter(2,0.001);
	}
	
	gJetResolutionAlpha -> Fit("fResolutionAlpha","QRB");	

      }
      

      resolution      = fResolutionAlpha -> GetParameter(0);
      resolutionError = fResolutionAlpha -> GetParError(0);
      
      cout<<"All = "<<endl;
      cout<<"Chi2 = "<<fResolutionAlpha -> GetChisquare()<<endl;
      cout<<"ndof = "<<fResolutionAlpha -> GetNDF()<<endl;
  
      cout<<"qprime = "<<fResolutionAlpha -> GetParameter(1)<<endl;
      cout<<"qprime Error = "<<fResolutionAlpha -> GetParError(1)<<endl;
      cout<<"cprime = "<<fResolutionAlpha -> GetParameter(0)<<endl;
      cout<<"cprime Error = "<<fResolutionAlpha -> GetParError(0)<<endl;
      cout<<"mprime = "<<fResolutionAlpha -> GetParameter(2)<<endl;
      cout<<"mprime Error = "<<fResolutionAlpha -> GetParError(2)<<endl<<endl<<endl;      
    }


    delete fScaleAlpha;
    delete fResolutionAlpha; 
  }
  
};

CScaleResAlpha::~CScaleResAlpha(){
    delete gJetScaleAlpha;
    delete gJetResolutionAlpha;  
};

//----
CScaleRes::CScaleRes(){
  
  scalenumber = 0;
  scalenumberError = 0;
  for(int i=0; i<pt_int; i++){
    
    pT[i] = 0;
    pTError[i] = 0;
    scale[i] = 0;
    scaleError[i] = 0;
    resolution[i] = 0;
    resolutionError[i] = 0;
  }
};

void CScaleRes::calculate(int length){


  gResolution = new TGraphErrors(length, pT , resolution, pTError, resolutionError);         
  gScale      = new TGraphErrors(length, pT , scale, pTError, scaleError);
  gq          = new TGraphErrors(length, pT , q, pTError, qError);
  gqprime     = new TGraphErrors(length, pT , qprime, pTError, qprimeError);
  
  
  fScale = new TF1("fScale","[0]",0,600);
  
  gScale -> Fit("fScale","QR");
  scalenumber      = fScale -> GetParameter(0); 
  scalenumberError = fScale -> GetParError(0);
  
  fResolution = new TF1("fResolution","TMath::Sqrt(TMath::Sign(1,[0])*TMath::Power([0]/x,2)+TMath::Power([1],2)*TMath::Power(x,[3]-1)+TMath::Power([2],2))", pT[0], 600);

  

  fResolution->SetParName(0,"N");
  fResolution->SetParName(1,"S");
  fResolution->SetParName(2,"C");
  fResolution->SetParName(3,"m");
  
  fResolution->FixParameter(2,0);
  fResolution->SetParameter("N",-1.13); 
  fResolution->SetParameter("S",0.611); 
  //fResolution->SetParameter("C",0.03121); 
  fResolution->SetParameter("m",0.143); 
  
  gResolution->Fit("fResolution","QR");
  
  N    = fResolution->GetParameter("N");
  NErr = fResolution->GetParError(0);
  //cout<<"N = "<<N<<endl;
  //cout<<"NErr = "<<NErr<<endl;
  S    = fResolution->GetParameter("S");
  SErr = fResolution->GetParError(1);
  //cout<<"S = "<<S<<endl;
  //cout<<"SErr = "<<SErr<<endl;
  C    = fResolution->GetParameter("C");
  CErr = fResolution->GetParError(2);
  //cout<<"C = "<<C<<endl;
  //cout<<"CErr = "<<CErr<<endl;
  m    = fResolution->GetParameter("m");
  mErr = fResolution->GetParError(3);
  //cout<<"m = "<<m<<endl;
  //cout<<"mErr = "<<mErr<<endl<<endl<<endl;
  

 
};
 

CScaleRes::~CScaleRes(){
  delete gScale;
  delete gResolution;
  delete gq;
  delete gqprime;
  delete fResolution;
  delete fScale;   
};



