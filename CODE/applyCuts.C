#include "myDeclarations.h"
#include "myClasses.h"
#include "TVector2.h"
#include "TMath.h"

bool applyCuts(){

 
    
  // Check if array is large enough
  if( nobjJet > NJETS ) {
    std::cerr << "\nERROR: 'nobjJet = " << nobjJet << " > NJETS = " << NJETS << "'" << std::endl;
    return 0;
  }
  if( nobjPhoton > NPHOTONS ) {
    std::cerr << "\nERROR: 'nobjPhoton = " << nobjJet << " > NPHOTONS = " << NPHOTONS << "'" << std::endl;
    return 0;
  }
  
  // Correct the photon Pt to evaluate the systematic uncertainty
  //if(photonEta[0]<1.479) photonPt[0]=photonPt[0]*(1.+0.006);
  //else photonPt[0]=photonPt[0]*(1.+0.014);   // But those are dismissed later on
   
  // Correct Jets and Order them 
  corrJets.clear();
  for(int j = 0; j < nobjJet; j++) {      
    corrJets.add(j,/*(1.-jetUncert[j])*/jetCorrL1[j]*jetCorrL2L3[j]*jetPt[j]);
  }
  
  corrJets.sort();
  
  // Definition of some variables
  
  double_t deta[nobjJet];
  double_t dR[nobjJet];
  //double_t photonpt=0.;
  double_t dphi[nobjJet];
  
  
  int j_jet=0;
  int j_photon = 0;
  lead_jet=-1;
  jet_2 = -1;
  photonidx = -1;
  
  
  // Find out to what index the photon belongs to in the jet sample
  for (int i=0 ; i<nobjJet ; i++) {
    
    if(j_jet==2 && j_photon==1) continue;
    dphi[i] = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(i)]-photonPhi[0]));
    //dphi[i] = std::abs(TVector2::Phi_mpi_pi(jetPhi[i]-photonPhi[0]));
    deta[i] = std::abs(jetEta[corrJets.idx(i)] - photonEta[0]);
    //deta[i] = std::abs(jetEta[i] - photonEta[0]);
    dR[i]   = TMath::Sqrt(dphi[i]*dphi[i] + deta[i]*deta[i]);
    //cout<<"dR["<<i<<"] = "<<dR[i]<<endl;      
    if(dR[i] > 0.5 && j_jet<2) {
      
      if(j_jet==0){
        lead_jet = i;
        j_jet=j_jet+1;
        //cout<<"lead_jet = "<<lead_jet<<endl;
      }       
      else if(j_jet==1){
	jet_2=i;
	j_jet=j_jet+1;
      }
      else cout<<"something wrong here"<<endl;
      
    }
    else if(dR[i]<=0.5){
      //photonpt = corrJets.pt(i);   
      photonidx = i; 
      //cout<<"photonidx = "<<photonidx<<endl;
      j_photon=j_photon+1;              
    }
  }

  //---------------------------------------------------------------------------------------------
  // Fill histogram for number of Photons per event
  hPhotonN->Fill(nobjPhoton,weight*PUWeight);
  
  //---------------------------------------(1. CUT)--
  // Select events with at least 1 Photon (1. CUT) 
  if( nobjPhoton < 1 ){
    cut1 = cut1 +1;
    //cout<<"Discarded because there is no photon in this event!"<<endl;
    return 0;}
  //-------------------------------------------------
  
  int num_tight=0;
  // Fill histogram with tight Photons per event
  for(int i=0; i<nobjPhoton;i++){
    if(tight[i]) num_tight=num_tight+1; 
  }
  
  hPhotonNtight->Fill(num_tight,weight*PUWeight);
  
   
  //---------------------------------------(2. CUT)--
  // Select only events with tight leading Photon (2. CUT)
  if(!tight[0]){
    //cut2 = cut2 +1;
    //cout<<"Discarded because leading Photon is not tight!"<<endl;
    return 0;
  }
  //-------------------------------------------------  
  
  
  // Fill Histogram for Eta plot
  hPhotonEta->Fill(photonEta[0],weight*PUWeight);    
  
  // Fill Pt histograms of first Photon (only tight Photons are used) 
  hPhotonT1Pt -> Fill(photonPt[0],weight*PUWeight);
  
  // Fill Histogram for phi plot
  hPhotonPhi->Fill(photonPhi[0],weight*PUWeight);
  
  //----------------------------JETS---------------------------------------------
  
  // Fill hsitogram with number of jets
  hJetN->Fill(nobjJet,weight*PUWeight);
  
  //---------------------------------------(3. CUT)--
  if(lead_jet==-1) {
    cut3 = cut3 +1;
    //cout<<"Discarded because there is no balancing Jet in the Jet Sample!"<<endl;
    return 0;   
  }
  //------------------------------------------------- 



  //---------------------------------------(2. CUT)--
  // Select only events with tight leading Photon (2. CUT)
  if(jet_2 != -1){

    jetPt2nd = corrJets.pt(jet_2);
    /*
    if(corrJets.pt(jet_2) < 11){
      
      jetPt2nd = 11;
      //cut2 = cut2 +1;
      //return 0;
    }
    else jetPt2nd = corrJets.pt(jet_2);
    */
  }
  else jetPt2nd = 0;
  
  //-------------------------------------------------    
  
  int detaGen = 0;
  genJetidx = -1;
  gen2ndJetidx = -1;
  
  float deltaR1stJet;
  float deltaPhi1stJet;
  float deltaEta1stJet;
  
  int genJet1stidx = -1;

  jetPt1stJet = 0.;
  float c = 1.;
  int test = 0;
  bool jet1 = true;
  bool jet2 = true;
  if(type==2){
    for (int i=0 ; i<nobjGenJet ; i++) {
      
      if(genJetColJetIdx[i] == corrJets.idx(lead_jet) && jet1){
      	genJetidx = i;
	test += 1;
	jet1 = false;
      }
      
      if(jet_2 != -1){
	if(genJetColJetIdx[i] == corrJets.idx(jet_2) && jet2){
	  gen2ndJetidx = i;
	  test += 1;
	  jet2 = false;
	}
      }
      
      if(test == 2) break;
      
    }


    // MC smearing Code piece
    /*
    if(std::abs(jetEta[corrJets.idx(lead_jet)])<0.5) c = 1.175;
    else if(std::abs(jetEta[corrJets.idx(lead_jet)])<1.1 && std::abs(jetEta[corrJets.idx(lead_jet)])>0.5) c = 1.150;
    else if(std::abs(jetEta[corrJets.idx(lead_jet)])<1.7 && std::abs(jetEta[corrJets.idx(lead_jet)])>1.1) c = 1.206;
    else if(std::abs(jetEta[corrJets.idx(lead_jet)])<2.3 && std::abs(jetEta[corrJets.idx(lead_jet)])>1.7) c = 1.430;
    else c = 1.000;
    */
  
    // Smear Second Jet
    /*
    if(jet_2 != -1 && gen2ndJetidx != -1){
      if(std::abs(jetEta[corrJets.idx(jet_2)])<0.5) c = 1.133;
      else if(std::abs(jetEta[corrJets.idx(jet_2)])<1.1 && std::abs(jetEta[corrJets.idx(jet_2)])>0.5) c = 1.083;
      else if(std::abs(jetEta[corrJets.idx(jet_2)])<1.7 && std::abs(jetEta[corrJets.idx(jet_2)])>1.1) c = 1.145;
      else if(std::abs(jetEta[corrJets.idx(jet_2)])<2.3 && std::abs(jetEta[corrJets.idx(jet_2)])>1.7) c = 1.288;
      else c = 1.000;
      jetPt2nd = genJetPt[gen2ndJetidx] + c*(jetPt2nd - genJetPt[gen2ndJetidx]);
    }

    // Smear first Jet
    if(std::abs(jetEta[corrJets.idx(lead_jet)])<0.5) c = 1.133;
    else if(std::abs(jetEta[corrJets.idx(lead_jet)])<1.1 && std::abs(jetEta[corrJets.idx(lead_jet)])>0.5) c = 1.083;
    else if(std::abs(jetEta[corrJets.idx(lead_jet)])<1.7 && std::abs(jetEta[corrJets.idx(lead_jet)])>1.1) c = 1.145;
    else if(std::abs(jetEta[corrJets.idx(lead_jet)])<2.3 && std::abs(jetEta[corrJets.idx(lead_jet)])>1.7) c = 1.288;
    else c = 1.000;
    */

    jetPt1stJet = genJetPt[genJetidx] + c*(corrJets.pt(lead_jet) - genJetPt[genJetidx]);
    intrinsic = jetPt1stJet/genJetPt[genJetidx];
    imbalance = genJetPt[genJetidx]/photonPt[0];
    
  }
  else if(type == 1) jetPt1stJet = corrJets.pt(lead_jet);
  
  // Calculate response    
  response  = jetPt1stJet/photonPt[0]; 
  
  //---------------------------------------(5. CUT)--
  // Take only tight Jets 
  if(!JetTight[corrJets.idx(lead_jet)]){
    cut5 = cut5 +1;
    // cout<<"Discarded because of tight Jet ID "<<endl;
    return 0;}
  //-------------------------------------------------
  
  
  // Take only events with vtxN>1 and plot eta spectrum again
  if(vtxN>1) hJetEta_1->Fill(jetEta[corrJets.idx(lead_jet)],weight*PUWeight);
  
  hJetEtaPt1->Fill(jetEta[corrJets.idx(lead_jet)],corrJets.pt(lead_jet),weight*PUWeight);
  
  
  //-------------------------------------(6. CUT)----
  if(std::abs(jetEta[corrJets.idx(lead_jet)])>5.0){  
    cut6 = cut6 +1;
    //cout<<"Discarded because of Jet eta!"<<endl; 
    return 0;
  }
  //-------------------------------------------------
  
  
  
  // Fill ptJet2/ptPhoton ratio
  if(jet_2!=-1) hRatioPt->Fill(corrJets.pt(jet_2)/photonPt[0],weight*PUWeight);
  
  
  //------------------------ VARIABLES from PHOTONS + JETS --------------------------------------
  
  // Fill dphi histogram
  float deltaphi=0;
  deltaphi=std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(lead_jet)]-photonPhi[0]));
  hDphi->Fill(deltaphi,weight*PUWeight);
  
  
  
  
  //-------------------------------------(8. CUT)----
  //if(corrJets.pt(lead_jet) < 10.){
  if(jetPt1stJet < 11.){
    cut8 = cut8 +1;
    //cout<<"Discarded because of too low jet pt (smaller than 10 GeV)!"<<endl; 
    return 0;}
  //-------------------------------------------------
  
  // Fill response Pt 2d-histogram
  hResponsePtBeforeEcal->Fill(corrJets.pt(lead_jet)/photonPt[0],photonPt[0],weight*PUWeight);
  
  //--------------------------------------(9. CUT)---
  if(!EcalBEFilter || !EcalTPFilter){
    cut9 = cut9+1;
    // cout<<"Discarded because of filters!"<<endl;
    return 0;}
  //-------------------------------------------------
  
  hResponsePhotonEta -> Fill(response, photonEta[0],weight*PUWeight);    
  
  //--------------------------------------(10. CUT)----
  if(std::abs(photonEta[0])>1.3){  
    cut10 = cut10+1;  
    //cout<<"Discarded because of photon etacut!"<<endl; 
    return 0;
  }
  //-------------------------------------------------

  
  
  
  // Fill response Pt 2d-histogram
  hResponsePt->Fill(corrJets.pt(lead_jet)/photonPt[0],photonPt[0],weight*PUWeight);
  
  if(std::abs(jetEta[corrJets.idx(lead_jet)]-photonEta[0])<=0.5) hResponsePt_0_5_Deta->Fill(corrJets.pt(lead_jet)/photonPt[0],photonPt[0],weight*PUWeight);
  if(std::abs(jetEta[corrJets.idx(lead_jet)]-photonEta[0])<=1.1) hResponsePt_1_1_Deta->Fill(corrJets.pt(lead_jet)/photonPt[0],photonPt[0],weight*PUWeight);
  if(std::abs(jetEta[corrJets.idx(lead_jet)]-photonEta[0])<=1.5) hResponsePt_1_5_Deta->Fill(corrJets.pt(lead_jet)/photonPt[0],photonPt[0],weight*PUWeight);

  
  float JetEnergySum=0.0;
  float VecJetEnergySumX=0.0;
  float VecJetEnergySumY=0.0;
  float VecJetEnergySum=0.0;
  int jetnumber=0;
  
  if(nobjJet==2) jetnumber=2;
  else if(nobjJet==3)  jetnumber=3;
  else if(nobjJet==4)  jetnumber=4;
  else if(nobjJet==5)  jetnumber=5;
  else if(nobjJet>=6)  jetnumber=6;
  
  // Calculate The scalar and vectorial Pt Sum
  for(int i=0; i<jetnumber; i++){
    JetEnergySum     = JetEnergySum+corrJets.pt(i);
    VecJetEnergySumX = VecJetEnergySumX+corrJets.pt(i)*sin(jetPhi[corrJets.idx(i)]);
    VecJetEnergySumY = VecJetEnergySumY+corrJets.pt(i)*cos(jetPhi[corrJets.idx(i)]);   
  }
  
  // Substract Energy from Photon (just in case there is a photon)
  if(photonidx!=-1){
    JetEnergySum=JetEnergySum-corrJets.pt(photonidx); 
    VecJetEnergySumX=VecJetEnergySumX-corrJets.pt(photonidx)*sin(jetPhi[corrJets.idx(photonidx)]);
    VecJetEnergySumY=VecJetEnergySumY-corrJets.pt(photonidx)*cos(jetPhi[corrJets.idx(photonidx)]);    
  }
  
  // Substract energy from leading Jet
  JetEnergySum=JetEnergySum-corrJets.pt(lead_jet);
  VecJetEnergySumX=VecJetEnergySumX-corrJets.pt(lead_jet)*sin(jetPhi[corrJets.idx(lead_jet)]);
  VecJetEnergySumY=VecJetEnergySumY-corrJets.pt(lead_jet)*cos(jetPhi[corrJets.idx(lead_jet)]);    
  VecJetEnergySum=TMath::Sqrt(TMath::Power(VecJetEnergySumX,2)+TMath::Power(VecJetEnergySumY,2));
  
  float JetEnergyFraction=0.0;
  float VecJetEnergyFraction=0.0;
  JetEnergyFraction=JetEnergySum/corrJets.pt(lead_jet);
  VecJetEnergyFraction=VecJetEnergySum/corrJets.pt(lead_jet);
  
  // Fill Diagram with Response and Fraction JetEnergySum/corrJets.pt(lead_jet)
  hResponseEnergyFraction->Fill(corrJets.pt(lead_jet)/photonPt[0],JetEnergyFraction,weight*PUWeight);
  hResponseVecEnergyFraction->Fill(corrJets.pt(lead_jet)/photonPt[0],VecJetEnergyFraction,weight*PUWeight);

  
  
  
  //------------------------------------------------------------------------------------------------------------------------------------------------------------------------       
  
  //Get rid of strange events (just for 2011 data)
  if(date == 2011 && photonidx!=-1 && corrJets.pt(photonidx)/photonPt[0]<0.5) return 0;
  
    
  alpha = jetPt2nd/photonPt[0] * 100.;
  //else alpha = 0;
  
  //Photon Pt against Number of Vertices
  // Fill response Pt 2d-histogram
  hPhotonPtVtx->Fill(photonPt[0],vtxN,weight*PUWeight);

  
  // Photon Pt against Number of Vertices for every PT Bin
  for(int i=0; i<numTrigger-2; i++){
    if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){ 
      hPhotonPtVtxBinned[i] ->Fill(photonPt[0],vtxN,weight*PUWeight);
      break;
    }
    else if(photonPt[0]>=bd[numTrigger-2]){
      hPhotonPtVtxBinned[numTrigger-2]->Fill(photonPt[0],vtxN,weight*PUWeight); 
      break;
    }
  }

  // Photon Pt for different pt bins
  for(int i=0; i<numTrigger-1; i++){
    if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){
      hPhotonPt[i] ->Fill(photonPt[0],weight*PUWeight);
      break;
    }
    else if(photonPt[0]>=bd[numTrigger-1]){
      hPhotonPt[numTrigger-1]->Fill(photonPt[0],weight*PUWeight); 
      break;
    }
  }

  for(int i=0; i<numTrigger-1; i++){
    if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){
      hTriggerEffBefore[i] ->Fill(vtxN,weight*PUWeight);
      break;
    }
    else if(photonPt[0]>=bd[numTrigger-1]){
      hTriggerEffBefore[numTrigger-1] ->Fill(vtxN,weight*PUWeight);
      break;
    }
  }

  for(int i=0; i<numTrigger-1; i++){
    if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){
      hRho[i] ->Fill(rho,weight*PUWeight);
      break;
    }
    else if(photonPt[0]>=bd[numTrigger-1]){
      hRho[numTrigger-1] ->Fill(rho,weight*PUWeight);
      break;
    }
  }

  if((photonPt[0] >= bd[0] && photonPt[0] < bd[1]) && !hltPhoton[0]){
    cut12[0] += 1; 
    //cout<<"Trigger 20!"<<endl;
    //cout<<"PhotonPt[0] = "<<photonPt[0]<<endl<<endl;
  }
  if((photonPt[0] >= bd[1] && photonPt[0] < bd[2]) && !hltPhoton[1]){
    cut12[1] += 1; 
    //cout<<"Trigger 30!"<<endl;
    //cout<<"PhotonPt[0] = "<<photonPt[0]<<endl<<endl;
  }
  if((photonPt[0] >= bd[2] && photonPt[0] < bd[3]) && !hltPhoton[2]){
    cut12[2] += 1; 
    //cout<<"Trigger 50!"<<endl;
    //cout<<"PhotonPt[0] = "<<photonPt[0]<<endl<<endl;
  }
  if((photonPt[0] >= bd[3] && photonPt[0] < bd[4]) && !hltPhoton[3]){
    cut12[3] += 1;
    //cout<<"Trigger 75!"<<endl;
    //cout<<"PhotonPt[0] = "<<photonPt[0]<<endl<<endl;
  } 
  if((photonPt[0] >= bd[4] && photonPt[0] < bd[5]) && !hltPhoton[4]){
    cut12[4] += 1;
    //cout<<"Trigger 90!"<<endl;
    //cout<<"PhotonPt[0] = "<<photonPt[0]<<endl<<endl;
  }   
  if(date == 2012){
    if((photonPt[0] >= bd[5] && photonPt[0] < bd[6]) && !hltPhoton[5]){
      cut12[5] += 1; 
      //cout<<"Trigger 135!"<<endl;
      //cout<<"PhotonPt[0] = "<<photonPt[0]<<endl<<endl;
    }
    if((photonPt[0] >= bd[6] && photonPt[0] < bd[7]) && !hltPhoton[6]){
      cut12[6] += 1;     
      //cout<<"Trigger 150!"<<endl;
      //cout<<"PhotonPt[0] = "<<photonPt[0]<<endl<<endl;
    }
    if((photonPt[0] >= bd[7]) && !hltPhoton[7]){
      cut12[7] += 1; 
      //cout<<"Trigger 160!"<<endl;
      //cout<<"PhotonPt[0] = "<<photonPt[0]<<endl<<endl;
    }
  }
  // --------------------------------------------------------------------------------------------------------------------------
  // Cuts on Trigger
  
  /* 
  if(type == 1 || type == 2){
    if((photonPt[0] >= bd[0] && photonPt[0] < bd[1]) && !hltPhoton[0]) return 0;
    if((photonPt[0] >= bd[1] && photonPt[0] < bd[2]) && !hltPhoton[1]) return 0;
    if((photonPt[0] >= bd[2] && photonPt[0] < bd[3]) && !hltPhoton[2]) return 0;
    if((photonPt[0] >= bd[3] && photonPt[0] < bd[4]) && !hltPhoton[3]) return 0;
    if((photonPt[0] >= bd[4] && photonPt[0] < bd[5]) && !hltPhoton[4]) return 0;
    
    if(date == 2012){
      if((photonPt[0] >= bd[5] && photonPt[0] < bd[6]) && !hltPhoton[5]) return 0;
      if((photonPt[0] >= bd[6] && photonPt[0] < bd[7]) && !hltPhoton[6]) return 0;    
      if((photonPt[0] >= bd[7]) && !hltPhoton[7]) return 0;
    }
  }
  */
  
  
  
  
  // --------------------------------------------------------------------------------------------------------------------------
  
  

  if(photonPt[0] >= bd[7] && jetEta[corrJets.idx(lead_jet)]) hAlpha->Fill(jetPt2nd/photonPt[0],weight*PUWeight);
  //-------------------------------------------(11. CUT)------
  if(jet_2!=-1){
    if(jetPt2nd > alphaBin[alpha_int]/100. * photonPt[0]){
      cut11=cut11+1;
      return 0;
      // cout<<"Discarded because of alpha constraint. alpha = "<<alpha <<endl;
    }
  }
  //---------------------------------------------------------

  if(photonPt[0] >= bd[4] && photonPt[0] < bd[5] && jetEta[corrJets.idx(lead_jet)] < etaBin[1]) hDeltaPhi->Fill(deltaphi,weight*PUWeight);
  //----------------------------------------(7. CUT)-
  // Select events with back-to-back photon and jet
  if( deltaphi < 2.7){
    cut7 = cut7 +1;
    //  cout<<"Discarded because of delta phi!"<<endl;
    return 0;}
  //-------------------------------------------------

  

  
  // massless particles: Et = Pt
  /*
  if(photonIsoEcal[0]> 5.0 + 0.012*photonPt[0]){
    //cout<<"because Ecal thrown out."<<endl;
    cut13[0] +=1; 
    return 0;
  }
  if(photonIsoHcal[0]> 3.0 + 0.005*photonPt[0]){
    cut13[1] +=1;
    //cout<<"because Hcal thrown out."<<endl;
    return 0;
  }
  if(photonIsoTrk[0]> 3.0 + 0.002*photonPt[0]){
    cut13[2] +=1;
    //cout<<"because Trk thrown out."<<endl;
    return 0;
  }
  */
  
  
  // For systematic PU uncertainty 
  // Temporary additional cuts
  
  //if(vtxN>10) return 0;
  //if(vtxN<=10 || vtxN>15) return 0;
  //if(vtxN<=15 || vtxN>20) return 0;
  //if(vtxN<=20) return 0;

  for(int i=0; i<numTrigger-1; i++){
    if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){
      hTriggerEffAfter[i] ->Fill(vtxN,weight*PUWeight);
      break;
    }
    else if(photonPt[0]>=bd[numTrigger-1]){
      hTriggerEffAfter[numTrigger-1] ->Fill(vtxN,weight*PUWeight);
      break;
    }
  }

  // Number of Vertices for different pt bins  
  for(int i=0; i<numTrigger-1; i++){
    if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){
      hVtxPtbinned[i] ->Fill(vtxN,weight*PUWeight);
      hRho[i] -> Fill(rho,weight*PUWeight);
      break;
    }
    else if(photonPt[0]>=bd[numTrigger-1]){
      hVtxPtbinned[numTrigger-1]->Fill(vtxN,weight*PUWeight); 
      hRho[numTrigger-1] -> Fill(rho,weight*PUWeight);
      break;
    }
  }
    
 
  if(vtxN<40) hRhoVtxBinned[vtxN] -> Fill(rho,weight*PUWeight);
 
  hPhotonIsoEcal -> Fill(photonIsoEcal[0],weight*PUWeight);
  hPhotonIsoHcal -> Fill(photonIsoHcal[0],weight*PUWeight);
  hPhotonIsoTrk  -> Fill(photonIsoTrk[0],weight*PUWeight);  
  /*
    if(photonPt[0] >= bd[pt_int-1] && photonPt[0] < bd[pt_int]){

    for(int k=0; k<alpha_int; k++){

    deltaphi=0;
    deltaphi=std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(jet_2)]-photonPhi[0]));

    if(alpha >= alphaBin[k] && alpha<alphaBin[k+1]){

	
    if( deltaphi > TMath::Pi()/2 ){
      
    h2ndJetPt1stJetHemisphere[k] -> Fill(jetPt[corrJets.idx(jet_2)],weight*PUWeight);
    h2ndJetEta1stJetHemisphere[k] -> Fill(jetEta[corrJets.idx(jet_2)],weight*PUWeight);
      
	  
    }
    else{
	  
    h2ndJetPtPhotonHemisphere[k] -> Fill(jetPt[corrJets.idx(jet_2)],weight*PUWeight);
    h2ndJetEtaPhotonHemisphere[k] -> Fill(jetEta[corrJets.idx(jet_2)],weight*PUWeight);
	  
    }

    break;
	
    }
      
    }
    }
  */
   
  return 1;
  
}
