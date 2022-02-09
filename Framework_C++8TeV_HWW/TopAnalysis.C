#define TopAnalysis_cxx
// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.

#include "TopAnalysis.h"
#include "Tophistograms.h"
#include <iostream>
#include <cstring>
#include <string>

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLorentzVector.h>

string name;

void TopAnalysis::Begin(TTree * )
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
}

void TopAnalysis::SlaveBegin(TTree * )
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  printf("Starting analysis with process option: %s \n", option.Data());

  name=option;

  define_histograms();

  FillOutputList();
}

Bool_t TopAnalysis::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fChain->GetTree()->GetEntry(entry);
  //  int cut1_mc = 0;

  if(fChain->GetTree()->GetEntries()>0){
    //Do analysis

    //SF
    Float_t scaleFactor = scaleFactor_ELE*scaleFactor_MUON*scaleFactor_TRIGGER;
    //EventW
    Float_t eventWeight = mcWeight*scaleFactor_PILEUP*scaleFactor_ZVERTEX;
    //weight = SF * EventW
    Double_t weight = scaleFactor*eventWeight;

    // Make difference between data and MC
    if (weight == 0. && channelNumber != 110090 && channelNumber != 110091) weight = 1.;

    // Missing Et of the event in GeV
    Float_t missingEt = met_et/1000.;


    // Preselection cut for electron/muon trigger, Good Run List, and good vertex
    if(trigE || trigM){
      if(passGRL){
        if(hasGoodVertex){
          //Find the good leptons
          int goodlep_index[5];
          int goodlep_n = 0;
          int lep_index = 0;

          for(unsigned int i=0; i<lep_n; i++){
            //Selection of goodleptons: pt>25GeV, isolated leptons:  
            if(lep_pt[i]>25000. && (lep_ptcone30[i]/lep_pt[i])<0.15 && (lep_etcone20[i]/lep_pt[i])<0.15){     
              // electron selection in fiducial region excluding candidates in the transition region between the barrel and endcap electromagnetic calorimeters
              if(lep_type[i]==11 && TMath::Abs(lep_eta[i])<2.47 && (TMath::Abs(lep_eta[i])<1.37 || TMath::Abs(lep_eta[i])>1.52)){
                //if( TMath::Abs(lep_trackd0pvunbiased->at(i))/lep_tracksigd0pvunbiased->at(i) < 5 && TMath::Abs(lep_z0->at(i)*TMath::Sin(leptemp.Theta())) < 0.5)
                goodlep_n++;
                goodlep_index[lep_index]=i;
                lep_index++;
              }  
              //muon selection:  
              if (lep_type[i]==13 && TMath::Abs(lep_eta[i])<2.5){
                goodlep_n++;
                goodlep_index[lep_index]=i;
                lep_index++;
              }
            }
          }

          //First cut: Exactly two good leptons
          if(goodlep_n==2){ 
             
            int goodlep1_index = goodlep_index[0];
            int goodlep2_index = goodlep_index[1];
            
            //Second cut: goodleptons of opposite charges
            if(lep_charge[goodlep1_index]*lep_charge[goodlep2_index]<0){
            
              //TLorentzVector definitions  
              TLorentzVector Lepton_1 = TLorentzVector();
              TLorentzVector Lepton_2 = TLorentzVector();
              TLorentzVector Lepton_12 = TLorentzVector();
              TLorentzVector MeT = TLorentzVector();              
                
              Lepton_1.SetPtEtaPhiE(lep_pt[goodlep1_index], lep_eta[goodlep1_index], lep_phi[goodlep1_index],lep_E[goodlep1_index]);
                
              Lepton_2.SetPtEtaPhiE(lep_pt[goodlep2_index], lep_eta[goodlep2_index], lep_phi[goodlep2_index],lep_E[goodlep2_index]);
                
              MeT.SetPtEtaPhiE(met_et, 0, met_phi , met_et);                
               
              Lepton_12 = Lepton_1 + Lepton_2; 
              
              float mLL = Lepton_12.Mag()/1000.; //dilepton system mass in GeV
              float ptLL = Lepton_12.Pt()/1000.; //dilepton system pt in GeV
              
              //Azimutal angle separation dphi between two goodleptons
              float dPhi_LL = TMath::Abs(lep_phi[goodlep1_index]-lep_phi[goodlep2_index]); 
              dPhi_LL = dPhi_LL<TMath::Pi() ? dPhi_LL : 2*TMath::Pi()-dPhi_LL;
                
              //Azimutal angle separation dphi between the goodlepton pair and MET 
              float dPhiLLmet = TMath::Abs(Lepton_12.Phi()-MeT.Phi());
              dPhiLLmet = dPhiLLmet<TMath::Pi() ? dPhiLLmet : 2*TMath::Pi()-dPhiLLmet;
                
              //Transverse mass in GeV
              float mt = TMath::Sqrt(2*Lepton_12.Pt()*MeT.Et()*(1-cos(Lepton_12.DeltaPhi(MeT))))/1000.;
                
              /*
                Third cut: pT lepton pair > 30 GeV;
                Four cut: Angular separation between lepton pair and Missing ET > π/2;
                Fifth cut: Reconstructed transverse mass W > 30 GeV;
                Sixth cut: Mass lepton pair < 55 GeV;
                Seventh cut: Angular separation between leptons < 1.8.
              */
                
              if(ptLL>30 && (dPhiLLmet>TMath::Pi()/2) && mt>30 && mLL<55 && dPhi_LL<1.8){
                /*
                Eighth cut: If leptons have same flavor:
                 -mass lepton pair > 12 GeV
                 -|reconstructed mass lepton pair - PDG mass Z| > 15 GeV
                 -Missing ET > 40 GeV  
                */
                //If goodleptons of same flavour:
                bool flavour_cut = false;
                if(lep_type[goodlep1_index]==lep_type[goodlep2_index]){
                  
                  float PDG_massZ=91.1876;
                  if(mLL>12){
                    if(TMath::Abs(mt-PDG_massZ)>15){
                      if(met_et>40000){ //met_et is in MeV
                          flavour_cut=true;
                      }
                    }
                  }
                }
                /*
                Ninth cut: Else:
                 -mass lepton pair > 10 GeV
                 -Missing ET > 20 GeV
                */
                //If goodleptons of different flavour:
                if(lep_type[goodlep1_index]!=lep_type[goodlep2_index]){                    
                    
                  if(mLL>10){
                    if(met_et>20000){ //met_et is in MeV
                      flavour_cut=true;    
                    }
                  }
                }
                
                if(flavour_cut==true){
                  //Variable used to count at least one jet with pt>25GeV  
                  int high_jet_pt_n = 0;
                  
                  for(unsigned int j=0; j<jet_n; j++){
                    if(jet_pt[j]>25000){
                        high_jet_pt_n++;
                    }
                  }
                  hist_jet_n->Fill(jet_n, weight);
                  //Tenth cut: No jets with pT > 25 GeV
                  if(high_jet_pt_n==0){    
  
                    float total_mass = (Lepton_1 + Lepton_2 + MeT + MeT).M()/1000.; //Sum of the four final objects in GeV
                    hist_total_mass->Fill(total_mass, weight); 
                    //Start to fill histograms : definitions of x-axis variables
                  
                    double names_of_global_variable[]={mLL, ptLL, dPhi_LL, dPhiLLmet, met_et/1000., mt};  
                
                    double names_of_leadlep_variable[]={Lepton_1.Pt()/1000., Lepton_1.Eta(), Lepton_1.Phi(), Lepton_1.E()/1000.,
                                                    (double)lep_charge[goodlep1_index], (double)lep_type[goodlep1_index], 
                                                    lep_ptcone30[goodlep1_index]/lep_pt[goodlep1_index], 
                                                    lep_etcone20[goodlep1_index]/lep_pt[goodlep1_index]};

                    double names_of_subleadlep_variable[]={Lepton_2.Pt()/1000., Lepton_2.Eta(), Lepton_2.Phi(), Lepton_2.E()/1000., 
                                                       (double)lep_charge[goodlep2_index], (double)lep_type[goodlep2_index],
                                                       lep_ptcone30[goodlep2_index]/lep_pt[goodlep2_index], 
                                                       lep_etcone20[goodlep2_index]/lep_pt[goodlep2_index]};

                    //Start to fill histograms : definitions of histogram names
  
                    TString histonames_of_global_variable[]={"hist_mLL", "hist_ptLL", "hist_dPhi_LL", "hist_dPhiLLmet", "hist_etmiss", 
                                                         "hist_mt"};
                
                    TString histonames_of_leadlep_variable[]={"hist_leadleptpt", "hist_leadlepteta", "hist_leadleptphi", "hist_leadleptE", 
                                                          "hist_leadleptch", "hist_leadleptID", "hist_leadlept_ptc", "hist_leadleptetc"};
          
                    TString histonames_of_subleadlep_variable[]={"hist_subleadleptpt", "hist_subleadlepteta", "hist_subleadleptphi", 
                                                             "hist_subleadleptE", "hist_subleadleptch", "hist_subleadleptID",
                                                             "hist_subleadlept_ptc", "hist_subleadleptetc"};

                    int length_global = sizeof(names_of_global_variable)/sizeof(names_of_global_variable[0]);
                    int length_leadlep = sizeof(names_of_leadlep_variable)/sizeof(names_of_leadlep_variable[0]);
                    int length_subleadlep = sizeof(names_of_subleadlep_variable)/sizeof(names_of_subleadlep_variable[0]);
                  
                    //Fill histograms
  
                    for(int i=0; i<length_global; i++){
                      FillHistogramsGlobal( names_of_global_variable[i], weight, histonames_of_global_variable[i]);
                    }
  
                    for(int i=0; i<length_leadlep; i++){
                      FillHistogramsLeadlept( names_of_leadlep_variable[i], weight, histonames_of_leadlep_variable[i]);
                    }
  
                    for(int i=0; i<length_subleadlep; i++){
                      FillHistogramsSubleadlept( names_of_subleadlep_variable[i], weight, histonames_of_subleadlep_variable[i]);
                    }

                  }
                
                }  
                  
              }
                
            }
              
          }
            
        }
          
      }
        
    }
  
  }
 
  return kTRUE;
}

void TopAnalysis::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void TopAnalysis::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  name="output_Top/"+name+".root";

  const char* filename = name.c_str();

  TFile physicsoutput_Top(filename,"recreate");
  WriteHistograms();
  physicsoutput_Top.Close();


}
