#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <iostream>
// Methods
////////////////////////////////////////////////////////////////////////////////
void TopAnalysis::define_histograms()
{
  // HISTOGRAMS
    
  // Global variables histograms
  hist_mLL          = new TH1F("hist_mLL",       "Mass of Dilepton System; m_{ll} [GeV];Events / bin", 30, 0, 110);
  hist_ptLL         = new TH1F("hist_ptLL",      "Transverse Momentum of Dilepton System; p_{T}^{ll} [GeV];Events / bin", 30, 0, 200);
  hist_dPhi_LL      = new TH1F("hist_dPhi_LL",   "dPhi_LL Dilepton System; #Delta#phi (ll);Events / bin", 20, 0, 3.2);
  hist_dPhiLLmet    = new TH1F("hist_dPhiLLmet", "dPhiLLmet Dilepton System; #Delta#phi (ll,E_{T}^{miss});Events / bin", 30, 1, 5.2);
  hist_etmiss       = new TH1F("hist_etmiss",    "Missing Transverse Momentum;E_{T}^{miss} [GeV];Events / bin", 20, 0,200);
  hist_mt           = new TH1F("hist_mt",        "Dilepton System Transverse Mass; m_{T} [GeV];Events / bin", 20, 30, 250);
  hist_total_mass   = new TH1F("hist_tmass",     "Total Mass of the four objects; m_{ll#nu#nu} [GeV];Events / bin", 20, 30, 250);    
  hist_jet_n        = new TH1F("hist_jet_n",     "Number of jets; N_{Jets};Events / bin", 7, 0, 7);
    
  // Leading Lepton histograms
  hist_leadleptpt   = new TH1F("hist_leadleptpt",  "Leading Lepton Transverse Momentum;p_{T}^{leadlep} [GeV];Events / bin", 15, 0, 150);
  hist_leadlepteta  = new TH1F("hist_leadlepteta", "Leading Lepton Pseudorapidity; #eta^{leadlep}; Events / bin", 10, -3, 3);
  hist_leadleptE    = new TH1F("hist_leadleptE",   "Leading Lepton Energy; E^{leadlep} [GeV]; Events / bin", 25, 10, 260);
  hist_leadleptphi  = new TH1F("hist_leadleptphi", "Leading Lepton Azimuthal Angle ; #phi^{leadlep}; Events / bin", 16, -3.2, 3.2);
  hist_leadleptch   = new TH1F("hist_leadleptch",  "Leading Lepton Charge; Q^{leadlep}; Events / bin", 7, -1.75, 1.75);
  hist_leadleptID   = new TH1F("hist_leadleptID",  "Leading Lepton Absolute PDG ID; |PDG ID|^{leadlep}; Events / bin", 15, 5.5, 20.5);
  hist_leadlept_ptc  = new TH1F("hist_leadlept_ptc", "Leading Lepton Relative Transverse Momentum Isolation; ptconerel30^{leadlep}; Events / bin", 40, -0.1, 0.4);
  hist_leadleptetc  = new TH1F("hist_leadleptetc", "Leading Lepton Relative Transverse Energy Isolation; etconerel20^{leadlep}; Events / bin", 40, -0.1, 0.4);

    
  // Subleading Lepton histograms
  hist_subleadleptpt  = new TH1F("hist_subleadleptpt", "Subleading Lepton Transverse Momentum;p_{T}^{traillep} [GeV];Events / bin", 10, 0, 100);
  hist_subleadlepteta = new TH1F("hist_subleadlepteta","Subleading Lepton Pseudorapidity; #eta^{traillep}; Events / bin", 10, -3, 3);
  hist_subleadleptE   = new TH1F("hist_subleadleptE",  "Subleading Lepton Energy; E^{traillep} [GeV]; Events / bin", 15, 0, 250);
  hist_subleadleptphi = new TH1F("hist_subleadleptphi","Subleading Lepton Azimuthal Angle ; #phi^{traillep}; Events / bin", 10, -3.2, 3.2);
  hist_subleadleptch  = new TH1F("hist_subleadleptch", "Subleading Lepton Charge; Q^{traillep}; Events / bin", 7, -1.75, 1.75);
  hist_subleadleptID  = new TH1F("hist_subleadleptID", "Subleading Lepton Absolute PDG ID; |PDG ID|^{traillep}; Events / bin",  15, 5.5, 20.5);
  hist_subleadlept_ptc = new TH1F("hist_subleadlept_ptc","Subleading Lepton Relative Transverse Momentum Isolation; ptconerel30^{traillep} [GeV]; Events / bin", 40, -0.1, 0.4);
  hist_subleadleptetc = new TH1F("hist_subleadleptetc","Subleading Lepton Relative Transverse Energy Isolation; etconerel20^{traillep} [GeV]; Events / bin", 40, -0.1, 0.4);

}

////////////////////////////////////////////////////////////////////////////////
void TopAnalysis::FillOutputList()
{
  // histograms

  // Add Global variables histograms
  GetOutputList()->Add(hist_etmiss);
  GetOutputList()->Add(hist_mLL);  
  GetOutputList()->Add(hist_ptLL);  
  GetOutputList()->Add(hist_dPhi_LL);
  GetOutputList()->Add(hist_dPhiLLmet);
  GetOutputList()->Add(hist_mt);
  GetOutputList()->Add(hist_total_mass);
  GetOutputList()->Add(hist_jet_n);  

  // Add Leading Lepton histograms
  GetOutputList()->Add(hist_leadleptpt);
  GetOutputList()->Add(hist_leadlepteta);
  GetOutputList()->Add(hist_leadleptE);
  GetOutputList()->Add(hist_leadleptphi);
  GetOutputList()->Add(hist_leadleptch);
  GetOutputList()->Add(hist_leadleptID);
  GetOutputList()->Add(hist_leadlept_ptc);
  GetOutputList()->Add(hist_leadleptetc);

  // Add Subleading Lepton histograms
  GetOutputList()->Add(hist_subleadleptpt);
  GetOutputList()->Add(hist_subleadlepteta);
  GetOutputList()->Add(hist_subleadleptE);
  GetOutputList()->Add(hist_subleadleptphi);
  GetOutputList()->Add(hist_subleadleptch);
  GetOutputList()->Add(hist_subleadleptID);
  GetOutputList()->Add(hist_subleadlept_ptc);
  GetOutputList()->Add(hist_subleadleptetc);
    
}

////////////////////////////////////////////////////////////////////////////////
void TopAnalysis::WriteHistograms()
{
  // histograms

  // Write Global histograms
  hist_etmiss->Write();
  hist_mLL->Write();  
  hist_ptLL->Write();  
  hist_dPhi_LL->Write();
  hist_dPhiLLmet->Write();
  hist_mt->Write();
  hist_total_mass->Write();   
  hist_jet_n->Write();  
    
  //Write Leading Lepton histograms
  hist_leadleptpt->Write();
  hist_leadlepteta->Write();
  hist_leadleptE->Write();
  hist_leadleptphi->Write();
  hist_leadleptch->Write();
  hist_leadleptID->Write();
  hist_leadlept_ptc->Write();
  hist_leadleptetc->Write();

  //Write Subleading Lepton histograms
  hist_subleadleptpt->Write();
  hist_subleadlepteta->Write();
  hist_subleadleptE->Write();
  hist_subleadleptphi->Write();
  hist_subleadleptch->Write();
  hist_subleadleptID->Write();
  hist_subleadlept_ptc->Write();
  hist_subleadleptetc->Write();    
    
}

void TopAnalysis::FillHistogramsGlobal( double m, float w , TString s)
{
    
  //Fill Global histograms
  if(s.Contains("hist_mLL")) hist_mLL->Fill(m,w);  
  if(s.Contains("hist_ptLL")) hist_ptLL->Fill(m,w);  
  if(s.Contains("hist_dPhi_LL")) hist_dPhi_LL->Fill(m,w);
  if(s.Contains("hist_etmiss")) hist_etmiss->Fill(m,w);
  if(s.Contains("hist_mt")) hist_mt->Fill(m,w);
  if(s.Contains("hist_dPhiLLmet")) hist_dPhiLLmet->Fill(m,w);
  if(s.Contains("hist_total_mass")) hist_total_mass->Fill(m,w);
  if(s.Contains("hist_jet_n")) hist_jet_n->Fill(m,w);

}

void TopAnalysis::FillHistogramsLeadlept( double m, float w , TString s)
{
    
  //Leading lepton histograms
  if(s.Contains("hist_leadleptpt")) hist_leadleptpt->Fill(m,w);
  if(s.Contains("hist_leadlepteta")) hist_leadlepteta->Fill(m,w);
  if(s.Contains("hist_leadleptE")) hist_leadleptE->Fill(m,w);
  if(s.Contains("hist_leadleptphi")) hist_leadleptphi->Fill(m,w);
  if(s.Contains("hist_leadleptch")) hist_leadleptch->Fill(m,w);
  if(s.Contains("hist_leadleptID")) hist_leadleptID->Fill(m,w);
  if(s.Contains("hist_leadlept_ptc")) hist_leadlept_ptc->Fill(m,w);
  if(s.Contains("hist_leadleptetc")) hist_leadleptetc->Fill(m,w);
    
}


void TopAnalysis::FillHistogramsSubleadlept( double m, float w , TString s)
{
  //Fill Subleading lepton histograms
  if (s.Contains("hist_subleadleptpt")) hist_subleadleptpt->Fill(m,w);
  if (s.Contains("hist_subleadlepteta")) hist_subleadlepteta->Fill(m,w);
  if (s.Contains("hist_subleadleptE")) hist_subleadleptE->Fill(m,w);
  if (s.Contains("hist_subleadleptphi")) hist_subleadleptphi->Fill(m,w);
  if (s.Contains("hist_subleadleptch")) hist_subleadleptch->Fill(m,w);
  if (s.Contains("hist_subleadleptID")) hist_subleadleptID->Fill(m,w);
  if (s.Contains("hist_subleadlept_ptc")) hist_subleadlept_ptc->Fill(m,w);
  if (s.Contains("hist_subleadleptetc")) hist_subleadleptetc->Fill(m,w);

}


////////////////////////////////////////////////////////////////////////////////
