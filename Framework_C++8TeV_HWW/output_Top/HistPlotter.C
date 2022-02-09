#include "include/AtlasStyle.C"
#include "include/AtlasLabels.C"
#include "include/AtlasUtils.C"


void groupSamples();
void makeGroupHist();

std::map<TString, float> odb;
std::map<TString, float> dy;
std::map<TString, float> wj;
std::map<TString, float> tq;
std::map<TString, float> ww;
std::map<TString, float> z;


TH1F *h_otherdb;
TH1F *h_drellyan;
TH1F *h_wjets;
TH1F *h_topqmc; //Single top-quark via Wtb vertex, ttbar?->supongo que leptónico
TH1F *h_wwmc;
TH1F *h_z;

//Global Histograms
//std::string histname = "hist_mLL";
//std::string histname = "hist_ptLL";
//std::string histname = "hist_dPhi_LL";
//std::string histname = "hist_dPhiLLmet";
//std::string histname = "hist_etmiss";
//std::string histname = "hist_mt";
//std::string histname = "hist_total_mass"; //NW
//std::string histname = "hist_jet_n";

//Leading lepton Histograms
//std::string histname = "hist_leadleptpt";
//std::string histname = "hist_leadleptE";
//std::string histname = "hist_leadleptch";
std::string histname = "hist_leadleptID";


//Sudleading lepton Histograms
//std::string histname = "hist_subleadleptpt"; //NW
//std::string histname = "hist_subleadlepteta"; //NW
//std::string histname = "hist_subleadleptE"; //NW
//std::string histname = "hist_subleadleptphi"; //NW
//std::string histname = "hist_subleadleptch"; //NW
//std::string histname = "hist_subleadleptID"; //NW
//std::string histname = "hist_subleadlept_ptc"; //NW
//std::string histname = "hist_subleadleptetc"; //NW


float lumi = 1000.;
//int rebin = 20; //20 - ZmInv //4 - wmt //1 - event selection
bool logy = true;
//bool logy = false;
//std::string signalmass = "2000";

void HistPlotter(){

  SetAtlasStyle();
  TH1::SetDefaultSumw2(true);
  gROOT->Reset();

  groupSamples();
  makeGroupHist();

  THStack *h_stack = new THStack("h_stack","");

  h_otherdb->SetFillColor(kRed);
  h_drellyan->SetFillColor(kGreen);
  h_wjets->SetFillColor(kMagenta);
  h_topqmc->SetFillColor(kCyan);
  h_wwmc->SetFillColor(kWhite);
  h_z->SetFillColor(kBlue);
    
  //h_ttbar->SetFillColor(kCyan);

  h_stack->Add(h_otherdb);
  h_stack->Add(h_drellyan);
  h_stack->Add(h_wjets);
  h_stack->Add(h_topqmc);
  h_stack->Add(h_wwmc);
  h_stack->Add(h_z);
  //h_stack->Add(h_ttbar);

  TH1F *h_err = (TH1F*) h_otherdb->Clone();
  h_err->Add(h_drellyan);
  h_err->Add(h_wjets);
  h_err->Add(h_topqmc);
  h_err->Add(h_wwmc);
  h_err->Add(h_z);
  //h_err->Add(h_ttbar);
  h_err->SetFillStyle(3004);
  h_err->SetFillColor(kBlack);
  h_err->SetLineColor(0);
  h_err->SetMarkerStyle(1);
    
  std::ostringstream Entries_topq;
  Entries_topq << setprecision (5) << h_topqmc->Integral(0,h_topqmc->GetNbinsX()+1);
  std::string sEntries_topq = Entries_topq.str();
  cout<<sEntries_topq<<"    "<<endl;
   
  std::ostringstream Entries_odb;
  Entries_odb << setprecision (5) << h_otherdb->Integral(0,h_otherdb->GetNbinsX()+1);
  std::string sEntries_odb = Entries_odb.str();
  cout<<sEntries_odb<<"    "<<endl;
   
  std::ostringstream Entries_drellyan;
  Entries_drellyan << setprecision (5) << h_drellyan->Integral(0,h_drellyan->GetNbinsX()+1);
  std::string sEntries_drellyan = Entries_drellyan.str();
  cout<<sEntries_drellyan<<"    "<<endl;
   
  std::ostringstream Entries_wjets;
  Entries_wjets << setprecision (5) << h_wjets->Integral(0,h_wjets->GetNbinsX()+1);
  std::string sEntries_wjets = Entries_wjets.str();
  cout<<sEntries_wjets<<"    "<<endl;
    
  std::ostringstream Entries_z;
  Entries_z << setprecision (5) << h_z->Integral(0,h_z->GetNbinsX()+1);
  std::string sEntries_z = Entries_z.str();
  cout<<sEntries_z<<"    "<<endl;

  TFile *f_data = TFile::Open("data.root");
  TH1F* h_data = (TH1F*) f_data->Get(histname.c_str());
 
  //  TH1F *h_data = (TH1F*) h_data_egamma->Clone();
  //h_data->Add(h_data_muon);
  h_data->SetMarkerSize(1);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerColor(kBlack);

  /// Chi2
  //double chi2_test = h_data->Chi2Test(h_err, "UW,P");
  //cout<<"chi2_test  "<<chi2_test<<endl;


  std::ostringstream Entries_data;
  Entries_data << setprecision (5) << h_data->Integral(0,h_data->GetNbinsX()+1);
  std::string sEntries_data = Entries_data.str();
  cout<<"Observed: "<<sEntries_data<<endl;
 
    
  std::ostringstream Entries_wwMC;
  Entries_wwMC << setprecision (5) << h_wwmc->Integral(0,h_wwmc->GetNbinsX()+1);
  std::string sEntries_wwMC = Entries_wwMC.str();
  cout<<"Signal: "<<sEntries_wwMC<<"    "<<endl;
 
  std::ostringstream Entries_bkg;
  Entries_bkg << setprecision (5) << h_topqmc->Integral(0,h_topqmc->GetNbinsX()+1) + h_otherdb->Integral(0,h_otherdb->GetNbinsX()+1) + h_wjets->Integral(0,h_wjets->GetNbinsX()+1) + h_drellyan->Integral(0,h_drellyan->GetNbinsX()+1) + h_z->Integral(0,h_z->GetNbinsX()+1);
  std::string sEntries_bkg = Entries_bkg.str();
  cout<<"Background:"<<sEntries_bkg<<endl;
    
    

  ///// Plotting
  TCanvas *can = new TCanvas("can","", 600, 600);
  can->cd();
  if(logy) can->SetLogy();


  h_stack->Draw("hist");
  h_err->Draw("e2same");
  h_data->Draw("epsame");
  //  h_signal->Draw("histsame");
   h_stack->SetMaximum(20*h_stack->GetMaximum());
  if (logy && histname.find("event_selection") != std::string::npos) h_stack->SetMinimum(10000);
  //else if(logy) h_stack->SetMinimum(0.5);

  h_stack->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
  h_stack->GetYaxis()->SetTitle("Events");
  h_stack->GetYaxis()->SetTitleOffset(1.60);
  if(histname.find("WtMass") != std::string::npos)h_stack->GetXaxis()->SetRangeUser(0, 500.);
  else if(histname.find("ZmInv") != std::string::npos)h_stack->GetXaxis()->SetRangeUser(0, 3775);
  else if(histname.find("ZmT") != std::string::npos)h_stack->GetXaxis()->SetRangeUser(0, 3575);
  else if(histname.find("lep_pt") != std::string::npos)h_stack->GetXaxis()->SetRangeUser(0, 1575);



  TLegend *leg = new TLegend(0.55, 0.55, 0.75, 0.9);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetHeader("     #sqrt{s} = 8 TeV, 1 fb^{-1}");
  leg->AddEntry(h_data,("Data: "+sEntries_data).c_str(),"p");
  leg->AddEntry(h_topqmc,("Top Quark MC: "+sEntries_topq).c_str(),"f");
  leg->AddEntry(h_wwmc,("WW MC: "+sEntries_wwMC).c_str(),"f");
  leg->AddEntry(h_otherdb,("Other diboson MC: "+sEntries_odb).c_str(),"f");
  leg->AddEntry(h_drellyan,("Drell-Yan: "+sEntries_drellyan).c_str(),"f");
  leg->AddEntry(h_wjets,("W+Jets: "+sEntries_wjets).c_str(),"f");
  leg->AddEntry(h_z,("Z: "+sEntries_z).c_str(),"f");
  leg->AddEntry(h_err,"MC stat. uncertainty","f");
  // leg->AddEntry(h_signal,("Z'(" + signalmass + " GeV)#rightarrow t#bar{t} x 10^{3}").c_str(),"l");
  leg->Draw();

  std::string outfile = histname.substr(0, std::string::npos);
  if(logy) can->SaveAs(("HWW-distribution-plots/"+outfile+"_log.pdf").c_str());
  else can->SaveAs(("HWW-distribution-plots/"+outfile+".pdf").c_str());
}

void groupSamples(){

  std::ifstream samples;
  samples.open("samples.txt");
  
  while(!samples.eof()){
    TString sample;
    float xsec;
    float red_eff;
    float sumw;
    float nevt;

    samples >> sample >> xsec >> red_eff >> sumw >> nevt;

    nevt=1;
    //sumw=1;
    if(sample.Contains("WZ")  || sample.Contains("ZZ")  ) odb[sample] = xsec*nevt/(red_eff*sumw);
      
    if(sample.Contains("WW") ) ww[sample] = xsec*nevt/(red_eff*sumw);      

    if(sample.Contains("DY")  ) dy[sample] = xsec*nevt/(red_eff*sumw);

    if(sample.Contains("Wenu")  || sample.Contains("Wmunu")  || sample.Contains("Wtaunu") ) wj[sample] = xsec*nevt/(red_eff*sumw); 

    if(sample.Contains("Zee")  || sample.Contains("Zmumu")  || sample.Contains("Ztautau")  ) z[sample] = xsec*nevt/(red_eff*sumw); 

    if(sample.Contains("stop") || sample.Contains("ttbar") ) tq[sample] = xsec*nevt/(red_eff*sumw);

  }


  return;
}

void makeGroupHist(){

  //for(unsigned int i = 0; i < db.size(); i++){
  int i = 0;
  for(std::map<TString, float>::iterator it = odb.begin(); it != odb.end(); ++it){
    TFile *f = TFile::Open(Form("%s.root",(it->first).Data()));
    TH1F *h = (TH1F*) f->Get(histname.c_str());    
    //h->Rebin(rebin);
    h->Scale(it->second*lumi);
    //h->Scale(it->second*lumi/h->Integral(0, h->GetNbinsX()+1));
    if(i==0) h_otherdb = (TH1F*) h->Clone();
    else h_otherdb->Add(h);
    i++;
  }

  i = 0;
  for(std::map<TString, float>::iterator it = dy.begin(); it != dy.end(); ++it){
    TFile *f = TFile::Open(Form("%s.root",(it->first).Data()));
    TH1F *h = (TH1F*) f->Get(histname.c_str());    
    //h->Rebin(rebin);
    h->Scale(it->second*lumi);
    //h->Scale(it->second*lumi/h->Integral(0, h->GetNbinsX()+1));
    if(i==0) h_drellyan = (TH1F*) h->Clone();
    else h_drellyan->Add(h);
    i++;
  }

  i = 0;
  for(std::map<TString, float>::iterator it = wj.begin(); it != wj.end(); ++it){
    TFile *f = TFile::Open(Form("%s.root",(it->first).Data()));
    TH1F *h = (TH1F*) f->Get(histname.c_str());    
    //h->Rebin(rebin);
    h->Scale(it->second*lumi);
    //h->Scale(it->second*lumi/h->Integral(0, h->GetNbinsX()+1));
    if(i==0) h_wjets = (TH1F*) h->Clone();
    else h_wjets->Add(h);
    i++;
  }

  i = 0;
  for(std::map<TString, float>::iterator it = tq.begin(); it != tq.end(); ++it){
    TFile *f = TFile::Open(Form("%s.root",(it->first).Data()));
    TH1F *h = (TH1F*) f->Get(histname.c_str());    
    //h->Rebin(rebin);
    h->Scale(it->second*lumi);
    //h->Scale(it->second*lumi/h->Integral(0, h->GetNbinsX()+1));
    if(i==0) h_topqmc = (TH1F*) h->Clone();
    else h_topqmc->Add(h);
    i++;
  }

  i = 0;
  for(std::map<TString, float>::iterator it = ww.begin(); it != ww.end(); ++it){
    TFile *f = TFile::Open(Form("%s.root",(it->first).Data()));
    TH1F *h = (TH1F*) f->Get(histname.c_str());    
    //h->Rebin(rebin);
    h->Scale(it->second*lumi);
    //h->Scale(it->second*lumi/h->Integral(0, h->GetNbinsX()+1));
    if(i==0) h_wwmc = (TH1F*) h->Clone();
    else h_wwmc->Add(h);
    i++;
  }
    
  i = 0;
  for(std::map<TString, float>::iterator it = z.begin(); it != z.end(); ++it){
    TFile *f = TFile::Open(Form("%s.root",(it->first).Data()));
    TH1F *h = (TH1F*) f->Get(histname.c_str());    
    //h->Rebin(rebin);
    h->Scale(it->second*lumi);
    //h->Scale(it->second*lumi/h->Integral(0, h->GetNbinsX()+1));
    if(i==0) h_z = (TH1F*) h->Clone();
    else h_z->Add(h);
    i++;
  }

}
