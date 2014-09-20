#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#endif

#include <iostream>
#include <fstream> 
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TLatex.h"


void MakeTimeResolutionPlot(string filename, string plotname, int energy, int run) {
  // Get the tree


  TFile *inputfile = new TFile(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("tree");

  // get the variables from the ntuple
  float t1gausroot = 0;
  float t2gausroot = 0;
  float t3gausroot = 0;
  float t4gausroot = 0;
  float t1RisingEdge = 0;
  float t2RisingEdge = 0;
  float t3RisingEdge = 0;
  float t4RisingEdge = 0;
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;
  unsigned int ch1QualityBit = 0;
  unsigned int ch2QualityBit = 0;
  unsigned int ch3QualityBit = 0;
  unsigned int ch4QualityBit = 0;

  tree->SetBranchAddress("t1gausroot",&t1gausroot);
  tree->SetBranchAddress("t2gausroot",&t2gausroot);
  tree->SetBranchAddress("t3gausroot",&t3gausroot);
  tree->SetBranchAddress("t4gausroot",&t4gausroot);
  tree->SetBranchAddress("ch1THM",&t1RisingEdge);
  tree->SetBranchAddress("ch2THM",&t2RisingEdge);
  tree->SetBranchAddress("ch3THM",&t3RisingEdge);
  tree->SetBranchAddress("ch4THM",&t4RisingEdge);
  tree->SetBranchAddress("ch1Amp",&ch1Amp);
  tree->SetBranchAddress("ch2Amp",&ch2Amp);
  tree->SetBranchAddress("ch3Amp",&ch3Amp);
  tree->SetBranchAddress("ch4Amp",&ch4Amp);
  tree->SetBranchAddress("ch1QualityBit",&ch1QualityBit);
  tree->SetBranchAddress("ch2QualityBit",&ch2QualityBit);
  tree->SetBranchAddress("ch3QualityBit",&ch3QualityBit);
  tree->SetBranchAddress("ch4QualityBit",&ch4QualityBit);

  //create histograms

  TH1F *dt = 0;
  if (energy == 32 || energy == 16) {
    dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 100, 3,5);
  }
  if (energy == 8
      || ( energy ==16 && run == 138)
    ) {
    dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 100, -1,5);
  }

  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;
  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
    tree->GetEntry(iEntry);
      

    // select electrons and good pulses
    if (energy == 32 ) {
      if (  ch1QualityBit == 0 && 
            ch4Amp > 0.01 &&
            ch3Amp > 0.02 && ch3Amp < 0.49
        ) {                       
        dt->Fill(-1*(t1gausroot-t3RisingEdge));
        //don't use ch2 here because it seems that there is something wrong with that channel in this run.
      }
    }
    if (energy == 16 && run == 140 ) {
      if (  ch1QualityBit == 0 && 
            ch4Amp > 0.01 &&
            ch2Amp > 0.02 && ch2Amp < 0.49 &&
            ch3Amp > 0.02 && ch3Amp < 0.49
        ) {                       
        dt->Fill(-1*(t1gausroot-t2RisingEdge+t1gausroot-t3RisingEdge)/2);
      }
    }
    if (energy == 16 && run == 138) {
      if( ch1QualityBit==0 &&       
          ch4Amp > 0.01 &&
          ch3Amp > 0.02 && ch3Amp < 0.49 
        ) {
        dt->Fill(-1*(t1gausroot-t3RisingEdge));
      }
    }
    if (energy == 8 ) {
      if( ch1QualityBit==0 &&       
          ch4Amp > 0.01 &&
          ch3Amp > 0.02 && ch3Amp < 0.49 
        ) {
        dt->Fill(-1*(t1gausroot-t3RisingEdge));
      }
    }
    

    }
  
  // Touch up and save
  TCanvas * c = new TCanvas("c","c",600,600);
  TF1 *f1 = new TF1("f1","gaus",dt->GetMean()-2*fabs(dt->GetRMS()),dt->GetMean()+2*fabs(dt->GetRMS()));
  dt->SetAxisRange(dt->GetMean()-3.0*fabs(dt->GetRMS()),dt->GetMean()+3.0*fabs(dt->GetRMS()),"X");
  dt->SetTitle("");
  dt->GetXaxis()->SetTitle("#Delta t [ns]");
  dt->GetYaxis()->SetTitle("Number of Events");
  dt->GetYaxis()->SetTitleOffset(1.25);
  dt->Draw();
  dt->SetStats(0);
  dt->Fit("f1","R");
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.040);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.13, 0.85, Form("%d GeV Electron Beam",energy));
  tex->DrawLatex(0.13, 0.80, Form("#sigma=%.0f +/- %.0f ps",1000*fitter->GetParameter(2), 1000*fitter->GetParError(2)));




  c->SaveAs( Form("%s.gif", plotname.c_str()) );
  c->SaveAs( Form("%s.pdf", plotname.c_str()) );



}


Double_t TOFResolutionFunction(Double_t *x,Double_t *par)
{
  Double_t arg = 0;
//if (par[2] != 0) arg = (x[0] - par[1])/par[2];
  arg = x[0];
  Double_t fitval = par[0]*(1.0 / sqrt(arg)) + par[1];
  return fitval;
}


void MakeTimeResolutionVsEnergyPlot() {

  float x[3] = {8,16,32};
  float xerr[3] = { 0.02*8, 0.02*16, 0.02*32 };
  
  //fit to full range
   float y[3] = {83,71,54 };
   float yerr[3] = { 3,4,5 };

  TGraphErrors *graph = new TGraphErrors(4,x,y,xerr,yerr);
  graph->SetLineWidth(2);
  TCanvas * c = new TCanvas("c","c",600,600);
  graph->Draw("AP");

  graph->SetTitle("");
  graph->GetXaxis()->SetTitle("Electron Beam Energy [GeV]");
  graph->GetYaxis()->SetTitle("Time Resolution [ps]");
  graph->GetYaxis()->SetTitleOffset(1.35);
  graph->GetYaxis()->SetRangeUser(0,150);

  TF1 *func = new TF1("TOFResolutionFunction",TOFResolutionFunction,0,100,2);

   // call TH1::Fit with the name of the TF1 object 
  func->SetParLimits(1,0,1000);
  graph->Fit("TOFResolutionFunction");


  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.040);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.30, 0.80, "TOF Resolution = #frac{A}{#sqrt{E}} #oplus C");
  tex->DrawLatex(0.30, 0.70, Form("A = %.0f +/- %.0f ps",func->GetParameter(0), func->GetParError(0)));
  tex->DrawLatex(0.30, 0.65, Form("C = %.0f +/- %.0f ps",func->GetParameter(1), func->GetParError(1)));

  TLatex *tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextSize(0.045);
  tex2->SetTextFont(42);
  tex2->SetTextColor(kBlack);
//   tex2->DrawLatex(0.15, 0.92, "Caltech Internal");


  c->SaveAs( "TimeResolutionVsEnergy_ShashlikSideReadout.gif" );
  c->SaveAs( "TimeResolutionVsEnergy_ShashlikSideReadout.pdf" );

}

void ShashlikSideReadoutAnalysis() {

  //This Run (140) has only side readout from one side. Run 140 has side readout from both sides
//   MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_may_run_138.ana.root","TOF_ShashlikSideReadout_Electron_16GeV",16,138);

  //*************************************
  //Best Results Here for Paper
  //*************************************
  // MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_may_run_137.ana.root","TOF_ShashlikSideReadout_Electron_8GeV",8,137);
  // MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_may_run_140.ana.root","TOF_ShashlikSideReadout_Electron_16GeV",16,140);  
  // MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_may_run_141.ana.root","TOF_ShashlikSideReadout_Electron_32GeV",32,141);
  MakeTimeResolutionVsEnergyPlot();


}
