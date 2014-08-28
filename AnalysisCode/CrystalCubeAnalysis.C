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


void MakeTimeResolutionPlot(string filename, string plotname, int energy, bool plotEnergyPeak = false) {
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
  float ch1Int = 0;
  float ch2Int = 0;
  float ch3Int = 0;
  float ch4Int = 0;
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
  tree->SetBranchAddress("ch1Int",&ch1Int);
  tree->SetBranchAddress("ch2Int",&ch2Int);
  tree->SetBranchAddress("ch3Int",&ch3Int);
  tree->SetBranchAddress("ch4Int",&ch4Int);
  tree->SetBranchAddress("ch1QualityBit",&ch1QualityBit);
  tree->SetBranchAddress("ch2QualityBit",&ch2QualityBit);
  tree->SetBranchAddress("ch3QualityBit",&ch3QualityBit);
  tree->SetBranchAddress("ch4QualityBit",&ch4QualityBit);

  //create histograms

  TH1F *dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 100, -1,1);
  TH1F *histIntegral;
  histIntegral = new TH1F("histIntegral","; Pulse Integral [V];Number of Events",25,0,200);

  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;
  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
      tree->GetEntry(iEntry);
      

      // select electrons and good pulses
      if (energy == 32) {
        if( ch1QualityBit==0 &&       
            ch3Amp > 0.2 && ch3Amp < 0.49 ) {
	  dt->Fill(-1*(t1gausroot-t3RisingEdge));
          histIntegral->Fill(ch3Int);
        }
      }
      if (energy == 16) {
        if( ch1QualityBit==0 &&       
            ch2Amp > 0.01 &&
            ch3Amp > 0.02 && ch3Amp < 0.49 ) {
	  dt->Fill(-1*(t1gausroot-t3RisingEdge));
          histIntegral->Fill(ch3Int);
        }
      }
      if (energy == 8 || energy == 4) {
        if( ch1QualityBit==0 &&       
            ch4Amp > 0.01 &&
            ch3Amp > 0.02 && ch3Amp < 0.49 ) {
	  dt->Fill(-1*(t1gausroot-t3RisingEdge));
          histIntegral->Fill(ch3Int);
        }
      }

  }
  
  // Touch up and save
  TCanvas * c = new TCanvas("c","c",600,600);

  dt->SetAxisRange(dt->GetMean()-1.0*fabs(dt->GetMean()),dt->GetMean()+1.0*fabs(dt->GetMean()),"X");
  dt->SetTitle("");
  dt->GetXaxis()->SetTitle("#Delta t [ns]");
  dt->GetYaxis()->SetTitle("Number of Events");
  dt->GetYaxis()->SetTitleOffset(1.25);
  dt->Draw();
  dt->SetStats(0);
  dt->Fit("gaus");
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.040);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.13, 0.80, Form("#sigma=%.1f +/- %.1f ps",1000*fitter->GetParameter(2), 1000*fitter->GetParError(2)));

  c->SaveAs( Form("%s.gif", plotname.c_str()) );
  c->SaveAs( Form("%s.pdf", plotname.c_str()) );


  if (plotEnergyPeak) {
    //Energy plot
    c = new TCanvas("c","c",600,600);
    
    histIntegral->SetAxisRange(histIntegral->GetMean()-2.0*fabs(histIntegral->GetMean()),histIntegral->GetMean()+2.0*fabs(histIntegral->GetMean()),"X");
    histIntegral->SetTitle("");
    histIntegral->GetXaxis()->SetTitle("Pulse Integral [V]");
    histIntegral->GetYaxis()->SetTitle("Number of Events");
    histIntegral->GetYaxis()->SetTitleOffset(1.25);
    histIntegral->Draw();
    histIntegral->SetStats(0);
    histIntegral->Fit("gaus");
    TVirtualFitter * fitter = TVirtualFitter::GetFitter();
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->SetTextColor(kBlack);
    tex->DrawLatex(0.13, 0.85, Form("%d GeV Electron Beam",energy));
    tex->DrawLatex(0.53, 0.80, Form("Resolution = %.0f %s",100*fitter->GetParameter(2)/fitter->GetParameter(1),"%"));
    
    c->SaveAs( Form("%s_energy.gif", plotname.c_str()) );
    c->SaveAs( Form("%s_energy.pdf", plotname.c_str()) );
  }
  


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

  //beam energy
  float x[4] = {4,8,16,32};
  float xerr[4] = { 0.02*4, 0.02*8, 0.02*16, 0.02*32 };

//   //amplitude
//   float x[4] = {170,300,600,1200};
//   float xerr[4] = { 0.02*150, 0.02*300, 0.02*600, 0.02*1200 };


  float y[4] = {70.4,57.4, 40.6, 33.5 };

  float yerr[4] = { 2.9,2.0,1.5,2.1 };

  TGraphErrors *graph = new TGraphErrors(4,x,y,xerr,yerr);
  graph->SetLineWidth(2);
  TCanvas * c = new TCanvas("c","c",600,600);
  graph->Draw("AP");

  graph->SetTitle("");
  graph->GetXaxis()->SetTitle("Electron Beam Energy [GeV]");
  //graph->GetXaxis()->SetTitle("Amplitude [mV]");
  graph->GetYaxis()->SetTitle("Time of Flight Resolution [ps]");
  graph->GetYaxis()->SetTitleOffset(1.1);


  TF1 *func = new TF1("TOFResolutionFunction",TOFResolutionFunction,0,100,2);
   // set the parameters to the mean and RMS of the histogram
   //func->SetParameters(0,1,1);
   // give the parameters meaningful names
   //func->SetParNames ("Constant","Mean_value","Sigma");

   // call TH1::Fit with the name of the TF1 object 
   graph->Fit("TOFResolutionFunction");

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.040);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.30, 0.80, "TOF Resolution = #frac{A}{#sqrt{E}} #oplus C");
  tex->DrawLatex(0.30, 0.70, Form("A = %.0f +/- %.0f ps",func->GetParameter(0), func->GetParError(0)));
  tex->DrawLatex(0.30, 0.65, Form("C = %.0f +/- %.0f ps",func->GetParameter(1), func->GetParError(1)));

  c->SaveAs( "TimeResolutionVsEnergy_CrystalCube.gif" );
  c->SaveAs( "TimeResolutionVsEnergy_CrystalCube.pdf" );

}

void CrystalCubeAnalysis() {

  MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_aug_run_070.ana.root","TOF_Electron_LYSOCube_4GeV",4,false);
   MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_aug_run_069.ana.root","TOF_Electron_LYSOCube_8GeV", 8,false);
   MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_aug_run_064.ana.root","TOF_Electron_LYSOCube_16GeV",16,false);
  MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_aug_run_065-068.ana.root","TOF_Electron_LYSOCube_32GeV",32,true);

  MakeTimeResolutionVsEnergyPlot();
}
