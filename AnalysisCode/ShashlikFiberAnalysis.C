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


void MakeTimeResolutionPlot(string filename, string plotname, int energy, int isDSBFiber, int run, bool plotEnergyPeak = false) {
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
  float ch1Risetime = 0;
  float ch2Risetime = 0;
  float ch3Risetime = 0;
  float ch4Risetime = 0;
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
  tree->SetBranchAddress("ch1Risetime",&ch1Risetime);
  tree->SetBranchAddress("ch2Risetime",&ch2Risetime);
  tree->SetBranchAddress("ch3Risetime",&ch3Risetime);
  tree->SetBranchAddress("ch4Risetime",&ch4Risetime);
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
  double histIntegralFitMin = 0;
  double histIntegralFitMax = 0;
  TH1F *dt;
  TH1F *histIntegral;
  if (energy == 4 | energy == 8 ) {
    histIntegral = new TH1F("histIntegral","; Pulse Integral [V*ns];Number of Events",150,0,300);
  } else if (energy == 16) {
    histIntegral = new TH1F("histIntegral","; Pulse Integral [V*ns];Number of Events",75,0,300);
  } else {
    histIntegral = new TH1F("histIntegral","; Pulse Integral [V*ns];Number of Events",50,0,200);
  }

  if (run <86) {
    dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 100, -5,5);
  } else if (run == 90) {
    dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 200, 0,5);
  } else if (run == 92) {
    dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 25, 0,5);
  } else {
    dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 50, 0,5);
  }

  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;
  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
      tree->GetEntry(iEntry);
      

      // select electrons and good pulses
      if (isDSBFiber) {
        //early runs
        if (run <= 78 && (energy == 32 || energy == 16 || energy == 8 || energy == 4 )) {
          if( ch4QualityBit==0 &&       
              ch2Amp > 0.01 &&
              ch3Amp > 0.02 && ch3Amp < 0.49 &&
              ch1Amp > 0.02 && ch1Amp < 0.49 
	    ) {
            dt->Fill(-1*(t4gausroot-t3RisingEdge+t4gausroot-t1RisingEdge)/2);
            histIntegral->Fill(ch3Int+ch1Int);
          }
        }
        //later runs : these are better runs
        else {
          //there are events with noise just before the main scintilation pulse
          //we remove them using a cut on the risetime
          if (energy == 4) {
            if( ch1QualityBit==0 &&       
                ch2Amp > 0.005 &&
                ch3Amp > 0.02 && ch3Amp < 0.49 &&
                ch4Amp > 0.02 && ch4Amp < 0.49 &&
                ch3Risetime < 5 && ch4Risetime < 5
              ) {
              dt->Fill(-1*(t1gausroot-0.5*(t3RisingEdge+t4RisingEdge)));
              histIntegral->Fill(ch3Int+ch4Int);
              histIntegralFitMin = 25;
              histIntegralFitMax = 50;
            }
          }
          if (energy == 8) {
            if( ch1QualityBit==0 &&       
                ch2Amp > 0.02 &&
                ch3Amp > 0.02 && ch3Amp < 0.49 &&
                ch4Amp > 0.02 && ch4Amp < 0.49 &&
                ch3Risetime < 4 && ch4Risetime < 4
              ) {
              dt->Fill(-1*(t1gausroot-0.5*(t3RisingEdge+t4RisingEdge)));
              histIntegral->Fill(ch3Int+ch4Int);
              histIntegralFitMin = 65;
              histIntegralFitMax = 100;
            }
          }
          if (energy == 16) {
            if( ch1QualityBit==0 &&       
                ch2Amp > 0.02 &&
                ch3Amp > 0.02 && ch3Amp < 0.49 &&
                ch4Amp > 0.02 && ch4Amp < 0.49 &&
                ch3Risetime < 4 && ch4Risetime < 4
              ) {
              dt->Fill(-1*(t1gausroot-0.5*(t3RisingEdge+t4RisingEdge)));
              histIntegral->Fill(ch3Int+ch4Int);
              histIntegralFitMin = 130;
              histIntegralFitMax = 200;
            }        
          }

          //cherenkov didn't work. we use the energy peak the do electron ID
          if (energy == 32) {
            if( ch1QualityBit==0 &&       
                ch3Amp > 0.02 && ch3Amp < 0.49 &&
                ch4Amp > 0.02 && ch4Amp < 0.49 
                //&& ch3Int+ch4Int > 140
              ) {
              histIntegral->Fill(ch3Int+ch4Int);
              histIntegralFitMin = 145;
              histIntegralFitMax = 200;
              if (ch3Int+ch4Int > 140) {
                dt->Fill(-1*(t1gausroot-0.5*(t3RisingEdge+t4RisingEdge)));
              }
            }   
          }
  
        }



      } else {
        if (run <= 78 && (energy == 32 || energy == 16 || energy == 8 || energy == 4) ) {
          if( ch3QualityBit==0 &&       
              ch4Amp > 0.01 &&
              ch2Amp > 0.02 && ch2Amp < 0.49 &&
              ch1Amp > 0.02 && ch1Amp < 0.49 
	    ) {
//             dt->Fill((t3gausroot-t2RisingEdge+t3gausroot-t1RisingEdge)/2);
//             dt->Fill((t3gausroot-t1RisingEdge));
             dt->Fill((t3gausroot-t2RisingEdge));
          }
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
  dt->SetMaximum(dt->GetMaximum()*1.2);
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



  if (plotEnergyPeak) {
    //Energy plot
    c = new TCanvas("c","c",600,600);

    histIntegral->SetAxisRange((histIntegralFitMin+histIntegralFitMax)*0.2,(histIntegralFitMin+histIntegralFitMax)*0.75,"X");
    histIntegral->SetTitle("");
    histIntegral->GetXaxis()->SetTitle("Pulse Integral [V*ns]");
    histIntegral->GetYaxis()->SetTitle("Number of Events");
    histIntegral->GetYaxis()->SetTitleOffset(1.25);
    histIntegral->Draw();
    histIntegral->SetStats(0);
    histIntegral->Fit("gaus","","",histIntegralFitMin,histIntegralFitMax);
    TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->SetTextColor(kBlack);
    tex->DrawLatex(0.13, 0.80, Form("Resolution = %.1f %s",100*fitter->GetParameter(2)/fitter->GetParameter(1), "%"));
    //tex->DrawLatex(0.53, 0.75, Form("RMS/Mean=%.1f %s",100*histIntegral->GetRMS()/histIntegral->GetMean(),"%"));

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


void MakeTimeResolutionVsEnergyPlot_DSB() {

  //use beam energy for xaxis
  float x[4] = {4,8,16,32};
  float xerr[4] = { 0.02*4, 0.02*8, 0.02*16, 0.02*32 };
  
//   // use amplitude for x axis 
//   float x[4] = {200,370,700,1400};
//   float xerr[4] = { 200*0.02, 370*0.02, 700*0.02, 1400*0.02 };
  
  // early runs with possible vertical misalignment with beam
//    float y[4] = {489,301,198, 144 };
//    float yerr[4] = { 32,23,12,9 };

//   // later runs with better beam alignment
   float y[4] = {296,222,149, 104 };
   float yerr[4] = { 61,14,10,5 };

  TGraphErrors *graph = new TGraphErrors(4,x,y,xerr,yerr);
  graph->SetLineWidth(2);
  TCanvas * c = new TCanvas("c","c",600,600);
  graph->Draw("AP");

  graph->SetTitle("");
  graph->GetXaxis()->SetTitle("Electron Beam Energy [GeV]");
//   graph->GetXaxis()->SetTitle("Amplitude");
  graph->GetYaxis()->SetTitle("Time Resolution [ps]");
  graph->GetYaxis()->SetTitleOffset(1.25);


  TF1 *func = new TF1("TOFResolutionFunction",TOFResolutionFunction,0,100,2);
   // set the parameters to the mean and RMS of the histogram
   //func->SetParameters(0,1,1);
   // give the parameters meaningful names
   //func->SetParNames ("Constant","Mean_value","Sigma");

   // call TH1::Fit with the name of the TF1 object 
  func->SetParLimits(1,0,1000);
  graph->Fit("TOFResolutionFunction");

//   TPaveText *pt = new TPaveText(.25,.55,.90,.85,"NDC");
//   pt->AddText(Form("TOF Resolution = (%.0f +/- %.0f ps) / sqrt(E) + (%.0f +/- %.0f ps)  ",func->GetParameter(0), func->GetParError(0),func->GetParameter(1), func->GetParError(10)));
//   pt->SetFillColor(0);
//   pt->SetShadowColor(0);
//   pt->SetBorderSize(0);
//   pt->Draw();

//   TLatex *tex = new TLatex();
//   tex->SetNDC();
//   tex->SetTextSize(0.040);
//   tex->SetTextFont(42);
//   tex->SetTextColor(kBlack);
//   tex->DrawLatex(0.30, 0.80, "TOF Resolution = ");
//   tex->DrawLatex(0.30, 0.70, Form("#frac{(%.0f +/- %.0f ps)}{#sqrt{E}} + (%.0f +/- %.0f ps)", func->GetParameter(0), func->GetParError(0),func->GetParameter(1), func->GetParError(1)));

//   TLatex *tex2 = new TLatex();
//   tex2->SetNDC();
//   tex2->SetTextSize(0.045);
//   tex2->SetTextFont(42);
//   tex2->SetTextColor(kBlack);
//   tex2->DrawLatex(0.15, 0.92, "Caltech Internal");




  c->SaveAs( "TimeResolutionVsEnergy_ShashlikDSB1Fiber.gif" );
  c->SaveAs( "TimeResolutionVsEnergy_ShashlikDSB1Fiber.pdf" );

}


void MakeTimeResolutionVsEnergyPlot_Y11() {

  float x[3] = {4,8,16};
  float xerr[3] = { 0.02*4, 0.02*8, 0.02*16 };
  
  //fit to full range
   float y[3] = {631,462,311 };
   float yerr[3] = { 84,17,11 };

  //fit to restricted range near peak
//   float y[4] = {443,268,192, 126 };
//   float yerr[4] = { 51,24,23,15 };

  TGraphErrors *graph = new TGraphErrors(4,x,y,xerr,yerr);
  graph->SetLineWidth(2);
  TCanvas * c = new TCanvas("c","c",600,600);
  graph->Draw("AP");

  graph->SetTitle("");
  graph->GetXaxis()->SetTitle("Electron Beam Energy [GeV]");
  graph->GetYaxis()->SetTitle("Time of Flight Resolution [ps]");
  graph->GetYaxis()->SetTitleOffset(1.35);
  graph->GetYaxis()->SetRangeUser(0,1000);

  TF1 *func = new TF1("TOFResolutionFunction",TOFResolutionFunction,0,100,2);
   // set the parameters to the mean and RMS of the histogram
   //func->SetParameters(0,1,1);
   // give the parameters meaningful names
   //func->SetParNames ("Constant","Mean_value","Sigma");

   // call TH1::Fit with the name of the TF1 object 
  func->SetParLimits(1,0,1000);
  graph->Fit("TOFResolutionFunction");

//   TPaveText *pt = new TPaveText(.25,.55,.90,.85,"NDC");
//   pt->AddText(Form("TOF Resolution = (%.0f +/- %.0f ps) / sqrt(E) + (%.0f +/- %.0f ps)  ",func->GetParameter(0), func->GetParError(0),func->GetParameter(1), func->GetParError(10)));
//   pt->SetFillColor(0);
//   pt->SetShadowColor(0);
//   pt->SetBorderSize(0);
//   pt->Draw();

//   TLatex *tex = new TLatex();
//   tex->SetNDC();
//   tex->SetTextSize(0.040);
//   tex->SetTextFont(42);
//   tex->SetTextColor(kBlack);
//   tex->DrawLatex(0.30, 0.80, "TOF Resolution = ");
//   tex->DrawLatex(0.30, 0.70, Form("#frac{(%.0f +/- %.0f ps)}{#sqrt{E}} + (%.0f +/- %.0f ps)", func->GetParameter(0), func->GetParError(0),func->GetParameter(1), func->GetParError(1)));

//   TLatex *tex2 = new TLatex();
//   tex2->SetNDC();
//   tex2->SetTextSize(0.045);
//   tex2->SetTextFont(42);
//   tex2->SetTextColor(kBlack);
// //   tex2->DrawLatex(0.15, 0.92, "Caltech Internal");




  c->SaveAs( "TimeResolutionVsEnergy_ShashlikY11Fiber.gif" );

}


void MakeTimeResolutionVsEnergyPlot_DSBAndCube() {

  //use beam energy for xaxis
//   // use amplitude for x axis 
   float x[4] = {200,370,700,1400};
   float xerr[4] = { 200*0.02, 370*0.02, 700*0.02, 1400*0.02 };
  
//   // later runs with better beam alignment
   float y[4] = {296,222,149, 104 };
   float yerr[4] = { 61,14,10,5 };

  TGraphErrors *graph = new TGraphErrors(4,x,y,xerr,yerr);
  graph->SetLineWidth(2);
  TCanvas * c = new TCanvas("c","c",600,600);
  c->DrawFrame(50,0,5000,700);
  c->SetLogx();

  graph->Draw("P");
  graph->SetMarkerColor(kRed);
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.25);
  graph->SetLineColor(kRed);


  graph->SetTitle("");
  graph->GetXaxis()->SetTitle("Amplitude");
  graph->GetYaxis()->SetTitle("Time of Flight Resolution [ps]");
  graph->GetYaxis()->SetTitleOffset(1.25);


  TF1 *func = new TF1("TOFResolutionFunction",TOFResolutionFunction,0,100,2);
   // set the parameters to the mean and RMS of the histogram
   //func->SetParameters(0,1,1);
   // give the parameters meaningful names
   //func->SetParNames ("Constant","Mean_value","Sigma");

   // call TH1::Fit with the name of the TF1 object 
  func->SetParLimits(1,0,1000);
   graph->Fit("TOFResolutionFunction");




  float x2[4] = {170*3.0,300*3.0,600*3.0,1200*3.0};
  float x2err[4] = { 0.02*150*3.0, 0.02*300*3.0, 0.02*600*3.0, 0.02*1200*3.0 };
  
//   // later runs with better beam alignment
  float y2[4] = {70.4,57.4, 40.6, 33.5 };
  float y2err[4] = { 2.9,2.0,1.5,2.1 };

  TGraphErrors *graph2 = new TGraphErrors(4,x2,y2,x2err,y2err);
  graph2->SetLineWidth(2);
  graph2->Draw("P");
  graph2->SetMarkerColor(kBlue);
  graph2->SetLineColor(kBlue);
  graph2->SetMarkerStyle(20);
  graph2->SetMarkerSize(1.25);

  graph2->SetTitle("");
//   graph2->GetXaxis()->SetTitle("Electron Beam Energy [GeV]");
  graph->GetXaxis()->SetTitle("Amplitude");
  graph2->GetYaxis()->SetTitle("Time of Flight Resolution [ps]");
  graph2->GetYaxis()->SetTitleOffset(1.25);


  TF1 *func2 = new TF1("TOFResolutionFunction2",TOFResolutionFunction,0,100,2);
  func2->SetLineColor(kBlue);
   // set the parameters to the mean and RMS of the histogram
   //func->SetParameters(0,1,1);
   // give the parameters meaningful names
   //func->SetParNames ("Constant","Mean_value","Sigma");

   // call TH1::Fit with the name of the TF1 object 
  func2->SetParLimits(1,0,1000);
   graph2->Fit("TOFResolutionFunction2");


  TLegend *legend = new TLegend (0.4,0.7,0.8,0.85);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(graph, "Shashlik cell w/DSB1 Fibers","LP");
  legend->AddEntry(graph2, "LYSO Cube","LP");
  legend->Draw();


  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.035);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.40, 0.65, "Shashlik Fiber TOF resolution  ");
  tex->DrawLatex(0.40, 0.60, "is about 2 times worse than ");
  tex->DrawLatex(0.40, 0.55, "LYSO cube at the same amplitude. ");
  
  tex->SetTextSize(0.030);
  tex->DrawLatex(0.40, 0.50, "Note: Rise time for DSB1 Fibers is");
  tex->DrawLatex(0.40, 0.45, "also 2 times slower.");



  c->SaveAs( "TimeResolutionVsEnergy_ShashlikDSB1FiberAndCube.gif" );
  c->SaveAs( "TimeResolutionVsEnergy_ShashlikDSB1FiberAndCube.pdf" );

}




void ShashlikFiberAnalysis() {

//   MakeTimeResolutionPlot("cpt_aug_run_074.ana.root","TOF_ShashlikDSB1Fiber_Electron_16GeV",16,true,74,false);
//        MakeTimeResolutionPlot("cpt_aug_run_075.ana.root","TOF_ShashlikDSB1Fiber_Electron_8GeV",8,true,75,false);
//         MakeTimeResolutionPlot("cpt_aug_run_079.ana.root","TOF_ShashlikDSB1Fiber_Electron_4GeV",4,true,79,false);
//        MakeTimeResolutionPlot("cpt_aug_run_078.ana.root","TOF_ShashlikDSB1Fiber_Electron_32GeV",32,true,78,false);

  //*************************************
  //Best Results Here for Paper
  //*************************************
   MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_aug_run_090And091.ana.root","TOF_ShashlikDSB1Fiber_Electron_32GeV",32,true,90, true);
//   MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_aug_run_088.ana.root","TOF_ShashlikDSB1Fiber_Electron_16GeV",16,true,88, false);
//   MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_aug_run_089.ana.root","TOF_ShashlikDSB1Fiber_Electron_8GeV",8,true,89, false);
//   MakeTimeResolutionPlot("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/Timing/cpt-aug-2014/cpt_aug_run_092.ana.root","TOF_ShashlikDSB1Fiber_Electron_4GeV",4,true,92, false);
  //MakeTimeResolutionVsEnergyPlot_DSB();
  // MakeTimeResolutionVsEnergyPlot_DSBAndCube();


  //*************************************
  //Results for Y11 Fibers from May Testbeam
  //*************************************
//   MakeTimeResolutionPlot("cpt_may_run_131.ana.root","TOF_ShashlikY11Fiber_Electron_4GeV",4,false);
//   MakeTimeResolutionPlot("cpt_may_run_132.ana.root","TOF_ShashlikY11Fiber_Electron_8GeV",8,false);
//   MakeTimeResolutionPlot("cpt_may_run_133.ana.root","TOF_ShashlikY11Fiber_Electron_16GeV",16,false);
//   MakeTimeResolutionVsEnergyPlot_Y11();




}
