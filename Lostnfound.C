#include <Riostream.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TString.h>
#include <TPaveText.h>
#include <TText.h>
#include <TLegend.h>
#include <TMath.h>

#include <iostream>
#include <ios>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <set>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

void makeStack (TString text_file, TString histName, TString stackName, int rebin, int sigAmpl);
void makeNormalised (TString text_file, TString histName, TString normName, int rebin);
void makeROC (TString text_file, TString histName, TString label, bool MinToX, bool XToMax);
void makeSignificance (TString text_file, TString histName, TString label, bool MinToX, bool XToMax);
void checkOverTraining (TString text_file, TString label, TString methodName, TString weightfile);  
void set_hstyle(TH1D* th, int icol, int itype, int iline, TString xlab, TString ylab, bool stack);
void read_lines(std::string tname, std::vector<std::string>& files);

void README(){
  std::cout<<"\n"
	   <<">>>MAIN Functions::\n"
	   <<"makeStack        (TString text_file, TString histName, TString stackName, int rebin, int sigAmpl)\n"
	   <<"makeNormalised   (TString text_file, TString histName, TString normName, int rebin)\n"
	   <<"makeROC          (TString text_file, TString histName, TString label, bool MinToX, bool XToMax)\n"
	   <<"makeSignificance (TString text_file, TString histName, TString label, bool MinToX, bool XToMax)\n"
	   <<"checkOverTraining (TString text_file, TString label, TString methodName, TString weightfile)\n"
	   <<"\n"
	   <<">>>Auxiliary Functions::\n"
	   <<"set_hstyle  (TH1D* th, int icol, int itype, int iline, TString xlab, TString ylab, bool stack)\n"
	   <<"read_lines  (std::string tname, std::vector<std::string>& files)\n"
	   <<"\n";
}

void makeStack (TString text_file, TString histName, TString stackName, int rebin, int sigAmpl) {
  gStyle->SetPalette(kOcean);
  THStack *hs = new THStack("hs","");
  TLegend *legend = new TLegend(0.7162726,0.4518072,0.9965229,0.9927711,NULL,"brNDC");
  legend->SetBorderSize(1);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(1001);
  legend->SetHeader("","");
  
  std::vector<std::string> lines;
  read_lines(text_file.Data(), lines);

  TFile* file_new = new TFile(stackName+".root", "RECREATE");

  size_t irow = 0;
  for (auto& ln: lines) {
    TString it = ln;
    it += "_hist.root";
    std::cout<<it<<"\n";
    TFile* file = TFile::Open(it);
    if (!file) continue;
    irow++;
    TH1D *hi = (TH1D*)(file->Get(histName));
    if (irow == 1) continue;
    set_hstyle(hi, 2*irow+3, 21, 1, stackName, "", 1);
    hi->Rebin(rebin);
    hs->Add(hi);
    legend->AddEntry (hi, TString(ln), "f");
  }
  TFile* file = TFile::Open(TString(lines.at(0))+"_hist.root");
  TH1D *hsig = (TH1D*)(file->Get(histName));
  hsig->Scale(sigAmpl);
  hsig->Rebin(rebin);
  hsig->SetLineStyle(7);
  hsig->SetLineWidth(3);
  hsig->SetLineColor(1);
  legend->AddEntry (hsig, "Signal", "l");
  TCanvas *cst = new TCanvas("cst","stacked hists",1600,1200);
  cst->cd(1);
  hs->SetTitle(stackName);
  hs->Draw("HIST");
  hsig->Draw("same HIST");
  legend->Draw();
  gPad->SetLogy();
  file_new->cd();
  hs->Write();
}

void makeNormalised (TString text_file, TString histName, TString normName, int rebin) {
  TLegend *legend = new TLegend(0.7162726,0.4518072,0.9965229,0.9927711,NULL,"brNDC");
  legend->SetBorderSize(1);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(1001);
  legend->SetHeader("","");
  
  std::vector<std::string> lines;
  read_lines(text_file.Data(), lines);

  size_t irow = 0;
  TCanvas *cst = new TCanvas("cst","morm hists",1600,1000);
  gStyle->SetOptStat(0);
  gPad->SetGrid();
  for (auto& ln: lines) {
    TString it = ln;
    it += "_hist.root";
    std::cout<<it<<"\n";
    TFile* file = TFile::Open(it);
    if (!file) continue;
    irow++;
    TH1D *hi = (TH1D*)(file->Get(histName));
    hi->Scale(1/hi->Integral());
    hi->Rebin(rebin);
    set_hstyle(hi, irow, 21, 1, normName, "", 0);
    legend->AddEntry (hi, TString(ln), "l");
    cst->cd(1);
    if (irow == 1) {
      hi->SetLineWidth(3);
      hi->GetXaxis()->SetNdivisions(512);
      hi->Draw("HIST");
      hi->GetYaxis()->SetRangeUser(0.0,0.5);
    }
    else{
      hi->GetXaxis()->SetNdivisions(512);
      hi->Draw("same HIST");
      hi->GetYaxis()->SetRangeUser(0.0,0.5);
    }
    legend->Draw();
  }
}

void makeROC (TString text_file, TString histName, TString label, bool MinToX, bool XToMax) {
  std::vector<std::string> lines;
  read_lines(text_file.Data(), lines);

  std::vector<TH1D*>hvec; 
  for (auto& ln: lines) {
    TString it = ln;
    it += "_hist.root";
    std::cout<<it<<"\n";
    TFile* file = TFile::Open(it);
    if (!file) continue;
    TH1D *hi = (TH1D*)(file->Get(histName));
    hvec.push_back(hi);
  }
  std::cout<<hvec.size()<<" Histograms are in histVec\n";

  int nBins = hvec[0]->GetNbinsX();
  float min = hvec[0]->GetBinCenter(1);
  float max = hvec[0]->GetBinCenter(nBins);
  std::cout<<nBins<<"   "<<min<<"   "<<max<<"\n";
  
  TH1F *sigEff_hist = new TH1F ("sigEff", "", nBins, min, max);
  TH1F *bkgEff_hist = new TH1F ("bkgEff", "", nBins, min, max);
  TH1F *bkgRej_hist = new TH1F ("bkgRej", "", nBins, min, max);
  TH1F *sigPur_hist = new TH1F ("sigPurity", "", nBins, min, max);
  TH1F *sigPurEff_hist = new TH1F ("sigPurity*sigEfficiency", "", nBins, min, max);
  TH1F *signf_hist  = new TH1F ("significance", "", nBins, min, max);
  TH1F *ROC_hist    = new TH1F ("ROC", "", nBins, 0.0, 1.0);

  std::cout<<"Cut: "<<setw(16)<<"nSignalEvt"<<setw(16)<<"SigEff"<<setw(16)<<"nBkgEvt"<<setw(16)<<"bkgEff"<<setw(16)<<"bkgRej"<<setw(16)<<"signficance"<<"\n";
  
  for (size_t ib = 1; ib <= nBins; ++ib) {
    int irow = 0;
    float S = 0.0;
    if (XToMax) S = hvec[0]->Integral(ib, nBins);
    else if (MinToX) S = hvec[0]->Integral(0, ib);
    float B = 0.0;
    float allBkgInt = 0.0;
    for (auto& h: hvec){
      irow++;
      if (irow == 1) continue; 
      allBkgInt += h->Integral();
      if (XToMax) B += h->Integral(ib, nBins);
      else if (MinToX) B += h->Integral(0, ib);
    }
    if ((S+B) == 0.0) continue;
    float signf  = S/TMath::Sqrt(S+B);
    float sigEff = S/hvec[0]->Integral();
    float bkgEff = B/allBkgInt;
    float bkgRej = 1.0 - bkgEff;
    float sigPur = S/(S+B);
    float sigPurEff = sigPur*sigEff;
    std::cout<<std::setprecision(5)
	     <<hvec[0]->GetBinCenter(ib)
	     <<setw(16)<<S
	     <<std::setprecision(5)
	     <<setw(16)<<sigEff
	     <<setw(16)<<B
	     <<std::setprecision(5)
	     <<setw(16)<<bkgEff
	     <<std::setprecision(5)
	     <<setw(16)<<bkgRej
	     <<std::setprecision(5)
	     <<setw(16)<<signf<<std::endl;
      
    sigEff_hist->SetBinContent(ib, sigEff);
    signf_hist ->SetBinContent(ib, signf);
    sigPur_hist->SetBinContent(ib, sigPur);
    sigPurEff_hist->SetBinContent(ib, sigPurEff);
    bkgEff_hist->SetBinContent(ib, bkgEff);
    bkgRej_hist->SetBinContent(ib, bkgRej);
    ROC_hist ->SetBinContent(ROC_hist->FindBin(sigEff), bkgRej);
  }
  
  
  gStyle->SetOptStat(0);
  TCanvas *cst = new TCanvas("cst","norm hists",1600,1200);
  cst->Divide(2,2);
  cst->cd(1);
  gPad->SetGrid();
  TLegend *leg = new TLegend(0.7904787,0.3876892,0.9817165,0.7948981,NULL,"brNDC");
  sigEff_hist->SetLineColor(kBlue);
  sigEff_hist->SetLineWidth(2);
  sigEff_hist->SetXTitle(label);
  sigEff_hist->SetYTitle("Efficiency (Purity)");
  sigEff_hist->SetTitle("Efficiency (Purity) Vs. "+label);
  sigEff_hist->GetXaxis()->SetNdivisions(512);
  leg -> AddEntry (sigEff_hist, "SigEff", "l");
  sigEff_hist->Draw("HIST");

  bkgEff_hist->SetLineColor(kMagenta);
  bkgEff_hist->SetLineWidth(2);
  bkgEff_hist->GetXaxis()->SetNdivisions(512);
  leg -> AddEntry (bkgEff_hist, "BkgEff", "l");
  bkgEff_hist ->Draw("same HIST");

  sigPur_hist->SetLineColor(kBlack);
  sigPur_hist->SetLineWidth(2);
  sigPur_hist->GetXaxis()->SetNdivisions(512);
  leg -> AddEntry (sigPur_hist, "SigPur", "l");
  sigPur_hist->Draw("same HIST");

  sigPurEff_hist->SetLineColor(kGreen+2);
  sigPurEff_hist->SetLineWidth(2);
  sigPurEff_hist->GetXaxis()->SetNdivisions(512);
  leg -> AddEntry (sigPurEff_hist, "Sig(Pur*Eff)", "l");
  sigPurEff_hist->Draw("same HIST");
  leg->Draw();
  
  cst->cd(2);
  gPad->SetGrid();
  signf_hist->SetLineColor(kRed+2);
  signf_hist->SetLineWidth(2);
  signf_hist->SetXTitle(label);
  signf_hist->SetYTitle("S/Sqrt(S+B)");
  signf_hist->SetTitle("Significance Vs. "+label);
  signf_hist->GetXaxis()->SetNdivisions(512);
  signf_hist ->Draw("HIST");

  cst->cd(3);
  gPad->SetGrid();
  TLegend *leg2 = new TLegend(0.825284,0.6790663,0.9817513,0.7945281,NULL,"brNDC");
  sigEff_hist->SetLineColor(kBlue);
  sigEff_hist->SetLineWidth(2);
  sigEff_hist->SetXTitle(label);
  sigEff_hist->SetYTitle("Efficiency");
  sigEff_hist->SetTitle("Efficiency Vs. "+label);
  sigEff_hist->GetXaxis()->SetNdivisions(512);
  leg2 -> AddEntry (sigEff_hist, "SigEff", "l");
  sigEff_hist->Draw("HIST");
  bkgRej_hist->SetLineColor(kRed);
  bkgRej_hist->SetLineWidth(2);
  bkgRej_hist->GetXaxis()->SetNdivisions(512);
  leg2 -> AddEntry (bkgRej_hist, "BkgRej", "l");
  bkgRej_hist ->Draw("same HIST");
  leg2->Draw();
  
  cst->cd(4);
  gPad->SetGrid();
  ROC_hist->SetMarkerColor(kCyan+2);
  ROC_hist->SetMarkerStyle(3);
  ROC_hist->SetMarkerSize(0.9);
  ROC_hist->SetXTitle("SignalEfficiency");
  ROC_hist->SetYTitle("BackgroundRejection");
  ROC_hist->SetTitle("ROC");
  ROC_hist->GetXaxis()->SetNdivisions(512);
  ROC_hist ->Draw("P");


  cst->SaveAs(label+"_ROC.png");
}

void makeSignificance (TString text_file, TString histName, TString label, bool MinToX, bool XToMax) {
  std::vector<std::string> lines;
  read_lines(text_file.Data(), lines);

  std::vector<TH1D*>hvec; 
  for (auto& ln: lines) {
    TString it = ln;
    it += "_hist.root";
    std::cout<<it<<"\n";
    TFile* file = TFile::Open(it);
    if (!file) continue;
    TH1D *hi = (TH1D*)(file->Get(histName));
    hvec.push_back(hi);
  }
  std::cout<<hvec.size()<<" Histograms are in histVec\n";

  int nBins = hvec[0]->GetNbinsX();
  float min = hvec[0]->GetBinCenter(1);
  float max = hvec[0]->GetBinCenter(nBins);
  
  TH1F *signf_hist  = new TH1F ("significance", "", nBins, min, max);
  
  for (size_t ib = 1; ib <= nBins; ++ib) {
    int irow = 0;
    float S = 0.0;
    if (XToMax) S = hvec[0]->Integral(ib, nBins);
    else if (MinToX) S = hvec[0]->Integral(0, ib);
    float B = 0.0;
    float allBkgInt = 0.0;
    for (auto& h: hvec){
      irow++;
      if (irow == 1) continue; 
      allBkgInt += h->Integral();
      if (XToMax) B += h->Integral(ib, nBins);
      else if (MinToX) B += h->Integral(0, ib);
    }
    if (B == 0.0) continue;
    float signf = S/TMath::Sqrt(S+B);
    signf_hist->SetBinContent(ib, signf);
  }
  
  gStyle->SetOptStat(0);
  TCanvas *cst = new TCanvas("cst","norm hists",1600,1200);
  cst->cd();
  gPad->SetGrid();
  signf_hist->SetLineColor(kRed);
  signf_hist->SetLineWidth(3);
  signf_hist->SetXTitle(label);
  signf_hist->SetYTitle("S/Sqrt(S+B)");
  signf_hist->SetTitle("Significance Vs. "+label);
  signf_hist->GetXaxis()->SetNdivisions(512);
  signf_hist ->Draw("HIST");

  cst->SaveAs(label+"_ROC.png");
}

void checkOverTraining (TString text_file, TString label, TString methodName, TString weightfile) {
  
  std::cout << "==> Start TMVAClassificationApplication" << std::endl;
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used

  //float met;
  //float HT;
  float HTvec;
  float HTMETfrac;
  float j1j2DR;
  //float j1j3DR;
  float minJJDR;
  float maxJJDR;
  //float minJJDPhi;
  float LT;
  //float LTMETvec;
  float LTMETfracInv;
  //float l1l2InvM;
  float l1l2DR;
  //   float l1l2DPhi;
  float l1MetDPhi;
  //   float l2MetDPhi;
  float j1MetDPhi;
  float j2MetDPhi;
  //float j3MetDPhi;
  float j1l1DPhi;
  //   float j1l2DPhi;
  float minJLDPhi;
  //float maxJLDPhi;
  float maxLepMetDPhi;
  float minJetMetDPhi;

  //reader->AddVariable("met", &met);
  //reader->AddVariable("HT", &HT);
  reader->AddVariable("HTvec", &HTvec);
  reader->AddVariable("HTMETfrac", &HTMETfrac);
  reader->AddVariable("j1j2DR", &j1j2DR);
  //reader->AddVariable("j1j3DR", &j1j3DR);
  reader->AddVariable("minJJDR", &minJJDR);
  reader->AddVariable("maxJJDR", &maxJJDR);
  //   reader->AddVariable("minJJDPhi", &minJJDPhi);
  reader->AddVariable("LT", &LT);
  //   reader->AddVariable("LTMETvec", &LTMETvec);
  reader->AddVariable("LTMETfracInv", &LTMETfracInv);
  //   reader->AddVariable("l1l2InvM", &l1l2InvM);
  reader->AddVariable("l1l2DR", &l1l2DR);
  //   reader->AddVariable("l1l2DPhi", &l1l2DPhi);
  reader->AddVariable("l1MetDPhi", &l1MetDPhi);
  //   reader->AddVariable("l2MetDPhi", &l2MetDPhi);
  reader->AddVariable("j1MetDPhi", &j1MetDPhi);
  reader->AddVariable("j2MetDPhi", &j2MetDPhi);
  //reader->AddVariable("j3MetDPhi", &j3MetDPhi);
  reader->AddVariable("j1l1DPhi", &j1l1DPhi);
  //   reader->AddVariable("j1l2DPhi", &j1l2DPhi);
  reader->AddVariable("minJLDPhi", &minJLDPhi);
  //   reader->AddVariable("maxJLDPhi", &maxJLDPhi);
  reader->AddVariable("maxLepMetDPhi", &maxLepMetDPhi);
  reader->AddVariable("minJetMetDPhi", &minJetMetDPhi);
  //   reader->AddVariable("jInvM", &jInvM);

  reader->BookMVA( methodName, weightfile );

  
  std::vector<std::string> lines;
  read_lines(text_file.Data(), lines);
  
  std::vector<TH1F*>htrain;
  std::vector<TH1F*>htest; 
  for (auto& ln: lines) {
    TString it = ln;
    it += ".root";
    std::cout<<it<<"\n";
    TFile* file = TFile::Open(it);
    if (!file) {
      std::cout<<">>>file not found!!!\n";
      continue;
    }
    TTree *tr = dynamic_cast<TTree*>(file->Get("train"));
    TTree *te = dynamic_cast<TTree*>(file->Get("test"));

    tr->SetBranchAddress("HTvec", &HTvec);
    tr->SetBranchAddress("HTMETfrac", &HTMETfrac);
    tr->SetBranchAddress("j1j2DR", &j1j2DR);
    //tr->SetBranchAddress("j1j3DR", &j1j3DR);
    tr->SetBranchAddress("minJJDR", &minJJDR);
    tr->SetBranchAddress("maxJJDR", &maxJJDR);
    //   tr->SetBranchAddress("minJJDPhi", &minJJDPhi);
    tr->SetBranchAddress("LT", &LT);
    //   tr->SetBranchAddress("LTMETvec", &LTMETvec);
    tr->SetBranchAddress("LTMETfracInv", &LTMETfracInv);
    //   tr->SetBranchAddress("l1l2InvM", &l1l2InvM);
    tr->SetBranchAddress("l1l2DR", &l1l2DR);
    //   tr->SetBranchAddress("l1l2DPhi", &l1l2DPhi);
    tr->SetBranchAddress("l1MetDPhi", &l1MetDPhi);
    //   tr->SetBranchAddress("l2MetDPhi", &l2MetDPhi);
    tr->SetBranchAddress("j1MetDPhi", &j1MetDPhi);
    tr->SetBranchAddress("j2MetDPhi", &j2MetDPhi);
    //tr->SetBranchAddress("j3MetDPhi", &j3MetDPhi);
    tr->SetBranchAddress("j1l1DPhi", &j1l1DPhi);
    //   tr->SetBranchAddress("j1l2DPhi", &j1l2DPhi);
    tr->SetBranchAddress("minJLDPhi", &minJLDPhi);
    //   tr->SetBranchAddress("maxJLDPhi", &maxJLDPhi);
    tr->SetBranchAddress("maxLepMetDPhi", &maxLepMetDPhi);
    tr->SetBranchAddress("minJetMetDPhi", &minJetMetDPhi);
    //   tr->SetBranchAddress("jInvM", &jInvM);

    te->SetBranchAddress("HTvec", &HTvec);
    te->SetBranchAddress("HTMETfrac", &HTMETfrac);
    te->SetBranchAddress("j1j2DR", &j1j2DR);
    //te->SetBranchAddress("j1j3DR", &j1j3DR);
    te->SetBranchAddress("minJJDR", &minJJDR);
    te->SetBranchAddress("maxJJDR", &maxJJDR);
    //   te->SetBranchAddress("minJJDPhi", &minJJDPhi);
    te->SetBranchAddress("LT", &LT);
    //   te->SetBranchAddress("LTMETvec", &LTMETvec);
    te->SetBranchAddress("LTMETfracInv", &LTMETfracInv);
    //   te->SetBranchAddress("l1l2InvM", &l1l2InvM);
    te->SetBranchAddress("l1l2DR", &l1l2DR);
    //   te->SetBranchAddress("l1l2DPhi", &l1l2DPhi);
    te->SetBranchAddress("l1MetDPhi", &l1MetDPhi);
    //   te->SetBranchAddress("l2MetDPhi", &l2MetDPhi);
    te->SetBranchAddress("j1MetDPhi", &j1MetDPhi);
    te->SetBranchAddress("j2MetDPhi", &j2MetDPhi);
    //te->SetBranchAddress("j3MetDPhi", &j3MetDPhi);
    te->SetBranchAddress("j1l1DPhi", &j1l1DPhi);
    //   te->SetBranchAddress("j1l2DPhi", &j1l2DPhi);
    te->SetBranchAddress("minJLDPhi", &minJLDPhi);
    //   te->SetBranchAddress("maxJLDPhi", &maxJLDPhi);
    te->SetBranchAddress("maxLepMetDPhi", &maxLepMetDPhi);
    te->SetBranchAddress("minJetMetDPhi", &minJetMetDPhi);
    //   te->SetBranchAddress("jInvM", &jInvM);


    TH1F *train_response = new TH1F ("train_response", "", 1000, -1.0, 1.0);
    TH1F *test_response  = new TH1F ("test_response", "", 1000, -1.0, 1.0);
    //Loop over training events
    for (size_t i = 0; i < tr->GetEntries(); ++i) {
      tr->GetEntry(i);
      train_response -> Fill(reader->EvaluateMVA(methodName));
    }
    htrain.push_back(train_response);

    //Loop over training events
    for (size_t i = 0; i < te->GetEntries(); ++i) {
      te->GetEntry(i);
      test_response -> Fill(reader->EvaluateMVA(methodName));
    }
    htest.push_back(test_response);
  }

  int nBins = htrain[0]->GetNbinsX();
  float min = htrain[0]->GetBinCenter(1);
  float max = htrain[0]->GetBinCenter(nBins);
  
  TH1F *ROC_train    = new TH1F ("ROC_train", "", nBins, 0.0, 1.0);

  for (size_t ib = 1; ib <= nBins; ++ib) {
    int irow = 0;
    float S = htrain[0]->Integral(ib, nBins);
    float B = 0.0;
    float allBkgInt = 0.0;
    for (auto& h: htrain){
      irow++;
      if (irow == 1) continue; 
      allBkgInt += h->Integral();
      B += h->Integral(ib, nBins);
    }
    if ((S+B) == 0.0) continue;
    float sigEff = S/htrain[0]->Integral();
    float bkgEff = B/allBkgInt;
    float bkgRej = 1.0 - bkgEff;

    ROC_train ->SetBinContent(ROC_train->FindBin(sigEff), bkgRej);
  }

  int nBins_ = htest[0]->GetNbinsX();
  float min_ = htest[0]->GetBinCenter(1);
  float max_ = htest[0]->GetBinCenter(nBins_);
  
  TH1F *ROC_test    = new TH1F ("ROC_test", "", nBins_, 0.0, 1.0);

  for (size_t ib = 1; ib <= nBins_; ++ib) {
    int irow = 0;
    float S = htest[0]->Integral(ib, nBins_);
    float B = 0.0;
    float allBkgInt = 0.0;
    for (auto& h: htest){
      irow++;
      if (irow == 1) continue; 
      allBkgInt += h->Integral();
      B += h->Integral(ib, nBins);
    }
    if ((S+B) == 0.0) continue;
    float sigEff = S/htest[0]->Integral();
    float bkgEff = B/allBkgInt;
    float bkgRej = 1.0 - bkgEff;

    ROC_test ->SetBinContent(ROC_test->FindBin(sigEff), bkgRej);
  }

  
  gStyle->SetOptStat(0);
  TCanvas *cst = new TCanvas("cst","norm hists",1600,1200);
  cst->cd(1);
  gPad->SetGrid();
  TLegend *leg = new TLegend(0.7904787,0.3876892,0.9817165,0.7948981,NULL,"brNDC");
  ROC_train->SetMarkerColor(kBlue);
  ROC_train->SetLineColor(kBlue);
  ROC_train->SetMarkerStyle(4);
  ROC_train->SetMarkerSize(0.5);
  ROC_train->SetXTitle("SignalEfficiency");
  ROC_train->SetYTitle("BackgroundRejection");
  ROC_train->SetTitle("ROC");
  ROC_train->GetXaxis()->SetNdivisions(512);
  leg -> AddEntry (ROC_train, "ROC_Train", "l");
  ROC_train ->Draw("P");
  ROC_test->SetMarkerColor(kRed);
  ROC_test->SetLineColor(kRed);
  ROC_test->SetMarkerStyle(8);
  ROC_test->SetMarkerSize(0.5);
  ROC_test->GetXaxis()->SetNdivisions(512);
  leg -> AddEntry (ROC_test, "ROC_Test", "l");
  ROC_test ->Draw("SAME P");
  leg->Draw();
}

void read_lines(std::string tname, std::vector<std::string>& files) {
  static const int BUF_SIZE = 512;

  ifstream myTxtFile;
  myTxtFile.open(tname.c_str(), ios::in);
  if (!myTxtFile) {
    std::cerr << "Input File: " << tname << " could not be opened!" << std::endl;
    return;
  }

  char buf[BUF_SIZE];
  if (myTxtFile) {
    while (myTxtFile.good()) {
      if (!myTxtFile.eof()) {
        myTxtFile.getline(buf, BUF_SIZE, '\n');
        string line(buf);
        if (line.empty()) continue;
        if (line.substr(0,2) == "//") continue;
        if (line.substr(0,1) == "#") continue;
        files.push_back(line);
      }
    }
  }
}

void set_hstyle(TH1D* th, int icol, int istyle, int iline, TString xlab, TString ylab, bool stack) {
  th->SetTitle(" ");

  if (icol == 10) icol = 46;

  th->SetMarkerStyle(istyle);
  th->SetLineStyle(iline);
  th->SetLineWidth(2);
  if (stack) th->SetFillColor(icol);
  th->SetMarkerSize(0.8);
  th->SetLineColor(icol);
  th->GetXaxis()->SetTitle(xlab);
  th->GetYaxis()->SetTitle(ylab);
}
