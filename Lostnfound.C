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

void makeStack (TString text_file, TString histName, TString stackName, int rebin, int sigAmpl);
void makeNormalised (TString text_file, TString histName, TString normName, int rebin);
void makeROC (TString text_file, TString histName, TString label, bool MinToX, bool XToMax);
void makeSignificance (TString text_file, TString histName, TString label, bool MinToX, bool XToMax);  
void set_hstyle(TH1D* th, int icol, int itype, int iline, TString xlab, TString ylab, bool stack);
void read_lines(std::string tname, std::vector<std::string>& files);

void README(){
  std::cout<<"\n"
	   <<">>>MAIN Functions::\n"
	   <<"makeStack        (TString text_file, TString histName, TString stackName, int rebin, int sigAmpl)\n"
	   <<"makeNormalised   (TString text_file, TString histName, TString normName, int rebin)\n"
	   <<"makeROC          (TString text_file, TString histName, TString label, bool MinToX, bool XToMax)\n"
	   <<"makeSignificance (TString text_file, TString histName, TString label, bool MinToX, bool XToMax)\n"
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
    if (B == 0.0) continue;
    float signf = S/TMath::Sqrt(S+B);
    float sigEff = S/hvec[0]->Integral();
    float bkgEff = B/allBkgInt;
    float bkgRej = 1.0 - bkgEff;
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
    bkgEff_hist->SetBinContent(ib, bkgEff);
    bkgRej_hist->SetBinContent(ib, bkgRej);
    ROC_hist ->SetBinContent(ROC_hist->FindBin(sigEff), bkgRej);
  }
  
  
  gStyle->SetOptStat(0);
  TCanvas *cst = new TCanvas("cst","norm hists",1600,1200);
  cst->Divide(2,2);
  cst->cd(1);
  gPad->SetGrid();
  gPad->SetFillColorAlpha(kGreen, 0.10);
  TLegend *leg = new TLegend(0.725284,0.5790663,0.8817513,0.6945281,NULL,"brNDC");
  sigEff_hist->SetLineColor(kBlue);
  sigEff_hist->SetLineWidth(3);
  sigEff_hist->SetXTitle(label);
  sigEff_hist->SetYTitle("SignalEfficiency & BkgRejection");
  sigEff_hist->SetTitle("SignalEfficiency & BkgRejection Vs. "+label);
  sigEff_hist->GetXaxis()->SetNdivisions(512);
  leg -> AddEntry (sigEff_hist, "SigEff", "l");
  sigEff_hist->Draw("HIST");
  bkgRej_hist->SetLineColor(kRed);
  bkgRej_hist->SetLineWidth(3);
  bkgRej_hist->SetXTitle(label);
  bkgEff_hist->GetXaxis()->SetNdivisions(512);
  leg -> AddEntry (bkgRej_hist, "BkgRej", "l");
  bkgRej_hist->Draw("same HIST");
  leg->Draw();
  
  cst->cd(2);
  gPad->SetGrid();
  gPad->SetFillColorAlpha(kPink, 0.10);
  signf_hist->SetLineColor(kBlack);
  signf_hist->SetLineWidth(3);
  signf_hist->SetXTitle(label);
  signf_hist->SetYTitle("S/Sqrt(S+B)");
  signf_hist->SetTitle("Significance Vs. "+label);
  signf_hist->GetXaxis()->SetNdivisions(512);
  signf_hist ->Draw("HIST");

  cst->cd(3);
  gPad->SetGrid();
  gPad->SetFillColorAlpha(kAzure, 0.10);
  bkgEff_hist->SetLineColor(kMagenta);
  bkgEff_hist->SetLineWidth(3);
  bkgEff_hist->SetXTitle(label);
  bkgEff_hist->SetYTitle("BackgroundEfficiency");
  bkgEff_hist->SetTitle("BkgEfficiency Vs. "+label);
  bkgEff_hist->GetXaxis()->SetNdivisions(512);
  bkgEff_hist ->Draw("HIST");

  cst->cd(4);
  gPad->SetGrid();
  gPad->SetFillColorAlpha(kYellow, 0.10);
  ROC_hist->SetMarkerColor(kGreen+2);
  ROC_hist->SetMarkerStyle(3);
  ROC_hist->SetMarkerSize(0.5);
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
