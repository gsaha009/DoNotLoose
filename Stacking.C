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

void set_hstyle(TH1D* th, int icol, int itype, int iline);
void read_lines(std::string tname, std::vector<std::string>& files);

void makeStack (TString text_file, TString histName, TString stackName, int rebin, int sigAmpl) {
  gStyle->SetPalette(kOcean);
  THStack *hs = new THStack("hs","");
  TLegend *legend = new TLegend(0.35,0.72,0.70,0.92);
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
    //    hi->SetFillColor(irow);
    set_hstyle(hi, (2*irow)+20, 21, 1);
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
  TCanvas *cst = new TCanvas("cst","stacked hists",10,10,700,700);
  cst->cd(1);
  hs->SetTitle(stackName);
  hs->Draw("HIST");
  hsig->Draw("same HIST");
  legend->Draw();
  gPad->SetLogy();
  hs->Write();
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

void set_hstyle(TH1D* th, int icol, int istyle, int iline) {
  //th->SetTitle(" ");

  if (icol == 10) icol = 46;

  th->SetMarkerStyle(istyle);
  th->SetLineStyle(iline);
  th->SetLineWidth(2);
  th->SetFillColor(icol);
  //th->SetMarkerSize(0.8);
  th->SetLineColor(icol);
  //  th->GetXaxis()->SetTitle(xlab);
  //  th->GetYaxis()->SetTitle(ylab);
}

