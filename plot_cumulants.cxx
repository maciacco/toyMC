#include <iostream>
#include <string>

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>

const std::string histos[]{"hk2sknet", "hk3k1net", "hk4k2net"};
const std::string labels[]{"#frac{#it{#kappa}_{2}(#Lambda - #bar{#Lambda})}{#LT#Lambda + #bar{#Lambda}#GT}", "#frac{#it{#kappa}_{3}(#Lambda - #bar{#Lambda})}{#it{#kappa}_{1}(#Lambda - #bar{#Lambda})}", "#frac{#it{#kappa}_{4}(#Lambda - #bar{#Lambda})}{#it{#kappa}_{2}(#Lambda - #bar{#Lambda})}"};
constexpr double nSample = 40.;

void plot_cumulants(std::string const in_file_name = "ex_lambda.root"){
  gStyle->SetOptStat(0);

  TFile* in_file = TFile::Open(in_file_name.data());
  if (!in_file || in_file->TestBit(TFile::kZombie)){
    std::cout << "No input file!" << std::endl;
    return;
  }

  TH1F* hCumulSubs[3];
  for (int iC{0}; iC < 3; ++iC) {
    hCumulSubs[iC] = static_cast<TH1F*>(in_file->Get(histos[iC].data()));
    if (!hCumulSubs[iC]) {
      std::cout << "No histos!" << std::endl;
      return;
    }
  }

  std::unique_ptr<TH1D> cumulants = std::make_unique<TH1D>("cumulants", ";;Cumulants ratio", 3, 0., 3.);
  for (int iB{1}; iB < cumulants->GetNbinsX() + 1; ++iB){
    cumulants->GetXaxis()->SetBinLabel(iB, labels[iB - 1].data());
    cumulants->SetBinContent(iB, hCumulSubs[iB - 1]->GetMean());
    cumulants->SetBinError(iB, hCumulSubs[iB - 1]->GetStdDev() / std::sqrt(nSample));
  }

  cumulants->SetMarkerStyle(20);
  cumulants->SetMarkerSize(1.);
  cumulants->SetLineColor(kRed);
  cumulants->SetMarkerColor(kRed);
  cumulants->SetLineWidth(2);
  cumulants->GetXaxis()->SetLabelFont(45);
  cumulants->GetXaxis()->SetLabelSize(40);
  cumulants->GetYaxis()->SetTitleOffset(1.3);
  cumulants->GetYaxis()->SetTitleFont(45);
  cumulants->GetYaxis()->SetTitleSize(30);
  cumulants->GetYaxis()->SetRangeUser(0., 1.1);

  TCanvas c("c", "c", 600, 600);
  c.SetTopMargin(0.02);
  c.SetRightMargin(0.02);
  c.SetBottomMargin(0.15);
  c.SetLeftMargin(0.13);

  TLatex txt;
  txt.SetTextFont(45);
  txt.SetTextSize(28);
  txt.SetNDC();

  c.cd();
  cumulants->Draw("pex0");
  txt.DrawLatex(0.19 /*0.45*/, 0.46/*0.92*/, "NA60+ Performance");
  txt.DrawLatex(0.19 /*0.45*/, 0.38/*0.84*/, "Pb#minusPb, #sqrt{#it{s}_{NN}} = 6.3 GeV, 0-20%");
//  txt.DrawLatex(0.22, 0.38, "0-10%, 10^{10} events");
  txt.DrawLatex(0.19, 0.3, "2 < #it{#eta} < 4");
  txt.DrawLatex(0.19, 0.22, "0.2 #leq #it{p}_{T} < 2.0 GeV/#it{c}");
  c.Print("final_plot.pdf");
}
