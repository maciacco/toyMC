std::unique_ptr<TH2D> divide(TH2D* const hNum, TH2D* const hDen){
  // TODO: check histo dims
  if (hNum->GetNbinsX() != hDen->GetNbinsX() || hNum->GetNbinsY() != hDen->GetNbinsY()){
    std::cout << "Axes do not match!" << std::endl;
    return std::move(nullptr);
  }

  auto res = std::make_unique<TH2D>("hEff2D", ";#it{p}_{T} (GeV/#it{c});#eta", 50, 0., 5, 12, 2., 5.);

  for (int iX{1}; iX <= hNum->GetNbinsX(); ++iX){
    for (int iY{1}; iY <= hNum->GetNbinsY(); ++iY){
      double eff = hNum->GetBinContent(iX, iY) / hDen->GetBinContent(iX, iY);
      res->SetBinContent(iX, iY, eff);
      res->SetBinError(iX, iY, std::sqrt(eff * (1. - eff) / hDen->GetBinContent(iX, iY)));
    }
  }
  return std::move(res);
}

void eff(){
  ROOT::EnableImplicitMT(4);
  ROOT::RDataFrame df("hyper", "/home/mciacco/Hypernuclei-Signal-ntuple.root");
  auto hGen = df.Filter("etaMC > 2 && etaMC < 5").Histo2D({"hGen", ";#it{p}_{T} (GeV/#it{c});#eta", 50, 0., 5, 12, 2., 5.}, "ptMC", "etaMC");
  auto hRec = df.Filter("etaMC > 2 && etaMC < 5 && eta > 2 && eta < 5 && reconstructed == 2 && cosPA > 0.99995").Histo2D({"hRec", ";#it{p}_{T} (GeV/#it{c});#eta", 50, 0., 5, 12, 2., 5.}, "pt", "eta");

  TH1D* proj_rec[2];
  for (int iC{0}; iC < 2; ++iC){
    TH1D* proj_gen = static_cast<TH1D*>(hGen->ProjectionX("gen"));
    proj_rec[iC] = static_cast<TH1D*>(hRec->ProjectionX(Form("rec_%d", iC)));
    proj_rec[iC]->Divide(proj_rec[iC], proj_gen, 1., 1., "B");
    proj_rec[iC]->SetName(Form("hEff_%d", iC));
  }

  auto hEff2D = divide(hRec.GetPtr(), hGen.GetPtr());
  TFile *fo = TFile::Open("fo_ef.root", "recreate");
  fo->cd();
  hGen->Write();
  hRec->Write();
  hEff2D->Write();
  for (int i{0}; i < 2; ++i) proj_rec[i]->Write();
  fo->Close();
}
