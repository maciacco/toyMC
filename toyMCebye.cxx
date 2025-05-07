#include <Riostream.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <string>
#include <TH2D.h>
#include <TH3D.h>

constexpr bool save_tree = false;
// constexpr double nEv = 6.e5; // for central (0-10%) --> 6e10 collisions
constexpr double pt_lim[] = {0.2, 2.};
constexpr double eta_lim[] = {2., 4.};

constexpr long int nS = 40;

constexpr double massLambda = 1.115683; // GeV/c^2
constexpr double brLambda = 0.639;

int main(int const argv, char* const argc[]){ //double const nEv = 6.e10, int const S_MIN = 0, int const S_MAX = 20){
  double nEv = std::stod(argc[1]);
  int S_MIN = std::stoi(argc[2]);
  int S_MAX = std::stoi(argc[3]);

  if (S_MIN >= S_MAX) {
    std::cout << "Fatal: S_MIN >= S_MAX!" << std::endl;
    return 1;
  }

  TFile *fDist = TFile::Open("fOut.root"); // distribution file
  TFile *fEff = TFile::Open("fo_ef.root"); // efficiency file
  TH3D *hCPtY = static_cast<TH3D*>(fDist->Get("hPtY"));
  TH2D *hnLaL = static_cast<TH2D*>(fDist->Get("hLAntiL"));
  if (!hCPtY || !hnLaL) {
    std::cout << "No input histos!" << std::endl;
    return -1;
  }
  TH2D *hPtY[2];
  TH2D *hEff[2];
  for (int iC{0}; iC < 2; ++iC){
    hCPtY->GetXaxis()->SetRangeUser(iC, iC + 1);
    hPtY[iC] = static_cast<TH2D*>(hCPtY->Project3D("zy"));
    hPtY[iC]->SetName(Form("hPtY_%d", iC));
    hEff[iC] = static_cast<TH2D*>(fEff->Get("hEff2D"));
    hEff[iC]->SetName(Form("hEff_%d", iC));

    if (!hEff[iC] || !hPtY[iC]) {
      std::cout << "No input histos (2)!" << std::endl;
      return -2;
    }
  }

  gRandom->SetSeed(10 * S_MIN);

  std::string fname{"ex_lambda_"};
  fname = fname + std::to_string(S_MIN) + std::string{".root"};
  TFile fo(fname.data(), "recreate");

  double n1_a[nS]{0};
  double n1_b[nS]{0};

  double n2_a[nS]{0};
  double n2_b[nS]{0};

  double ne2_a[nS]{0};
  double ne2_b[nS]{0};

  double n11_ab[nS]{0};

  double q1_1_1[nS]{0.};
  double q2_1_1[nS]{0.};
  double q1_2_1[nS]{0.};
  double q1_2_2[nS]{0.};
  double q3_1_1[nS]{0.};
  double q1_1_1_x_q1_2_1[nS]{0.};
  double q1_1_1_x_q1_2_2[nS]{0.};
  double q1_3_1[nS]{0.};
  double q1_3_2[nS]{0.};
  double q1_3_3[nS]{0.};
  double q4_1_1[nS]{0.};
  double q2_1_1_x_q1_2_1[nS]{0.};
  double q2_1_1_x_q1_2_2[nS]{0.};
  double q1_1_1_x_q1_3_1[nS]{0.};
  double q2_2_1[nS]{0.};
  double q2_2_2[nS]{0.};
  double q1_1_1_x_q1_3_2[nS]{0.};
  double q1_1_1_x_q1_3_3[nS]{0.};
  double q1_2_1_x_q1_2_2[nS]{0.};
  double q1_4_1[nS]{0.};
  double q1_4_2[nS]{0.};
  double q1_4_3[nS]{0.};
  double q1_4_4[nS]{0.};

  for (long int i{0}; i < nEv; ++i){
    if (i%1000000 == 0) std::cout << i << "\n";
    double k1[2]{0.};
    hnLaL->GetRandom2(k1[0], k1[1]);

    double q_[2][4]{0.};

    auto qa_b_c = [&](int const a, int const b, int const c) -> double
    {
      double sgn = (b % 2) == 0 ? 1. : -1.;
      return std::pow(q_[0][c - 1] + sgn * q_[1][c - 1], a);
    };

//    std::cout << k1[0] << "\t" << k1[1] << std::endl;

    for (int iC{0}; iC < 2; ++iC) {
      for (int ik{0}; ik < k1[iC]; ++ik){
        double pt_rndm = -999., y_rndm = -999.;
        hPtY[iC]->GetRandom2(pt_rndm, y_rndm);
        y_rndm += 1.8; // beam rapidity
        double p = pt_rndm * std::cosh(y_rndm);
        double pz = std::sqrt(p * p - pt_rndm * pt_rndm);
        double eta = 0.5 * std::log((p + pz) / (p - pz));
        if (pt_rndm < pt_lim[0] || pt_rndm > pt_lim[1] || eta < eta_lim[0] || eta > eta_lim[1] || gRandom->Rndm() > brLambda)
          continue;
        double eff_tmp = 1.; // hEff[iC]->GetBinContent(hEff[iC]->FindBin(pt_rndm, eta));
//        std::cout << eff_tmp << std::endl;
        if (gRandom->Rndm() < eff_tmp) {
          for (int i{0}; i < 4; ++i)
            q_[iC][i] = q_[iC][i] + (1. / std::pow(eff_tmp, i + 1));
        }
      }
    }

    long int iS = static_cast<long int>(gRandom->Uniform(S_MIN, S_MAX));

    n1_a[iS] += (q_[0][0]);
    n1_b[iS] += (q_[1][0]);

    n2_a[iS] += (q_[0][0] * q_[0][0]);
    n2_b[iS] += (q_[1][0] * q_[1][0]);

    ne2_a[iS] += (q_[0][1]);
    ne2_b[iS] += (q_[1][1]);

    n11_ab[iS] += (q_[0][0] * q_[1][0]);

    // full formula
    // 1st order
    q1_1_1[iS] += qa_b_c(1, 1, 1);

    // 2nd order
    q2_1_1[iS] += qa_b_c(2, 1, 1);
    q1_2_1[iS] += qa_b_c(1, 2, 1);
    q1_2_2[iS] += qa_b_c(1, 2, 2);

    // 3rd order
    q3_1_1[iS] += qa_b_c(3, 1, 1);
    q1_1_1_x_q1_2_1[iS] += qa_b_c(1, 1, 1) * qa_b_c(1, 2, 1);
    q1_1_1_x_q1_2_2[iS] += qa_b_c(1, 1, 1) * qa_b_c(1, 2, 2);
    q1_3_1[iS] += qa_b_c(1, 3, 1);
    q1_3_2[iS] += qa_b_c(1, 3, 2);
    q1_3_3[iS] += qa_b_c(1, 3, 3);

    // 4th order
    q4_1_1[iS] += qa_b_c(4, 1, 1);
    q2_1_1_x_q1_2_1[iS] += qa_b_c(2, 1, 1) * qa_b_c(1, 2, 1);
    q2_1_1_x_q1_2_2[iS] += qa_b_c(2, 1, 1) * qa_b_c(1, 2, 2);
    q1_1_1_x_q1_3_1[iS] += qa_b_c(1, 1, 1) * qa_b_c(1, 3, 1);
    q2_2_1[iS] += qa_b_c(2, 2, 1);
    q2_2_2[iS] += qa_b_c(2, 2, 2);
    q1_1_1_x_q1_3_2[iS] += qa_b_c(1, 1, 1) * qa_b_c(1, 3, 2);
    q1_1_1_x_q1_3_3[iS] += qa_b_c(1, 1, 1) * qa_b_c(1, 3, 3);
    q1_2_1_x_q1_2_2[iS] += qa_b_c(1, 2, 1) * qa_b_c(1, 2, 2);
    q1_4_1[iS] += qa_b_c(1, 4, 1);
    q1_4_2[iS] += qa_b_c(1, 4, 2);
    q1_4_3[iS] += qa_b_c(1, 4, 3);
    q1_4_4[iS] += qa_b_c(1, 4, 4);
  }

  const long int nEvS = nEv / (S_MAX - S_MIN);
  for (long int i{S_MIN}; i < S_MAX; ++i){
    n1_a[i] /= nEvS;
    n1_b[i] /= nEvS;

    n2_a[i] /= nEvS;
    n2_b[i] /= nEvS;

    ne2_a[i] /= nEvS;
    ne2_b[i] /= nEvS;

    n11_ab[i] /= nEvS;

    q1_1_1[i] /= nEvS;

    q2_1_1[i] /= nEvS;
    q1_2_1[i] /= nEvS;
    q1_2_2[i] /= nEvS;

    q3_1_1[i] /= nEvS;
    q1_1_1_x_q1_2_1[i] /= nEvS;
    q1_1_1_x_q1_2_2[i] /= nEvS;
    q1_3_1[i] /= nEvS;
    q1_3_2[i] /= nEvS;
    q1_3_3[i] /= nEvS;

    q4_1_1[i] /= nEvS;
    q2_1_1_x_q1_2_1[i] /= nEvS;
    q2_1_1_x_q1_2_2[i] /= nEvS;
    q1_1_1_x_q1_3_1[i] /= nEvS;
    q2_2_1[i] /= nEvS;
    q2_2_2[i] /= nEvS;
    q1_1_1_x_q1_3_2[i] /= nEvS;
    q1_1_1_x_q1_3_3[i] /= nEvS;
    q1_2_1_x_q1_2_2[i] /= nEvS;
    q1_4_1[i] /= nEvS;
    q1_4_2[i] /= nEvS;
    q1_4_3[i] /= nEvS;
    q1_4_4[i] /= nEvS;
  }

  double c2_a[nS];
  double c2_b[nS];
  double c2_net[nS];
  double c11[nS];
  double rho[nS];
  double k1_small[nS];
  double k2sk_small[nS];
  double k2_small[nS];
  double k3_small[nS];
  double k4_small[nS];
  TH1F ha("ha", "ha", 10000, .95, 1.05);
  TH1F hb("hb", "hb", 10000, .95, 1.05);
  TH1F hnet("hnet", "hnet", 10000, .95, 1.05);

  TH1F hk2sknet("hk2sknet", ";#kappa_{2}/Skellam;Entries", 10000, .0, 2.);
  TH1F hk4k2net("hk4k2net", ";#kappa_{4}/#kappa_{2};Entries", 20000, -20., 20.);
  TH1F hk3k1net("hk3k1net", ";#kappa_{3}/#kappa_{1};Entries", 20000, 0., 2.);

  for (long int i{S_MIN}; i < S_MAX; ++i){
    // std::cout << n1_a[i] << "\n";
    c2_a[i] = n2_a[i] - n1_a[i] * n1_a[i] + n1_a[i] - ne2_a[i];
    c2_b[i] = n2_b[i] - n1_b[i] * n1_b[i] + n1_b[i] - ne2_b[i];
    c11[i] = n11_ab[i] - n1_a[i] * n1_b[i];
    c2_net[i] = c2_a[i] + c2_b[i] - 2. * c11[i];

    double Q1_small = q1_1_1[i];
    double Q2_small = q2_1_1[i] + q1_2_1[i] - q1_2_2[i];
    double Q3_small = q3_1_1[i] + (3. * q1_1_1_x_q1_2_1[i]) - (3. * q1_1_1_x_q1_2_2[i]) + q1_3_1[i] - (3. * q1_3_2[i]) + (2. * q1_3_3[i]);
    double Q4_small = q4_1_1[i] + (6. * q2_1_1_x_q1_2_1[i]) - (6. * q2_1_1_x_q1_2_2[i]) + (4. * q1_1_1_x_q1_3_1[i]) + (3. * q2_2_1[i])
                      + (3. * q2_2_2[i]) - (12. * (q1_1_1_x_q1_3_2[i])) + (8. * q1_1_1_x_q1_3_3[i]) - (6. * q1_2_1_x_q1_2_2[i])
                      + (q1_4_1[i]) - (7. * q1_4_2[i]) + (12. * q1_4_3[i]) - (6. * q1_4_4[i]);
    // std::cout << Q3_small << std::endl;

    k1_small[i] = q1_1_1[i];
    k2sk_small[i] = q1_2_1[i];
    k2_small[i] = Q2_small - std::pow(Q1_small, 2);
    k3_small[i] = Q3_small - (3. * Q2_small * Q1_small) + (2. * std::pow(Q1_small, 3));
    k4_small[i] = Q4_small - (4. * Q3_small * Q1_small) - (3. * std::pow(Q2_small, 2)) + (12. * Q2_small * std::pow(Q1_small, 2)) - (6. * std::pow(Q1_small, 4));

    ha.Fill(c2_a[i] / n1_a[i]);
    hb.Fill(c2_b[i] / n1_b[i]);
    hnet.Fill(c2_net[i] / (n1_a[i] + n1_b[i]));

    std::cout << k2_small[i] / k2sk_small[i] << std::endl;

    hk3k1net.Fill(k3_small[i] / k1_small[i]);
    hk4k2net.Fill(k4_small[i] / k2_small[i]);
    hk2sknet.Fill(k2_small[i] / k2sk_small[i]);
  }
  fo.cd();
  ha.Write();
  hb.Write();
  hnet.Write();
  hk2sknet.Write();
  hk4k2net.Write();
  hk3k1net.Write();

  fo.Close();

  return 0;
}
