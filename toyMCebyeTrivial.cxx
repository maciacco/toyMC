#include <Riostream.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <string>

constexpr bool save_tree = false;
// constexpr double nEv = 6.e5; // for central (0-10%) --> 6e10 collisions
constexpr double pt_a[] = {0.4, 3.};
constexpr double pt_b[] = {0.4, 3.};

constexpr double pt_bins_a[] = {0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3.0};
constexpr double pt_bins_b[] = {0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3.0};

constexpr double eff_pt_bins_a[] = {0.009, 0.012, 0.014, 0.016, 0.0175, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018};
constexpr double eff_pt_bins_b[] = {0.009, 0.012, 0.014, 0.016, 0.0175, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018};

constexpr long int nS = 40;

constexpr double massLambda = 1.115683; // GeV/c^2
constexpr double T[]{0.244, 0.339}; // L, Lbar
constexpr double brLambda = 0.639;

constexpr double yield_a = 13.4;
constexpr double yield_b = 0.10;

int main(int const argv, char* const argc[]){ //double const nEv = 6.e10, int const S_MIN = 0, int const S_MAX = 20){
  double nEv = std::stod(argc[1]);
  int S_MIN = std::stoi(argc[2]);
  int S_MAX = std::stoi(argc[3]);

  if (S_MIN >= S_MAX) {
    std::cout << "Fatal: S_MIN >= S_MAX!" << std::endl;
    return 1;
  }

  gRandom->SetSeed(10 * S_MIN);

  TF1* bw_a = new TF1("bw_a", "x * TMath::Exp(-TMath::Sqrt([0] * [0] + x * x) / [1])", 0., 10.);
  TF1* bw_b = new TF1("bw_a", "x * TMath::Exp(-TMath::Sqrt([0] * [0] + x * x) / [1])", 0., 10.);
  bw_a->SetParameter(0, massLambda);
  bw_a->SetParameter(1, T[0]); // L
  bw_b->SetParameter(0, massLambda);
  bw_b->SetParameter(1, T[1]); // Lbar
  double norm_a = yield_a * brLambda / bw_a->Integral(0., 10., 1.e-5);
  double norm_b = yield_b * brLambda / bw_b->Integral(0., 10., 1.e-5);

  double k1_a = norm_a * bw_a->Integral(pt_a[0], pt_a[1], 1.e-5);
  double k1_b = norm_b * bw_b->Integral(pt_b[0], pt_b[1], 1.e-5);
  std::cout << "k1_a = " << k1_a << ", k1_b = " << k1_b << "\n";

  std::string fname{"ex_lambda_"};
  fname = fname + std::to_string(S_MIN) + std::string{".root"};
  TFile fo(fname.data(), "recreate");

  TH1D heff_a("heff_a", ";#it{p}_{T} (GeV/#it{c});Efficiency", 13, pt_bins_a);
  TH1D heff_b("heff_b", ";#it{p}_{T} (GeV/#it{c});Efficiency", 13, pt_bins_b);
  for (int iB{1}; iB < 14; ++iB){
    heff_a.SetBinContent(iB, eff_pt_bins_a[iB - 1]);
    heff_b.SetBinContent(iB, eff_pt_bins_b[iB - 1]);
  }

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
    double a = gRandom->Poisson(k1_a);
    double b = gRandom->Poisson(k1_b);

    double q_p[4]{0.};
    double q_n[4]{0.};

    auto qa_b_c = [&](int const a, int const b, int const c) -> double
    {
      double sgn = (b % 2) == 0 ? 1. : -1.;
      return std::pow(q_p[c - 1] + sgn * q_n[c - 1], a);
    };

    for (long int ia{0}; ia < a; ++ia) {
      double pt_rndm = bw_a->GetRandom(pt_a[0], pt_a[1]);
      double eff_a = 1.; // heff_a.GetBinContent(heff_a.FindBin(pt_rndm));
//       std::cout << eff_a << std::endl;
      if (gRandom->Rndm() < eff_a) {
        for (int i{0}; i < 4; ++i)
          q_p[i] = q_p[i] + (1. / std::pow(eff_a, i + 1));
      }
    }

    for (long int ib{0}; ib < b; ++ib) {
      double pt_rndm = bw_b->GetRandom(pt_b[0], pt_b[1]);
      double eff_b = 1.; // heff_b.GetBinContent(heff_b.FindBin(pt_rndm));
      if (gRandom->Rndm() < eff_b) {
        for (int i{0}; i < 4; ++i)
          q_n[i] = q_n[i] + (1. / std::pow(eff_b, i + 1));
      }
    }

    long int iS = static_cast<long int>(gRandom->Uniform(S_MIN, S_MAX));

    n1_a[iS] += (q_p[0]);
    n1_b[iS] += (q_n[0]);

    n2_a[iS] += (q_p[0] * q_p[0]);
    n2_b[iS] += (q_n[0] * q_n[0]);

    ne2_a[iS] += (q_p[1]);
    ne2_b[iS] += (q_n[1]);

    n11_ab[iS] += (q_p[0] * q_n[0]);

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

  TH1F hk2sknet("hk2sknet", ";#kappa_{2}/Skellam;Entries", 10000, .95, 1.05);
  TH1F hk4k2net("hk4k2net", ";#kappa_{4}/#kappa_{2};Entries", 20000, -300., 300.);
  TH1F hk3k1net("hk3k1net", ";#kappa_{3}/#kappa_{1};Entries", 20000, -10., 10.);

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
