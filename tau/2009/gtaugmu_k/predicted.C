{
// http://arxiv.org/pdf/0709.0282v1 Eqn 4.1

  const double M_Mu = 0.1056584;

  const double M_K = 493.677e-3; // +- 0.013e-3 GeV
  const double Tau_K = 1.2379e-8; // +- 0.0021e-8 s
  const double BR_KToMuNu = 63.60e-2;// +- 0.16e-2

  const double M_Tau = 1.77677; // +- 0.15e-3 GeV
  const double Tau_Tau = 290.6e-15; // +- 1.0e-15 s

  const double Delta_K = 1.0090; // +-0.0022

  double BR_TauToKNu = BR_KToMuNu 
                     * Delta_K
                     * pow(M_Tau,3) * Tau_Tau * 1./(2*M_K * pow(M_Mu,2)* Tau_K) 
                     * pow((1-pow(M_K/M_Tau,2))/(1-pow(M_Mu/M_K,2)),2);

  cout << "BR_TauToKNu = " << BR_TauToKNu*100 << " %" << endl;

}
