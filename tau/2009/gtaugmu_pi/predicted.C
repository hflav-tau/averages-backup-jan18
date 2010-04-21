{
// http://arxiv.org/pdf/0709.0282v1 Eqn 4.1

  const double M_Mu = 0.1056584;

  const double M_Pi = 0.1395702;
  const double Tau_Pi = 2.6033e-8; // +- 0.0005e-8 s
  const double BR_PiToMuNu = 99.98770e-2; // +- 0.00004e-2

  const double M_Tau = 1.77677; // +- 0.15e-3 GeV
  const double Tau_Tau = 290.6e-15; // +- 1.0e-15 s

  const double Delta_Pi = 1.0016; // +-0.0014

  double BR_TauToPiNu = BR_PiToMuNu 
                      * Delta_Pi 
                      * pow(M_Tau,3) * Tau_Tau * 1./(2*M_Pi * pow(M_Mu,2)* Tau_Pi) 
                      * pow((1-pow(M_Pi/M_Tau,2))/(1-pow(M_Mu/M_Pi,2)),2);

  cout << "BR_TauToPiNu = " << BR_TauToPiNu*100 << " %" << endl;

}
