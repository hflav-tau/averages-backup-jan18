{
  // gtauge^2  = B_tautomu  * (tau_mu/tau_tau) * (m_mu/m_tau)^5 * (f(m^2_e/m^2_mu)  / f(m^2_mu/m^2_tau)) * (delta^mu / delta^tau)

  // B_tautomu  = gtauge^2  * (tau_tau/tau_mu) * (m_tau/m_mu)^5 * (f(m^2_mu/m^2_tau))/ (f(m^2_e/m^2_mu)) * (delta^tau/ delta_mu )
  //
  //--- from PDG 2010
  //
  //    	MOHR	 08 RVUE	 	2006 CODATA value
  //  MOHR 2008	    	Reviews of Modern Physics 80 (2008) 633 
  // Published in Rev.Mod.Phys.80:633-730,2008.  e-Print: arXiv:0801.0028 [physics.atom-ph]
  // CODATA Recomended Values of the Fundamental Physical Constants: 2006
  // http://www.slac.stanford.edu/spires/find/hep/www?eprint=arXiv:0801.0028
  //
  double m_e = 	0.510998910; double e_m_e = 0.000000013 ; // MeV
  double m_mu = 105.6583668 ; double e_m_mu = 0.0000038 ; // MeV
  //
  double tau_tau = 290.6e-15; double e_tau_tau = 1.0e-15;
  
  //--- from HFAG 2009
  double m_tau =  1776.7673082; double e_m_tau = 0.1507259; // MeV

  
  //--- from PDG 2010
  double m_W = 80.399*1e3 ; double e_m_W =  0.023*1e3; // MeV

  double tau_mu = 2.197034e-6; double e_tau_mu = 0.000021e-6;

  
  // rad. corrections from to get Be from tau lifetime
  // values from 10.1103/RevModPhys.78.1043 p.1047, arXiv:hep-ph/0507078v2 p.7, could be recomputed
  // - delta^L_gamma = 1 + alpha(mL)/2pi * (25/4 - pi^2)
  // - delta^L_W = 1 + 3/5* m_L^2/M_W^2
  //
  double delta_tau_gamma = 1 - 43.2e-4;
  double delta_mu_gamma = 1 - 42.4e-4;
  double delta_tau_W = 1 + 2.9e-4;
  double delta_mu_W = 1 + 1.0e-6;
  
  double delta_tau = delta_tau_W * delta_tau_gamma;
  double delta_mu = delta_mu_W * delta_mu_gamma;

  //
  // Bmu from tau lifetime
  // Bmu = tau_tau / tau_mu (m_tau/m_mu)^5 f(m^2_mu/m^2_tau)/f(m^2_e/m^2_mu) (delta^tau_gamma delta^tau_W)/(delta^mu_gamma delta^mu_W)
  //

  double mmubymtau = (m_mu * m_mu) / (m_tau * m_tau); 
  double phspf_mmubymtau = 
    1 
    - 8 * mmubymtau 
    + 8 * mmubymtau * mmubymtau * mmubymtau 
    - mmubymtau * mmubymtau * mmubymtau * mmubymtau
    - 12 * mmubymtau * mmubymtau * log(mmubymtau);
  
  double mebymmu = (m_e * m_e) / (m_mu * m_mu); 
  double phspf_mebymmu = 
    1 
    - 8 * mebymmu
    + 8 * mebymmu * mebymmu * mebymmu
    - mebymmu * mebymmu * mebymmu * mebymmu
    - 12 * mebymmu * mebymmu * log(mebymmu);

  double Bmu_from_taulife = (tau_tau/tau_mu) * TMath::Power((m_tau/m_mu),5) * (phspf_mmubymtau/phspf_mebymmu) * (delta_tau/delta_mu);
  cout << "Bmu_from_taulife = " << Bmu_from_taulife*100 << "%" << endl;

  double Bmu_measured = 1.7407e-01; double e_Bmu_measured = 3.8506e-04;
  double gtauge_sq = Bmu_measured * (tau_mu/tau_tau) * TMath::Power((m_mu/m_tau),5) * (phspf_mebymmu/phspf_mmubymtau) * (delta_mu/delta_tau);
  double gtauge = TMath::Sqrt(gtauge_sq);
  cout << "gtauge = " << gtauge << endl;
}
