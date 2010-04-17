{

  const double elecmass=0.000510999;
  const double muonmass=0.1056584;
  const double taumass=1.77675;
 
  double x1 = (elecmass * elecmass) / (taumass * taumass);
  double fx1 = 1 
             - 8 * x1 
             + 8 * x1 * x1 * x1 
             - x1 * x1 * x1 * x1
             - 12 * x1 * x1 * log(x1);

  cout << x1 << " " << fx1 <<  endl;

  double x2 = (muonmass * muonmass) / (taumass * taumass);
  double fx2 = 1 
             - 8 * x2 
             + 8 * x2 * x2 * x2 
             - x2 * x2 * x2 * x2
             - 12 * x2 * x2 * log(x2);
  
    cout << x2 << " " << fx2 <<  endl;

    cout << fx1/fx2 << " " << fx2/fx1 << endl;

}
