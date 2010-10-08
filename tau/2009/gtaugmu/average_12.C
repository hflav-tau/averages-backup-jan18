{
  double x1 = 0.996564 ; double e_x1 = 0.00299738; // pi
  double x2 = 0.986033 ; double e_x2 = 0.007233;   // K
  
  double cov_11 = e_x1 * e_x1;
  double cov_22 = e_x2 * e_x2;


  double cov_12, cov_21; 

  cov_12 = cov_21 = (-0.004707 * 0.00236167 * 0.00680697) ; // pi - K

  TMatrixD covmat(2,2);
  covmat[0][0] = cov_11; covmat[0][1] = cov_12 ; 
  covmat[1][0] = cov_21; covmat[1][1] = cov_22 ; 


  TMatrixD inverse = covmat;
  double determinant;
  inverse.Invert(&determinant);
  
  double w_1 = inverse[0][0] + inverse[0][1] ;
  double w_2 = inverse[1][0] + inverse[1][1] ;

  double denominator = w_1 + w_2 ;
  
  w_1/=denominator;
  w_2/=denominator;


  double average = w_1 * x1 + w_2 * x2 ;
  double error = sqrt(1/denominator);

  double pull[2] = { x1 - average, x2 - average};
  
  double chi2 = 0;
  for (int i=0;i<2;++i) {
    for (int j=0;j<2;++j) {
      chi2+= pull[i] * inverse[i][j] * pull[j];
    }
  }

  cout << " average = " << average << " error = " << error << " chi2 = " << chi2 << " prob =  " << TMath::Prob(chi2,1) << endl;
  
}
