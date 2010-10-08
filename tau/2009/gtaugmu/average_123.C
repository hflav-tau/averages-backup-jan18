{
  double x1 = 0.996564 ; double e_x1 = 0.00299738; // pi
  double x2 = 0.986033 ; double e_x2 = 0.007233;   // K
  double x3 = 1.00105  ; double e_x3 = 0.00207702; // e
  
  double cov_11 = e_x1 * e_x1;
  double cov_22 = e_x2 * e_x2;
  double cov_33 = e_x3 * e_x3;

  double cov_12, cov_21, cov_13, cov_31, cov_23, cov_32;

  cov_12 = cov_21 = (-0.004707 * 0.00236167 * 0.00680697) ; // pi - K
  cov_23 = cov_32 = (0.041556 * 0.00680697 * 0.00112125);  //  e - K
  cov_13 = cov_31 = (0.000961 * 0.00236167 * 0.00112125);  // e - pi

  TMatrixD covmat(3,3);
  covmat[0][0] = cov_11; covmat[0][1] = cov_12 ; covmat[0][2] = cov_13;
  covmat[1][0] = cov_21; covmat[1][1] = cov_22 ; covmat[1][2] = cov_23;
  covmat[2][0] = cov_31; covmat[2][1] = cov_32 ; covmat[2][2] = cov_33;

  TMatrixD inverse = covmat;
  double determinant;
  inverse.Invert(&determinant);
  
  double w_1 = inverse[0][0] + inverse[0][1] + inverse[0][2];
  double w_2 = inverse[1][0] + inverse[1][1] + inverse[1][2];
  double w_3 = inverse[2][0] + inverse[2][1] + inverse[2][2];

  double denominator = w_1 + w_2 + w_3;
  
  w_1/=denominator;
  w_2/=denominator;
  w_3/=denominator;

  double average = w_1 * x1 + w_2 * x2 + w_3 * x3;
  double error = sqrt(1/denominator);

  double pull[3] = { x1 - average, x2 - average, x3 - average};
  
  double chi2 = 0;
  for (int i=0;i<3;++i) {
    for (int j=0;j<3;++j) {
      chi2+= pull[i] * inverse[i][j] * pull[j];
    }
  }

  cout << " average = " << average << " error = " << error << " chi2 = " << chi2 << " prob =  " << TMath::Prob(chi2,2) << endl;
  
}
