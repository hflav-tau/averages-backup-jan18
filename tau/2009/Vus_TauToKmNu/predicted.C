{
 const double Gmu = 1.16637e-5 ; // +- 0.00001e-5 GeV-2
 const double fpi = 132e-3 ; // +- 2e-3 GeV // arXiv:0706.1726 [hep-lat]
 const double fK  = 157e-3 ; // +- 2e-3 GeV // arXiv:0706.1726 [hep-lat]
 const double Vud = 0.97425; // +- 0.00022 // arXiv:0812.1202 [nucl-ex]
 const double Vus = 0.2255; // +- 0.0010 // Unitarity 
 const double mtau = 1.77677; // +- 0.15e-3 GeV
 const double mrho = 775.49e-3; // +- 0.34e-3 GeV
 const double mpi  = 139.57018e-3; // +- 0.00035e-3 GeV
 const double mK   = 493.677e-3; // +- 0.013e-3 GeV
 const double mZ   = 91.1876; // +- 0.0021 GeV
 const double alpha = 1./133.3;
 const double pi    = acos(-1.);
 const double TauLifeTime = 290.6e-15; // +-1e-15 s
 const double hbar =  6.58211899e-25; // +- 0.00000016e-25; # GeV s 
 //
 double SEW_Old = (1 + ((2*alpha)/pi)*(log(mZ/mtau))); 
 // cout << "SEW_Old = " << SEW_Old << endl; // SEW_Old = 1.01881
 double SEW = 1.0201;// +- 0.0003
 //
 double BR_TauToKNu = (1./(16.*pi*hbar)) * pow(Gmu,2) * pow(fK,2) * pow(Vus,2) * pow(mtau,3) * TauLifeTime*  pow(1 - (pow(mK,2))/(pow(mtau,2)),2) * SEW;
 cout << "Vus = 0.2255 => BR_TauToKNu = " << BR_TauToKNu*100. << " % " << endl;
 //
 BR_TauToKNu = 0.6969439e-2; // +- 0.00975197e-2; # measured average
 double Vus = sqrt(
		   BR_TauToKNu * (16*pi*hbar)* 1./( pow(Gmu,2) * pow(fK,2) * pow(mtau,3) * TauLifeTime*  pow(1 - (pow(mK,2))/(pow(mtau,2)),2) * SEW)
                    );
 cout << "BR_TauToKNu = " << BR_TauToKNu*100. << " % =>  Vus = " <<  Vus << endl;
}
