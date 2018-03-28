#include "junkInterface.h"
#include "mclimit_csm.C"

// Class for individually create templates that
// can be added to a csm_model
// Diego Martinez Santos. Universidade de Santiago de Compostela
// Diego.Martinez.Santos@cern.ch
//
// Memory Leak fix by Swagato Banerjee 
//
Int_t k;

//creator
csm_template::csm_template(Int_t np)
{
  lowshape = new TH1*[np];
  highshape = new TH1*[np];
  lowsigma = new Double_t[np];
  highsigma = new Double_t[np];
  nps_low = new Double_t[np];
  nps_high = new Double_t[np];
  current_np = 0;
  npnames = new char*[np];
  nnp = np;

}

//destructor
csm_template::~csm_template()
{
  /*
  for(k=0;k<=nnp;k++){
    delete lowshape[k];
    delete highshape[k];
    delete npnames[k];
  }
  */
  delete [] lowshape;
  delete [] highshape;
  delete lowsigma;
  delete highsigma;
  delete nps_low;
  delete nps_high;
  delete [] npnames;
  //  delete template_hist;
  delete chaname;
}
// Addition of one nuisance parameter and its properties...
void csm_template::set_np(char *npname, Double_t low, Double_t high, TH1 *h_high, Double_t sigma_high, TH1 *h_low, Double_t sigma_low)
//  void csm_template::set_np(char *npname, Double_t high, Double_t low, TH1 *h_high, Double_t sigma_high, TH1 *h_low, Double_t sigma_low)
{
  char *s;
  lowshape[current_np] = h_low;
  highshape[current_np] = h_high;
  lowsigma[current_np] = sigma_low;
  highsigma[current_np] = sigma_high;
  nps_low[current_np] = low;
  nps_high[current_np] = high;
  s = new char[strlen(npname)+1];
  strcpy(s,npname);

  npnames[current_np] = s;

  current_np++;
}
// Setting the template hisogram, the normalization, poisson and signal flags and the channel name
void csm_template::make(TH1 *templ, Double_t sf_, Int_t pflag_, Int_t sflag_, char *chaname_)
{
  char *s2;
  sf = sf_;
  pflag = pflag_;
  sflag = sflag_;
  s2 = new char[strlen(chaname_)+1];
  strcpy(s2,chaname_);
  chaname = s2;
  template_hist = templ;
}
// Adding the template to a given csm_model
void csm_template::add_to(csm_model *themodel)
{
  cout << "sf = " << sf << " nnp = " << nnp << " pflag = " << pflag << " sflag = " << sflag << endl;
  if (nnp==1) {
    cout << npnames[0] << " " << nps_low[0] << " " << nps_high[0] << " "
	 << lowsigma[0] << " " << highsigma[0] << endl;
  }
  themodel->add_template(template_hist, sf, nnp, npnames, nps_low, nps_high, lowshape, lowsigma, highshape, highsigma, pflag, sflag, chaname);
}
