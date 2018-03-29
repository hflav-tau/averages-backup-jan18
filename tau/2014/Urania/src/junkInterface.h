# include "TH1.h"
# include "mclimit_csm.h"
//# include "mclimit_csm.C"

// Class for individually create templates that
// can be added to a csm_model
// Diego Martinez Santos. Universidade de Santiago de Compostela
// Diego.Martinez.Santos@cern.ch

class csm_template
{
 public:
  csm_template(Int_t);
  ~csm_template();
  void set_np(char *, Double_t, Double_t, TH1 *, Double_t, TH1 *, Double_t);
  // set_np: Addition of one nuisance parameter and its properties..

  void make(TH1 *, Double_t, Int_t, Int_t, char *);
  //make: Sets the template hisogram, the scale fraction, poisson and signal flags and the channel name

  void add_to(csm_model *);
  //add_to:  Adds itself to a given csm_model

 private:
  Int_t current_np;
  char **npnames;
  TH1 **lowshape;
  TH1 **highshape;
  TH1 *template_hist; 
  Double_t *lowsigma;
  Double_t *highsigma;
  Double_t *nps_low;
  Double_t *nps_high;
  Double_t sf;
  Int_t sflag;
  Int_t pflag;
  Int_t nnp;
  char *chaname;
};
