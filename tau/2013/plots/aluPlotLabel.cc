//
// plot label with title and subtitle
//
void aluPlotLabel(const TString& title="",
		  const TString& subtitle="",
		  Double_t xpos=0, Double_t ypos=0,
		  Double_t scale=1);

void aluPlotLabel(const TString& title,
		  const TString& subtitle,
		  Double_t xpos, Double_t ypos,
		  Double_t scale)
{
  TVirtualPad* thePad;
  
  if ((thePad = TVirtualPad::Pad()) == 0) return;
  
  const Float_t font_size_px(20*scale);
  const Float_t label_height_px(1.25 * font_size_px);
  const Float_t fsratio(1.2);

  Float_t x1, x2, y1, y2;

  //--- get pad dimensions
  UInt_t pad_width(thePad->XtoPixel(thePad->GetX2()));
  UInt_t pad_height(thePad->YtoPixel(thePad->GetY1()));
  
  //--- title label height in NDC coordinates
  Float_t label_height(Double_t(label_height_px)/Double_t(pad_height));

  Float_t label_width(max(Double_t(title.Length())*1.2, Double_t(subtitle.Length()))
		      /Double_t(2.5)
		      *label_height
		      *Double_t(pad_height)/Double_t(pad_width)
		      );

  if (xpos >= 0) {
    x1 = xpos;
    x2 = xpos + label_width;
  } else {
    x1 = 1 + xpos - label_width;
    x2 = 1 + xpos;
  }
  
  if (ypos >= 0) {
    y1 = ypos+0.9*label_height;
    y2 = ypos+0.9*label_height + label_height;
  } else {
    y1 = 1 + ypos - label_height;
    y2 = 1 + ypos;
  }
  
  if (title != "") {
  TPaveText *tbox1 = new TPaveText(x1, y1, x2, y2, "BRNDC");
  // tbox1->SetLineColor(1);
  // tbox1->SetLineStyle(1);
  // tbox1->SetLineWidth(2);
  tbox1->SetFillColor(kBlack);
  tbox1->SetFillStyle(1001);
  // tbox1->SetBorderSize(1);
  tbox1->SetShadowColor(kWhite);
  tbox1->SetTextColor(kWhite);
  tbox1->SetTextFont(76);
  tbox1->SetTextSize(20*scale);
  tbox1->SetTextAlign(22); //center-adjusted and vertically centered
  tbox1->AddText(title);
  tbox1->Draw();
  }

  if (subtitle != "") {
  TPaveText *tbox2 = new TPaveText(x1, y1-0.9*label_height, x2, y2-label_height, "BRNDC");
  // tbox2->SetLineColor(1);
  // tbox2->SetLineStyle(1);
  // tbox2->SetLineWidth(2);
  tbox2->SetFillColor(kWhite);
  tbox2->SetFillStyle(1001);
  // tbox2->SetBorderSize(1);
  tbox2->SetShadowColor(kWhite);
  tbox2->SetTextColor(kBlack);
  tbox2->SetTextFont(76);
  tbox2->SetTextSize(18*scale);
  tbox2->SetTextAlign(22); //center-adjusted and vertically centered
  tbox2->AddText(subtitle);
  tbox2->Draw();
  }
  
  return;
}
