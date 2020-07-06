void inc_mlb_mc()
{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Thu Mar  5 20:31:01 2020) by ROOT version 6.14/09
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",514,1154,538,323);
   Canvas_1->Range(-31.25,-0.000341009,281.25,0.003069081);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   TH1F *inc_mlb__1 = new TH1F("inc_mlb__1","",20,0,250);
   inc_mlb__1->SetBinContent(2,3.018237e-05);
   inc_mlb__1->SetBinContent(3,0.0002789073);
   inc_mlb__1->SetBinContent(4,0.0006226255);
   inc_mlb__1->SetBinContent(5,0.001006529);
   inc_mlb__1->SetBinContent(6,0.001526482);
   inc_mlb__1->SetBinContent(7,0.002084685);
   inc_mlb__1->SetBinContent(8,0.00247323);
   inc_mlb__1->SetBinContent(9,0.002585545);
   inc_mlb__1->SetBinContent(10,0.002458502);
   inc_mlb__1->SetBinContent(11,0.001956044);
   inc_mlb__1->SetBinContent(12,0.001373706);
   inc_mlb__1->SetBinContent(13,0.0009230652);
   inc_mlb__1->SetBinContent(14,0.0006954003);
   inc_mlb__1->SetBinContent(15,0.0005750366);
   inc_mlb__1->SetBinContent(16,0.0004803006);
   inc_mlb__1->SetBinContent(17,0.0004101323);
   inc_mlb__1->SetBinContent(18,0.0003382284);
   inc_mlb__1->SetBinContent(19,0.0002856746);
   inc_mlb__1->SetBinContent(20,0.0002331041);
   inc_mlb__1->SetBinContent(21,0.00128284);
   inc_mlb__1->SetBinError(2,1.352739e-06);
   inc_mlb__1->SetBinError(3,4.142026e-06);
   inc_mlb__1->SetBinError(4,6.120605e-06);
   inc_mlb__1->SetBinError(5,7.823249e-06);
   inc_mlb__1->SetBinError(6,9.636042e-06);
   inc_mlb__1->SetBinError(7,1.122243e-05);
   inc_mlb__1->SetBinError(8,1.228758e-05);
   inc_mlb__1->SetBinError(9,1.261874e-05);
   inc_mlb__1->SetBinError(10,1.218588e-05);
   inc_mlb__1->SetBinError(11,1.097791e-05);
   inc_mlb__1->SetBinError(12,9.245799e-06);
   inc_mlb__1->SetBinError(13,7.599891e-06);
   inc_mlb__1->SetBinError(14,6.615801e-06);
   inc_mlb__1->SetBinError(15,5.987244e-06);
   inc_mlb__1->SetBinError(16,5.492647e-06);
   inc_mlb__1->SetBinError(17,5.07204e-06);
   inc_mlb__1->SetBinError(18,4.618878e-06);
   inc_mlb__1->SetBinError(19,4.267263e-06);
   inc_mlb__1->SetBinError(20,3.891169e-06);
   inc_mlb__1->SetBinError(21,9.196499e-06);
   inc_mlb__1->SetEntries(1832727);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_LaTex = ptstats->AddText("inc_mlb");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 1832727");
   ptstats_LaTex = ptstats->AddText("Mean  =  116.1");
   ptstats_LaTex = ptstats->AddText("Std Dev   =  44.98");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   inc_mlb__1->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(inc_mlb__1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   inc_mlb__1->SetLineColor(ci);
   inc_mlb__1->GetXaxis()->SetTitle("m(l,b) [GeV]");
   inc_mlb__1->GetXaxis()->SetLabelFont(42);
   inc_mlb__1->GetXaxis()->SetLabelSize(0.035);
   inc_mlb__1->GetXaxis()->SetTitleSize(0.035);
   inc_mlb__1->GetXaxis()->SetTitleFont(42);
   inc_mlb__1->GetYaxis()->SetTitle("Events");
   inc_mlb__1->GetYaxis()->SetLabelFont(42);
   inc_mlb__1->GetYaxis()->SetLabelSize(0.035);
   inc_mlb__1->GetYaxis()->SetTitleSize(0.035);
   inc_mlb__1->GetYaxis()->SetTitleOffset(0);
   inc_mlb__1->GetYaxis()->SetTitleFont(42);
   inc_mlb__1->GetZaxis()->SetLabelFont(42);
   inc_mlb__1->GetZaxis()->SetLabelSize(0.035);
   inc_mlb__1->GetZaxis()->SetTitleSize(0.035);
   inc_mlb__1->GetZaxis()->SetTitleFont(42);
   inc_mlb__1->Draw("");
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
