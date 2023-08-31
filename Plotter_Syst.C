#include "TMath.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TF1.h"
#include "THashList.h"
#include "TList.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDatime.h"

using namespace std;

void Plotter_Syst(const char *filename = "AE_pp13TeV_Charmonia_Ideal_07_New.root")
{
//______________________
    //Variables declaration
    int entries_number = 0;
    int entries_number2 = 0;
    Double_t ptbins[19] = {0.,0.5,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,18.,20.,22.,26.,30.};
    //Double_t ptbins[16] = {0.,0.5,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,18.,20.};
    Double_t ptbins2[14] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,20.};
    Double_t ptbins3[9] = {0.,1.,2.,3.,4.,5.,7.,10.,20.};
    
    vector<double> ptCenter = {0.5,1.5,2.5,3.5,4.5,6.,8.5,15};
    vector<double> eptCenter = {0.5,0.5,0.5,0.5,0.5,1.,1.5,5};
    vector<double> yCenter = {2.625,2.875,3.125,3.375,3.625,3.875};
    vector<double> eyCenter = {0.125,0.125,0.125,0.125,0.125,0.125};
    //Double_t ybins[7] = {2.5,2.75,3.,3.25,3.5,3.75,4.};
    Double_t ybins[7] = {-4.,-3.75,-3.5,-3.25,-3.,-2.75,-2.5};
    
    TH1F *hTime   = new TH1F("hTime","Time distribution",700000,-200000,500000);
    TH1F *hTimeMC   = new TH1F("hTimeMC","MC time distribution",700000,-200000,500000);
    TH1F *hTimeEmbed   = new TH1F("hTimeEmbed","Embed time distribution",700000,-200000,500000);
    
    TH1F *hTimou, *hTimou1, *hTimouAll, *hTimouAll_, *hRec, *hSim, *hTimouY, *hRecY, *hSimY, *hTimou2, *hRec2, *hSim2, *hTimouY2, *hRecY2, *hSimY2, *hTimouAllY2, *hTimouAllY2_;
    
    
    vector<vector<double>> AxEVecJPsi(8, vector<double>());
    vector<vector<double>> AxEVecPsi2S(8, vector<double>());
    vector<vector<double>> AxEVecJPsiError(8, vector<double>());
    vector<vector<double>> AxEVecPsi2SError(8, vector<double>());
    vector<vector<double>> AxEVecRatio(8, vector<double>());
    vector<vector<double>> AxEVecRatioError(8, vector<double>());
    
    vector<vector<double>> AxEVecJPsiY(6, vector<double>());
    vector<vector<double>> AxEVecPsi2SY(6 , vector<double>());
    vector<vector<double>> AxEVecJPsiErrorY(6, vector<double>());
    vector<vector<double>> AxEVecPsi2SErrorY(6 , vector<double>());
    vector<vector<double>> AxEVecRatioY(6, vector<double>());
    vector<vector<double>> AxEVecRatioErrorY(6, vector<double>());
    
    vector<vector<vector<double>>> AxEJPsiTen(8, vector<vector<double>>(4, vector<double>(14, 0.0)));
    vector<vector<vector<double>>> AxEPsi2STen(8, vector<vector<double>>(4, vector<double>(12, 0.0)));
    vector<vector<vector<double>>> AxEJPsiErrorTen(8, vector<vector<double>>(4, vector<double>(14, 0.0)));
    vector<vector<vector<double>>> AxEPsi2SErrorTen(8, vector<vector<double>>(4, vector<double>(12, 0.0)));
    
    vector<double> SystJPsiVec, SystJPsiVecY, SystPsi2SVec, SystPsi2SVecY, SystRatioVec, SystRatioVecY, diffJPsiVec, diffJPsiVecY, diffPsi2SVec, diffPsi2SVecY, diffRatioVec, diffRatioVecY, xVec, xVecY, xVecRatio, xVecRatioY, AxeIntJPsi, AxeIntPsi2S, AxeIntRatio, diffJPsiIntVec, diffPsi2SIntVec, diffRatioIntVec, SystJPsiIntVec, SystPsi2SIntVec, SystRatioIntVec, AxeIntJPsiError, AxeIntPsi2SError, AxeIntRatioError;
    
  TGraphErrors* scatterPlotPtJPsi = new TGraphErrors();
  TGraphErrors* scatterPlotPtPsi2S = new TGraphErrors();
  TGraphErrors* scatterPlotYJPsi = new TGraphErrors();
  TGraphErrors* scatterPlotYPsi2S = new TGraphErrors();
  scatterPlotPtJPsi->GetXaxis()->SetBit(TAxis::kLabelsHori);
  scatterPlotPtPsi2S->GetXaxis()->SetBit(TAxis::kLabelsHori);
  scatterPlotYJPsi->GetXaxis()->SetBit(TAxis::kLabelsHori);
  scatterPlotYPsi2S->GetXaxis()->SetBit(TAxis::kLabelsHori);

TCanvas * cJPsiAlone = new TCanvas("cJPsiAlone","Acceptance efficiency vs Pt (J/Psi)");
TCanvas * cPsi2SAlone = new TCanvas("cPsi2SAlone","Acceptance efficiency vs Pt (Psi2S)");
TCanvas * cJPsiAloneY = new TCanvas("cJPsiAloneY","Acceptance efficiency vs Y (J/Psi)");
TCanvas * cPsi2SAloneY = new TCanvas("cPsi2SAloneY","Acceptance efficiency vs Y (Psi2S)");
//TCanvas * cRatioAlone = new TCanvas("cRatioAlone","Ratio of acceptance efficiency vs Pt ");
    TLegend* legendJPsi = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend* legendPsi2S = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend* legendJPsiY = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend* legendPsi2SY = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend* legendRatioAlone = new TLegend(0.7, 0.7, 0.9, 0.9);
for (int i=0;i<4;i++){
	for (int i2=0;i2<15;i2++){
	//int i2=0;

    // J/Psi extraction
    TFile *file = TFile::Open(filename,"read");
    if(!file){cout << "Error :: The data file is not found, please check"<< endl; return;}
    hTimouAll = (TH1F*) file->Get(("Jpsi/fhRecPt"+to_string(i)+"-"+to_string(i2)).c_str());
    hTimouAll_ = (TH1F*) file->Get(("Jpsi/fhRecPt"+to_string(i)+"-"+to_string(i2)).c_str());
    hTimou = (TH1F*) file->Get(("Jpsi/fhRecPt"+to_string(i)+"-"+to_string(i2)).c_str());
    hRec = (TH1F*) file->Get(("Jpsi/fhRecPt"+to_string(i)+"-"+to_string(i2)).c_str());
    hSim = (TH1F*) file->Get(("Jpsi/fhSimPt"+to_string(i)+"-"+to_string(i2)).c_str());
    hTimouY = (TH1F*) file->Get(("Jpsi/fhRecY"+to_string(i)+"-"+to_string(i2)).c_str());
    hRecY = (TH1F*) file->Get(("Jpsi/fhRecY"+to_string(i)+"-"+to_string(i2)).c_str());
    hSimY = (TH1F*) file->Get(("Jpsi/fhSimY"+to_string(i)+"-"+to_string(i2)).c_str());
    if(!hTimou){cout << "Error :: The Histogram 1 is not found, please check"<< endl; return;}
    
    
    for (int j=0;j<4;j++){
    for (int j2=0;j2<13;j2++){
    //int j2 = 1;
    // Psi(2S) extraction
    hTimou2 = (TH1F*) file->Get(("PsiPrime/fhRecPt"+to_string(j)+"-"+to_string(j2)).c_str());
    hRec2 = (TH1F*) file->Get(("PsiPrime/fhRecPt"+to_string(j)+"-"+to_string(j2)).c_str());
    hSim2 = (TH1F*) file->Get(("PsiPrime/fhSimPt"+to_string(j)+"-"+to_string(j2)).c_str());
    hTimouAllY2 = (TH1F*) file->Get(("PsiPrime/fhRecY"+to_string(j)+"-"+to_string(j2)).c_str());
    hTimouAllY2_ = (TH1F*) file->Get(("PsiPrime/fhRecY"+to_string(j)+"-"+to_string(j2)).c_str());
    hTimouY2 = (TH1F*) file->Get(("PsiPrime/fhRecY"+to_string(j)+"-"+to_string(j2)).c_str());
    hRecY2 = (TH1F*) file->Get(("PsiPrime/fhRecY"+to_string(j)+"-"+to_string(j2)).c_str());
    hSimY2 = (TH1F*) file->Get(("PsiPrime/fhSimY"+to_string(j)+"-"+to_string(j2)).c_str());
    if(!hTimou2){cout << "Error :: The Histogram 2 is not found, please check"<< endl; return;}
    
    
    TH1F *hRecnew = dynamic_cast<TH1F*>(hRec->Rebin(8,"",ptbins3));
    TH1F *hSimnew = dynamic_cast<TH1F*>(hSim->Rebin(8,"",ptbins3));
    TH1F *hRecnew2 = dynamic_cast<TH1F*>(hRec2->Rebin(8,"",ptbins3));
    TH1F *hSimnew2 = dynamic_cast<TH1F*>(hSim2->Rebin(8,"",ptbins3));
    hTimou = hRecnew;
    hTimou2 = hRecnew2;
    hTimou->Divide(hRecnew,hSimnew, 1.0, 1.0, "B");
    hTimou2->Divide(hRecnew2,hSimnew2, 1.0, 1.0, "B");
    hTimou->SetLineColor(kBlack);
    hTimou->SetMarkerColor(kBlack);
    hTimou2->SetLineColor(kBlack);
    hTimou2->SetMarkerColor(kBlack);
    hTimou->SetTitle("Acceptance efficiency vs transverse momentum;p_{T} (GeV/c); acc #times eff");
    hTimou2->SetTitle("Acceptance efficiency vs transverse momentum;p_{T} (GeV/c); acc #times eff");
    
	cJPsiAlone->cd();   
    hTimou->SetLineWidth(2);
    if (j==0 and j2==0){
    Double_t sum_Rec = hRec->Integral();
    Double_t sum_Sim = hSim->Integral();
    Double_t err_Rec;
    Double_t err_Sim;
    hRec->IntegralAndError(1, -1, err_Rec);
    hSim->IntegralAndError(1, -1, err_Sim);
    AxeIntJPsi.push_back(sum_Rec/sum_Sim);
    AxeIntJPsiError.push_back(sum_Rec/sum_Sim*sqrt((err_Rec/sum_Rec)*(err_Rec/sum_Rec)+(err_Sim/sum_Sim)*(err_Sim/sum_Sim)));
    //cout << i << " , " << i2 << " : AxE integrated J/Psi : " << sum_Rec/sum_Sim << endl;
    
    if (i==0) {hTimou->SetLineColor(kBlack);
    		if (i2==0){legendJPsi->AddEntry(hTimou, "Original");}
	hTimou->GetYaxis()->SetRangeUser(0.2, 1.0);
    hTimou->Draw();}
    else {
    if (i==1) {hTimou->SetLineColor(kBlue);
    		if (i2==0){legendJPsi->AddEntry(hTimou, "2.5 < y < 3.0");}}
    if (i==2) {hTimou->SetLineColor(kRed);
    		if (i2==0){legendJPsi->AddEntry(hTimou, "3.0 < y < 3.5");}}
    if (i==3) {hTimou->SetLineColor(kOrange);
    		if (i2==0){legendJPsi->AddEntry(hTimou, "3.5 < y < 4.0");}
    		legendJPsi->Draw();}
    hTimou->Draw("same");}
    }
    
	cPsi2SAlone->cd();   
    hTimou2->SetLineWidth(2);
    if (i==0 and i2==0){
    Double_t sum_Rec = hRec2->Integral();
    Double_t sum_Sim = hSim2->Integral();
    Double_t err_Rec;
    Double_t err_Sim;
    hRec2->IntegralAndError(1, -1, err_Rec);
    hSim2->IntegralAndError(1, -1, err_Sim);
    AxeIntPsi2S.push_back(sum_Rec/sum_Sim);
    AxeIntPsi2SError.push_back(sum_Rec/sum_Sim*sqrt((err_Rec/sum_Rec)*(err_Rec/sum_Rec)+(err_Sim/sum_Sim)*(err_Sim/sum_Sim)));
    cout << (err_Rec/sum_Rec) << " , " << (err_Sim/sum_Sim) << " : " << sum_Rec/sum_Sim*sqrt((err_Rec/sum_Rec)*(err_Rec/sum_Rec)+(err_Sim/sum_Sim)*(err_Sim/sum_Sim)) << endl;
    //cout << j << " , " << j2 << " : AxE integrated Psi(2S) : " << sum_Rec/sum_Sim << endl;
    if (j==0) {hTimou2->SetLineColor(kBlack);
    		if (j2==0){legendPsi2S->AddEntry(hTimou2, "Original");}
	hTimou2->GetYaxis()->SetRangeUser(0.2, 1.0);
    hTimou2->Draw();}
    else {
    if (j==1) {hTimou2->SetLineColor(kBlue);
    		if (j2==0){legendPsi2S->AddEntry(hTimou2, "2.5 < y < 3.0");}}
    if (j==2) {hTimou2->SetLineColor(kRed);
    		if (j2==0){legendPsi2S->AddEntry(hTimou2, "3.0 < y < 3.5");}}
    if (j==3) {hTimou2->SetLineColor(kOrange);
    		if (j2==0){legendPsi2S->AddEntry(hTimou2, "3.5 < y < 4.0");}
    		legendPsi2S->Draw();}
    hTimou2->Draw("same");}
    }
    
    
    TH1F *hRecnewY = dynamic_cast<TH1F*>(hRecY->Rebin(6,"",ybins));
    TH1F *hSimnewY = dynamic_cast<TH1F*>(hSimY->Rebin(6,"",ybins));
    TH1F *hRecnewY2 = dynamic_cast<TH1F*>(hRecY2->Rebin(6,"",ybins));
    TH1F *hSimnewY2 = dynamic_cast<TH1F*>(hSimY2->Rebin(6,"",ybins));
    hTimouY = hRecnewY;
    hTimouY2 = hRecnewY2;
    hTimouY->Divide(hRecnewY,hSimnewY, 1.0, 1.0, "B");
    hTimouY2->Divide(hRecnewY2,hSimnewY2, 1.0, 1.0, "B");
    hTimouY->SetLineColor(kBlack);
    hTimouY->SetMarkerColor(kBlack);
    hTimouY2->SetLineColor(kBlack);
    hTimouY2->SetMarkerColor(kBlack);
    hTimouY->SetTitle("Acceptance efficiency vs rapidity;y; acc #times eff");
    hTimouY2->SetTitle("Acceptance efficiency vs rapidity;y; acc #times eff");
    
    
    
    cJPsiAloneY->cd();   
    hTimouY->SetLineWidth(2);
    if (j==0 and j2==0){
    if (i==0) {hTimouY->SetLineColor(kBlack);
    		if (i2==0){legendJPsiY->AddEntry(hTimouY, "Original");}
	hTimouY->GetYaxis()->SetRangeUser(0., 0.6);
    hTimouY->Draw();}
    else {
    if (i==1) {hTimouY->SetLineColor(kBlue);
    		if (i2==0){legendJPsiY->AddEntry(hTimouY, "2.5 < y < 3.0");}}
    if (i==2) {hTimouY->SetLineColor(kRed);
    		if (i2==0){legendJPsiY->AddEntry(hTimouY, "3.0 < y < 3.5");}}
    if (i==3) {hTimouY->SetLineColor(kOrange);
    		if (i2==0){legendJPsiY->AddEntry(hTimouY, "3.5 < y < 4.0");}
    		legendJPsiY->Draw();}
    hTimouY->Draw("same");}
    }
    
	cPsi2SAloneY->cd();   
    hTimouY2->SetLineWidth(2);
    if (i==0 and i2==0){
    if (j==0) {hTimouY2->SetLineColor(kBlack);
    		if (j2==0){legendPsi2SY->AddEntry(hTimouY2, "Original");}
	hTimouY2->GetYaxis()->SetRangeUser(0., 0.6);
    hTimouY2->Draw();}
    else {
    if (j==1) {hTimouY2->SetLineColor(kBlue);
    		if (j2==0){legendPsi2SY->AddEntry(hTimouY2, "2.5 < y < 3.0");}}
    if (j==2) {hTimouY2->SetLineColor(kRed);
    		if (j2==0){legendPsi2SY->AddEntry(hTimouY2, "3.0 < y < 3.5");}}
    if (j==3) {hTimouY2->SetLineColor(kOrange);
    		if (j2==0){legendPsi2SY->AddEntry(hTimouY2, "3.5 < y < 4.0");}
    		legendPsi2SY->Draw();}
    hTimouY2->Draw("same");}
    }
    
    TDatime date;
    //TFile *fout = new TFile(Form("Ratio_%d.root",date.GetDate()),"recreate");
    //fout->cd();
    //hTimou->Write("hAxEvsPt_J/Psi");
    //hTimou2->Write("hAxEvsPt_Psi(2S)");
    
    // Get the number of bins in the histogram
    int nBins = hTimou->GetNbinsX();
    int nBins2 = hTimou2->GetNbinsX();
    

    // Loop over each bin and print its content
    //cout << "--------- J/Psi y " << i << " and Psi(2S) y " << j << " : -----------" << endl; 
    for (int bin = 1; bin <= nBins; ++bin) {
        double binContent1 = hTimou->GetBinContent(bin);
        double binError1 = hTimou->GetBinError(bin);
        double binLowEdge = hTimou->GetXaxis()->GetBinLowEdge(bin);
    	double binUpEdge = hTimou->GetXaxis()->GetBinUpEdge(bin);
        double binContent2 = hTimou2->GetBinContent(bin);
        double binError2 = hTimou2->GetBinError(bin);
        
        if (j == 0 and j2==0){
    		AxEVecJPsi[bin-1].push_back(binContent1);
    		AxEVecJPsiError[bin-1].push_back(binError1);
    		AxEJPsiTen[bin-1][i].push_back(binContent1);
    		AxEJPsiErrorTen[bin-1][i].push_back(binError1);
    	}
    	if (i == 0 and i2==0){
    		AxEVecPsi2S[bin-1].push_back(binContent2);
    		AxEVecPsi2SError[bin-1].push_back(binError2);
    		AxEPsi2STen[bin-1][i].push_back(binContent2);
    		AxEPsi2SErrorTen[bin-1][i].push_back(binError2);
    		
    		//cout << j << " : " << binContent2 << " +/- " << binError2 << endl;
    	}
        
        if(i==j and (i2 == j2+2 or (i2 < 2 and j2==0)) ) {
    	AxEVecRatio[bin-1].push_back(binContent2/binContent1);
    	AxEVecRatioError[bin-1].push_back(binContent2/binContent1*sqrt((binError1/binContent1)*(binError1/binContent1)+(binError2/binContent2)*(binError2/binContent2)));
    	}
    }
    
    // Get the number of bins in the histogram
    int nBinsY = hTimouY->GetNbinsX();
    int nBinsY2 = hTimouY2->GetNbinsX();

    // Loop over each bin and print its content
    //cout << "--------- J/Psi y " << i << " and Psi(2S) y " << j << " : -----------" << endl; 
    for (int bin = 1; bin <= nBinsY; ++bin) {
        double binContent1 = hTimouY->GetBinContent(bin);
        double binError1 = hTimouY->GetBinError(bin);
        double binLowEdge = hTimouY->GetXaxis()->GetBinLowEdge(bin);
    	double binUpEdge = hTimouY->GetXaxis()->GetBinUpEdge(bin);
        double binContent2 = hTimouY2->GetBinContent(bin);
        double binError2 = hTimouY2->GetBinError(bin);
        
        if (j == 0 and j2==0){
    		AxEVecJPsiY[bin-1].push_back(binContent1);
    		AxEVecJPsiErrorY[bin-1].push_back(binError1);
    	}
    	if (i == 0 and i2==0){
    		AxEVecPsi2SY[bin-1].push_back(binContent2);
    		AxEVecPsi2SErrorY[bin-1].push_back(binError2);
    	}
        if(i==j and (i2 == j2+2 or (i2 < 2 and j2==0)) ) {
    	AxEVecRatioY[bin-1].push_back(binContent2/binContent1);
    	AxEVecRatioErrorY[bin-1].push_back(binContent2/binContent1*sqrt((binError1/binContent1)*(binError1/binContent1)+(binError2/binContent2)*(binError2/binContent2)));
    	}
    }
    }
    }
    }
    }
    
    double u = 0;
    
    cout << "JPsi integrated  : " << endl; 
    	u = 0;
    for (int i = 0; i < AxeIntJPsi.size(); i++) {
    		diffJPsiIntVec.push_back(abs(AxeIntJPsi[i]-AxeIntJPsi[0])/AxeIntJPsi[0]*100);
    		if ( abs(AxeIntJPsi[i] - AxeIntJPsi[0]) > u ) {
    			u = abs(AxeIntJPsi[i] - AxeIntJPsi[0]);
    		}
	    }
    	std::cout << "Integrated : " << AxeIntJPsi[0] << " +/- " << AxeIntJPsiError[0] << " +/- " << u << endl;
	    SystJPsiIntVec.push_back(u);
	    
    cout << "Psi2S integrated  : " << endl; 
    	u = 0;
    for (int i = 0; i < AxeIntPsi2S.size(); i++) {
    		diffPsi2SIntVec.push_back(abs(AxeIntPsi2S[i]-AxeIntPsi2S[0])/AxeIntPsi2S[0]*100);
    		if ( abs(AxeIntPsi2S[i] - AxeIntPsi2S[0]) > u ) {
    			u = abs(AxeIntPsi2S[i] - AxeIntPsi2S[0]);
    		}
    		AxeIntRatio.push_back(AxeIntPsi2S[i]/AxeIntJPsi[i]);
    		AxeIntRatioError.push_back(AxeIntPsi2S[i]/AxeIntJPsi[i]*sqrt( (AxeIntPsi2SError[i]/AxeIntPsi2S[i])*(AxeIntPsi2SError[i]/AxeIntPsi2S[i]) + (AxeIntJPsiError[i]/AxeIntJPsi[i])*(AxeIntJPsiError[i]/AxeIntJPsi[i]) ));
	    }
	    
    	cout << "Integrated : " << AxeIntPsi2S[0] << " +/- " << AxeIntPsi2SError[0] << " +/- " << u << endl;
	SystPsi2SIntVec.push_back(u);
	
    cout << "Ratio integrated  : " << endl; 
    	u = 0;
    for (int i = 0; i < AxeIntRatio.size(); i++) {
    		diffRatioIntVec.push_back(abs(AxeIntRatio[i]-AxeIntRatio[0])/AxeIntRatio[0]*100);
    		if ( abs(AxeIntRatio[i] - AxeIntRatio[0]) > u ) {
    			u = abs(AxeIntRatio[i] - AxeIntRatio[0]);
    		}
	    }
    	cout << "Integrated : " << AxeIntRatio[0] << " +/- " << AxeIntRatioError[0] << " +/- " << u << endl;
	SystRatioIntVec.push_back(u);
	
	
    cout << "JPsi vs Pt  : " << endl; 
    for (int j = 0; j < AxEVecJPsi.size(); j++) {
    	u = 0;
    for (int i = 0; i < AxEVecJPsi[j].size(); i++) {
    		//xVec.push_back(i+j*AxEVecJPsi[j].size());
    		xVec.push_back(ptCenter[j]);
    		diffJPsiVec.push_back(abs(AxEVecJPsi[j][i]-AxEVecJPsi[j][0])/AxEVecJPsi[j][0]*100);
    		if ( abs(AxEVecJPsi[j][i] - AxEVecJPsi[j][0]) > u ) {
    			u = abs(AxEVecJPsi[j][i] - AxEVecJPsi[j][0]);
    		}
	    }
	    double binLowEdge = hTimou->GetXaxis()->GetBinLowEdge(j+1);
    	double binUpEdge = hTimou->GetXaxis()->GetBinUpEdge(j+1);
    	std::cout << "Pt = [" << binLowEdge << "," << binUpEdge << "] : " << AxEVecJPsi[j][0] << " +/- " << AxEVecJPsiError[j][0] << " +/- " << u << endl;
	    SystJPsiVec.push_back(u);
	    }
	    
    cout << "JPsi vs Y  : " << endl; 
    for (int j = 0; j < AxEVecJPsiY.size(); j++) {
    	u = 0;
    for (int i = 0; i < AxEVecJPsiY[j].size(); i++) {
    		//xVecY.push_back(i+j*AxEVecJPsiY[j].size());
    		xVecY.push_back(yCenter[j]);
    		diffJPsiVecY.push_back(abs(AxEVecJPsiY[j][i]-AxEVecJPsiY[j][0])/AxEVecJPsiY[j][0]*100);
    		if ( abs(AxEVecJPsiY[j][i] - AxEVecJPsiY[j][0]) > u ) {
    			u = abs(AxEVecJPsiY[j][i] - AxEVecJPsiY[j][0]);
    		}
	    }
	    double binLowEdge = hTimouY->GetXaxis()->GetBinLowEdge(j+1);
    	double binUpEdge = hTimouY->GetXaxis()->GetBinUpEdge(j+1);
    	std::cout << "Y = [" << binLowEdge << "," << binUpEdge << "] : " << AxEVecJPsiY[j][0] << " +/- " << AxEVecJPsiErrorY[j][0] << " +/- " << u << endl;
	    SystJPsiVecY.push_back(u);
	    }
	    
    cout << "Psi(2S) vs Pt  : " << endl; 
    for (int j = 0; j < AxEVecPsi2S.size(); j++) {
    	u = 0;
    for (int i = 0; i < AxEVecPsi2S[j].size(); i++) {
    		diffPsi2SVec.push_back(abs(AxEVecPsi2S[j][i]-AxEVecPsi2S[j][0])/AxEVecPsi2S[j][0]*100);
    		if ( abs(AxEVecPsi2S[j][i] - AxEVecPsi2S[j][0]) > u ) {
    			u = abs(AxEVecPsi2S[j][i] - AxEVecPsi2S[j][0]);
    		}
	    }
	    double binLowEdge = hTimou2->GetXaxis()->GetBinLowEdge(j+1);
    	double binUpEdge = hTimou2->GetXaxis()->GetBinUpEdge(j+1);
    	std::cout << "Pt = [" << binLowEdge << "," << binUpEdge << "] : " << AxEVecPsi2S[j][0] << " +/- " << AxEVecPsi2SError[j][0] << " +/- " << u << endl;
	    SystPsi2SVec.push_back(u);
	    }
	    
    cout << "Psi(2S) vs Y  : " << endl; 
    for (int j = 0; j < AxEVecPsi2SY.size(); j++) {
    	u = 0;
    for (int i = 0; i < AxEVecPsi2SY[j].size(); i++) {
    		diffPsi2SVecY.push_back(abs(AxEVecPsi2SY[j][i]-AxEVecPsi2SY[j][0])/AxEVecPsi2SY[j][0]*100);
    		if ( abs(AxEVecPsi2SY[j][i] - AxEVecPsi2SY[j][0]) > u ) {
    			u = abs(AxEVecPsi2SY[j][i] - AxEVecPsi2SY[j][0]);
    		}
	    }
	    double binLowEdge = hTimouY2->GetXaxis()->GetBinLowEdge(j+1);
    	double binUpEdge = hTimouY2->GetXaxis()->GetBinUpEdge(j+1);
    	std::cout << "Y = [" << binLowEdge << "," << binUpEdge << "] : " << AxEVecPsi2SY[j][0] << " +/- " << AxEVecPsi2SErrorY[j][0] << " +/- " << u << endl;
	    SystPsi2SVecY.push_back(u);
	    }
    cout << "Ratio vs Pt  : " << endl; 
    for (int j = 0; j < AxEVecRatio.size(); j++) {
    	u = 0;
    for (int i = 0; i < AxEVecRatio[j].size(); i++) {
    		//xVecRatio.push_back(i+j*AxEVecRatio[j].size());
    		xVecRatio.push_back(ptCenter[j]);
    		diffRatioVec.push_back(abs(AxEVecRatio[j][i]-AxEVecRatio[j][0])/AxEVecRatio[j][0]*100);
    		if ( abs(AxEVecRatio[j][i] - AxEVecRatio[j][0]) > u ) {
    			u = abs(AxEVecRatio[j][i] - AxEVecRatio[j][0]);
    		}
	    }
	    double binLowEdge = hTimou->GetXaxis()->GetBinLowEdge(j+1);
    	double binUpEdge = hTimou->GetXaxis()->GetBinUpEdge(j+1);
    	std::cout << "Pt = [" << binLowEdge << "," << binUpEdge << "] : " << AxEVecRatio[j][0] << " +/- " << AxEVecRatioError[j][0] << " +/- " << u << endl;
    	//cout << u << ", ";
	    SystRatioVec.push_back(u);
	    }
	    
    cout << "Ratio vs Y  : " << endl; 
    for (int j = 0; j < AxEVecRatioY.size(); j++) {
    	u = 0;
    for (int i = 0; i < AxEVecRatioY[j].size(); i++) {
    		//xVecRatioY.push_back(i+j*AxEVecRatio[j].size());
    		xVecRatioY.push_back(yCenter[j]);
    		diffRatioVecY.push_back(abs(AxEVecRatioY[j][i]-AxEVecRatioY[j][0])/AxEVecRatioY[j][0]*100);
    		if ( abs(AxEVecRatioY[j][i] - AxEVecRatioY[j][0]) > u ) {
    			u = abs(AxEVecRatioY[j][i] - AxEVecRatioY[j][0]);
    		}
	    }
	    double binLowEdge = hTimouY->GetXaxis()->GetBinLowEdge(j+1);
    	double binUpEdge = hTimouY->GetXaxis()->GetBinUpEdge(j+1);
    	std::cout << "Y = [" << binLowEdge << "," << binUpEdge << "] : " << AxEVecRatioY[j][0] << " +/- " << AxEVecRatioErrorY[j][0] << " +/- " << u << endl;
    	//cout << u << ", ";
	    SystRatioVecY.push_back(u);
	    }
	    cout << endl;
	    
	    vector<double> AxEVecJPsi0, AxEVecJPsiError0, AxEVecPsi2S0, AxEVecPsi2SError0, AxEVecJPsiY0, AxEVecJPsiErrorY0, AxEVecPsi2SY0, AxEVecPsi2SErrorY0, AxEVecRatio0, AxEVecRatioError0, AxEVecRatioY0, AxEVecRatioErrorY0;
	    
	    for (size_t i = 0; i < AxEVecJPsi.size(); i++) {
		AxEVecJPsi0.push_back(AxEVecJPsi[i][0]);
		AxEVecJPsiError0.push_back(AxEVecJPsiError[i][0]);
		AxEVecPsi2S0.push_back(AxEVecPsi2S[i][0]);
		AxEVecPsi2SError0.push_back(AxEVecPsi2SError[i][0]);
		AxEVecRatio0.push_back(AxEVecRatio[i][0]);
		AxEVecRatioError0.push_back(AxEVecRatioError[i][0]);
		}
		
	    for (size_t i = 0; i < AxEVecJPsiY.size(); i++) {
		AxEVecJPsiY0.push_back(AxEVecJPsiY[i][0]);
		AxEVecJPsiErrorY0.push_back(AxEVecJPsiErrorY[i][0]);
		AxEVecPsi2SY0.push_back(AxEVecPsi2SY[i][0]);
		AxEVecPsi2SErrorY0.push_back(AxEVecPsi2SErrorY[i][0]);
		AxEVecRatioY0.push_back(AxEVecRatioY[i][0]);
		AxEVecRatioErrorY0.push_back(AxEVecRatioErrorY[i][0]);
	    }
	    
TCanvas * cJPsi = new TCanvas("cJPsi","Acceptance efficiency vs Pt (J/Psi)");
	    cJPsi->cd();    
    TGraphMultiErrors* gme = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), AxEVecJPsi0.data(), eptCenter.data(), eptCenter.data(), AxEVecJPsiError0.data(), AxEVecJPsiError0.data());
gme->AddYError(SystJPsiVec.size(), SystJPsiVec.data(), SystJPsiVec.data());
gme->SetTitle("Acceptance efficiency vs pT (JPsi)");
    gme->GetXaxis()->SetTitle("pT");
    gme->GetYaxis()->SetTitle("AxE");
gme->SetMarkerStyle(20);
gme->SetMarkerColor(kBlue);
gme->SetLineColor(kBlue);
gme->GetAttLine(0)->SetLineColor(kBlue);
gme->GetAttLine(1)->SetLineColor(kBlue);
gme->GetAttFill(1)->SetFillStyle(0);
gme->GetYaxis()->SetRangeUser(0., 0.6);
    gme->Draw("a p s ; ; 5 s=1");


TCanvas * cPsi2S = new TCanvas("cPsi2S","Acceptance efficiency vs Pt (Psi2S)");
cPsi2S->cd();
TGraphMultiErrors* gme2 = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), AxEVecPsi2S0.data(), eptCenter.data(), eptCenter.data(), AxEVecPsi2SError0.data(), AxEVecPsi2SError0.data());
gme2->AddYError(SystPsi2SVec.size(), SystPsi2SVec.data(), SystPsi2SVec.data());
gme2->SetTitle("Acceptance efficiency vs pT (Psi2S)");
gme2->GetXaxis()->SetTitle("pT");
gme2->GetYaxis()->SetTitle("AxE");
gme2->SetMarkerStyle(20);
gme2->SetMarkerColor(kRed);
gme2->SetLineColor(kRed);
gme2->GetAttLine(0)->SetLineColor(kRed);
gme2->GetAttLine(1)->SetLineColor(kRed);
gme2->GetAttFill(1)->SetFillStyle(0);
gme2->GetYaxis()->SetRangeUser(0., 0.6);
    gme2->Draw("a p s ; ; 5 s=1");
    
TCanvas * cRatio = new TCanvas("cRatio","Ratio of acceptance efficiency vs Pt");
	    cRatio->cd();    
    TGraphMultiErrors* gmeRatio = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), AxEVecRatio0.data(), eptCenter.data(), eptCenter.data(), AxEVecRatioError0.data(), AxEVecRatioError0.data());
gmeRatio->AddYError(SystRatioVec.size(), SystRatioVec.data(), SystRatioVec.data());
gmeRatio->SetTitle("Ratio of acceptance efficiency vs pT");
    gmeRatio->GetXaxis()->SetTitle("pT");
    gmeRatio->GetYaxis()->SetTitle("AxE");
gmeRatio->SetMarkerStyle(20);
gmeRatio->SetMarkerColor(kOrange);
gmeRatio->SetLineColor(kOrange);
gmeRatio->GetAttLine(0)->SetLineColor(kOrange);
gmeRatio->GetAttLine(1)->SetLineColor(kOrange);
gmeRatio->GetAttFill(1)->SetFillStyle(0);
gmeRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
    gmeRatio->Draw("a p s ; ; 5 s=1");
    
TCanvas * cJPsiY = new TCanvas("cJPsiY","Acceptance efficiency vs Y (J/Psi)");
	    cJPsiY->cd();    
    TGraphMultiErrors* gmeY = new TGraphMultiErrors(yCenter.size(), yCenter.data(), AxEVecJPsiY0.data(), eyCenter.data(), eyCenter.data(), AxEVecJPsiErrorY0.data(), AxEVecJPsiErrorY0.data());
gmeY->AddYError(SystJPsiVecY.size(), SystJPsiVecY.data(), SystJPsiVecY.data());
gmeY->SetTitle("Acceptance efficiency vs y (JPsi)");
    gmeY->GetXaxis()->SetTitle("y");
    gmeY->GetYaxis()->SetTitle("AxE");
gmeY->SetMarkerStyle(20);
gmeY->SetMarkerColor(kBlue);
gmeY->SetLineColor(kBlue);
gmeY->GetAttLine(0)->SetLineColor(kBlue);
gmeY->GetAttLine(1)->SetLineColor(kBlue);
gmeY->GetAttFill(1)->SetFillStyle(0);
gmeY->GetYaxis()->SetRangeUser(0., 0.5);
    gmeY->Draw("a p s ; ; 5 s=1");


TCanvas * cPsi2SY = new TCanvas("cPsi2SY","Acceptance efficiency vs Y (Psi2S)");
cPsi2SY->cd();
TGraphMultiErrors* gme2Y = new TGraphMultiErrors(yCenter.size(), yCenter.data(), AxEVecPsi2SY0.data(), eyCenter.data(), eyCenter.data(), AxEVecPsi2SErrorY0.data(), AxEVecPsi2SErrorY0.data());
gme2Y->AddYError(SystPsi2SVecY.size(), SystPsi2SVecY.data(), SystPsi2SVecY.data());
gme2Y->SetTitle("Acceptance efficiency vs y (Psi2S)");
gme2Y->GetXaxis()->SetTitle("y");
gme2Y->GetYaxis()->SetTitle("AxE");
gme2Y->SetMarkerStyle(20);
gme2Y->SetMarkerColor(kRed);
gme2Y->SetLineColor(kRed);
gme2Y->GetAttLine(0)->SetLineColor(kRed);
gme2Y->GetAttLine(1)->SetLineColor(kRed);
gme2Y->GetAttFill(1)->SetFillStyle(0);
gme2Y->GetYaxis()->SetRangeUser(0., 0.5);
    gme2Y->Draw("a p s ; ; 5 s=1");

TCanvas * cRatioY = new TCanvas("cRatioY","Ratio of acceptance efficiency vs Y");
cRatioY->cd();
TGraphMultiErrors* gmeRatioY = new TGraphMultiErrors(yCenter.size(), yCenter.data(), AxEVecRatioY0.data(), eyCenter.data(), eyCenter.data(), AxEVecRatioErrorY0.data(), AxEVecRatioErrorY0.data());
gmeRatioY->AddYError(SystRatioVecY.size(), SystRatioVecY.data(), SystRatioVecY.data());
gmeRatioY->SetTitle("Ratio of acceptance efficiency vs y");
gmeRatioY->GetXaxis()->SetTitle("y");
gmeRatioY->GetYaxis()->SetTitle("AxE");
gmeRatioY->SetMarkerStyle(20);
gmeRatioY->SetMarkerColor(kOrange);
gmeRatioY->SetLineColor(kOrange);
gmeRatioY->GetAttLine(0)->SetLineColor(kOrange);
gmeRatioY->GetAttLine(1)->SetLineColor(kOrange);
gmeRatioY->GetAttFill(1)->SetFillStyle(0);
gmeRatioY->GetYaxis()->SetRangeUser(0.9, 1.3);
    gmeRatioY->Draw("a p s ; ; 5 s=1");
    
TCanvas * cPT = new TCanvas("cPT","MC Input uncertainty vs Pt");
cPT->cd();    
TGraph* MCInputRatioJPsi = new TGraph(diffJPsiVec.size(), xVec.data(), diffJPsiVec.data());
TGraph* MCInputRatioPsi2S = new TGraph(diffPsi2SVec.size(), xVec.data(), diffPsi2SVec.data());
MCInputRatioJPsi->SetTitle("MC Input uncertainty vs Pt");
MCInputRatioJPsi->GetXaxis()->SetTitle("Pt (GeV/c)");
MCInputRatioJPsi->GetYaxis()->SetTitle("|AxE-AxE0|/AxE0 (%)");
MCInputRatioJPsi->SetMarkerStyle(20);
MCInputRatioJPsi->SetMarkerColor(kBlue);
MCInputRatioPsi2S->SetMarkerStyle(20);
MCInputRatioPsi2S->SetMarkerColor(kRed);
MCInputRatioJPsi->GetYaxis()->SetRangeUser(0., 10);
    TLegend* legendRatio = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendRatio->AddEntry(MCInputRatioJPsi, "J/Psi");
    legendRatio->AddEntry(MCInputRatioPsi2S, "Psi(2S)");

MCInputRatioJPsi->Draw("a p s");
MCInputRatioPsi2S->Draw("p s");
legendRatio->Draw();


TCanvas * cY = new TCanvas("cY","MC Input uncertainty vs Y");
cY->cd();    
TGraph* MCInputRatioJPsiY = new TGraph(diffJPsiVecY.size(), xVecY.data(), diffJPsiVecY.data());
TGraph* MCInputRatioPsi2SY = new TGraph(diffPsi2SVecY.size(), xVecY.data(), diffPsi2SVecY.data());
MCInputRatioJPsiY->SetTitle("MC Input uncertainty vs Y");
MCInputRatioJPsiY->GetXaxis()->SetTitle("y");
MCInputRatioJPsiY->GetYaxis()->SetTitle("|AxE-AxE0|/AxE0 (%)");
MCInputRatioJPsiY->SetMarkerStyle(20);
MCInputRatioJPsiY->SetMarkerColor(kBlue);
MCInputRatioPsi2SY->SetMarkerStyle(20);
MCInputRatioPsi2SY->SetMarkerColor(kRed);
MCInputRatioJPsiY->GetYaxis()->SetRangeUser(0, 10);
    TLegend* legendRatioY = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendRatioY->AddEntry(MCInputRatioJPsiY, "J/Psi");
    legendRatioY->AddEntry(MCInputRatioPsi2SY, "Psi(2S)");

MCInputRatioJPsiY->Draw("a p s");
MCInputRatioPsi2SY->Draw("p s");
legendRatioY->Draw();



TCanvas * cPtRatio = new TCanvas("cPtRatio","MC Input ratio uncertainty vs Pt");
cPtRatio->cd();    
TGraph* MCInputRatio = new TGraph(diffRatioVec.size(), xVecRatio.data(), diffRatioVec.data());
MCInputRatio->SetTitle("MC Input ratio uncertainty vs Pt");
MCInputRatio->GetXaxis()->SetTitle("Pt (GeV/c)");
MCInputRatio->GetYaxis()->SetTitle("|AxE-AxE0|/AxE0 (%)");
MCInputRatio->SetMarkerStyle(20);
MCInputRatio->SetMarkerColor(kOrange);
MCInputRatio->GetYaxis()->SetRangeUser(0., 10.);

MCInputRatio->Draw("a p s");

TCanvas * cYRatio = new TCanvas("cYRatio","MC Input ratio uncertainty vs Y");
cYRatio->cd();    
TGraph* MCInputRatioY = new TGraph(diffRatioVecY.size(), xVecRatioY.data(), diffRatioVecY.data());
MCInputRatioY->SetTitle("MC Input ratio uncertainty vs Y");
MCInputRatioY->GetXaxis()->SetTitle("y");
MCInputRatioY->GetYaxis()->SetTitle("|AxE-AxE0|/AxE0 (%)");
MCInputRatioY->SetMarkerStyle(20);
MCInputRatioY->SetMarkerColor(kOrange);
MCInputRatioY->GetYaxis()->SetRangeUser(0., 10.);

MCInputRatioY->Draw("a p s");

	   return;
	   
    
    TCanvas * cTime = new TCanvas("cTime","Acceptance efficiency vs Pt");
    cTime->Divide(1,2);
    cTime->cd(1);
    hTimou->Draw();
    hTimou2->SetLineColor(kRed);
    hTimou2->Draw("same");
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hTimou, "J/Psi");
    TLegend* legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2->AddEntry(hTimou2, "Psi(2S)");
    // Draw the legend
    legend->Draw();
    legend2->Draw();
     //gPad->BuildLegend();
    
    //TCanvas * cTime2 = new TCanvas("cTime2","Real data time Distribution 2");
    TH1F *hRecnew1 = dynamic_cast<TH1F*>(hRec->Rebin(8,"",ptbins3));
    TH1F *hSimnew1 = dynamic_cast<TH1F*>(hSim->Rebin(8,"",ptbins3));
    TH1F *hTimounew2 = dynamic_cast<TH1F*>(hTimou2->Rebin(8,"",ptbins3));
    
    hTimouAll_ = (TH1F*) hTimouAll_->Rebin(8,"",ptbins3);
    hTimouAll = (TH1F*) hTimouAll->Rebin(8,"",ptbins3);
    hTimouAll->SetTitle("Ratio of acceptance efficiency vs transverse momentum;p_{T} (GeV/c);Ratio of acc #times eff");
    hTimouAll_->Divide(hRecnew1,hSimnew1, 1.0, 1.0, "B");
    hTimouAll->Divide(hTimounew2,hTimouAll_, 1.0, 1.0, "B");
    
    cTime->cd(2);
    hTimouAll->Draw();
    TLegend* legendR = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendR->AddEntry(hTimouAll, "Psi(2S) acc #times eff / J/Psi acc #times eff");
    legendR->Draw("same");
    
    TCanvas * cTimeY = new TCanvas("cTimeY","Acceptance efficiency vs Y");
    cTimeY->Divide(1,2);
    cTimeY->cd(1);
    hTimouY->SetMaximum(0.45);
    hTimouY->SetLineColor(kBlack);
    hTimouY->SetTitle("Acceptance efficiency vs rapidity;y;acc #times eff");
    hTimouY2->SetTitle("Acceptance efficiency vs rapidity;y;acc #times eff");
    //fout->cd();
    //hTimouY->Write("hAxEvsY_J/Psi");
    //hTimouY2->Write("hAxEvsY_Psi(2S)");
    
    
    //hTimouY->SetMarkerStyle(6);
    hTimouY->Draw();
    hTimouY2->SetLineColor(kRed);
    hTimouY2->SetMarkerColor(kRed);
    //hTimouY2->SetMarkerStyle(6);
    hTimouY2->Draw("same"); 
    TLegend* legendY = new TLegend(0.7, 0.7, 0.9, 0.9);
    //legendY->AddEntry(hTimouY, "J/Psi");
    TLegend* legendY2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    //legendY2->AddEntry(hTimouY2, "Psi(2S)");
    // Draw the legend
    //legendY->Draw();
    //legendY2->Draw();
    //gPad->BuildLegend();
    
    cTimeY->cd(2);
    hTimouAllY2->SetTitle("Ratio of acceptance efficiency vs rapidity;Y;Ratio of acc #times eff");
    hTimouAllY2_->Divide(hRecY2,hSimY2, 1.0, 1.0, "B");
    hTimouAllY2->Divide(hTimouAllY2_,hTimouY, 1.0, 1.0, "B");
    hTimouAllY2->SetMarkerColor(kBlue);
    //hTimouAllY2->SetMarkerStyle(6);
    hTimouAllY2->Draw();
    TLegend* legendRY = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendRY->AddEntry(hTimouAllY2, "Psi(2S) acc #times eff / J/Psi acc #times eff");
    legendRY->Draw();
    
    return;
    
    //fout->Close();
    
    return;
}
