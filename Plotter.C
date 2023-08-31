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

void Plotter(const char *filename = "AE_pp13TeV_Charmonia_Real_0.root", const char *filename2 = "AE_pp13TeV_Charmonia_20230710_Psi2S.root")
{
//______________________
    //Variables declaration
    int entries_number = 0;
    int entries_number2 = 0;
    Double_t ptbins[19] = {0.,0.5,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,18.,20.,22.,26.,30.};
    //Double_t ptbins[16] = {0.,0.5,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,18.,20.};
    Double_t ptbins2[14] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,20.};
    Double_t ptbins3[9] = {0.,1.,2.,3.,4.,5.,7.,10.,20.};
    //Double_t ybins[7] = {2.5,2.75,3.,3.25,3.5,3.75,4.};
    Double_t ybins[7] = {-4.,-3.75,-3.5,-3.25,-3.,-2.75,-2.5};
    
    TH1F *hTime   = new TH1F("hTime","Time distribution",700000,-200000,500000);
    //TH1F *hTime1   = new TH1F("hTime1","Time distribution for rof 1",10000,-200000,500000);
    //TH1F *hTime2   = new TH1F("hTime2","Time distribution for rof 2",10000,-200000,500000);
    TH1F *hTimeMC   = new TH1F("hTimeMC","MC time distribution",700000,-200000,500000);
    TH1F *hTimeEmbed   = new TH1F("hTimeEmbed","Embed time distribution",700000,-200000,500000);
    //TH1F *hTimeDif   = new TH1F("hTimeDif","Time difference between first and last digit distribution",10100,-10000,100);
    //TH1F *hTimeDifRof   = new TH1F("hTimeDifRof","Time difference between consecutive rofs distribution",10100,-100,10000);
    //TH1F *hTimeDifDig   = new TH1F("hTimeDifDig","Time difference between consecutive digits distribution",2000,-1000,1000);
    

    // J/Psi extraction
    TFile *file = TFile::Open(filename,"read");
    if(!file){cout << "Error :: The data file is not found, please check"<< endl; return;}
    TH1F *hTimou = (TH1F*) file->Get("Jpsi/fhEffPt0-0");
    TH1F *hTimou1 = (TH1F*) file->Get("Jpsi/fhEffPt0-0");
    TH1F *hTimouAll = (TH1F*) file->Get("Jpsi/fhEffPt0-0");
    TH1F *hRec = (TH1F*) file->Get("Jpsi/fhRecPt0-0");
    TH1F *hSim = (TH1F*) file->Get("Jpsi/fhSimPt0-0");
    TH1F *hTimouY = (TH1F*) file->Get("Jpsi/fhEffY0-0");
    TH1F *hRecY = (TH1F*) file->Get("Jpsi/fhRecY0-0");
    TH1F *hSimY = (TH1F*) file->Get("Jpsi/fhSimY0-0");
    if(!hTimou){cout << "Error :: The Histogram 1 is not found, please check"<< endl; return;}
    
    
    // Psi(2S) extraction
    TFile *file2 = TFile::Open(filename,"read");
    if(!file2){cout << "Error :: The data file is not found, please check"<< endl; return;}
    TH1F *hTimou2 = (TH1F*) file->Get("PsiPrime/fhEffPt0-0");
    TH1F *hRec2 = (TH1F*) file->Get("PsiPrime/fhRecPt0-0");
    TH1F *hSim2 = (TH1F*) file->Get("PsiPrime/fhSimPt0-0");
    TH1F *hTimouY2 = (TH1F*) file->Get("PsiPrime/fhEffY0-0");
    TH1F *hTimouAllY2 = (TH1F*) file->Get("PsiPrime/fhRecY0-0");
    TH1F *hTimouAllY2_ = (TH1F*) file->Get("PsiPrime/fhRecY0-0");
    TH1F *hRecY2 = (TH1F*) file->Get("PsiPrime/fhRecY0-0");
    TH1F *hSimY2 = (TH1F*) file->Get("PsiPrime/fhSimY0-0");
    if(!hTimou2){cout << "Error :: The Histogram 2 is not found, please check"<< endl; return;}
    
  
  double ErrorJPsiRec = 0;
  double ErrorPsi2SRec = 0;
  double ErrorJPsiSim = 0;
  double ErrorPsi2SSim = 0;
  for (int i=0;i<hRec->GetNbinsX();i++) {
  	ErrorJPsiRec += hRec->GetBinError(i)*hRec->GetBinError(i);
  	ErrorPsi2SRec += hRec2->GetBinError(i)*hRec2->GetBinError(i);
  	ErrorJPsiSim += hSim->GetBinError(i)*hSim->GetBinError(i);
  	ErrorPsi2SSim += hSim2->GetBinError(i)*hSim2->GetBinError(i);
  }
  double ErrorJPsi = hRec->GetEntries()/hSim->GetEntries()*sqrt(ErrorJPsiRec/(hRec->GetEntries()*hRec->GetEntries()) + ErrorJPsiSim/(hSim->GetEntries()*hSim->GetEntries()));
  double ErrorPsi2S = hRec2->GetEntries()/hSim2->GetEntries()*sqrt(ErrorPsi2SRec/(hRec2->GetEntries()*hRec2->GetEntries()) + ErrorPsi2SSim/(hSim2->GetEntries()*hSim2->GetEntries()));
  
  
  cout << "Global ratio JPsi : " << hRec->GetEntries()/hSim->GetEntries() << " +/- " << ErrorJPsi << endl;
  cout << "Global ratio Psi(2S) : " << hRec2->GetEntries()/hSim2->GetEntries() << " +/- " << ErrorPsi2S << endl;
    
    TCanvas * cTime = new TCanvas("cTime","Acceptance efficiency vs Pt");
    cTime->Divide(1,2);
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
    
   
    
    TDatime date;
    TFile *fout = new TFile(Form("Ratio_%d.root",date.GetDate()),"recreate");
    fout->cd();
    hTimou->Write("hAxEvsPt_J/Psi");
    hTimou2->Write("hAxEvsPt_Psi(2S)");
    
    // Get the number of bins in the histogram
    int nBins = hTimou->GetNbinsX();
    int nBins2 = hTimou2->GetNbinsX();

    // Loop over each bin and print its content
    cout << " J/Psi :" << endl; 
    for (int bin = 1; bin <= nBins; ++bin) {
        double binContent = hTimou->GetBinContent(bin);
        double binError = hTimou->GetBinError(bin);
        double binLowEdge = hTimou->GetXaxis()->GetBinLowEdge(bin);
    	double binUpEdge = hTimou->GetXaxis()->GetBinUpEdge(bin);

    	std::cout << "pT = [" << binLowEdge << "," << binUpEdge << "] GeV/c : " << binContent << " +/- " << binError << std::endl;
    }
    
    // Loop over each bin and print its content
    cout << " Psi(2S) :" << endl; 
    for (int bin = 1; bin <= nBins2; ++bin) {
        double binContent = hTimou2->GetBinContent(bin);
        double binError = hTimou2->GetBinError(bin);
        double binLowEdge = hTimou2->GetXaxis()->GetBinLowEdge(bin);
    	double binUpEdge = hTimou2->GetXaxis()->GetBinUpEdge(bin);

    	std::cout << "pT = [" << binLowEdge << "," << binUpEdge << "] GeV/c : " << binContent << " +/- " << binError << std::endl;
    }
    
    // Get the number of bins in the histogram
    int nBinsY = hTimouY->GetNbinsX();
    int nBinsY2 = hTimouY2->GetNbinsX();

    // Loop over each bin and print its content
    cout << " J/Psi :" << endl; 
    for (int bin = 1; bin <= nBinsY; ++bin) {
        double binContent = hTimouY->GetBinContent(bin);
        double binError = hTimouY->GetBinError(bin);
        double binLowEdge = hTimouY->GetXaxis()->GetBinLowEdge(bin);
    	double binUpEdge = hTimouY->GetXaxis()->GetBinUpEdge(bin);

    	std::cout << "y = [" << binLowEdge << "," << binUpEdge << "] : " << binContent << " +/- " << binError << std::endl;
    }
    
    // Loop over each bin and print its content
    cout << " Psi(2S) :" << endl; 
    for (int bin = 1; bin <= nBinsY2; ++bin) {
        double binContent = hTimouY2->GetBinContent(bin);
        double binError = hTimouY2->GetBinError(bin);
        double binLowEdge = hTimouY2->GetXaxis()->GetBinLowEdge(bin);
    	double binUpEdge = hTimouY2->GetXaxis()->GetBinUpEdge(bin);

    	std::cout << "y = [" << binLowEdge << "," << binUpEdge << "] : " << binContent << " +/- " << binError << std::endl;
    }
    
    cTime->cd(1);
    hTimou->Draw();
    hTimou2->SetLineColor(kRed);
    hTimou2->Draw("same");
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hTimou, "J/Psi");
    legend->AddEntry(hTimou2, "Psi(2S)");
    // Draw the legend
    legend->Draw();
     //gPad->BuildLegend();
    
    //TCanvas * cTime2 = new TCanvas("cTime2","Real data time Distribution 2");
    TH1F *hRecnew1 = dynamic_cast<TH1F*>(hRec->Rebin(13,"",ptbins2));
    TH1F *hSimnew1 = dynamic_cast<TH1F*>(hSim->Rebin(13,"",ptbins2));
    hTimou1 = hRecnew1;
    hTimouAll = hRecnew1;
    
    cTime->cd(2);
    hTimouAll->SetTitle("Ratio of acceptance efficiency vs transverse momentum;p_{T} (GeV/c);Ratio of acc #times eff");
    hTimou1->Divide(hRecnew1,hSimnew1, 1.0, 1.0, "B");
    hTimouAll->Divide(hTimou2,hTimou1, 1.0, 1.0, "B");
    hTimouAll->Draw();
    //TLegend* legendR = new TLegend(0.7, 0.7, 0.9, 0.9);
    //legendR->AddEntry(hTimouAll, "Psi(2S) acc #times eff / J/Psi acc #times eff");
    //legendR->Draw();
    
    TCanvas * cTimeY = new TCanvas("cTimeY","Acceptance efficiency vs Y");
    cTimeY->Divide(1,2);
    cTimeY->cd(1);
    hTimouY->SetMaximum(0.5);
    hTimouY->SetLineColor(kBlack);
    hTimouY->SetTitle("Acceptance efficiency vs rapidity;y;acc #times eff");
    hTimouY->SetTitle("Acceptance efficiency vs rapidity;y;acc #times eff");
    
    
    //fout->cd();
    //hTimouY->Write("hAxEvsY_J/Psi");
    //hTimouY2->Write("hAxEvsY_Psi(2S)");
    //fout->Close();
    
    //hTimouY->SetMarkerStyle(6);
    hTimouY->Draw();
    hTimouY2->SetLineColor(kRed);
    hTimouY2->SetMarkerColor(kRed);
    //hTimouY2->SetMarkerStyle(6);
    hTimouY2->Draw("same"); 
    TLegend* legendY = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendY->AddEntry(hTimouY, "J/Psi");
    legendY->AddEntry(hTimouY2, "Psi(2S)");
    // Draw the legend
    legendY->Draw();
    //gPad->BuildLegend();
    
    cTimeY->cd(2);
    hTimouAllY2->SetTitle("Ratio of acceptance efficiency vs rapidity;Y;Ratio of acc #times eff");
    hTimouAllY2_->Divide(hRecY2,hSimY2, 1.0, 1.0, "B");
    hTimouAllY2->Divide(hTimouAllY2_,hTimouY, 1.0, 1.0, "B");
    hTimouAllY2->SetMarkerColor(kBlue);
    //hTimouAllY2->SetMarkerStyle(6);
    hTimouAllY2->Draw();
    //TLegend* legendRY = new TLegend(0.7, 0.7, 0.9, 0.9);
    //legendRY->AddEntry(hTimouAllY2, "Psi(2S) acc #times eff / J/Psi acc #times eff");
    //legendRY->Draw();
    
    return;
}
