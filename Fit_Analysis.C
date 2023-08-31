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
#include <stdio.h>
#include <dirent.h>
#include <string.h>

Double_t getBinContent(TH1* histogram, Int_t bin) {
    if (!histogram) {
        std::cerr << "Invalid histogram pointer!" << std::endl;
        return 0.0;
    }

    if (bin < 1 || bin > histogram->GetNbinsX()) {
        std::cerr << "Invalid bin number!" << std::endl;
        return 0.0;
    }

    return histogram->GetBinContent(bin);
}

pair<double,double> printFileContent(const char* file_path) {
    // Open the ROOT file
    TFile* file = TFile::Open(file_path);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return make_pair(0,0);
    }

    // Get the list of keys in the file
    TList* keyList = file->GetListOfKeys();
    if (!keyList || keyList->GetSize() == 0) {
        std::cerr << "No objects found in the file: " << file_path << std::endl;
        file->Close();
        return make_pair(0,0);
    }

    // Get the first object from the list of keys (assumed to be a TH1)
    TKey* key = dynamic_cast<TKey*>(keyList->At(0));
    if (!key) {
        std::cerr << "No valid object found in the file: " << file_path << std::endl;
        file->Close();
        return make_pair(0,0);
    }

    // Cast the TKey to TH1
    TH1* histogram = dynamic_cast<TH1*>(key->ReadObj());
    if (!histogram) {
        std::cerr << "Failed to cast to TH1!" << std::endl;
        file->Close();
        return make_pair(0,0);
    }

    // Access the bin content in the histogram
    TString label_1 = "sig_Jpsi";
    TString label_2 = "sig_Psi2s";
    Int_t binNumber_JPsi = histogram->GetXaxis()->FindBin(label_1);
    Int_t binNumber_Psi2S = histogram->GetXaxis()->FindBin(label_2);
    Double_t binContent_JPsi = histogram->GetBinContent(binNumber_JPsi);
    Double_t binContent_Psi2S = histogram->GetBinContent(binNumber_Psi2S);
    Double_t binError_JPsi = histogram->GetBinError(binNumber_JPsi);
    Double_t binError_Psi2S = histogram->GetBinError(binNumber_Psi2S);
    //std::cout << "J/Psi : " << binContent_JPsi << " - Psi(2S) : " << binContent_Psi2S << " - Ratio : " << binContent_Psi2S/binContent_JPsi << std::endl;
    

    // Close the file
    file->Close();
    
    
    //cout << "J/Psi : " << (binError_JPsi/binContent_JPsi)*(binError_JPsi/binContent_JPsi) << ", Psi(2S) : " << (binError_Psi2S/binContent_Psi2S)*(binError_Psi2S/binContent_Psi2S) << endl;
    
    return make_pair(binContent_Psi2S/binContent_JPsi, binContent_Psi2S/binContent_JPsi*sqrt( (binError_JPsi/binContent_JPsi)*(binError_JPsi/binContent_JPsi) + (binError_Psi2S/binContent_Psi2S)*(binError_Psi2S/binContent_Psi2S)));
    
}

pair<double,double> printFileContentJPsi(const char* file_path) {
    // Open the ROOT file
    TFile* file = TFile::Open(file_path);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return make_pair(0,0);
    }

    // Get the list of keys in the file
    TList* keyList = file->GetListOfKeys();
    if (!keyList || keyList->GetSize() == 0) {
        std::cerr << "No objects found in the file: " << file_path << std::endl;
        file->Close();
        return make_pair(0,0);
    }

    // Get the first object from the list of keys (assumed to be a TH1)
    TKey* key = dynamic_cast<TKey*>(keyList->At(0));
    if (!key) {
        std::cerr << "No valid object found in the file: " << file_path << std::endl;
        file->Close();
        return make_pair(0,0);
    }

    // Cast the TKey to TH1
    TH1* histogram = dynamic_cast<TH1*>(key->ReadObj());
    if (!histogram) {
        std::cerr << "Failed to cast to TH1!" << std::endl;
        file->Close();
        return make_pair(0,0);
    }

    // Access the bin content in the histogram
    TString label_1 = "sig_Jpsi";
    TString label_2 = "sig_Psi2s";
    Int_t binNumber_JPsi = histogram->GetXaxis()->FindBin(label_1);
    Int_t binNumber_Psi2S = histogram->GetXaxis()->FindBin(label_2);
    Double_t binContent_JPsi = histogram->GetBinContent(binNumber_JPsi);
    Double_t binContent_Psi2S = histogram->GetBinContent(binNumber_Psi2S);
    Double_t binError_JPsi = histogram->GetBinError(binNumber_JPsi);
    Double_t binError_Psi2S = histogram->GetBinError(binNumber_Psi2S);
    

    // Close the file
    file->Close();
    
    return make_pair(binContent_JPsi, binError_JPsi);
}

pair<double,double> printFileContentPsi2S(const char* file_path) {
    // Open the ROOT file
    TFile* file = TFile::Open(file_path);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return make_pair(0,0);
    }

    // Get the list of keys in the file
    TList* keyList = file->GetListOfKeys();
    if (!keyList || keyList->GetSize() == 0) {
        std::cerr << "No objects found in the file: " << file_path << std::endl;
        file->Close();
        return make_pair(0,0);
    }

    // Get the first object from the list of keys (assumed to be a TH1)
    TKey* key = dynamic_cast<TKey*>(keyList->At(0));
    if (!key) {
        std::cerr << "No valid object found in the file: " << file_path << std::endl;
        file->Close();
        return make_pair(0,0);
    }

    // Cast the TKey to TH1
    TH1* histogram = dynamic_cast<TH1*>(key->ReadObj());
    if (!histogram) {
        std::cerr << "Failed to cast to TH1!" << std::endl;
        file->Close();
        return make_pair(0,0);
    }

    // Access the bin content in the histogram
    TString label_1 = "sig_Jpsi";
    TString label_2 = "sig_Psi2s";
    Int_t binNumber_JPsi = histogram->GetXaxis()->FindBin(label_1);
    Int_t binNumber_Psi2S = histogram->GetXaxis()->FindBin(label_2);
    Double_t binContent_JPsi = histogram->GetBinContent(binNumber_JPsi);
    Double_t binContent_Psi2S = histogram->GetBinContent(binNumber_Psi2S);
    Double_t binError_JPsi = histogram->GetBinError(binNumber_JPsi);
    Double_t binError_Psi2S = histogram->GetBinError(binNumber_Psi2S);
    

    // Close the file
    file->Close();
    
    return make_pair(binContent_Psi2S, binError_Psi2S);
}

void Fit_Analysis() {
    // Directory path where the files are located
    const char* folder_path = "Fit_Histo_New_2";
    
    vector<double> YLucaPtPsi2S = {3383.,7043.,6837.,7268.,2941.,3645.,1871.,1087.};
    vector<double> YLucaPtJPsi = {196221.,380385.,324689.,187944.,132629.,144829.,70523.,27298. };
    vector<double> YLucaPt = {};
    
    
    vector<double> YLucaPtStatJPsi = { 933,1275,1082,728,757,771,449,285} ;
    vector<double> YLucaPtStatPsi2S = {299,404,470,342,227,267,160,103} ;
    vector<double> YLucaPtSystJPsi = {9403,21244,16066,7489,6167,8528,4245,1599} ;
    vector<double> YLucaPtSystPsi2S = {612,1421,1928,1049,942,893,450,207} ;
    
   vector<double> YLucaYPsi2S = {2048.,5224.,6039.,5145.,3245.,663.};
    vector<double> YLucaYJPsi = {71106,217468,286154,255322,173291,55822};
    vector<double> YLucaY = {};
    
    vector<double> YLucaYStatJPsi = {741,990,1451,1006,1215,462} ;
    vector<double> YLucaYStatPsi2S = {372,451,469,403,338,159} ;
    vector<double> YLucaYSystJPsi = {6562,19529,25199,25786,17386,5409} ;
    vector<double> YLucaYSystPsi2S = {693,1019,1280,1688,872,372} ;
    
    vector<double> YLucaRatio = {0.01673, 0.01904, 0.01918, 0.02826, 0.02748, 0.03019, 0.03047, 0.04112};
    vector<double> YLucaRatioStat = {0.00165, 0.00131, 0.00137, 0.00184, 0.00219, 0.00195, 0.00322, 0.00490};
    vector<double> YLucaRatioSyst = {0.00455, 0.00482, 0.00498, 0.00662, 0.00666, 0.00764, 0.00920, 0.00959};
    vector<double> YLucaRatioY = {0.02089, 0.02346, 0.02058, 0.02140, 0.01964, 0.01444};
    vector<double> YLucaRatioStatY = {0.00332, 0.00176, 0.00123, 0.00121, 0.00121, 0.00241};
    vector<double> YLucaRatioSystY = {0.00555, 0.00558, 0.00518, 0.00516, 0.00603, 0.00665};
    vector<double> YLucaPtStatRatio, YLucaPtSystRatio,YLucaYStatRatio, YLucaYSystRatio;
    
    for (int j = 0; j < YLucaPtStatJPsi.size() ; ++j) {
        YLucaPt.push_back(YLucaPtPsi2S[j]/YLucaPtJPsi[j]);
    	YLucaPtStatRatio.push_back(YLucaPt[j]*sqrt((YLucaPtStatJPsi[j]/YLucaPtJPsi[j])*(YLucaPtStatJPsi[j]/YLucaPtJPsi[j])+(YLucaPtStatPsi2S[j]/YLucaPtPsi2S[j])*(YLucaPtStatPsi2S[j]/YLucaPtPsi2S[j])));
    	YLucaPtSystRatio.push_back(YLucaPt[j]*sqrt((YLucaPtSystJPsi[j]/YLucaPtJPsi[j])*(YLucaPtSystJPsi[j]/YLucaPtJPsi[j])+(YLucaPtSystPsi2S[j]/YLucaPtPsi2S[j])*(YLucaPtSystPsi2S[j]/YLucaPtPsi2S[j])));
    }
    
     for (int j = 0; j < YLucaYStatJPsi.size() ; ++j) {
        YLucaY.push_back(YLucaYPsi2S[j]/YLucaYJPsi[j]);
    	YLucaYStatRatio.push_back(YLucaY[j]*sqrt((YLucaYStatJPsi[j]/YLucaYJPsi[j])*(YLucaYStatJPsi[j]/YLucaYJPsi[j])+(YLucaYStatPsi2S[j]/YLucaYPsi2S[j])*(YLucaYStatPsi2S[j]/YLucaYPsi2S[j])));
    	YLucaYSystRatio.push_back(YLucaY[j]*sqrt((YLucaYSystJPsi[j]/YLucaYJPsi[j])*(YLucaYSystJPsi[j]/YLucaYJPsi[j])+(YLucaYSystPsi2S[j]/YLucaYPsi2S[j])*(YLucaYSystPsi2S[j]/YLucaYPsi2S[j])));
    	}
    
    vector<double> XLucaPt = {0.5,1.5,2.5,3.5,4.5,6.,8.5,15};
    vector<double> eXLucaPt = {0.5,0.5,0.5,0.5,0.5,1.,1.5,5};
    vector<double> XLucaY = {2.625,2.875,3.125,3.375,3.625,3.875};
    vector<double> eXLucaY = {0.125,0.125,0.125,0.125,0.125,0.125};

    // Loop through all files in the folder
    TSystemDirectory dir("dir", folder_path);
    TList* files = dir.GetListOfFiles();
    files->Sort();
    vector<TString> Fit_name_Vec;
    vector<double> XVector,XVector2,XVector2Y,YVector,eYVector,YJPsi,eYJPsi,YPsi2S,eYPsi2S, MeanVector,eMeanVector,eRMSMeanVector,MeanVectorY,eMeanVectorY,eRMSMeanVectorY, MeanVector_Luca, eMeanVector_Luca;
    int i = 0;
    int j = 0;
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                TString file_path = TString::Format("%s/%s", folder_path, fname.Data());
                pair<double,double> result = printFileContent(file_path.Data());
                pair<double,double> resultJPsi = printFileContentJPsi(file_path.Data());
                pair<double,double> resultPsi2S = printFileContentPsi2S(file_path.Data());
                double ratiou = result.first;
                double errorou = result.second; 
                double sigJPsi = resultJPsi.first;
                double errJPsi = resultJPsi.second;
                double sigPsi2S = resultPsi2S.first;
                double errPsi2S = resultPsi2S.second;
                XVector.push_back(i);
                YJPsi.push_back(sigJPsi);
                eYJPsi.push_back(errJPsi);
                YPsi2S.push_back(sigPsi2S);
                eYPsi2S.push_back(errPsi2S);
                YVector.push_back(ratiou);
                eYVector.push_back(errorou);
                string modifiedString = fname.Data();
                modifiedString = TString(modifiedString).ReplaceAll(".root", "").Data();
    		modifiedString = TString(modifiedString).ReplaceAll("hist_mass_", "").Data();
    		modifiedString = TString(modifiedString).ReplaceAll("all_histo_PairsMuonSEPM_matchedMchMid_", "").Data();
		cout << i << " : " << modifiedString << endl;
                Fit_name_Vec.push_back(modifiedString);
                //cout << ratiou << " , " << errorou << endl;
            	i++;
            	}
            j++;
        }
    }
    
    int startIdx = 12; // Index of the 3rd element (0-based indexing)

    // Calculate the sum of the elements in the section
    for (int j = -1; j < 14; ++j) {
    	double sum = 0.0;
    	double sum_rat = 0.0;
    	double sum_JPsi = 0.0;
    	double sum_Psi2S = 0.0;
    	double err_JPsi = 0.0;
    	double err_Psi2S = 0.0;
    	vector<double> XVector_Temp,YVector_Temp,eYVector_Temp;
    	vector<TString> nameVector_Temp;
    	TString canvasName = TString::Format("canvas_%d", j);
        TCanvas* canvas = new TCanvas(canvasName, canvasName, 800, 600);
   canvas->SetBottomMargin(0.3);
   	if (j==-1) {
   		for (int i = startIdx + j*12; i < startIdx + (j+1)*12 ; ++i) {
        	sum += YVector[i];
        	sum_rat += eYVector[i];
        	sum_JPsi += YJPsi[i];
        	sum_Psi2S += YPsi2S[i];
        	err_JPsi += eYJPsi[i]*eYJPsi[i];
        	err_Psi2S += eYPsi2S[i]*eYPsi2S[i];
        	
        	XVector_Temp.push_back(i);
        	YVector_Temp.push_back(YVector[i]);
        	eYVector_Temp.push_back(eYVector[i]);
        	nameVector_Temp.push_back(Fit_name_Vec[i]);
        	}
        	TGraphErrors* graph_Temp = new TGraphErrors(XVector_Temp.size(),XVector_Temp.data(),YVector_Temp.data(),nullptr,eYVector_Temp.data());
	graph_Temp->SetMarkerStyle(20);
    graph_Temp->SetMarkerSize(1.5);
    for (int i = 0; i < nameVector_Temp.size(); i++) {
        graph_Temp->GetXaxis()->SetBinLabel(7.5*(i+1)+2, nameVector_Temp[i]);
    }
     graph_Temp->GetXaxis()->SetLabelSize(0.02);
     
        	int numElementsInSection = 11 + 1; // Number of elements in the section
        	double mean = sum / numElementsInSection;
	    	double mean_err = sum_rat/numElementsInSection;
	    	double rms_temp = 0;
		 for (int i = startIdx + j*12; i < startIdx + (j+1)*12 ; ++i) {
		 	rms_temp += (YVector[i]-mean)*(YVector[i]-mean);
		 }
		 double rms = sqrt(rms_temp/12);
		 
		 TLine* lineMean = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean,
                                       graph_Temp->GetXaxis()->GetXmax(), mean);
    lineMean->SetLineColor(kBlack);
    lineMean->SetLineWidth(1);
	 
	 TLine* lineMeanPlusRMS = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean + rms,
                                       graph_Temp->GetXaxis()->GetXmax(), mean + rms);
    lineMeanPlusRMS->SetLineColor(kGray);
    lineMeanPlusRMS->SetLineWidth(2);

    // Create a TLine for the mean - rms value
    TLine* lineMeanMinusRMS = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean - rms,
                                        graph_Temp->GetXaxis()->GetXmax(), mean - rms);
    lineMeanMinusRMS->SetLineColor(kGray);
    lineMeanMinusRMS->SetLineWidth(2);
    
    TLine* lineMeanPlusErr = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean + mean_err,
                                       graph_Temp->GetXaxis()->GetXmax(), mean + mean_err);
    lineMeanPlusErr->SetLineColor(kRed);
    lineMeanPlusErr->SetLineWidth(2);

    // Create a TLine for the mean - rms value
    TLine* lineMeanMinusErr = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean - mean_err,
                                        graph_Temp->GetXaxis()->GetXmax(), mean - mean_err);
    lineMeanMinusErr->SetLineColor(kRed);
    lineMeanMinusErr->SetLineWidth(2);

    // Draw the zone on the canvas
    graph_Temp->SetTitle("Ratio per fit type");
    graph_Temp->GetXaxis()->SetTitle("Fit type");
    graph_Temp->GetYaxis()->SetTitle("Ratio");
    graph_Temp->Draw("AP");
    lineMeanPlusRMS->Draw();
    lineMeanMinusRMS->Draw();
    lineMeanPlusErr->Draw();
    lineMeanMinusErr->Draw();
    lineMean->Draw();
		 
		 	
	    // Print the result
	    std::cout << "All integrated in section [" << startIdx + j*12 << " to " << startIdx + (j+1)*12-1 << "]: " << mean << " $\pm$ " << mean_err << " (" << mean_err/mean*100 << "\%)" << " $\pm$ " << rms << " (" << rms/mean*100 << "\%) " << std::endl;
	    
    		}
	    	else if (j<8) {
    	for (int i = startIdx + j*12; i < startIdx + (j+1)*12 ; ++i) {
        	sum += YVector[i];
        	sum_rat += eYVector[i];
        	sum_JPsi += YJPsi[i];
        	sum_Psi2S += YPsi2S[i];
        	err_JPsi += eYJPsi[i]*eYJPsi[i];
        	err_Psi2S += eYPsi2S[i]*eYPsi2S[i];
        	
        	XVector_Temp.push_back(i);
        	YVector_Temp.push_back(YVector[i]);
        	eYVector_Temp.push_back(eYVector[i]);
        	nameVector_Temp.push_back(Fit_name_Vec[i]);
    	}
    	TGraphErrors* graph_Temp = new TGraphErrors(XVector_Temp.size(),XVector_Temp.data(),YVector_Temp.data(),nullptr,eYVector_Temp.data());
	graph_Temp->SetMarkerStyle(20);
    graph_Temp->SetMarkerSize(1.5);
    for (int i = 0; i < nameVector_Temp.size(); i++) {
        graph_Temp->GetXaxis()->SetBinLabel(7.5*(i+1)+2, nameVector_Temp[i]);
    }
     graph_Temp->GetXaxis()->SetLabelSize(0.02);
	
	
	    // Calculate the mean
	    int numElementsInSection = 11 + 1; // Number of elements in the section
	    double mean = sum / numElementsInSection;
	    double mean_err = sum_rat/numElementsInSection;
	    
	    double mean_Luca = sum_Psi2S / sum_JPsi;
	    double mean_err_Luca = sum_Psi2S / sum_JPsi * sqrt( err_JPsi/(sum_JPsi*sum_JPsi) + err_Psi2S/(sum_Psi2S*sum_Psi2S) );
	    
	    double rms_temp = 0;
	 for (int i = startIdx + j*12; i < startIdx + (j+1)*12 ; ++i) {
	 	rms_temp += (YVector[i]-mean)*(YVector[i]-mean);
	 }
	 double rms = sqrt(rms_temp/12);
	 
	 TLine* lineMean = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean,
                                       graph_Temp->GetXaxis()->GetXmax(), mean);
    lineMean->SetLineColor(kBlack);
    lineMean->SetLineWidth(1);
	 
	 TLine* lineMeanPlusRMS = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean + rms,
                                       graph_Temp->GetXaxis()->GetXmax(), mean + rms);
    lineMeanPlusRMS->SetLineColor(kGray);
    lineMeanPlusRMS->SetLineWidth(2);

    // Create a TLine for the mean - rms value
    TLine* lineMeanMinusRMS = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean - rms,
                                        graph_Temp->GetXaxis()->GetXmax(), mean - rms);
    lineMeanMinusRMS->SetLineColor(kGray);
    lineMeanMinusRMS->SetLineWidth(2);
    
    TLine* lineMeanPlusErr = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean + mean_err,
                                       graph_Temp->GetXaxis()->GetXmax(), mean + mean_err);
    lineMeanPlusErr->SetLineColor(kRed);
    lineMeanPlusErr->SetLineWidth(2);

    // Create a TLine for the mean - rms value
    TLine* lineMeanMinusErr = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean - mean_err,
                                        graph_Temp->GetXaxis()->GetXmax(), mean - mean_err);
    lineMeanMinusErr->SetLineColor(kRed);
    lineMeanMinusErr->SetLineWidth(2);

    // Draw the zone on the canvas
    graph_Temp->SetTitle("Ratio per fit type");
    graph_Temp->GetXaxis()->SetTitle("Fit type");
    graph_Temp->GetYaxis()->SetTitle("Ratio");
    graph_Temp->Draw("AP");
    lineMeanPlusRMS->Draw();
    lineMeanMinusRMS->Draw();
    lineMeanPlusErr->Draw();
    lineMeanMinusErr->Draw();
    lineMean->Draw();
	 
			XVector2.push_back(j);
			MeanVector.push_back(mean);
			eMeanVector.push_back(mean_err);
			eRMSMeanVector.push_back(rms);
			
			MeanVector_Luca.push_back(mean_Luca);
			eMeanVector_Luca.push_back(mean_err_Luca);
			
			
	    // Print the result
	    std::cout << "Mean of section [" << startIdx + j*12 << " to " << startIdx + (j+1)*12-1 << "]: " << mean << " $\pm$ " << mean_err << " ( " << mean_err/mean*100 << "\%)" << " $\pm$ " << rms << " (" << rms/mean*100 << "\%) " << std::endl;
	    
	    //std::cout << "Mean of section [" << startIdx + j*12 << " to " << startIdx + (j+1)*12-1 << "] (Luca's method) : " << mean_Luca << " +/- " << mean_err_Luca << std::endl;
		}
		
		else {
		startIdx = 108;
		
	for (int i = startIdx + (j-8)*12; i < startIdx + (j-7)*12 ; ++i) {
        	sum += YVector[i];
        	sum_rat += eYVector[i];
        	XVector_Temp.push_back(i);
        	YVector_Temp.push_back(YVector[i]);
        	eYVector_Temp.push_back(eYVector[i]);
        	nameVector_Temp.push_back(Fit_name_Vec[i]);
    	}
    	TGraphErrors* graph_Temp = new TGraphErrors(XVector_Temp.size(),XVector_Temp.data(),YVector_Temp.data(),nullptr,eYVector_Temp.data());
	graph_Temp->SetMarkerStyle(20);
    graph_Temp->SetMarkerSize(1.5);
    for (int i = 0; i < nameVector_Temp.size(); i++) {
        graph_Temp->GetXaxis()->SetBinLabel(7.5*(i+1)+2, nameVector_Temp[i]);
    }
     graph_Temp->GetXaxis()->SetLabelSize(0.02);
	
	
	
	
	    // Calculate the mean
	    int numElementsInSection = 11 + 1; // Number of elements in the section
	    double mean = sum / numElementsInSection;
	    double mean_err = sum_rat/numElementsInSection;
	    double rms_temp = 0;
	    //cout << "Mean : " << mean << endl; 
	 for (int i = startIdx + (j-8)*12; i < startIdx + (j-7)*12 ; ++i) {
	 	rms_temp += (YVector[i]-mean)*(YVector[i]-mean);
	 	//cout << "( " << YVector[i] << " - " << mean << " )**2 = " << rms_temp << endl;
	 }
	 double rms = sqrt(rms_temp/12);
	 
	 TLine* lineMean = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean,
                                       graph_Temp->GetXaxis()->GetXmax(), mean);
    lineMean->SetLineColor(kBlack);
    lineMean->SetLineWidth(1);
	 
    TLine* lineMeanPlusRMS = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean + rms,
                                       graph_Temp->GetXaxis()->GetXmax(), mean + rms);
    lineMeanPlusRMS->SetLineColor(kGray);
    lineMeanPlusRMS->SetLineWidth(2);

    // Create a TLine for the mean - rms value
    TLine* lineMeanMinusRMS = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean - rms,
                                        graph_Temp->GetXaxis()->GetXmax(), mean - rms);
    lineMeanMinusRMS->SetLineColor(kGray);
    lineMeanMinusRMS->SetLineWidth(2);
    
    TLine* lineMeanPlusErr = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean + mean_err,
                                       graph_Temp->GetXaxis()->GetXmax(), mean + mean_err);
    lineMeanPlusErr->SetLineColor(kRed);
    lineMeanPlusErr->SetLineWidth(2);

    // Create a TLine for the mean - rms value
    TLine* lineMeanMinusErr = new TLine(graph_Temp->GetXaxis()->GetXmin(), mean - mean_err,
                                        graph_Temp->GetXaxis()->GetXmax(), mean - mean_err);
    lineMeanMinusErr->SetLineColor(kRed);
    lineMeanMinusErr->SetLineWidth(2);

    // Draw the zone on the canvas
    graph_Temp->SetTitle("Ratio per fit type");
    graph_Temp->GetXaxis()->SetTitle("Fit type");
    graph_Temp->GetYaxis()->SetTitle("Ratio");
    graph_Temp->Draw("AP");
    lineMeanPlusRMS->Draw();
    lineMeanMinusRMS->Draw();
    lineMeanPlusErr->Draw();
    lineMeanMinusErr->Draw();
    lineMean->Draw();
    
			MeanVectorY.push_back(mean);
			XVector2Y.push_back(j);
			eMeanVectorY.push_back(mean_err);
			eRMSMeanVectorY.push_back(rms);
			
			
	    // Print the result
	    std::cout << "Mean of section [" << startIdx + (j-8)*12 << " to " << startIdx + (j-7)*12 -1 << "]: " << mean << " $\pm$ " << mean_err << " (" << mean_err/mean*100 << "\%) " << " $\pm$ " << rms << " (" << rms/mean*100 << "\%) " << std::endl;
		}
	 }
	double elementToMove = MeanVector[1];
	double elementToMoveErr = eMeanVector[1];
	double elementToMoveErrRMS = eRMSMeanVector[1];
	
	double elementToMove_Luca = MeanVector_Luca[1];
	double elementToMoveErr_Luca = eMeanVector_Luca[1];
	
	double sum_tot = std::accumulate(MeanVector.begin(), MeanVector.end(), 0.0);
    	double mean_tot = sum_tot / MeanVector.size();
    	double sum_tot_err = std::accumulate(eMeanVector.begin(), eMeanVector.end(), 0.0);
    	double mean_tot_err = sum_tot_err / eMeanVector.size();
    	double sum_tot_err_sys = std::accumulate(eRMSMeanVector.begin(), eRMSMeanVector.end(), 0.0);
    	double mean_tot_err_sys = sum_tot_err_sys / eRMSMeanVector.size();
    	
    	//cout << "Average of all ratio : " << mean_tot << " +/- " << mean_tot_err << " ( " << mean_tot_err/mean_tot*100 << " % )" << " +/- " << mean_tot_err_sys << " ( " << mean_tot_err_sys/mean_tot*100 << " % ) " << endl;

        // Erase the element from its current position
        MeanVector.erase(MeanVector.begin() + 1);
        eMeanVector.erase(eMeanVector.begin() + 1);
        eRMSMeanVector.erase(eRMSMeanVector.begin() + 1);
        
        MeanVector_Luca.erase(MeanVector_Luca.begin() + 1);
        eMeanVector_Luca.erase(eMeanVector_Luca.begin() + 1);

        // Insert the element at the new position
        MeanVector.push_back(elementToMove);
        eMeanVector.push_back(elementToMoveErr);
        eRMSMeanVector.push_back(elementToMoveErrRMS);
        
        MeanVector_Luca.push_back(elementToMove_Luca);
        eMeanVector_Luca.push_back(elementToMoveErr_Luca);
    
    //TGraph* graph = new TGraph(XVector.size(),XVector.data(),YVector.data());
    TGraphErrors* graph = new TGraphErrors(XVector.size(),XVector.data(),YVector.data(),nullptr,eYVector.data());
     
     // Customize the appearance of the graph
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.5);
    graph->SetTitle("Ratio per fit type");
    graph->GetXaxis()->SetTitle("Fit function");
    graph->GetYaxis()->SetTitle("Ratio");
    TCanvas * cRatio = new TCanvas("cRatio","Ratio per fit type");
    graph->Draw("AP");
    
    TGraphMultiErrors* gme = new TGraphMultiErrors(XLucaPt.size(), XLucaPt.data(), MeanVector.data(), eXLucaPt.data(), eXLucaPt.data(), eMeanVector.data(), eMeanVector.data());
gme->AddYError(XLucaPt.size(), eRMSMeanVector.data(), eRMSMeanVector.data());
gme->SetTitle("Ratio vs pT");
    gme->GetXaxis()->SetTitle("Pt (GeV/c)");
    gme->GetYaxis()->SetTitle("Ratio");
gme->SetMarkerStyle(20);
gme->GetAttFill(1)->SetFillStyle(0);
gme->GetYaxis()->SetRangeUser(0.001, 0.1);

    TCanvas * cRatio2 = new TCanvas("cRatio2","Ratio vs pT");
    cRatio2->SetLogy();
    
    //TGraphMultiErrors* gme2 = new TGraphMultiErrors(XLucaPt.size(), XLucaPt.data(), YLucaPt.data(), eXLucaPt.data(), eXLucaPt.data(), YLucaPtStatRatio.data(), YLucaPtStatRatio.data());
    TGraphMultiErrors* gme2 = new TGraphMultiErrors(XLucaPt.size(), XLucaPt.data(), YLucaRatio.data(), eXLucaPt.data(), eXLucaPt.data(), YLucaRatioStat.data(), YLucaRatioStat.data());
gme2->AddYError(XLucaPt.size(), YLucaRatioSyst.data(), YLucaRatioSyst.data());
//TGraphMultiErrors* gme2 = new TGraphMultiErrors(XLucaPt.size(), XLucaPt.data(), MeanVector_Luca.data(), eXLucaPt.data(), eXLucaPt.data(), eMeanVector_Luca.data(), eMeanVector_Luca.data());
gme2->SetTitle("Ratio vs pT");
    gme2->GetXaxis()->SetTitle("Pt (GeV/c)");
    gme2->GetYaxis()->SetTitle("Ratio");
gme2->SetMarkerStyle(20);
gme2->SetMarkerColor(kRed);
gme2->SetLineColor(kRed);
gme2->GetAttLine(0)->SetLineColor(kRed);
gme2->GetAttLine(1)->SetLineColor(kRed);
gme2->GetAttFill(1)->SetFillStyle(0);
gme2->GetYaxis()->SetRangeUser(0.01, 0.1);
TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
legend->AddEntry(gme2, "Method 1 (mean then ratio)");
legend->AddEntry(gme, "Method 2 (ratio then mean)");
    gme2->Draw("a p s ; ; 5 s=1");
    gme->Draw("p s ; ; 5 s=1");
    legend->Draw();
    
    TGraphMultiErrors* gmeY = new TGraphMultiErrors(XLucaY.size(), XLucaY.data(), MeanVectorY.data(), eXLucaY.data(), eXLucaY.data(), eMeanVectorY.data(), eMeanVectorY.data());
gmeY->AddYError(XLucaY.size(), eRMSMeanVectorY.data(), eRMSMeanVectorY.data());
gmeY->SetTitle("Ratio vs y");
    gmeY->GetXaxis()->SetTitle("y");
    gmeY->GetYaxis()->SetTitle("Ratio");
gmeY->SetMarkerStyle(20);
gmeY->GetAttFill(1)->SetFillStyle(0);
gmeY->GetYaxis()->SetRangeUser(0.001, 0.1);
    
    TCanvas * cRatio2Y = new TCanvas("cRatio2Y","Ratio vs y");
    cRatio2Y->SetLogy();
    
    //TGraphMultiErrors* gmeY2 = new TGraphMultiErrors(XLucaY.size(), XLucaY.data(), YLucaY.data(), eXLucaY.data(), eXLucaY.data(), YLucaYStatRatio.data(), YLucaYStatRatio.data());
//gmeY2->AddYError(XLucaY.size(), YLucaYSystRatio.data(), YLucaYSystRatio.data());
    TGraphMultiErrors* gmeY2 = new TGraphMultiErrors(XLucaY.size(), XLucaY.data(), YLucaRatioY.data(), eXLucaY.data(), eXLucaY.data(), YLucaRatioStatY.data(), YLucaRatioStatY.data());
gmeY2->AddYError(XLucaY.size(), YLucaRatioSystY.data(), YLucaRatioSystY.data());
gmeY2->SetTitle("Ratio vs y");
    gmeY2->GetXaxis()->SetTitle("y");
    gmeY2->GetYaxis()->SetTitle("Ratio");
gmeY2->SetMarkerStyle(20);
gmeY2->SetMarkerColor(kRed);
gmeY2->SetLineColor(kRed);
gmeY2->GetAttLine(0)->SetLineColor(kRed);
gmeY2->GetAttLine(1)->SetLineColor(kRed);
gmeY2->GetAttFill(1)->SetFillStyle(0);
gmeY2->GetYaxis()->SetRangeUser(0.001, 0.1);

TLegend* legendY = new TLegend(0.7, 0.7, 0.9, 0.9);
legendY->AddEntry(gmeY2, "Method 1 (mean then ratio)");
legendY->AddEntry(gmeY, "Method 2 (ratio then mean)");
    gmeY2->Draw("a p s ; ; 5 s=1");
    gmeY->Draw("p s ; ; 5 s=1");
    legendY->Draw();

    return;
}
