#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include <iostream>
#include <TVirtualFitter.h>
#include "TChain.h"
#include "TDatime.h"
// #include "TForm.h"
#include "TString.h"
#include <fstream>

Double_t PtJpsiPara(Double_t x, Double_t p0, Double_t p1, Double_t p2, Double_t p3)
{
    return p0 * x / TMath::Power(p1 + TMath::Power(x, p2), p3);
}

Double_t YJpsiPara(Double_t x, Double_t p0, Double_t p1)
{ 
  return p0 + p1*x;
}

// Fit function to wrap the PtJpsiPara function with the specified parameters
Double_t fitFunction(Double_t* x, Double_t* par)
{
    Double_t p0 = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];

    return PtJpsiPara(x[0], p0, p1, p2, p3);
}

// Fit function to wrap the YJpsiPara function with the specified parameters
Double_t fitFunctionY(Double_t* x, Double_t* par)
{
    Double_t p0 = par[0];
    Double_t p1 = par[1];

    return YJpsiPara(x[0], p0, p1);
}

using namespace std;
Bool_t isPsiPrime=true;
Bool_t isScale = false;

Double_t scaleX( double x= 0)
{

  double pt_average = (isPsiPrime) ? 3.2 : 2.7;
  x = x/pt_average;
  return x ;
}
void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t))
{
  a->Set( a->GetNbins(),
            Scale(a->GetXmin()), // new Xmin
            Scale(a->GetXmax()) ); // new Xmax
  return;
}
void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t))
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetXaxis(), Scale);
  return;
}

///////////////////////////////////////////////////////////////////////////

void computeEff(TChain &trec, TChain &tsim, TF1 ratioFunc, TF1 ratioFuncPt, TH1F *&hEffY, TH1F *&hEffPt, TH1F *&hSimY, TH1F *&hSimPt, TH1F *&hRecY, TH1F *&hRecPt,Bool_t isweight=false,  int iter=0)

{
  // double ptrec=0;


  double  ptrec, massrec, yrec, ptrecsim, massrecsim, yrecsim;
  // cout << "Theraa"<< ptrec << endl;
  trec.SetBranchAddress("ptRec",&ptrec);
  trec.SetBranchAddress("massRec",&massrec);
  trec.SetBranchAddress("yRec",&yrec);
  trec.SetBranchAddress("ptMC",&ptrecsim);
  trec.SetBranchAddress("massMC",&massrecsim);
  trec.SetBranchAddress("yMC",&yrecsim);

  //Address to the simulated branches
  double ptsim, ysim, masssim;
  tsim.SetBranchAddress("ptMC",&ptsim);
  tsim.SetBranchAddress("yMC",&ysim);
  tsim.SetBranchAddress("massMC",&masssim);

  for (Int_t i = 0; i<trec.GetEntries(); i++)
  {

      trec.GetEntry(i);
      double weight = 1.0;
      if(isweight && iter>0)
      {

          for(int jter=0; jter<iter; jter++)
          {
          	weight *= ratioFunc.Eval(ptrecsim)*ratioFuncPt.Eval(yrecsim);
              //cout << "ERR: you are entring loop" << jter << endl;
              if (i==0) {
              	cout << " Pt : " << ptrec << " , " << ptrecsim << endl; 
              	cout << " Y : " << yrec << " , " << yrecsim << endl; 
              }
          }
      }
      hRecY->Fill( yrec, weight );
      hRecPt->Fill( ptrec, weight );
  }
  // Loop on simulated tree
  for (Int_t i = 0; i<tsim.GetEntries(); i++) {

      tsim.GetEntry(i);
      double weight = 1.0;
      if(isweight && iter>0){

          for(int jter=0; jter<iter; jter++) {
          	weight *= ratioFunc.Eval(ptsim)*ratioFuncPt.Eval(ysim);
            //cout << "ERR: you are entring loop" << jter << endl;
              // weight *= weightY[iMult][jter]->Eval( ysim );
              // weight *= weightPt[iMult][jter]->Eval( ptsim );
          }
      }
      hSimY->Fill( ysim, weight );
      hSimPt->Fill( ptsim, weight );
  hEffY->Reset();
  hEffY->Divide( hRecY, hSimY , 1.0, 1.0, "B");
  hEffPt->Reset();
  hEffPt->Divide( hRecPt, hSimPt, 1.0, 1.0, "B" );
  }
  //cout << "Global ratio : " << hRecY->GetEntries()/hSimY->GetEntries() << endl;
  // checks canvases
  if (isScale) {ScaleXaxis(hEffPt,scaleX);}
  TCanvas *c = new TCanvas(Form("c_%d",iter),Form("checks for tree iteration %d",iter));
  c->Divide(2,3);
  c->cd(1);
  hSimPt->Draw();
  c->cd(3);
  hRecPt->Draw();
  c->cd(2);
  hSimY->Draw();
  c->cd(4);
  hRecY->Draw();
  c->cd(5);
  hEffPt->Draw();
  c->cd(6);
  hEffY->Draw();


} // end of compute eff

///////////////////////////////////////////////////////////////////////////


void accEffpp()
{
  TDatime date;
  TString path = "../AE_CharmAnalysisSept2022/";
  //TString file = "AnalysisResults_AllRuns.root";
  TString file = "AnalysisResults-Ideal-0.root";
  TString periods[13] = {"16h","16k","16p","17h","17k","17l","17o","17r","18c","18d","18e","18l","18p"};
  
vector<double> FuncParamJPsiY0 = { 21534.9 , 14.3287 , 1.75964 , 4.16559 };
vector<double> FuncParamJPsiY1 = { 3539.78 , 16.3812 , 1.96168 , 3.38207 };
vector<double> FuncParamJPsiY2 = { 45044.4 , 14.9258 , 1.72341 , 4.38368 };
vector<double> FuncParamJPsiY3 = { 451990 , 16.0551 , 1.6672 , 5.07906 };

vector<double> FuncParamJPsiPtN = { 0.11638 , -0.0105578 };
vector<double> FuncParamJPsiPt0 = { 0.567244 , -0.0514596 };
vector<double> FuncParamJPsiPt1 = { 0.63685 , -0.072877 };
vector<double> FuncParamJPsiPt2 = { 0.731716 , -0.102066 };
vector<double> FuncParamJPsiPt3 = { 0.755882 , -0.109502 };
vector<double> FuncParamJPsiPt4 = { 0.837982 , -0.134764 };
vector<double> FuncParamJPsiPt5 = { 0.910163 , -0.156973 };
vector<double> FuncParamJPsiPt6 = { 0.957658 , -0.171587 };
vector<double> FuncParamJPsiPt7 = { 0.970288 , -0.175473 };
vector<double> FuncParamJPsiPt8 = { 1.0181 , -0.190184 };
vector<double> FuncParamJPsiPt9 = { 1.02733 , -0.193025 };
vector<double> FuncParamJPsiPt10 = { 1.08792 , -0.211669 };
vector<double> FuncParamJPsiPt11 = { 1.16547 , -0.235529 };
vector<double> FuncParamJPsiPt12 = { 1.10368 , -0.216518 };
vector<double> FuncParamJPsiPt13 = { 1.24902 , -0.261237 };

//vector<double> FuncParamPsi2SY0 = { 164.118 , 37.1172 , 2.82922 , 1.88823 };
vector<double> FuncParamPsi2SY0 = { 163.315, 37.1172, 2.82922, 1.88823};
vector<double> FuncParamPsi2SY1 = { 8221.79 , 19.9294 , 1.87069 , 3.53309 };
vector<double> FuncParamPsi2SY2 = { 42332.2 , 21.2148 , 1.81359 , 3.98966 };
vector<double> FuncParamPsi2SY3 = { 31226.2 , 22.8099 , 1.89369 , 3.80658 };

vector<double> FuncParamPsi2SPtN = { 0.00861902 , -0.00135947 };
vector<double> FuncParamPsi2SPt0 = { 0.820712 , -0.12945 };
vector<double> FuncParamPsi2SPt1 = { 0.912407 , -0.157664 };
vector<double> FuncParamPsi2SPt2 = { 0.90423 , -0.155148 };
vector<double> FuncParamPsi2SPt3 = { 0.933543 , -0.164167 };
vector<double> FuncParamPsi2SPt4 = { 0.95915 , -0.172046 };
vector<double> FuncParamPsi2SPt5 = { 1.00046 , -0.184757 };
vector<double> FuncParamPsi2SPt6 = { 1.07252 , -0.20693 };
vector<double> FuncParamPsi2SPt7 = { 1.08323 , -0.210224 };
vector<double> FuncParamPsi2SPt8 = { 1.12709 , -0.22372 };
vector<double> FuncParamPsi2SPt9 = { 1.07998 , -0.209224 };
vector<double> FuncParamPsi2SPt10 = { 1.22119 , -0.252674 };
vector<double> FuncParamPsi2SPt11 = { 1.12589 , -0.223352 };
  
  TF1* fitFuncJPsi_N = new TF1("fitFuncJPsi_N", fitFunction, 0, 20, 4);
  fitFuncJPsi_N->SetParameters(FuncParamJPsiY0.data());
  TF1* fitFuncJPsi_1 = new TF1("fitFuncJPsi_1", fitFunction, 0, 20, 4);
  fitFuncJPsi_1->SetParameters(FuncParamJPsiY1.data());
  TF1* fitFuncJPsi_2 = new TF1("fitFuncJPsi_2", fitFunction, 0, 20, 4);
  fitFuncJPsi_2->SetParameters(FuncParamJPsiY2.data());
  TF1* fitFuncJPsi_3 = new TF1("fitFuncJPsi_3", fitFunction, 0, 20, 4);
  fitFuncJPsi_3->SetParameters(FuncParamJPsiY3.data());
  
  TF1* fitFuncJPsi_Pt_N = new TF1("fitFuncJPsi_Pt_N", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_N->SetParameters(FuncParamJPsiPtN.data());
  TF1* fitFuncJPsi_Pt_0 = new TF1("fitFuncJPsi_Pt_0", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_0->SetParameters(FuncParamJPsiPt0.data());
  TF1* fitFuncJPsi_Pt_1 = new TF1("fitFuncJPsi_Pt_1", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_1->SetParameters(FuncParamJPsiPt1.data());
  TF1* fitFuncJPsi_Pt_2 = new TF1("fitFuncJPsi_Pt_2", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_2->SetParameters(FuncParamJPsiPt2.data());
  TF1* fitFuncJPsi_Pt_3 = new TF1("fitFuncJPsi_Pt_3", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_3->SetParameters(FuncParamJPsiPt3.data());
  TF1* fitFuncJPsi_Pt_4 = new TF1("fitFuncJPsi_Pt_4", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_4->SetParameters(FuncParamJPsiPt4.data());
  TF1* fitFuncJPsi_Pt_5 = new TF1("fitFuncJPsi_Pt_5", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_5->SetParameters(FuncParamJPsiPt5.data());
  TF1* fitFuncJPsi_Pt_6 = new TF1("fitFuncJPsi_Pt_6", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_6->SetParameters(FuncParamJPsiPt6.data());
  TF1* fitFuncJPsi_Pt_7 = new TF1("fitFuncJPsi_Pt_7", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_7->SetParameters(FuncParamJPsiPt7.data());
  TF1* fitFuncJPsi_Pt_8 = new TF1("fitFuncJPsi_Pt_8", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_8->SetParameters(FuncParamJPsiPt8.data());
  TF1* fitFuncJPsi_Pt_9 = new TF1("fitFuncJPsi_Pt_9", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_9->SetParameters(FuncParamJPsiPt9.data());
  TF1* fitFuncJPsi_Pt_10 = new TF1("fitFuncJPsi_Pt_10", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_10->SetParameters(FuncParamJPsiPt10.data());
  TF1* fitFuncJPsi_Pt_11 = new TF1("fitFuncJPsi_Pt_11", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_11->SetParameters(FuncParamJPsiPt11.data());
  TF1* fitFuncJPsi_Pt_12 = new TF1("fitFuncJPsi_Pt_12", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_12->SetParameters(FuncParamJPsiPt12.data());
  TF1* fitFuncJPsi_Pt_13 = new TF1("fitFuncJPsi_Pt_13", fitFunctionY, 2.0, 4.5, 2);
  fitFuncJPsi_Pt_13->SetParameters(FuncParamJPsiPt13.data());
  
  TF1* fitFuncPsi2S_N = new TF1("fitFuncPsi2S_N", fitFunction, 0, 20, 4);
  fitFuncPsi2S_N->SetParameters(FuncParamPsi2SY0.data());
  TF1* fitFuncPsi2S_1 = new TF1("fitFuncPsi2S_1", fitFunction, 0, 20, 4);
  fitFuncPsi2S_1->SetParameters(FuncParamPsi2SY1.data());
  TF1* fitFuncPsi2S_2 = new TF1("fitFuncPsi2S_2", fitFunction, 0, 20, 4);
  fitFuncPsi2S_2->SetParameters(FuncParamPsi2SY2.data());
  TF1* fitFuncPsi2S_3 = new TF1("fitFuncPsi2S_3", fitFunction, 0, 20, 4);
  fitFuncPsi2S_3->SetParameters(FuncParamPsi2SY3.data());
  
  TF1* fitFuncPsi2S_Pt_N = new TF1("fitFuncPsi2S_Pt_N", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_N->SetParameters(FuncParamPsi2SPtN.data());
  TF1* fitFuncPsi2S_Pt_0 = new TF1("fitFuncPsi2S_Pt_0", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_0->SetParameters(FuncParamPsi2SPt0.data());
  TF1* fitFuncPsi2S_Pt_1 = new TF1("fitFuncPsi2S_Pt_1", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_1->SetParameters(FuncParamPsi2SPt1.data());
  TF1* fitFuncPsi2S_Pt_2 = new TF1("fitFuncPsi2S_Pt_2", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_2->SetParameters(FuncParamPsi2SPt2.data());
  TF1* fitFuncPsi2S_Pt_3 = new TF1("fitFuncPsi2S_Pt_3", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_3->SetParameters(FuncParamPsi2SPt3.data());
  TF1* fitFuncPsi2S_Pt_4 = new TF1("fitFuncPsi2S_Pt_4", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_4->SetParameters(FuncParamPsi2SPt4.data());
  TF1* fitFuncPsi2S_Pt_5 = new TF1("fitFuncPsi2S_Pt_5", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_5->SetParameters(FuncParamPsi2SPt5.data());
  TF1* fitFuncPsi2S_Pt_6 = new TF1("fitFuncPsi2S_Pt_6", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_6->SetParameters(FuncParamPsi2SPt6.data());
  TF1* fitFuncPsi2S_Pt_7 = new TF1("fitFuncPsi2S_Pt_7", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_7->SetParameters(FuncParamPsi2SPt7.data());
  TF1* fitFuncPsi2S_Pt_8 = new TF1("fitFuncPsi2S_Pt_8", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_8->SetParameters(FuncParamPsi2SPt8.data());
  TF1* fitFuncPsi2S_Pt_9 = new TF1("fitFuncPsi2S_Pt_9", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_9->SetParameters(FuncParamPsi2SPt9.data());
  TF1* fitFuncPsi2S_Pt_10 = new TF1("fitFuncPsi2S_Pt_10", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_10->SetParameters(FuncParamPsi2SPt10.data());
  TF1* fitFuncPsi2S_Pt_11 = new TF1("fitFuncPsi2S_Pt_11", fitFunctionY, 2.0, 4.5, 2);
  fitFuncPsi2S_Pt_11->SetParameters(FuncParamPsi2SPt11.data());
  
  TF1 ratioFunc, ratioFuncPt;
	    
  vector<TF1*> fitFuncJPsiVec = {fitFuncJPsi_N, fitFuncJPsi_1, fitFuncJPsi_2, fitFuncJPsi_3};
  vector<TF1*> fitFuncPsi2SVec = {fitFuncPsi2S_N, fitFuncPsi2S_1, fitFuncPsi2S_2, fitFuncPsi2S_3};
  
  vector<TF1*> fitFuncJPsiPtVec = {fitFuncJPsi_Pt_N, fitFuncJPsi_Pt_0, fitFuncJPsi_Pt_1, fitFuncJPsi_Pt_2, fitFuncJPsi_Pt_3, fitFuncJPsi_Pt_4, fitFuncJPsi_Pt_5, fitFuncJPsi_Pt_6, fitFuncJPsi_Pt_7, fitFuncJPsi_Pt_8, fitFuncJPsi_Pt_9, fitFuncJPsi_Pt_10,  fitFuncJPsi_Pt_11, fitFuncJPsi_Pt_12, fitFuncJPsi_Pt_13 };
  vector<TF1*> fitFuncPsi2SPtVec = {fitFuncPsi2S_Pt_N, fitFuncPsi2S_Pt_0, fitFuncPsi2S_Pt_1, fitFuncPsi2S_Pt_2, fitFuncPsi2S_Pt_3, fitFuncPsi2S_Pt_4, fitFuncPsi2S_Pt_5, fitFuncPsi2S_Pt_6, fitFuncPsi2S_Pt_7, fitFuncPsi2S_Pt_8, fitFuncPsi2S_Pt_9, fitFuncPsi2S_Pt_10, fitFuncPsi2S_Pt_11 };
	   
  TString rfileName = (isScale) ? "_PtScaled": "";
  TFile *fout = new TFile("AE_pp13TeV_Charmonia_Ideal_0_New"+rfileName+".root","update");
  TString particle = (isPsiPrime)? "PsiPrime":"Jpsi";
  TDirectoryFile *particleDire = (TDirectoryFile*)fout->mkdir(particle);
  
  TH1F *fhEffPt, *fhEffPtSame, *fhEffY;
  TH1F *fhRecPt, *fhRecPtSame, *fhRecY;
  TH1F *fhSimPt, *fhSimPtSame, *fhSimY;
  
  int PtIter = 0;
  
  for (int i=0; i<fitFuncJPsiVec.size();i++){
    if (!isPsiPrime) {
    	// Define the ratio function using a lambda function
    	PtIter = fitFuncJPsiPtVec.size();
    	auto ratioFunctionJPsi = [&fitFuncJPsiVec,i](Double_t* x, Double_t* par) {
        	return fitFuncJPsiVec[i]->EvalPar(x) / fitFuncJPsiVec[0]->EvalPar(x);
    	}; 
	ratioFunc = TF1("ratioFuncJPsi", ratioFunctionJPsi, 0, 14, 0);
    }
    else {
    	// Define the ratio function using a lambda function
    	PtIter = fitFuncPsi2SPtVec.size();
    	auto ratioFunctionPsi2S = [&fitFuncPsi2SVec,i](Double_t* x, Double_t* par) {
        	return fitFuncPsi2SVec[i]->EvalPar(x) / fitFuncPsi2SVec[0]->EvalPar(x);
    	}; 
	ratioFunc = TF1("ratioFuncPsi2S", ratioFunctionPsi2S, 0, 14, 0);
    }
    
    for (int j=0; j<PtIter; j++){
    
    if (!isPsiPrime) {
    	auto ratioFunctionJPsiPt = [&fitFuncJPsiPtVec,j](Double_t* x, Double_t* par) {
        	return fitFuncJPsiPtVec[j]->EvalPar(x) / fitFuncJPsiPtVec[0]->EvalPar(x);
    	}; 
	ratioFuncPt = TF1("ratioFuncJPsiPt", ratioFunctionJPsiPt, 2.0, 4.5, 0);
    }
    else {
    	auto ratioFunctionPsi2SPt = [&fitFuncPsi2SPtVec,j](Double_t* x, Double_t* par) {
        	return fitFuncPsi2SPtVec[j]->EvalPar(x) / fitFuncPsi2SPtVec[0]->EvalPar(x);
    	}; 
	ratioFuncPt = TF1("ratioFuncPsi2SPt", ratioFunctionPsi2SPt, 2.0, 4.5, 0);
    }

  // cout << particle << endl;
  // return ;

  TChain chainReco("Charmonia"+particle+"/RecoCandTree"+particle);
  TChain chainSim("Charmonia"+particle+"/SimCandTree"+particle);

  for( int per=0; per< 13; per++)
    {
    // cout << path+periods[per]+file << endl;
    //chainReco.Add(path+periods[per]+file);
    //chainSim.Add(path+periods[per]+file);
    //chainReco.Add(file);
    //chainSim.Add(file);
    }
  chainReco.Add(file);
  chainSim.Add(file);

  TString x_axis_title = (isScale) ? "/<p_{T}>" : "";

  Double_t ptbins = {0};
  //Double_t ptbins[19] = {0.,0.5,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,18.,20.,22.,26.,30.};
  //Double_t ptbins2[14] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,20.};
  Double_t ybins[7] = {2.5,2.75,3.,3.25,3.5,3.75,4.};
  Double_t ptbinou = 0;
  int nbinou = 0;
  if (!isPsiPrime) {
  cout << "J/Psi : Y = " << i << ", Pt = " << j << endl;
  	//ptbinou = ptbins;
  	//Double_t ptbins[14] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,20.};
  	//nbinou = 13;
  	//Double_t ptbins[19] = {0.,0.5,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,18.,20.,22.,26.,30.};
  	//nbinou = 18;
  	//fhEffPt = new TH1F("fhEffPt","fhEffPt"+particle+";p_{T}"+ x_axis_title+";Eff",nbinou,ptbins);
  	//fhRecPt = new TH1F("fhRecPt","fhRecPt"+particle+";pt;Rec",nbinou,ptbins);
	//fhSimPt = new TH1F("fhSimPt","fhSimPt"+particle+";pt;Sim",nbinou,ptbins);
	fhEffPt = new TH1F(TString::Format("fhEffPt%d-%d",i,j),"fhEffPt"+particle+";p_{T}"+ x_axis_title+";Eff"+i,300,0,30);
  	fhRecPt = new TH1F(TString::Format("fhRecPt%d-%d",i,j),"fhRecPt"+particle+";pt;Rec"+i,300,0,30);
  	fhSimPt = new TH1F(TString::Format("fhSimPt%d-%d",i,j),"fhSimPt"+particle+";pt;Sim"+i,300,0,30);
	
	//fhEffPtSame = new TH1F("fhEffPtSame","fhEffPt"+particle+";p_{T}"+ x_axis_title+";Eff;Same",13,ptbinssame);
  	//fhRecPtSame = new TH1F("fhRecPtSame","fhRecPt"+particle+";pt;Rec;Same",13,ptbinssame);
	//fhSimPtSame = new TH1F("fhSimPtSame","fhSimPt"+particle+";pt;Sim;Same",13,ptbinssame);
  }
  else {
  cout << "Psi(2S) : Y = " << i << ", Pt = " << j << endl;
  	//ptbinou = ptbins2;
  	//Double_t ptbins[14] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,20.};
  	//nbinou = 13;
  	//fhEffPt = new TH1F("fhEffPt","fhEffPt"+particle+";p_{T}"+ x_axis_title+";Eff",nbinou,ptbins);
  	//fhRecPt = new TH1F("fhRecPt","fhRecPt"+particle+";pt;Rec",nbinou,ptbins);
	//fhSimPt = new TH1F("fhSimPt","fhSimPt"+particle+";pt;Sim",nbinou,ptbins);
	fhEffPt = new TH1F(TString::Format("fhEffPt%d-%d",i,j),"fhEffPt"+particle+";p_{T}"+ x_axis_title+";Eff"+i,200,0,20);
  	fhRecPt = new TH1F(TString::Format("fhRecPt%d-%d",i,j),"fhRecPt"+particle+";pt;Rec"+i,200,0,20);
  	fhSimPt = new TH1F(TString::Format("fhSimPt%d-%d",i,j),"fhSimPt"+particle+";pt;Sim"+i,200,0,20);
  }
  //fhEffPt = new TH1F("fhEffPt","fhEffPt"+particle+";p_{T}"+ x_axis_title+";Eff",200,0,20);
  //fhRecPt = new TH1F("fhRecPt","fhRecPt"+particle+";pt;Rec",200,0,20);
  //fhSimPt = new TH1F("fhSimPt","fhSimPt"+particle+";pt;Sim",200,0,20);
  fhEffY = new TH1F(TString::Format("fhEffY%d-%d",i,j),"fhEffY"+particle+";Y;Eff"+i,6,-4,-2.5);
  //fhEffY = new TH1F("fhEffY","fhEffY"+particle+";Y;Eff",7,ybins);
  fhRecY = new TH1F(TString::Format("fhRecY%d-%d",i,j),"fhRecY"+particle+";Y;Reco"+i,6,-4,-2.5);
  //fhRecY = new TH1F("fhRecY","fhRecY"+particle+";Y;Reco",7,ybins);
  fhSimY = new TH1F(TString::Format("fhSimY%d-%d",i,j),"fhSimY"+particle+";Y;Sim"+i,6,-4,-2.5);
  //fhSimY = new TH1F("fhSimY","fhSimY"+particle+";Y;Sim",7,ybins);

  computeEff(chainReco,chainSim,ratioFunc,ratioFuncPt,fhEffY,fhEffPt,fhSimY,fhSimPt,fhRecY,fhRecPt,true,1);
  //computeEff(chainReco,chainSim,fhEffY,fhEffPtSame,fhSimY,fhSimPtSame,fhRecY,fhRecPtSame,false,1);
  
  fout->cd();
  particleDire->cd();
  fhEffPt->Write();
  fhRecPt->Write();
  fhSimPt->Write();
  //if (!isPsiPrime) {
  //	fhEffPtSame->Write();
  //	fhRecPtSame->Write();
  //	fhSimPtSame->Write();
  //	}
  fhEffY->Write();
  fhRecY->Write();
  fhSimY->Write();
  
  }
  }

  fout->Close();
} // end of accEff
///////////////////////////////////////////////////////////////////////////
