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

void Comparator_Run2() {

    vector<double> ptCenter = {0.5,1.5,2.5,3.5,4.5,6.,8.5,15};
    vector<double> eptCenter = {0.5,0.5,0.5,0.5,0.5,1.,1.5,5};
    vector<double> yCenter = {2.625,2.875,3.125,3.375,3.625,3.875};
    vector<double> eyCenter = {0.125,0.125,0.125,0.125,0.125,0.125};
    
vector<double> pt1 = {0.261996 ,0.251658 ,0.251606 ,0.268117 , 0.307243 , 0.371151 , 0.463477 , 0.545844 };
vector<double> pt2 = {0.283543, 0.275227, 0.270678, 0.275373, 0.295734, 0.34517, 0.431867, 0.545437 };
vector<double> y1 = {0.115446 ,0.316881,0.424771,0.424496,0.305863,0.0875632};
vector<double> y2 = {0.121253,0.3351,0.449883,0.448982,0.328012,0.101495};

vector<double> ept1 = {0.00083838,0.000593484,0.000643208,0.000818693,0.00112836,0.00127405,0.00212171,0.00373881};
vector<double> ept2 = {0.00109694,0.000679303,0.000644784,0.000759304,0.000994407,0.00110067,0.00172884,0.00260469};
vector<double> ept1_Sys = {0.000644863,0.000646681,0.000687003,0.00068751,0.000671387,0.000980556,0.00118741,0.00241387};
vector<double> ept2_Sys = {0.00036481,0.000348151,0.000262529,0.000365138,0.000508517,0.00129342,0.000587553,0.00395638
};
vector<double> ey1 = {0.000627714,0.000871793,0.000887063, 0.000851115,0.000765296,0.000453738};
vector<double> ey2 = {
0.000694534,0.000921461,0.000904581,0.000846805,0.000756878,0.00046241};
vector<double> ey1_Sys = {0.00185333,0.00282231,0.00272366,0.00317413,0.00379628,0.00204325};
vector<double> ey2_Sys = {0.00195061,0.00447041,0.00639349,0.00661352,0.00497282,0.00212514};
vector<double> ptRat = {1.0824, 1.09402, 1.07568, 1.02731, 0.963102, 0.930407, 0.931724, 0.999548};
vector<double> eptRat = {0.00541304, 0.00371928, 0.00374408, 0.00420953, 0.00477695, 0.00434068, 0.00564073, 0.00831697};
//vector<double> eptRat_Sys = {0.00284124, 0.00280433, 0.00298496, 0.00276543, 0.00228919, 0.00430644, 0.0035023, 0.0106967};
vector<double> eptRat_Sys = {0.0101158,0.0046995,0.0018817,0.00390058,0.00440724,0.00464119,0.00265542,0.0128326};
vector<double> yRat = {1.0502, 1.05792, 1.05941, 1.05773, 1.07331, 1.15907};
vector<double> eyRat = {0.00826077, 0.00409921, 0.00305841, 0.0028996, 0.00363834, 0.00796757};
//vector<double> eyRat_Sys = {0.0261997, 0.0200692, 0.0204547, 0.0205111, 0.0229221, 0.0361182};
vector<double> eyRat_Sys = {0.011667,0.0110391,0.0121205,0.0120937,0.0107666,0.0133846};

vector<double> pt1_07 = {0.261455,0.250202,0.244132 , 0.253834 ,0.294758,0.360442,0.455924,0.541018};
vector<double> pt2_07 = {0.283197,0.274765,0.269787,0.271976,0.286951,0.335671,0.424423,0.539172};
vector<double> y1_07 = {0.111725,0.310393,0.417912,0.415556,0.295949,0.084689 };
vector<double> y2_07 = {0.119066,0.332235,0.446764,0.444921,0.323107,0.0995201};
vector<double> ept1_07 = {0.0008344,0.000589941,0.000634173,0.000801098,0.00111077,0.00126144,0.00211078,0.0037277};
vector<double> ept2_07 = {0.00109227,0.000676186,0.000641533,0.000753375,0.000981606,0.00108857,0.0017179,0.00259708};
vector<double> ept1_07_Sys = {0.000641912,0.0006513,0.000661686,0.000668883,0.000670582,0.000970006,0.00122747};
vector<double> ept2_07_Sys = {0.000366271,0.000349194,0.000263304,0.000335157,0.000466406,0.00129926,0.000597656,0.00410563};
vector<double> ey1_07 = {0.000616365,0.000863417,0.000881594,0.00084516,0.000755124,0.000445144};
vector<double> ey2_07 = {0.000686339,0.000915771,0.00090032,0.000842651,0.000750966,0.000456484};
vector<double> ey1_07_Sys = {0.00165317,0.00252172,0.0023993,0.00278705,0.00338683,0.00186267};
vector<double> ey2_07_Sys = {0.0017591,0.00424889,0.00608671,0.00629592,0.00470012,0.00192892};
vector<double> ptRat_07 = {1.08316, 1.09817, 1.10509, 1.07147, 0.973513, 0.931277, 0.930907, 0.996586};
vector<double> eptRat_07 = {0.00542235, 0.00374279, 0.00389179, 0.00449932, 0.00495468, 0.00444334, 0.00572467, 0.00837819};
//vector<double> eptRat_07_Sys = {0.00283633, 0.00285122, 0.00304358, 0.00296241, 0.00240605, 0.00444795, 0.00366058, 0.0110392};
vector<double> eptRat_07_Sys = {0.00977466,0.00468989,0.00206231,0.00389188,0.00440403,0.00476381,0.00271247,0.0130308};
vector<double> yRat_07 = {1.06571, 1.07037, 1.06904, 1.07066, 1.09176, 1.17512};
vector<double> eyRat_07 = {0.00850314, 0.00419162, 0.00311879, 0.00297548, 0.00376813, 0.00819788};
//vector<double> eyRat_07_Sys = {0.0252971, 0.0198075, 0.0200562, 0.0202127, 0.0228082, 0.0349148};
vector<double> eyRat_07_Sys = {0.011667,0.0110391,0.0121205,0.0120937,0.0107666,0.0133846};

vector<double> pt1_real = {0.235604,0.225407,0.224796,0.24298,0.279811,0.339583,0.42148,0.506932};
vector<double> pt2_real = {0.257256,0.250147,0.244888,0.250062,0.269573,0.315926,0.397601,0.494803};
vector<double> ept1_real = {0.000190169,0.000134322,0.000145495,0.000186559,0.000257783,0.000293325,0.00049079,0.000885531};
vector<double> ept2_real = {0.000250147,0.000154894,0.000146907,0.000172677,0.00022756,0.000251975,0.000401683,0.00061411};
vector<double> ept1_real_Sys = {0.00183755,0.00284187,0.00257169,0.00128973,0.000810444,0.000763088,0.000951529,0.00225586};
vector<double> ept2_real_Sys = {0.00112823,0.00162362,0.00158735,0.00172585,0.00129369,0.00112996,0.000556678,0.00510213};
vector<double> y1_real = {0.107963,0.284784,0.376185,0.378999,0.278639,0.0837987};
vector<double> y2_real = {0.115078,0.302324,0.398504,0.403105,0.304001,0.0982866};
vector<double> ey1_real = {0.000143035,0.000198507,0.000204283,0.000196547,0.000175109,0.00010465};
vector<double> ey2_real = {0.000159634,0.000210754,0.000209018,0.000196507,0.000174164,0.000107139};
vector<double> ey1_real_Sys = {0.00539406,0.0191261,0.0287115,0.0334479,0.0296458,0.00964544};
vector<double> ey2_real_Sys = {0.00714706,0.0216293,0.0295715,0.0301064,0.0260516,0.00971732};
vector<double> ptRat_real = {1.0919, 1.10975, 1.08938, 1.02915, 0.96341, 0.930333, 0.943345, 0.976075};
vector<double> eptRat_real = {0.00137986, 0.000953696, 0.000961358, 0.00106274, 0.00120381, 0.00109378, 0.00145427, 0.00209159};
//vector<double> eptRat_real_Sys = {0.0132017, 0.019695, 0.0191587, 0.0124992, 0.00739244, 0.00435867, 0.00343746, 0.0132805};
vector<double> eptRat_real_Sys = {0.0105948,0.0133405,0.0123216,0.00801935,0.00567421,0.00488565,0.00291796,0.0139345};
vector<double> yRat_real = {1.0659, 1.06159, 1.05933, 1.0636, 1.09102, 1.17289};
vector<double> eyRat_real = {0.00204461, 0.00104653, 0.000799777, 0.000757015, 0.00092779, 0.00194425};
//vector<double> eyRat_real_Sys = {0.125735, 0.157847, 0.160287, 0.159249, 0.189421, 0.206757};
vector<double> eyRat_real_Sys = {0.102088,0.118235,0.136256,0.132213,0.174871,0.177413};

vector<double> pt1_real_07 = {0.238256,0.226886,0.221263,0.233359,0.271588,0.334074,0.418918,0.508534};
vector<double> pt2_real_07 = {0.260022,0.252617,0.246817,0.250277,0.264722,0.311318,0.394913,0.495348};
vector<double> ept1_real_07 = {0.000188343,0.000132836,0.000142696,0.000181452,0.000251889,0.00028803,0.000484359, 0.000873737};
vector<double> ept2_real_07 = {0.000247741,0.000153381,0.000145286,0.00017043,0.000223256,0.000247461,0.000395965,0.00060659};
vector<double> ept1_real_07_Sys = {0.00194296,0.00292867,0.00249518,0.00126211,0.0008021,0.000777483,0.00094834,0.00223631};
vector<double> ept2_real_07_Sys = {0.00117686,0.00159219,0.00170399,0.00172755,0.00130445,0.00114912,0.000582159,0.00525516};
vector<double> y1_real_07 = {0.105295,0.281575,0.3734,0.376721,0.274417,0.0827011};
vector<double> y2_real_07 = {0.113718,0.302448,0.400228,0.404472,0.303732,0.0980973};
vector<double> ey1_real_07 = {0.000139554,0.000195215,0.000201243,0.000193582,0.000171926,0.000102604};
vector<double> ey2_real_07 = {0.000156717,0.00020812,0.000206398,0.000193991,0.000171773,0.000105575};
vector<double> ey1_real_07_Sys = {0.00533754,0.019334,0.0299695,0.0347045,0.030133,0.00978906};
vector<double> ey2_real_07_Sys = {0.00723777,0.0220253,0.0302652,0.0292021,0.0263622,0.0104931};
vector<double> ptRat_real_07 = {1.09135, 1.11341, 1.11549, 1.0725, 0.974721, 0.931886, 0.942698, 0.97407};
vector<double> eptRat_real_07 = {0.00135111, 0.00093912, 0.000974009, 0.00110853, 0.00122189, 0.0010928, 0.00144272, 0.00205517};
//vector<double> eptRat_real_07_Sys = {0.0137275, 0.0205572, 0.020512, 0.0131325, 0.00765914, 0.00452786, 0.00351578, 0.0135612};
vector<double> eptRat_real_07_Sys = {0.0103876,0.013741,0.0123661,0.00838845,0.00583551,0.00505432,0.00290593,0.014071};
vector<double> yRat_real_07 = {1.07999, 1.07413, 1.07185, 1.07367, 1.10683, 1.18617};
vector<double> eyRat_real_07 = {0.00206496, 0.00104922, 0.000799525, 0.00075469, 0.000934177, 0.00194817};
//vector<double> eyRat_real_07_Sys = {0.130078, 0.163179, 0.157311, 0.161544, 0.196074, 0.220399};
vector<double> eyRat_real_07_Sys = {0.125027,0.133109,0.145226,0.181754,0.184137};


vector<double> ptRatio, ptRatio_07, ptRatio_real, ptRatio_real_07, eptRatio, eptRatio_Sys, eptRatio_07, eptRatio_07_Sys, eptRatio_real, eptRatio_real_Sys, eptRatio_real_07, eptRatio_real_07_Sys, yRatio, yRatio_07, yRatio_real, yRatio_real_07, eyRatio, eyRatio_Sys, eyRatio_07, eyRatio_07_Sys, eyRatio_real, eyRatio_real_Sys, eyRatio_real_07, eyRatio_real_07_Sys, ptDoubleRatio, eptDoubleRatio, eptDoubleRatio_Sys, ptDoubleRatio_07, eptDoubleRatio_07, eptDoubleRatio_07_Sys, yDoubleRatio, eyDoubleRatio, eyDoubleRatio_Sys, yDoubleRatio_07, eyDoubleRatio_07, eyDoubleRatio_07_Sys;

for (int i=0; i<pt1.size();i++){
	ptRatio.push_back(pt2[i]/pt1[i]);
	ptRatio_07.push_back(pt2_07[i]/pt1_07[i]);
	ptRatio_real.push_back(pt2_real[i]/pt1_real[i]);
	ptRatio_real_07.push_back(pt2_real_07[i]/pt1_real_07[i]);
	
	ptDoubleRatio.push_back(pt1[i]/pt2[i]*pt2_real[i]/pt1_real[i]);
	ptDoubleRatio_07.push_back(pt2_real_07[i]/pt1_real_07[i]*pt1_07[i]/pt2_07[i]);
	
	eptRatio.push_back(pt2[i]/pt1[i]*sqrt( (ept2[i]/pt2[i])*(ept2[i]/pt2[i]) + (ept1[i]/pt1[i])*(ept1[i]/pt1[i]) ));
	eptRatio_Sys.push_back(pt2[i]/pt1[i]*sqrt( (ept2_Sys[i]/pt2[i])*(ept2_Sys[i]/pt2[i]) + (ept1_Sys[i]/pt1[i])*(ept1_Sys[i]/pt1[i]) ));
	eptRatio_07.push_back(pt2_07[i]/pt1_07[i]*sqrt( (ept2_07[i]/pt2_07[i])*(ept2_07[i]/pt2_07[i]) + (ept1_07[i]/pt1_07[i])*(ept1_07[i]/pt1_07[i]) ));
	eptRatio_07_Sys.push_back(pt2_07[i]/pt1_07[i]*sqrt( (ept2_07_Sys[i]/pt2_07[i])*(ept2_07_Sys[i]/pt2_07[i]) + (ept1_07_Sys[i]/pt1_07[i])*(ept1_07_Sys[i]/pt1_07[i]) ));
	eptRatio_real.push_back(pt2_real[i]/pt1_real[i]*sqrt( (ept2_real[i]/pt2_real[i])*(ept2_real[i]/pt2_real[i]) + (ept1_real[i]/pt1_real[i])*(ept1_real[i]/pt1_real[i]) ));
	eptRatio_real_Sys.push_back(pt2_real[i]/pt1_real[i]*sqrt( (ept2_real_Sys[i]/pt2_real[i])*(ept2_real_Sys[i]/pt2_real[i]) + (ept1_real_Sys[i]/pt1_real[i])*(ept1_real_Sys[i]/pt1_real[i]) ));
	eptRatio_real_07.push_back(pt2_real_07[i]/pt1_real_07[i]*sqrt( (ept2_real_07[i]/pt2_real_07[i])*(ept2_real_07[i]/pt2_real_07[i]) + (ept1_real_07[i]/pt1_real_07[i])*(ept1_real_07[i]/pt1_real_07[i]) ));
	eptRatio_real_07_Sys.push_back(pt2_real_07[i]/pt1_real_07[i]*sqrt( (ept2_real_07_Sys[i]/pt2_real_07[i])*(ept2_real_07_Sys[i]/pt2_real_07[i]) + (ept1_real_07_Sys[i]/pt1_real_07[i])*(ept1_real_07_Sys[i]/pt1_real_07[i]) ));
	
	eptDoubleRatio.push_back(pt1[i]/pt2[i]*pt2_real[i]/pt1_real[i]*sqrt((ept2[i]/pt2[i])*(ept2[i]/pt2[i]) + (ept1[i]/pt1[i])*(ept1[i]/pt1[i]) + (ept2_real[i]/pt2_real[i])*(ept2_real[i]/pt2_real[i]) + (ept1_real[i]/pt1_real[i])*(ept1_real[i]/pt1_real[i])) );
	eptDoubleRatio_07.push_back(pt2_real_07[i]/pt1_real_07[i]*pt1_07[i]/pt2_07[i]*sqrt((ept2_07[i]/pt2_07[i])*(ept2_07[i]/pt2_07[i]) + (ept1_07[i]/pt1_07[i])*(ept1_07[i]/pt1_07[i]) + (ept2_real_07[i]/pt2_real_07[i])*(ept2_real_07[i]/pt2_real_07[i]) + (ept1_real_07[i]/pt1_real_07[i])*(ept1_real_07[i]/pt1_real_07[i])) );
	eptDoubleRatio_Sys.push_back(pt1[i]/pt2[i]*pt2_real[i]/pt1_real[i]*sqrt((ept2_Sys[i]/pt2[i])*(ept2_Sys[i]/pt2[i]) + (ept1_Sys[i]/pt1[i])*(ept1_Sys[i]/pt1[i]) + (ept2_real_Sys[i]/pt2_real[i])*(ept2_real_Sys[i]/pt2_real[i]) + (ept1_real_Sys[i]/pt1_real[i])*(ept1_real_Sys[i]/pt1_real[i])) );
	eptDoubleRatio_07_Sys.push_back(pt2_real_07[i]/pt1_real_07[i]*pt1_07[i]/pt2_07[i]*sqrt((ept2_07_Sys[i]/pt2_07[i])*(ept2_07_Sys[i]/pt2_07[i]) + (ept1_07_Sys[i]/pt1_07[i])*(ept1_07_Sys[i]/pt1_07[i]) + (ept2_real_07_Sys[i]/pt2_real_07[i])*(ept2_real_07_Sys[i]/pt2_real_07[i]) + (ept1_real_07_Sys[i]/pt1_real_07[i])*(ept1_real_07_Sys[i]/pt1_real_07[i])) );
}

for (int i=0; i<y1.size();i++){
	yRatio.push_back(y2[i]/y1[i]);
	yRatio_07.push_back(y2_07[i]/y1_07[i]);
	yRatio_real.push_back(y2_real[i]/y1_real[i]);
	yRatio_real_07.push_back(y2_real_07[i]/y1_real_07[i]);
	
	yDoubleRatio.push_back(y1[i]/y2[i]*y2_real[i]/y1_real[i]);
	yDoubleRatio_07.push_back(y2_real_07[i]/y1_real_07[i]*y1_07[i]/y2_07[i]);
	
	eyRatio.push_back(y2[i]/y1[i]*sqrt( (ey2[i]/y2[i])*(ey2[i]/y2[i]) + (ey1[i]/y1[i])*(ey1[i]/y1[i]) ));
	eyRatio_Sys.push_back(y2[i]/y1[i]*sqrt( (ey2_Sys[i]/y2[i])*(ey2_Sys[i]/y2[i]) + (ey1_Sys[i]/y1[i])*(ey1_Sys[i]/y1[i]) ));
	eyRatio_07.push_back(y2_07[i]/y1_07[i]*sqrt( (ey2_07[i]/y2_07[i])*(ey2_07[i]/y2_07[i]) + (ey1_07[i]/y1_07[i])*(ey1_07[i]/y1_07[i]) ));
	eyRatio_07_Sys.push_back(y2_07[i]/y1_07[i]*sqrt( (ey2_07_Sys[i]/y2_07[i])*(ey2_07_Sys[i]/y2_07[i]) + (ey1_07_Sys[i]/y1_07[i])*(ey1_07_Sys[i]/y1_07[i]) ));
	eyRatio_real.push_back(y2_real[i]/y1_real[i]*sqrt( (ey2_real[i]/y2_real[i])*(ey2_real[i]/y2_real[i]) + (ey1_real[i]/y1_real[i])*(ey1_real[i]/y1_real[i]) ));
	eyRatio_real_Sys.push_back(y2_real[i]/y1_real[i]*sqrt( (ey2_real_Sys[i]/y2_real[i])*(ey2_real_Sys[i]/y2_real[i]) + (ey1_real_Sys[i]/y1_real[i])*(ey1_real_Sys[i]/y1_real[i]) ));
	eyRatio_real_07.push_back(y2_real_07[i]/y1_real_07[i]*sqrt( (ey2_real_07[i]/y2_real_07[i])*(ey2_real_07[i]/y2_real_07[i]) + (ey1_real_07[i]/y1_real_07[i])*(ey1_real_07[i]/y1_real_07[i]) ));
	eyRatio_real_07_Sys.push_back(y2_real_07[i]/y1_real_07[i]*sqrt( (ey2_real_07_Sys[i]/y2_real_07[i])*(ey2_real_07_Sys[i]/y2_real_07[i]) + (ey1_real_07_Sys[i]/y1_real_07[i])*(ey1_real_07_Sys[i]/y1_real_07[i]) ));
	
	eyDoubleRatio.push_back(y1[i]/y2[i]*y2_real[i]/y1_real[i]*sqrt((ey2[i]/y2[i])*(ey2[i]/y2[i]) + (ey1[i]/y1[i])*(ey1[i]/y1[i]) + (ey2_real[i]/y2_real[i])*(ey2_real[i]/y2_real[i]) + (ey1_real[i]/y1_real[i])*(ey1_real[i]/y1_real[i])) );
	eyDoubleRatio_07.push_back(y2_real_07[i]/y1_real_07[i]*y1_07[i]/y2_07[i]*sqrt((ey2_07[i]/y2_07[i])*(ey2_07[i]/y2_07[i]) + (ey1_07[i]/y1_07[i])*(ey1_07[i]/y1_07[i]) + (ey2_real_07[i]/y2_real_07[i])*(ey2_real_07[i]/y2_real_07[i]) + (ey1_real_07[i]/y1_real_07[i])*(ey1_real_07[i]/y1_real_07[i])) );
	eyDoubleRatio_Sys.push_back(y1[i]/y2[i]*y2_real[i]/y1_real[i]*sqrt((ey2_Sys[i]/y2[i])*(ey2_Sys[i]/y2[i]) + (ey1_Sys[i]/y1[i])*(ey1_Sys[i]/y1[i]) + (ey2_real_Sys[i]/y2_real[i])*(ey2_real_Sys[i]/y2_real[i]) + (ey1_real_Sys[i]/y1_real[i])*(ey1_real_Sys[i]/y1_real[i])) );
	eyDoubleRatio_07_Sys.push_back(y2_real_07[i]/y1_real_07[i]*y1_07[i]/y2_07[i]*sqrt((ey2_07_Sys[i]/y2_07[i])*(ey2_07_Sys[i]/y2_07[i]) + (ey1_07_Sys[i]/y1_07[i])*(ey1_07_Sys[i]/y1_07[i]) + (ey2_real_07_Sys[i]/y2_real_07[i])*(ey2_real_07_Sys[i]/y2_real_07[i]) + (ey1_real_07_Sys[i]/y1_real_07[i])*(ey1_real_07_Sys[i]/y1_real_07[i])) );
}

TCanvas * cJPsi = new TCanvas("cJPsi","Acceptance efficiency vs Pt (J/Psi)");
	    cJPsi->cd();    
	    cJPsi->SetLogy();
    TGraphMultiErrors* gme = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), pt1.data(), eptCenter.data(), eptCenter.data(), ept1.data(), ept1.data());
gme->AddYError(ept1_Sys.size(), ept1_Sys.data(), ept1_Sys.data());
gme->SetTitle("Acceptance efficiency vs pT (JPsi)");
    gme->GetXaxis()->SetTitle("pT");
    gme->GetYaxis()->SetTitle("AxE");
gme->SetMarkerStyle(20);
gme->SetMarkerColor(kBlue);
gme->SetLineColor(kBlue);
gme->GetAttLine(0)->SetLineColor(kBlue);
gme->GetAttLine(1)->SetLineColor(kBlue);
gme->GetAttFill(1)->SetFillStyle(0);
gme->GetYaxis()->SetRangeUser(0.2, 0.7);

    TGraphMultiErrors* gme_07 = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), pt1_07.data(), eptCenter.data(), eptCenter.data(), ept1_07.data(), ept1_07.data());
gme_07->AddYError(ept1_07_Sys.size(), ept1_07_Sys.data(), ept1_07_Sys.data());
gme_07->SetMarkerStyle(20);
gme_07->SetMarkerColor(kBlack);
gme_07->SetLineColor(kBlack);
gme_07->GetAttLine(0)->SetLineColor(kBlack);
gme_07->GetAttLine(1)->SetLineColor(kBlack);
gme_07->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors* gme_real = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), pt1_real.data(), eptCenter.data(), eptCenter.data(), ept1_real.data(),ept1_real.data());
gme_real->AddYError(ept1_real_Sys.size(), ept1_real_Sys.data(), ept1_real_Sys.data());
gme_real->SetMarkerStyle(20);
gme_real->SetMarkerColor(kRed);
gme_real->SetLineColor(kRed);
gme_real->GetAttLine(0)->SetLineColor(kRed);
gme_real->GetAttLine(1)->SetLineColor(kRed);
gme_real->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors* gme_real_07 = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), pt1_real_07.data(), eptCenter.data(), eptCenter.data(), ept1_real_07.data(), ept1_real_07.data());
gme_real_07->AddYError(ept1_real_07_Sys.size(), ept1_real_07_Sys.data(), ept1_real_07_Sys.data());
gme_real_07->SetMarkerStyle(20);
gme_real_07->SetMarkerColor(kOrange);
gme_real_07->SetLineColor(kOrange);
gme_real_07->GetAttLine(0)->SetLineColor(kOrange);
gme_real_07->GetAttLine(1)->SetLineColor(kOrange);
gme_real_07->GetAttFill(1)->SetFillStyle(0);


    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(gme, "Ideal, no cut");
    legend->AddEntry(gme_07, "Ideal, pt > 0.7");
    legend->AddEntry(gme_real, "Real, no cut");
    legend->AddEntry(gme_real_07, "Real, pt > 0.7");

    gme->Draw("a p s ; ; 5 s=1");
    gme_07->Draw("p s ; ; 5 s=1");
    gme_real->Draw("p s ; ; 5 s=1");
    gme_real_07->Draw("p s ; ; 5 s=1");
    legend->Draw();
    
TCanvas * cPsi2S = new TCanvas("cPsi2S","Acceptance efficiency vs Pt (Psi2S)");
	    cPsi2S->cd();    
	    cPsi2S->SetLogy();
    TGraphMultiErrors* gme2 = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), pt2.data(), eptCenter.data(), eptCenter.data(), ept2.data(), ept2.data());
gme2->AddYError(ept2_Sys.size(), ept2_Sys.data(), ept2_Sys.data());
gme2->SetTitle("Acceptance efficiency vs pT (Psi2S)");
    gme2->GetXaxis()->SetTitle("pT");
    gme2->GetYaxis()->SetTitle("AxE");
gme2->SetMarkerStyle(20);
gme2->SetMarkerColor(kBlue);
gme2->SetLineColor(kBlue);
gme2->GetAttLine(0)->SetLineColor(kBlue);
gme2->GetAttLine(1)->SetLineColor(kBlue);
gme2->GetAttFill(1)->SetFillStyle(0);
gme2->GetYaxis()->SetRangeUser(0.2, 0.7);

    TGraphMultiErrors* gme2_07 = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), pt2_07.data(), eptCenter.data(), eptCenter.data(), ept2_07.data(), ept2_07.data());
gme2_07->AddYError(ept2_07_Sys.size(), ept2_07_Sys.data(), ept2_07_Sys.data());
gme2_07->SetMarkerStyle(20);
gme2_07->SetMarkerColor(kBlack);
gme2_07->SetLineColor(kBlack);
gme2_07->GetAttLine(0)->SetLineColor(kBlack);
gme2_07->GetAttLine(1)->SetLineColor(kBlack);
gme2_07->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors* gme2_real = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), pt2_real.data(), eptCenter.data(), eptCenter.data(), ept2_real.data(),ept2_real.data());
gme2_real->AddYError(ept2_real_Sys.size(), ept2_real_Sys.data(), ept2_real_Sys.data());
gme2_real->SetMarkerStyle(20);
gme2_real->SetMarkerColor(kRed);
gme2_real->SetLineColor(kRed);
gme2_real->GetAttLine(0)->SetLineColor(kRed);
gme2_real->GetAttLine(1)->SetLineColor(kRed);
gme2_real->GetAttFill(1)->SetFillStyle(0);

TGraphMultiErrors* gme2_real_07 = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), pt2_real_07.data(), eptCenter.data(), eptCenter.data(), ept2_real_07.data(), ept2_real_07.data());
gme2_real_07->AddYError(ept2_real_07_Sys.size(), ept2_real_07_Sys.data(), ept2_real_07_Sys.data());
gme2_real_07->SetMarkerStyle(20);
gme2_real_07->SetMarkerColor(kOrange);
gme2_real_07->SetLineColor(kOrange);
gme2_real_07->GetAttLine(0)->SetLineColor(kOrange);
gme2_real_07->GetAttLine(1)->SetLineColor(kOrange);
gme2_real_07->GetAttFill(1)->SetFillStyle(0);


    TLegend* legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2->AddEntry(gme2, "Ideal, no cut");
    legend2->AddEntry(gme2_07, "Ideal, pt > 0.7");
    legend2->AddEntry(gme2_real, "Real, no cut");
    legend2->AddEntry(gme2_real_07, "Real, pt > 0.7");

    gme2->Draw("a p s ; ; 5 s=1");
    gme2_07->Draw("p s ; ; 5 s=1");
    gme2_real->Draw("p s ; ; 5 s=1");
    gme2_real_07->Draw("p s ; ; 5 s=1");
    legend2->Draw();
    
TCanvas * cJPsiY = new TCanvas("cJPsiY","Acceptance efficiency vs Y (J/Psi)");
	    cJPsiY->cd();    
	    cJPsiY->SetLogy();
    TGraphMultiErrors* gmeY = new TGraphMultiErrors(yCenter.size(), yCenter.data(), y1.data(), eyCenter.data(), eyCenter.data(), ey1.data(), ey1.data());
gmeY->AddYError(ey1_Sys.size(), ey1_Sys.data(), ey1_Sys.data());
gmeY->SetTitle("Acceptance efficiency vs Y (JPsi)");
    gmeY->GetXaxis()->SetTitle("y");
    gmeY->GetYaxis()->SetTitle("AxE");
gmeY->SetMarkerStyle(20);
gmeY->SetMarkerColor(kBlue);
gmeY->SetLineColor(kBlue);
gmeY->GetAttLine(0)->SetLineColor(kBlue);
gmeY->GetAttLine(1)->SetLineColor(kBlue);
gmeY->GetAttFill(1)->SetFillStyle(0);
gmeY->GetYaxis()->SetRangeUser(0.05, 1);

    TGraphMultiErrors* gmeY_07 = new TGraphMultiErrors(yCenter.size(), yCenter.data(), y1_07.data(), eyCenter.data(), eyCenter.data(), ey1_07.data(), ey1_07.data());
gmeY_07->AddYError(ey1_07_Sys.size(), ey1_07_Sys.data(), ey1_07_Sys.data());
gmeY_07->SetMarkerStyle(20);
gmeY_07->SetMarkerColor(kBlack);
gmeY_07->SetLineColor(kBlack);
gmeY_07->GetAttLine(0)->SetLineColor(kBlack);
gmeY_07->GetAttLine(1)->SetLineColor(kBlack);
gmeY_07->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors* gmeY_real = new TGraphMultiErrors(yCenter.size(), yCenter.data(), y1_real.data(), eyCenter.data(), eyCenter.data(), ey1_real.data(),ey1_real.data());
gmeY_real->AddYError(ey1_real_Sys.size(), ey1_real_Sys.data(), ey1_real_Sys.data());
gmeY_real->SetMarkerStyle(20);
gmeY_real->SetMarkerColor(kRed);
gmeY_real->SetLineColor(kRed);
gmeY_real->GetAttLine(0)->SetLineColor(kRed);
gmeY_real->GetAttLine(1)->SetLineColor(kRed);
gmeY_real->GetAttFill(1)->SetFillStyle(0);


TGraphMultiErrors* gmeY_real_07 = new TGraphMultiErrors(yCenter.size(), yCenter.data(), y1_real_07.data(), eyCenter.data(), eyCenter.data(), ey1_real_07.data(), ey1_real_07.data());
gmeY_real_07->AddYError(ey1_real_07_Sys.size(), ey1_real_07_Sys.data(), ey1_real_07_Sys.data());
gmeY_real_07->SetMarkerStyle(20);
gmeY_real_07->SetMarkerColor(kOrange);
gmeY_real_07->SetLineColor(kOrange);
gmeY_real_07->GetAttLine(0)->SetLineColor(kOrange);
gmeY_real_07->GetAttLine(1)->SetLineColor(kOrange);
gmeY_real_07->GetAttFill(1)->SetFillStyle(0);


    TLegend* legendY = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendY->AddEntry(gmeY, "Ideal, no cut");
    legendY->AddEntry(gmeY_07, "Ideal, pt > 0.7");
    legendY->AddEntry(gmeY_real, "Real, no cut");
    legendY->AddEntry(gmeY_real_07, "Real, pt > 0.7");

    gmeY->Draw("a p s ; ; 5 s=1");
    gmeY_07->Draw("p s ; ; 5 s=1");
    gmeY_real->Draw("p s ; ; 5 s=1");
    gmeY_real_07->Draw("p s ; ; 5 s=1");
    legendY->Draw();
    
    
TCanvas * cJPsiY2 = new TCanvas("cJPsiY2","Acceptance efficiency vs Y (Psi2S)");
	    cJPsiY2->cd();    
	    cJPsiY2->SetLogy();
    TGraphMultiErrors* gmeY2 = new TGraphMultiErrors(yCenter.size(), yCenter.data(), y2.data(), eyCenter.data(), eyCenter.data(), ey2.data(), ey2.data());
gmeY2->AddYError(ey2_Sys.size(), ey2_Sys.data(), ey2_Sys.data());
gmeY2->SetTitle("Acceptance efficiency vs Y (Psi2S)");
    gmeY2->GetXaxis()->SetTitle("y");
    gmeY2->GetYaxis()->SetTitle("AxE");
gmeY2->SetMarkerStyle(20);
gmeY2->SetMarkerColor(kBlue);
gmeY2->SetLineColor(kBlue);
gmeY2->GetAttLine(0)->SetLineColor(kBlue);
gmeY2->GetAttLine(1)->SetLineColor(kBlue);
gmeY2->GetAttFill(1)->SetFillStyle(0);
gmeY2->GetYaxis()->SetRangeUser(0.05, 1);

    TGraphMultiErrors* gmeY2_07 = new TGraphMultiErrors(yCenter.size(), yCenter.data(), y2_07.data(), eyCenter.data(), eyCenter.data(), ey2_07.data(), ey2_07.data());
gmeY2_07->AddYError(ey2_07_Sys.size(), ey2_07_Sys.data(), ey2_07_Sys.data());
gmeY2_07->SetMarkerStyle(20);
gmeY2_07->SetMarkerColor(kBlack);
gmeY2_07->SetLineColor(kBlack);
gmeY2_07->GetAttLine(0)->SetLineColor(kBlack);
gmeY2_07->GetAttLine(1)->SetLineColor(kBlack);
gmeY2_07->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors* gmeY2_real = new TGraphMultiErrors(yCenter.size(), yCenter.data(), y2_real.data(), eyCenter.data(), eyCenter.data(), ey2_real.data(),ey2_real.data());
gmeY2_real->AddYError(ey2_real_Sys.size(), ey2_real_Sys.data(), ey2_real_Sys.data());
gmeY2_real->SetMarkerStyle(20);
gmeY2_real->SetMarkerColor(kRed);
gmeY2_real->SetLineColor(kRed);
gmeY2_real->GetAttLine(0)->SetLineColor(kRed);
gmeY2_real->GetAttLine(1)->SetLineColor(kRed);
gmeY2_real->GetAttFill(1)->SetFillStyle(0);


TGraphMultiErrors* gmeY2_real_07 = new TGraphMultiErrors(yCenter.size(), yCenter.data(), y2_real_07.data(), eyCenter.data(), eyCenter.data(), ey2_real_07.data(), ey2_real_07.data());
gmeY2_real_07->AddYError(ey2_real_07_Sys.size(), ey2_real_07_Sys.data(), ey2_real_07_Sys.data());
gmeY2_real_07->SetMarkerStyle(20);
gmeY2_real_07->SetMarkerColor(kOrange);
gmeY2_real_07->SetLineColor(kOrange);
gmeY2_real_07->GetAttLine(0)->SetLineColor(kOrange);
gmeY2_real_07->GetAttLine(1)->SetLineColor(kOrange);
gmeY2_real_07->GetAttFill(1)->SetFillStyle(0);


    TLegend* legendY2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendY2->AddEntry(gmeY2, "Ideal, no cut");
    legendY2->AddEntry(gmeY2_07, "Ideal, pt > 0.7");
    legendY2->AddEntry(gmeY2_real, "Real, no cut");
    legendY2->AddEntry(gmeY2_real_07, "Real, pt > 0.7");

    gmeY2->Draw("a p s ; ; 5 s=1");
    gmeY2_07->Draw("p s ; ; 5 s=1");
    gmeY2_real->Draw("p s ; ; 5 s=1");
    gmeY2_real_07->Draw("p s ; ; 5 s=1");
    legendY2->Draw();
    
TCanvas * cRatio = new TCanvas("cRatio","Ratio of acceptance efficiency vs pT");
	    cRatio->cd();    
	   // cRatio->SetLogy();
    TGraphMultiErrors* gmeRatio = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), ptRat.data(), eptCenter.data(), eptCenter.data(), eptRat.data(), eptRat.data());
gmeRatio->AddYError(eptRat_Sys.size(), eptRat_Sys.data(), eptRat_Sys.data());
gmeRatio->SetTitle("Ratio of Acceptance efficiency vs Pt ");
    gmeRatio->GetXaxis()->SetTitle("pT");
    gmeRatio->GetYaxis()->SetTitle("AxE(Psi2S)/AxE(J/Psi)");
gmeRatio->SetMarkerStyle(20);
gmeRatio->SetMarkerColor(kBlue);
gmeRatio->SetLineColor(kBlue);
gmeRatio->GetAttLine(0)->SetLineColor(kBlue);
gmeRatio->GetAttLine(1)->SetLineColor(kBlue);
gmeRatio->GetAttFill(1)->SetFillStyle(0);
gmeRatio->GetYaxis()->SetRangeUser(0.9, 1.15);

    TGraphMultiErrors* gmeRatio_07 = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), ptRat_07.data(), eptCenter.data(), eptCenter.data(), eptRat_07.data(), eptRat_07.data());
gmeRatio_07->AddYError(eptRat_07_Sys.size(), eptRat_07_Sys.data(), eptRat_07_Sys.data());
gmeRatio_07->SetMarkerStyle(20);
gmeRatio_07->SetMarkerColor(kBlack);
gmeRatio_07->SetLineColor(kBlack);
gmeRatio_07->GetAttLine(0)->SetLineColor(kBlack);
gmeRatio_07->GetAttLine(1)->SetLineColor(kBlack);
gmeRatio_07->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors* gmeRatio_real = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), ptRat_real.data(), eptCenter.data(), eptCenter.data(), eptRat_real.data(), eptRat_real.data());
gmeRatio_real->AddYError(eptRat_real_Sys.size(), eptRat_real_Sys.data(), eptRat_real_Sys.data());
gmeRatio_real->SetMarkerStyle(20);
gmeRatio_real->SetMarkerColor(kRed);
gmeRatio_real->SetLineColor(kRed);
gmeRatio_real->GetAttLine(0)->SetLineColor(kRed);
gmeRatio_real->GetAttLine(1)->SetLineColor(kRed);
gmeRatio_real->GetAttFill(1)->SetFillStyle(0);


TGraphMultiErrors* gmeRatio_real_07 = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), ptRat_real_07.data(), eptCenter.data(), eptCenter.data(), eptRat_real_07.data(), eptRat_real_07.data());
gmeRatio_real_07->AddYError(eptRat_real_07_Sys.size(), eptRat_real_07_Sys.data(), eptRat_real_07_Sys.data());
gmeRatio_real_07->SetMarkerStyle(20);
gmeRatio_real_07->SetMarkerColor(kOrange);
gmeRatio_real_07->SetLineColor(kOrange);
gmeRatio_real_07->GetAttLine(0)->SetLineColor(kOrange);
gmeRatio_real_07->GetAttLine(1)->SetLineColor(kOrange);
gmeRatio_real_07->GetAttFill(1)->SetFillStyle(0);


    TLegend* legendRatio = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendRatio->AddEntry(gmeRatio, "Ideal, no cut");
    legendRatio->AddEntry(gmeRatio_07, "Ideal, pt > 0.7");
    legendRatio->AddEntry(gmeRatio_real, "Real, no cut");
    legendRatio->AddEntry(gmeRatio_real_07, "Real, pt > 0.7");

    gmeRatio->Draw("a p s ; ; 5 s=1");
    gmeRatio_07->Draw(" p s ; ; 5 s=1");
    gmeRatio_real->Draw("p s ; ; 5 s=1");
    gmeRatio_real_07->Draw("p s ; ; 5 s=1");
    legendRatio->Draw();
    
    
TCanvas * cRatioY = new TCanvas("cRatioY","Ratio of acceptance efficiency vs Y");
	    cRatioY->cd();    
	   // cRatio->SetLogy();
    TGraphMultiErrors* gmeRatioY = new TGraphMultiErrors(yCenter.size(), yCenter.data(), yRat.data(), eyCenter.data(), eyCenter.data(), eyRat.data(), eyRat.data());
gmeRatioY->AddYError(eyRat_Sys.size(), eyRat_Sys.data(), eyRat_Sys.data());
gmeRatioY->SetTitle("Ratio of Acceptance efficiency vs y ");
    gmeRatioY->GetXaxis()->SetTitle("y");
    gmeRatioY->GetYaxis()->SetTitle("AxE(Psi2S)/AxE(J/Psi)");
gmeRatioY->SetMarkerStyle(20);
gmeRatioY->SetMarkerColor(kBlue);
gmeRatioY->SetLineColor(kBlue);
gmeRatioY->GetAttLine(0)->SetLineColor(kBlue);
gmeRatioY->GetAttLine(1)->SetLineColor(kBlue);
gmeRatioY->GetAttFill(1)->SetFillStyle(0);
gmeRatioY->GetYaxis()->SetRangeUser(0.8, 1.4);

    TGraphMultiErrors* gmeRatioY_07 = new TGraphMultiErrors(yCenter.size(), yCenter.data(), yRat_07.data(), eyCenter.data(), eyCenter.data(), eyRat_07.data(), eyRat_07.data());
gmeRatioY_07->AddYError(eyRat_07_Sys.size(), eyRat_07_Sys.data(), eyRat_07_Sys.data());
gmeRatioY_07->SetMarkerStyle(20);
gmeRatioY_07->SetMarkerColor(kBlack);
gmeRatioY_07->SetLineColor(kBlack);
gmeRatioY_07->GetAttLine(0)->SetLineColor(kBlack);
gmeRatioY_07->GetAttLine(1)->SetLineColor(kBlack);
gmeRatioY_07->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors* gmeRatioY_real = new TGraphMultiErrors(yCenter.size(), yCenter.data(), yRat_real.data(), eyCenter.data(), eyCenter.data(), eyRat_real.data(),eyRat_real.data());
gmeRatioY_real->AddYError(eyRatio_real_Sys.size(), eyRat_real_Sys.data(), eyRat_real_Sys.data());
gmeRatioY_real->SetMarkerStyle(20);
gmeRatioY_real->SetMarkerColor(kRed);
gmeRatioY_real->SetLineColor(kRed);
gmeRatioY_real->GetAttLine(0)->SetLineColor(kRed);
gmeRatioY_real->GetAttLine(1)->SetLineColor(kRed);
gmeRatioY_real->GetAttFill(1)->SetFillStyle(0);


TGraphMultiErrors* gmeRatioY_real_07 = new TGraphMultiErrors(yCenter.size(), yCenter.data(), yRat_real_07.data(), eyCenter.data(), eyCenter.data(), eyRat_real_07.data(), eyRat_real_07.data());
gmeRatioY_real_07->AddYError(eyRat_real_07_Sys.size(), eyRat_real_07_Sys.data(), eyRat_real_07_Sys.data());
gmeRatioY_real_07->SetMarkerStyle(20);
gmeRatioY_real_07->SetMarkerColor(kOrange);
gmeRatioY_real_07->SetLineColor(kOrange);
gmeRatioY_real_07->GetAttLine(0)->SetLineColor(kOrange);
gmeRatioY_real_07->GetAttLine(1)->SetLineColor(kOrange);
gmeRatioY_real_07->GetAttFill(1)->SetFillStyle(0);


    TLegend* legendRatioY = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendRatioY->AddEntry(gmeRatioY, "Ideal, no cut");
    legendRatioY->AddEntry(gmeRatioY_07, "Ideal, pT > 0.7");
    legendRatioY->AddEntry(gmeRatioY_real, "Real, no cut");
    legendRatioY->AddEntry(gmeRatioY_real_07, "Real, pT > 0.7");

    gmeRatioY->Draw("a p s ; ; 5 s=1");
    gmeRatioY_07->Draw(" p s ; ; 5 s=1");
    gmeRatioY_real->Draw("p s ; ; 5 s=1");
    gmeRatioY_real_07->Draw("p s ; ; 5 s=1");
    legendRatioY->Draw();
    
    TCanvas * cDoubleRatio = new TCanvas("cDoubleRatio","Double ratio of acceptance efficiency vs pT");
	    cDoubleRatio->cd();    
	   // cRatio->SetLogy();
    TGraphMultiErrors* gmeDoubleRatio = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), ptDoubleRatio.data(), eptCenter.data(), eptCenter.data(), eptDoubleRatio.data(),eptDoubleRatio.data());
gmeDoubleRatio->AddYError(eptDoubleRatio_Sys.size(), eptDoubleRatio_Sys.data(), eptDoubleRatio_Sys.data());
gmeDoubleRatio->SetTitle("Double ratio of Acceptance efficiency vs pT");
    gmeDoubleRatio->GetXaxis()->SetTitle("pT");
    gmeDoubleRatio->GetYaxis()->SetTitle("AxE Real / AxE Ideal");
gmeDoubleRatio->SetMarkerStyle(20);
gmeDoubleRatio->SetMarkerColor(kBlue);
gmeDoubleRatio->SetLineColor(kBlue);
gmeDoubleRatio->GetAttLine(0)->SetLineColor(kBlue);
gmeDoubleRatio->GetAttLine(1)->SetLineColor(kBlue);
gmeDoubleRatio->GetAttFill(1)->SetFillStyle(0);
gmeDoubleRatio->GetYaxis()->SetRangeUser(0.94, 1.04);


    TGraphMultiErrors* gmeDoubleRatio_07 = new TGraphMultiErrors(ptCenter.size(), ptCenter.data(), ptDoubleRatio_07.data(), eptCenter.data(), eptCenter.data(), eptDoubleRatio_07.data(), eptDoubleRatio_07.data());
gmeDoubleRatio_07->AddYError(eptDoubleRatio_07_Sys.size(), eptDoubleRatio_07_Sys.data(), eptDoubleRatio_07_Sys.data());
gmeDoubleRatio_07->SetMarkerStyle(20);
gmeDoubleRatio_07->SetMarkerColor(kBlack);
gmeDoubleRatio_07->SetLineColor(kBlack);
gmeDoubleRatio_07->GetAttLine(0)->SetLineColor(kBlack);
gmeDoubleRatio_07->GetAttLine(1)->SetLineColor(kBlack);
gmeDoubleRatio_07->GetAttFill(1)->SetFillStyle(0);


    TLegend* legendDoubleRatio = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendDoubleRatio->AddEntry(gmeDoubleRatio, "No cut");
    legendDoubleRatio->AddEntry(gmeDoubleRatio_07, "pT > 0.7 GeV");

    gmeDoubleRatio->Draw("a p s ; ; 5 s=1");
    gmeDoubleRatio_07->Draw("p s ; ; 5 s=1");
    legendDoubleRatio->Draw();
    
    
TCanvas * cDoubleRatioY = new TCanvas("cDoubleRatioY","Double ratio of acceptance efficiency vs y ");
	    cDoubleRatioY->cd();    
	   // cRatio->SetLogy();
    TGraphMultiErrors* gmeDoubleRatioY = new TGraphMultiErrors(yCenter.size(), yCenter.data(), yDoubleRatio.data(), eyCenter.data(), eyCenter.data(), eyDoubleRatio.data(),eyDoubleRatio.data());
gmeDoubleRatioY->AddYError(eyDoubleRatio_Sys.size(), eyDoubleRatio_Sys.data(), eyDoubleRatio_Sys.data());
gmeDoubleRatioY->SetTitle("Double ratio of Acceptance efficiency vs y ");
    gmeDoubleRatioY->GetXaxis()->SetTitle("y");
    gmeDoubleRatioY->GetYaxis()->SetTitle("AxE Real / AxE Ideal");
gmeDoubleRatioY->SetMarkerStyle(20);
gmeDoubleRatioY->SetMarkerColor(kBlue);
gmeDoubleRatioY->SetLineColor(kBlue);
gmeDoubleRatioY->GetAttLine(0)->SetLineColor(kBlue);
gmeDoubleRatioY->GetAttLine(1)->SetLineColor(kBlue);
gmeDoubleRatioY->GetAttFill(1)->SetFillStyle(0);
gmeDoubleRatioY->GetYaxis()->SetRangeUser(0.8, 1.2);


    TGraphMultiErrors* gmeDoubleRatioY_07 = new TGraphMultiErrors(yCenter.size(), yCenter.data(), yDoubleRatio_07.data(), eyCenter.data(), eyCenter.data(), eyDoubleRatio_07.data(), eyDoubleRatio_07.data());
gmeDoubleRatioY_07->AddYError(eyDoubleRatio_07_Sys.size(), eyDoubleRatio_07_Sys.data(), eyDoubleRatio_07_Sys.data());
gmeDoubleRatioY_07->SetMarkerStyle(20);
gmeDoubleRatioY_07->SetMarkerColor(kBlack);
gmeDoubleRatioY_07->SetLineColor(kBlack);
gmeDoubleRatioY_07->GetAttLine(0)->SetLineColor(kBlack);
gmeDoubleRatioY_07->GetAttLine(1)->SetLineColor(kBlack);
gmeDoubleRatioY_07->GetAttFill(1)->SetFillStyle(0);


    TLegend* legendDoubleRatioY = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendDoubleRatioY->AddEntry(gmeDoubleRatioY, "No cut");
    legendDoubleRatioY->AddEntry(gmeDoubleRatioY_07, "pT > 0.7 GeV");

    gmeDoubleRatioY->Draw("a p s ; ; 5 s=1");
    gmeDoubleRatioY_07->Draw("p s ; ; 5 s=1");
    legendDoubleRatioY->Draw();

return;
}
