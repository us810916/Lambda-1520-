/*******************************************************************************
	This program is for analysis of g12 data for the K+ pi+ pi- skim.
	U. Shrestha, May. 2019	

	For the new standard g12 skimmed data
*******************************************************************************/

#include "TRint.h"
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <TVector3.h>
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include <stdlib.h>
#include "g12_corrections.hpp"

#include "/home/shrestha/project/utsav-g12/kppippim/kppippim_g12data/commonToData&MC/TreeFunctions_timingCut.C"	//to make PID timing cuts
#include "/home/shrestha/project/utsav-g12/kppippim/kppippim_g12data/commonToData&MC/TreeFunctions_timingCut.h"

//const Float_t vlight=30.0;

#define PI 3.14159265358979323846
#define WBINS 12
#define WMIN 2.05	//minimum CM enegy or W value or sqrt(s_value) considered
#define WMAX 3.25	//maximum CM energy or W value or sqrt(s_value) considered
#define WBINSIZE ((WMAX-WMIN)/WBINS) 	// size of each W bin

#define CBINS 18
#define CMIN -0.9
#define CMAX 0.9
#define CBINSIZE ((CMAX-CMIN)/CBINS)

#define WCMIN 1.44
#define WCMAX 1.6
#define WCBIN 100 

#define minsigp 1.15
#define maxsigp 1.25

#define minsigm 1.15
#define maxsigm 1.25

#define sigm 1.197
#define sigp 1.189

using namespace std;

/*
//MomCut function for mom cut using mc
float MomCut(int sect, float p, int part)
{
	float y = -9.999;
	if(part == 0)
	{
		if(sect == 1)	y = 1.487*exp(-1.33*p);
		if(sect == 2)	y = 1.4849*exp(-1.32*p);
		if(sect == 3)	y = 1.4697*exp(-1.29*p);
		if(sect == 4)	y = 1.5137*exp(-1.35*p);
		if(sect == 5)	y = 1.4889*exp(-1.29*p);
		if(sect == 6)	y = 1.4999*exp(-1.34*p);
		
	}
	if(part == 1)
	{				
		if(sect == 1)	y = 1.4822*exp(-1.22*p);
		if(sect == 2)	y = 1.4765*exp(-1.22*p);
		if(sect == 3)	y = 1.4851*exp(-1.23*p);
		if(sect == 4)	y = 1.4713*exp(-1.20*p);
		if(sect == 5)	y = 1.4864*exp(-1.23*p);
		if(sect == 6)	y = 1.462*exp(-1.19*p);
	}
	if(part == 2)
	{
		if(sect == 1)	y = 1.5501*exp(-0.59*p);
		if(sect == 2)	y = 1.4937*exp(-0.54*p);
		if(sect == 3)	y = 1.5516*exp(-0.58*p);
		if(sect == 4)	y = 1.5373*exp(-0.57*p);
		if(sect == 5)	y = 1.5463*exp(-0.58*p);
		if(sect == 6)	y = 1.5365*exp(-0.57*p);
	
	}
	return y;
}	
*/


int main(int __argc,char *__argv[])
{
	int start_s = clock();
	cout << "Start time: " << start_s/double(CLOCKS_PER_SEC) << " s" << endl;

	char *outFileName = "g12Kpipi.root";// output filename
 	extern int optind;   
 	TROOT troot();  
	TFile outFile(outFileName,"recreate");  


	//Int_t q[3], id[3];
	//Int_t sec[3], npd[3];

	Float_t run, evt;
	Float_t egam;
	Float_t pxkp1pcor, pykp1pcor, pzkp1pcor, pxpip1pcor, pypip1pcor, pzpip1pcor, pxpim1pcor, pypim1pcor, pzpim1pcor;
	Float_t pxpip2pcor, pypip2pcor, pzpip2pcor, pxpim2pcor, pypim2pcor, pzpim2pcor;
	Float_t vxkp1, vykp1, vzkp1, vxpip1, vypip1, vzpip1, vxpim1, vypim1, vzpim1;
	Float_t vxpip2, vypip2, vzpip2, vxpim2, vypim2, vzpim2;
	Float_t nn, n0, np;
	Float_t nkp, npip, npim;
	Float_t vtime, stvtime, delvt;
	Float_t kp1toflen, kp1beta, kp1Beta, pip1toflen, pip1beta, pip1Beta, pim1toflen, pim1beta, pim1Beta;
	Float_t pip2toflen, pip2beta, pip2Beta, pim2toflen, pim2beta, pim2Beta;
	Float_t seckp1, secpip1, secpim1, padkp1, padpip1, padpim1;
	Float_t secpip2, secpim2, padpip2, padpim2;
	Float_t tpho, helicity;

	Float_t xmmall, xmm2pi, xmmkp, xm2pi, xmall, xmmkppim, xmmkppip;
	Float_t cosTkpCM, cmeg, s_value;

	TLorentzVector P4pim1, P4pip1, P4kp1, P4miss, P4target, P4pho, P4inc, P4all, P4pipi, P4kpCM, P4phoCM; 
	TVector3 V3pim1, V3pip1, V3kp, Bcm;
	TVector3 P3pim, P3pip, P3kp;

	Float_t kp1Efficiency, pip1Efficiency, pim1Efficiency, pip2Efficiency, pim2Efficiency;

	Float_t dtkp1, dtpip1, dtpim1;

	Int_t sec[3], npd[3];

	Float_t Efficiency[5];

	Int_t totsigmsigp;
	
	//Int_t evnt = 0;
	//Int_t evnt1 = 0;
	Int_t ngamma = 0;
	Int_t total = 0;
	Int_t npho1 = 0, npho12 = 0;
	Int_t npho1sigp = 0, npho12sigp = 0;
	Int_t npho1sigm = 0, npho12sigm = 0;
 	const Float_t xmp =0.938272;
 	const Float_t xmpi = 0.139570;
 	const Float_t xmkp = 0.493677;
 	const Float_t vlight = 30.0;

 	TH1F *Hist_ngamma = new TH1F("Ngamma","Ngamma",100, 0., 3.);
 	TH1F *Hist_ngam = new TH1F("Ngam","Ngam",100, -1.0, 3.);
 	//TH2F *Histdtkp1 = new TH2F("dtkp1","dtkp1", 500, 0., 3., 500, 0. , 1.2);

 	TH2F *HistbetaVspkp1 = new TH2F("betaVspkp1","betaVspkp1", 500, 0., 3., 500, 0. , 1.2);
 	TH2F *HistbetaVsppip1 = new TH2F("betaVsppip1","betaVsppip1", 500, 0., 3., 500, 0. , 1.2);
 	TH2F *HistbetaVsppim1 = new TH2F("betaVsppim1","betaVsppim1", 500, 0., 3., 500, 0. , 1.2);
 	TH2F *HistBetaVspkp1 = new TH2F("BetaVspkp1","BetaVspkp1", 500, 0., 3., 500, 0. , 1.2);
 	TH2F *HistBetaVsppip1 = new TH2F("BetaVsppip1","BetaVsppip1", 500, 0., 3., 500, 0. , 1.2);
 	TH2F *HistBetaVsppim1 = new TH2F("BetaVsppim1","BetaVsppim1", 500, 0., 3., 500, 0. , 1.2);

 	TH2F *Histdtkp1 = new TH2F("dtkp1","dtkp1", 500, 0., 4., 500, -4. , 4.);
 	TH2F *Histdtpip1 = new TH2F("dtpip1","dtpip1", 500, 0., 4., 500, -4. , 4.);
 	TH2F *Histdtpim1 = new TH2F("dtpim1","dtpim1", 500, 0., 4., 500, -4. , 4.);

 	TH2F *Histdtkp1timcut = new TH2F("dtkp1timcut","dtkp1timcut", 500, 0., 4., 500, -4. , 4.);
 	TH2F *Histdtpip1timcut = new TH2F("dtpip1timcut","dtpip1timcut", 500, 0., 4., 500, -4. , 4.);
 	TH2F *Histdtpim1timcut = new TH2F("dtpim1timcut","dtpim1timcut", 500, 0., 4., 500, -4. , 4.);

	TH1F *HistEgam = new TH1F("Egam","E_{#gamma}",150,0.8,3.8);
	TH1F *HistCMenergy = new TH1F("CMenergy", "CM energy",100, 1.0,4.0);
	TH1F *HistcosTkpCM = new TH1F("cosTkpCM", "Boost K+ Pz/P", 200, -1.0, 1.0);
	TH2F *HistWVsCos = new TH2F("WVsCos", "WVsCos", 300, -1.05, 1.05, 300, 1.5, 3.5);

    //*************************histos with no cuts************************//
	TH1F *HistMMsq = new TH1F("MMsq","MM^{2}(#pi #pi K+)",800, -1.,4.);
	TH1F *HistMM = new TH1F("MM","MM(#pi #pi K+)",200, 0.0,2.);
	TH1F *HistM2pi = new TH1F("M2pi","M(#pi^{+} #pi^{-})",200, 0.,2.);
	TH1F *HistM2piK = new TH1F("M2piK","M(#pi #pi K+)",200, 0.,3.);
        TH2F *HistMMvM2pi = new TH2F("MMvM2pi","MM(#pi #pi) v. M(#pi^{+}#pi^{-})",200,0.8,2.8,200,0.46,0.52);
        TH1F *HistMMKp = new TH1F("MMKp","MM(K+)",300,0.0,2.5);
        TH2F *HistMMvM2piK = new TH2F("MMvM2piK","MM(#pi #pi K+) v. M(#pi #pi K+)",200,0.0,2.0,200,0.5,2.5);
    //*************************histos with no cuts************************//

    //*************** histos with cut 1 (exclusive, MM = neutron)*****************//
	TH1F *Hist1MMsq = new TH1F("1MMsq","MM^{2}(#pi #pi K+)",200,0.,4.);
	TH1F *Hist1M2pi = new TH1F("1M2pi","M(#pi^{+} #pi^{-})",300, 0.,1.5);
	TH1F *Hist1MM2pi = new TH1F("1MM2pi","MM(#pi #pi)",200, 1.4,2.4);
    	TH2F *Hist1MMvM2pi = new TH2F("1MMvM2pi","MM(#pi #pi) v. M(#pi^{+}#pi^{-})",200,0.8,2.8,200,0.46,0.52);
	TH1F *Hist1MMkpim = new TH1F("1MMkpim","MM(#pi- K+)",200, 1.,2.);
	TH1F *Hist1MMkpip = new TH1F("1MMkpip","MM(#pi+ K+)",200, 1.,2.);
        TH1F *Hist12MMKp = new TH1F("12MMKp","MM(K+) (#Sigma^{+} -K^{0})",300,1.0,2.5);
        TH1F *Hist13MMKp = new TH1F("13MMKp","MM(K+) (#Sigma^{-} -K^{0})",300,1.0,2.5);
        //Global Spectrum
		TH1F *HistGlb12MMKp = new TH1F("Glb12MMKp","MM(K+) (#Sigma^{+} -K^{0})",WCBIN, WCMIN, WCMAX);
		TH1F *HistGlb13MMKp = new TH1F("Glb13MMKp","MM(K+) (#Sigma^{-} -K^{0})",WCBIN, WCMIN, WCMAX);
	//*************** histos with cut 1 (exclusive, MM = neutron)*****************//

    //**************** histos with cut 2 (exclusive, MM = Lambda)*********************//
	TH1F *Hist2MMsq = new TH1F("2MMsq","MM^{2}(#pi #pi K+)",200,0.,4.);
	TH1F *Hist2M2pi = new TH1F("2M2pi","M(#pi^{+} #pi^{-})",300, 0.,1.5);
	TH1F *Hist2MM2pi = new TH1F("2MM2pi","MM(#pi #pi)",200, 1.4,2.4);
        TH2F *Hist2MMvM2pi = new TH2F("2MMvM2pi","MM(#pi #pi) v. M(#pi^{+}#pi^{-})",200,0.8,2.8,200,0.46,0.52);
	TH1F *Hist2MMkpim = new TH1F("2MMkpim","MM(#pi- K+)",200, 1.,2.);
	TH1F *Hist2MMkpip = new TH1F("2MMkpip","MM(#pi+ K+)",200, 1.,2.);
        TH2F *Hist2MMvMM = new TH2F("2MMvMM","MM(#pi- K+) v. MM(#pi+ K+)",200,1.2,2.0,200,1.2,2.0);
        TH1F *Hist22MMKp = new TH1F("22MMKp","MM(K+) (#Sigma^{*+} -K^{0})",300,1.4,2.9);
        TH1F *Hist23MMKp = new TH1F("23MMKp","MM(K+) (#Sigma^{*-} -K^{0})",300,1.4,2.9);
    //**************** histos with cut 2 (exclusive, MM = Lambda)*********************//

    //***************** histos with cut 3 (exclusive, MM = Sigma)**********************//
	TH1F *Hist3MMsq = new TH1F("3MMsq","MM^{2}(#pi #pi K+) (3)",200,0.,4.);
	TH1F *Hist3M2pi = new TH1F("3M2pi","M(#pi^{+} #pi^{-})",300, 0.,1.5);
	TH1F *Hist3MM2pi = new TH1F("3MM2pi","MM(#pi #pi)",200, 1.4,2.4);
        TH2F *Hist3MMvM2pi = new TH2F("3MMvM2pi","MM(#pi #pi) v. M(#pi^{+}#pi^{-})",200,0.8,2.8,200,0.46,0.52);
	TH1F *Hist3MMkpim = new TH1F("3MMkpim","MM(#pi- K+) (-K^{0})",200, 1.,2.);
	TH1F *Hist3MMkpip = new TH1F("3MMkpip","MM(#pi+ K+) (-K^{0})",200, 1.,2.);
        TH2F *Hist3MMvMM = new TH2F("3MMvMM","MM(#pi- K+) v. MM(#pi+ K+)",200,1.2,2.0,200,1.2,2.0);
        TH1F *Hist32MMKp = new TH1F("32MMKp","MM(K+) (#Sigma^{**+} -K^{0})",300,1.4,2.9);
        TH1F *Hist33MMKp = new TH1F("33MMKp","MM(K+) (#Sigma^{**-} -K^{0})",300,1.4,2.9);
    //***************** histos with cut 3 (exclusive, MM = Sigma)**********************//


	//**********************Fiducial and Paddle Cut Plots***********************//
		//theta vs phi
	TH2F *Hist_fidkp = new TH2F("Fiducial_kp","Data: K^{+}", 500, 0., PI, 500, -PI, PI);	//k+
	TH2F *Hist_fidpip = new TH2F("Fiducial_pip","Data: #pi^{+}", 500, 0., PI, 500, -PI, PI);	//pi+
	TH2F *Hist_fidpim = new TH2F("Fiducial_pim","Data: #pi^{-}", 500, 0., PI, 500, -PI, PI);	//pi-
	TH2F *Hist_fidcutkp = new TH2F("Fiducialcut_kp","Data: K^{+}", 500, 0., PI, 500, -PI, PI);	//k+
	TH2F *Hist_fidcutpip = new TH2F("Fiducialcut_pip","Data: #pi^{+}", 500, 0., PI, 500, -PI, PI);	//pi+
	TH2F *Hist_fidcutpim = new TH2F("Fiducialcut_pim","Data: #pi^{-}", 500, 0., PI, 500, -PI, PI);	//pi-
	
	//mom vs phi
	TH2F *Hist_fidkpmomphi = new TH2F("Fiducial_kpmomphi","Data: K^{+}", 500, 0., 5., 500, -PI, PI);	//k+
	TH2F *Hist_fidpipmomphi = new TH2F("Fiducial_pipmomphi","Data: #pi^{+}", 500, 0., 5., 500, -PI, PI);	//pi+
	TH2F *Hist_fidpimmomphi = new TH2F("Fiducial_pimmomphi","Data: #pi^{-}", 500, 0., 5., 500, -PI, PI);	//pi-
	TH2F *Hist_fidcutkpmomphi= new TH2F("Fiducialcut_kpmomphi","Data: K^{+}", 500, 0., 5., 500, -PI, PI);	//k+
	TH2F *Hist_fidcutpipmomphi = new TH2F("Fiducialcut_pipmomphi","Data: #pi^{+}", 500, 0., 5., 500, -PI, PI);	//pi+
	TH2F *Hist_fidcutpimmomphi = new TH2F("Fiducialcut_pimmomphi","Data: #pi^{-}", 500, 0., 5., 500, -PI, PI);	//pi-
	
	//theta(cosphi) vs theta(sinphi)
	TH2F *Hist_fidkpang = new TH2F("Fiducial_kpang","Data: K^{+}", 500, -2/3*PI, 2/3*PI, 500, -2/3*PI, 2/3*PI);	//k+
	TH2F *Hist_fidpipang = new TH2F("Fiducial_pipang","Data: #pi^{+}", 500, -2/3*PI, 2/3*PI, 500, -2/3*PI, 2/3*PI);	//pi+
	TH2F *Hist_fidpimang = new TH2F("Fiducial_pimang","Data: #pi^{-}", 500, -2/3*PI, 2/3*PI, 500, -2/3*PI, 2/3*PI);	//pi-
	TH2F *Hist_fidcutkpang= new TH2F("Fiducialcut_kpang","Data: K^{+}", 500, -2/3*PI, 2/3*PI, 500, -2/3*PI, 2/3*PI);	//k+
	TH2F *Hist_fidcutpipang = new TH2F("Fiducialcut_pipang","Data: #pi^{+}", 500, -2/3*PI, 2/3*PI, 500, -2/3*PI, 2/3*PI);	//pi+
	TH2F *Hist_fidcutpimang = new TH2F("Fiducialcut_pimang","Data: #pi^{-}",500, -2/3*PI, 2/3*PI, 500, -2/3*PI, 2/3*PI);	//pi-
	
	//histos for mom vs theta
	int sect;
	char parsect[30];
	TH2F *Histfidpipmomthetasec[6];
	TH2F *Histfidcutpipmomthetasec[6];
	TH2F *Histfidpimmomthetasec[6];
	TH2F *Histfidcutpimmomthetasec[6];
	TH2F *Histfidkpmomthetasec[6];
	TH2F *Histfidcutkpmomthetasec[6];
	
	TH2F *Histpadpipmomsec[6];
	TH2F *Histpadcutpipmomsec[6];
	TH2F *Histpadpimmomsec[6];
	TH2F *Histpadcutpimmomsec[6];
	TH2F *Histpadkpmomsec[6];
	TH2F *Histpadcutkpmomsec[6];
	
	for(sect = 1; sect < 7; ++sect)
	{
		sprintf(parsect, "Fiducial_pipmomthetasec%d",sect);		//without fiducal cut pi+ all sector
		Histfidpipmomthetasec[sect-1] = new TH2F(parsect, "Data: #pi^{+}",500, 0., 5., 500, 0.0, PI);
		
		sprintf(parsect, "Fiducialcut_pipmomthetasec%d",sect);		//with fiducal cut pi+ all sector
		Histfidcutpipmomthetasec[sect-1] = new TH2F(parsect, "Data: #pi^{+}", 500, 0., 5., 500, 0.0, PI);
		
		sprintf(parsect, "Fiducial_pimmomthetasec%d",sect);
		Histfidpimmomthetasec[sect-1] = new TH2F(parsect, "Data: #pi^{-}", 500, 0., 5., 500, 0.0, PI);
		
		sprintf(parsect, "Fiducialcut_pimmomthetasec%d",sect);
		Histfidcutpimmomthetasec[sect-1] = new TH2F(parsect, "Data: #pi^{-}", 500, 0., 5., 500, 0.0, PI);
		
		sprintf(parsect, "Fiducial_kpmomthetasec%d",sect);
		Histfidkpmomthetasec[sect-1] = new TH2F(parsect, "Data: K^{+}", 500, 0., 5., 500, 0.0, PI);
		
		sprintf(parsect, "Fiducialcut_kpmomthetasec%d",sect);
		Histfidcutkpmomthetasec[sect-1] = new TH2F(parsect, "Data: K^{+}", 500, 0., 5., 500, 0.0, PI);
		
		
		sprintf(parsect, "Paddle_pipmomsec%d",sect);		//without paddle cut pi+ all sector
		Histpadpipmomsec[sect-1] = new TH2F(parsect, "Data: #pi^{+}", 50,1.,50.,500, 0., 5.);
		
		sprintf(parsect, "Paddlecut_pipmomsec%d",sect);		//with paddle cut pi+ all sector
		Histpadcutpipmomsec[sect-1] = new TH2F(parsect, "Data: #pi^{+}", 50,1.,50.,500, 0., 5.);
		
		sprintf(parsect, "Paddle_pimmomsec%d",sect);
		Histpadpimmomsec[sect-1] = new TH2F(parsect, "Data: #pi^{-}", 50,1.,50.,500, 0., 5.);
		
		sprintf(parsect, "Paddlecut_pimmomsec%d",sect);
		Histpadcutpimmomsec[sect-1] = new TH2F(parsect, "Data: #pi^{-}", 50,1.,50.,500, 0., 5.);
		
		sprintf(parsect, "Paddle_kpmomsec%d",sect);
		Histpadkpmomsec[sect-1] = new TH2F(parsect, "Data: K^{+}", 50,1.,50.,500, 0., 5.);
		
		sprintf(parsect, "Paddlecut_kpmomsec%d",sect);
		Histpadcutkpmomsec[sect-1] = new TH2F(parsect, "Data: K^{+}", 50,1.,50.,500, 0., 5.);
	}
	//**********************Fiducial and Paddle CutPlots***********************//

	//*********************** Track Efficiency Values Histos ******************************//
	TH1F *Histkp1Eff = new TH1F("kp1Efficiency","Data: K^{+} Efficiency", 200, 0., 4.);
	TH1F *Histpip1Eff = new TH1F("pip1Efficiency","Data: #pi^{+} Efficiency", 200, 0., 4.);
	TH1F *Histpim1Eff = new TH1F("pim1Efficiency","Data: #pi^{-} Efficiency", 200, 0., 4.);
	TH1F *Histpip2Eff = new TH1F("pip2Efficiency","Data: #pi^{+} Efficiency", 200, 0., 4.);
	TH1F *Histpim2Eff = new TH1F("pim2Efficiency","Data: #pi^{-} Efficiency", 200, 0., 4.);

	TH1F *Histkp1Effallcut = new TH1F("kp1Efficiencyallcut","Data: K^{+} Efficiency allcut", 200, 0., 4.);
	TH1F *Histpip1Effallcut = new TH1F("pip1Efficiencyallcut","Data: #pi^{+} Efficiency allcut", 200, 0., 4.);
	TH1F *Histpim1Effallcut = new TH1F("pim1Efficiencyallcut","Data: #pi^{-} Efficiency allcut", 200, 0., 4.);
	// TH1F *Histpip2Effallcut = new TH1F("pip2Efficiencyallcut","Data: #pi^{+} Efficiency", 200, 0., 4.);
	// TH1F *Histpim2Effallcut = new TH1F("pim2Efficiencyallcut","Data: #pi^{-} Efficiency", 200, 0., 4.);
	//*********************** Track Efficiency Values Histos ******************************//


	//******************histos for W and Cos#theta****************//
	int w, c;
	char wc[20];
	TH1F *HistWCKPM[WBINS*CBINS];
	TH1F *HistWCKPP[WBINS*CBINS];
	for(w = 0; w < WBINS; ++w)
	{
		for(c = 0; c < CBINS; ++c)
		{
			sprintf(wc, "W%dC%02d(Kppim #Sigma^{+})", w, c);
			HistWCKPM[w*CBINS+c] = new TH1F(wc, "Data: Cos#theta binned in W(CMenergy)", WCBIN, WCMIN, WCMAX);

			sprintf(wc, "W%dC%02d(Kppip #Sigma^{-})", w, c);
			HistWCKPP[w*CBINS+c] = new TH1F(wc, "Data: Cos#theta binned in W(CMenergy)", WCBIN, WCMIN, WCMAX);
		}
	} 
	//******************histos for W and Cos#theta****************//     



//******************************* MAIN PROGRAM ********************************//
	for(int n_arg = optind; n_arg < __argc; n_arg++)
	{	//loop over data files
   
    	TFile inFile(__argv[n_arg]);

		// open each input file        
		if( TTree *p1 = (TTree*)inFile.Get("kpipi"))
		{
			p1->SetBranchAddress("run", &run);
			p1->SetBranchAddress("evt", &evt);
			p1->SetBranchAddress("ebeamcor", &egam);
			p1->SetBranchAddress("pxkp1pcor", &pxkp1pcor);
			p1->SetBranchAddress("pykp1pcor", &pykp1pcor);
			p1->SetBranchAddress("pzkp1pcor", &pzkp1pcor);
			p1->SetBranchAddress("pxpip1pcor", &pxpip1pcor);
			p1->SetBranchAddress("pypip1pcor", &pypip1pcor);
			p1->SetBranchAddress("pzpip1pcor", &pzpip1pcor);
			p1->SetBranchAddress("pxpim1pcor", &pxpim1pcor);
			p1->SetBranchAddress("pypim1pcor", &pypim1pcor);
			p1->SetBranchAddress("pzpim1pcor", &pzpim1pcor);
			p1->SetBranchAddress("vxkp1", &vxkp1);
			p1->SetBranchAddress("vykp1", &vykp1);
			p1->SetBranchAddress("vzkp1", &vzkp1);
			p1->SetBranchAddress("vxpip1", &vxpip1);
			p1->SetBranchAddress("vypip1", &vypip1);
			p1->SetBranchAddress("vzpip1", &vzpip1);
			p1->SetBranchAddress("vxpim1", &vxpim1);
			p1->SetBranchAddress("vypim1", &vypim1);
			p1->SetBranchAddress("vzpim1", &vzpim1);
			p1->SetBranchAddress("nn", &nn);
			p1->SetBranchAddress("n0", &n0);
			p1->SetBranchAddress("np", &np);
			p1->SetBranchAddress("nkp", &nkp);
			p1->SetBranchAddress("npip", &npip);
			p1->SetBranchAddress("npim", &npim);
			p1->SetBranchAddress("vtime", &vtime);
			p1->SetBranchAddress("stvtime", &stvtime);
			p1->SetBranchAddress("delvt", &delvt);
			p1->SetBranchAddress("kp1toflen", &kp1toflen);
			p1->SetBranchAddress("kp1beta", &kp1beta);
			p1->SetBranchAddress("kp1Beta", &kp1Beta);
			p1->SetBranchAddress("pip1toflen", &pip1toflen);
			p1->SetBranchAddress("pip1beta", &pip1beta);
			p1->SetBranchAddress("pip1Beta", &pip1Beta);
			p1->SetBranchAddress("pim1toflen", &pim1toflen);
			p1->SetBranchAddress("pim1beta", &pim1beta);
			p1->SetBranchAddress("pim1Beta", &pim1Beta);
			p1->SetBranchAddress("seckp1", &seckp1);
			p1->SetBranchAddress("secpip1", &secpip1);
			p1->SetBranchAddress("secpim1", &secpim1);
			p1->SetBranchAddress("padkp1", &padkp1);
			p1->SetBranchAddress("padpip1", &padpip1);
			p1->SetBranchAddress("padpim1", &padpim1);
			p1->SetBranchAddress("tpho", &tpho);
			p1->SetBranchAddress("pxpip2pcor", &pxpip2pcor);
			p1->SetBranchAddress("pypip2pcor", &pypip2pcor);
			p1->SetBranchAddress("pzpip2pcor", &pzpip2pcor);
			p1->SetBranchAddress("pxpim2pcor", &pxpim2pcor);
			p1->SetBranchAddress("pypim2pcor", &pypim2pcor);
			p1->SetBranchAddress("pzpim2pcor", &pzpim2pcor);
			p1->SetBranchAddress("kp1Efficiency", &kp1Efficiency);
			p1->SetBranchAddress("pip1Efficiency", &pip1Efficiency);
			p1->SetBranchAddress("pim1Efficiency", &pim1Efficiency);
			p1->SetBranchAddress("pip2Efficiency", &pip2Efficiency);
			p1->SetBranchAddress("pim2Efficiency", &pim2Efficiency);

			Int_t nentries = (Int_t)p1->GetEntries();
			//do{
			for (Int_t j=0;j<nentries;j++)
			{ 
				p1->GetEntry(j);	//Loop over entries
				
				if (evt == evt)
				{
					ngamma = 1;
					++npho1;
				}
				if (evt == evt -1)
				{
					ngamma = 2;
					++npho12;
				}
				


				for ( Int_t ngam = 0; ngam < ngamma; ngam++)
				{

					++total;
					// {
					// 	++evnt;
					// 	//cout << evt << endl;
					// 	//cout << evnt << endl;
					// }
					// 	//cout << evt << endl;
					// 	++evnt1;
					// 	cout << "****" << evnt1 << "****" << endl;
					// 	cout << "-------------" << evnt << "--------------" << endl;
					//for (Int_t ievnt = 0; ievnt < evt+1; e)
					//cout << evt << endl;
					// do{
					// 	p1->GetEntry(j);	//Loop over entries

					for(Int_t ipip = 1; ipip <= npip; ipip++)	// Loop over 1 and 2 piplus
					{

					for(Int_t ipim = 1; ipim <= npim; ipim++)	// Loop over 1 and 2 piminus
					{	
						// Set the target and photon values
						P4target.SetPxPyPzE(0.0,0.0,0.0,xmp);		
						P4pho.SetPxPyPzE(0.0,0.0,egam,egam);
						P4inc =  P4target + P4pho;

						// Four-vector of the particles

						P4kp1.SetPxPyPzE(pxkp1pcor,pykp1pcor,pzkp1pcor,sqrt(pxkp1pcor*pxkp1pcor+pykp1pcor*pykp1pcor+pzkp1pcor*pzkp1pcor+xmkp*xmkp));
						sec[0] = seckp1;
						npd[0] = padkp1;
						Efficiency[0] = kp1Efficiency;


						if(ipip == 1)
						{
							P4pip1.SetPxPyPzE(pxpip1pcor,pypip1pcor,pzpip1pcor,sqrt(pxpip1pcor*pxpip1pcor+pypip1pcor*pypip1pcor+pzpip1pcor*pzpip1pcor+xmpi*xmpi));
							sec[1] = secpip1;
							npd[1] = padpip1;
							Efficiency[1] = pip1Efficiency;
							
						}
						if(ipip == 2)
						{	
							P4pip1.SetPxPyPzE(pxpip2pcor,pypip2pcor,pzpip2pcor,sqrt(pxpip2pcor*pxpip2pcor+pypip2pcor*pypip2pcor+pzpip2pcor*pzpip2pcor+xmpi*xmpi)); 
							sec[1] = secpip2;
							npd[1] = padpip2;
							Efficiency[1] = pip2Efficiency;
							Efficiency[3] = pip2Efficiency;
							pip1beta = pip2beta;
							pip1toflen = pip2toflen;
							pip1Beta = pip2Beta;

						}
						if(ipim == 1)
						{
							P4pim1.SetPxPyPzE(pxpim1pcor,pypim1pcor,pzpim1pcor,sqrt(pxpim1pcor*pxpim1pcor+pypim1pcor*pypim1pcor+pzpim1pcor*pzpim1pcor+xmpi*xmpi));
							sec[2] = secpim1;
							npd[2] = padpim1;
							Efficiency[2] = pim1Efficiency;
						}
						if(ipim == 2)
						{	
							P4pim1.SetPxPyPzE(pxpim2pcor,pypim2pcor,pzpim2pcor,sqrt(pxpim2pcor*pxpim2pcor+pypim2pcor*pypim2pcor+pzpim2pcor*pzpim2pcor+xmpi*xmpi));
							sec[2] = secpim2;
							npd[2] = padpim2;
							Efficiency[2] = pim2Efficiency;
							Efficiency[4] = pim2Efficiency;
							pim1beta = pim2beta;
							pim1toflen = pim2toflen;
							pim1Beta = pim2Beta;  
						}
						  


						// Vertices 
						V3pip1.SetXYZ(vxpip1,vypip1,vzpip1);  
						V3pim1.SetXYZ(vxpim1,vypim1,vzpim1);  
						V3kp.SetXYZ(vxkp1,vykp1,vzkp1);

						//momentum vectors
						P3pip = P4pip1.Vect();
						P3pim = P4pim1.Vect();
						P3kp = P4kp1.Vect();

						// Calculate various 4-vectors
						P4pipi = P4pip1 + P4pim1;
						P4all = P4pip1 + P4pim1 + P4kp1;
						P4miss = P4inc - P4all;

						// Get invariant masses
						xmmall = P4miss.M();
						xmm2pi = (P4inc - P4pipi).M();
						xmmkp = (P4inc - P4kp1).M();
						xmmkppim = (P4inc - P4kp1 - P4pim1).M();
						xmmkppip = (P4inc - P4kp1 - P4pip1).M();
						xm2pi = P4pipi.M();
						xmall = P4all.M();

						//Boost into CM frame
						Bcm.SetZ((-egam/(egam+xmp)));	//eq46.4 pdg
						P4kpCM = P4kp1;
						P4phoCM = P4pho;
						P4kpCM.Boost(Bcm);					//Boost the k+ to CM frame
						P4phoCM.Boost(Bcm);					//Boost the photon to CM frame
						cosTkpCM = P4kpCM.CosTheta();		//cosTkpCM.Boost(Bcm);
						s_value = (P4inc).M2();				//Mandelstam s-variable
						cmeg = sqrt(s_value);

						//**************** pid Histos without pidcuts***********//
						//Float_t dtkp1 = stvtime - kp1toflen/(vlight*kp1beta);
						//Float_t dtkp1 = kp1beta;
						//Histdtkp1->Fill(P3kp.Mag(),dtkp1);
						HistbetaVspkp1->Fill(P3kp.Mag(),kp1beta);
						HistbetaVsppip1->Fill(P3pip.Mag(),pip1beta);
						HistbetaVsppim1->Fill(P3pim.Mag(),pim1beta);

						HistBetaVspkp1->Fill(P3kp.Mag(),kp1Beta);
						HistBetaVsppip1->Fill(P3pip.Mag(),pip1Beta);
						HistBetaVsppim1->Fill(P3pim.Mag(),pim1Beta);

						dtkp1 = kp1toflen*(1/kp1Beta -1/kp1beta)/vlight;
						dtpip1 = pip1toflen*(1/pip1Beta -1/pip1beta)/vlight;
						dtpim1 = pim1toflen*(1/pim1Beta -1/pim1beta)/vlight;
						//**************** pid Histos without pidcuts***********//


						//****************** pidcut ****************************//
						bool pidcut = false;
						if(fabs(dtkp1)< 1.0 && fabs(dtpip1) < 1.0 && fabs(dtpim1) < 1.0)
						{
							Histdtkp1-> Fill(P3kp.Mag(), dtkp1);
							Histdtpip1-> Fill(P3pip.Mag(), dtpip1);
							Histdtpim1-> Fill(P3pim.Mag(), dtpim1);

							//pidcut = true;

							if(CutTimingPID(0, 0, P3kp.Mag(), dtkp1, 3.0) && CutTimingPID(1, 0, P3pip.Mag(), dtpip1, 3.0) && CutTimingPID(2, 0, P3pim.Mag(), dtpim1, 3.0))
							{
								Histdtkp1timcut-> Fill(P3kp.Mag(), dtkp1);
								Histdtpip1timcut-> Fill(P3pip.Mag(), dtpip1);
								Histdtpim1timcut-> Fill(P3pim.Mag(), dtpim1);
								pidcut = true;
							}
						}
						//****************** pidcut ****************************//


						//***************Fiducial and Paddle Cuts******************//
						//Fiducial Cuts
						bool fiducialcut = false;
						if(clas::g12::g12_PosParticle_fiducial_cuts(P3kp.Mag(),	P4kp1.Theta()*180/PI, P4kp1.Phi()*180/PI, "nominal")== true && 
							clas::g12::g12_PosParticle_fiducial_cuts(P3pip.Mag(), P4pip1.Theta()*180/PI, P4pip1.Phi()*180/PI, "nominal")== true && 							clas::g12::g12_NegParticle_fiducial_cuts(P3pim.Mag(), P4pim1.Theta()*180/PI, P4pim1.Phi()*180/PI, "nominal")== true)
						{
							fiducialcut = true;
						}
						
						//Bad paddle knock out
						//if(sec[1]==2)	Histpaddle10->Fill(npd[1],P3pip.Mag());
						bool padcut1 = false;	// paddle knock off based on bad occupancy
						if(clas::g12::pass_g12_TOFKO(sec[0],npd[0], 0) == true && 
							clas::g12::pass_g12_TOFKO(sec[1],npd[1], 0) == true &&
								clas::g12::pass_g12_TOFKO(sec[2],npd[2], 0) == true)
						{
							padcut1 = true;
						}
						
						bool padcut2 = false;	// paddle knock off based on bad resolution
						if(clas::g12::pass_g12_TOFKO(sec[0],npd[0], 1) == true && 
							clas::g12::pass_g12_TOFKO(sec[1],npd[1], 1) == true &&
								clas::g12::pass_g12_TOFKO(sec[2],npd[2], 1) == true)
						{
							padcut2 = true;
						}
						
						bool padcut12 = false;
						if(padcut1 == true && padcut2 == true)
						{
							padcut12 = true;
						}
						//if(padcut12 == true && sec[1] == 2)	Histpaddle11->Fill(npd[1],P3pip.Mag());
						//***************Fiducial and Paddle Cuts******************//
						
						/*
						// MomCut functin

						//cout << P4pip1.Theta() << "\t" << MomCut(sec[0],P3pip.Mag(),0) << endl;
						bool cutM = false;
						if(P4pip1.Theta() <= MomCut(sec[0],P3pip.Mag(),0))
							cutM = true;
						else
							cutM = false;
			
						if((P4pim1.Theta() <= MomCut(sec[1],P3pim.Mag(),1)) && cutM == true)
							cutM = true;
						else
							cutM = false;
			
						if((P4kp1.Theta() <= MomCut(sec[2],P3kp.Mag(),2)) && cutM == true)
							cutM = true;
						else
							cutM = false;	
						*/


						//**********Fiducal and Paddle Histos without cuts**********//
						Hist_fidkp->Fill(P4kp1.Theta(), P4kp1.Phi());
						Hist_fidpip->Fill(P4pip1.Theta(), P4pip1.Phi());
						Hist_fidpim->Fill(P4pim1.Theta(), P4pim1.Phi());

						Hist_fidkpmomphi->Fill(P3kp.Mag(), P4kp1.Phi());
						Hist_fidpipmomphi->Fill(P3pip.Mag(), P4pip1.Phi());
						Hist_fidpimmomphi->Fill(P3pim.Mag(), P4pim1.Phi());

						Hist_fidkpang->Fill(P4kp1.Theta()*sin(P4kp1.Phi()), P4kp1.Theta()*cos(P4kp1.Phi()));
						Hist_fidpipang->Fill(P4pip1.Theta()* sin(P4pip1.Phi()), P4pip1.Theta()*cos(P4pip1.Phi()));
						Hist_fidpimang->Fill(P4pim1.Theta()*sin(P4pim1.Phi()), P4pim1.Theta()* cos(P4pim1.Phi()));
						
						for(sect = 1; sect < 7 ; ++sect)
						{
							if(sec[0] == sect)	Histfidpipmomthetasec[sect-1]->Fill(P3pip.Mag(), P4pip1.Theta());
							if(sec[1] == sect)	Histfidpimmomthetasec[sect-1]->Fill(P3pim.Mag(), P4pim1.Theta());
							if(sec[2] == sect)	Histfidkpmomthetasec[sect-1]->Fill(P3kp.Mag(), P4kp1.Theta());
						}
															
						//Paddle histos without cut
						for(sect = 1; sect < 7 ; ++sect)
						{
							if(sec[0] == sect)	Histpadpipmomsec[sect-1]->Fill(npd[0], P3pip.Mag());
							if(sec[1] == sect)	Histpadpimmomsec[sect-1]->Fill(npd[1], P3pim.Mag());
							if(sec[2] == sect)	Histpadkpmomsec[sect-1]->Fill(npd[2], P3kp.Mag());
						}
						//**********Fiducal and Paddle Histos without cuts**********//


						//******************** Track Efficiency Histos without any cuts************//
						Histkp1Eff->Fill(Efficiency[0]);
						Histpip1Eff->Fill(Efficiency[1]);
						Histpim1Eff->Fill(Efficiency[2]);
						Histpip2Eff->Fill(Efficiency[3]);
						Histpim2Eff->Fill(Efficiency[4]);
						//******************** Track Efficiency Histos without any cuts************//


						//************************* all cuts ************************//
						bool allcuts = false;
						if(pidcut == true && fiducialcut == true & padcut12 == true)
						{
							allcuts = true;
						}
						//************************* all cuts ************************//


						if(allcuts == true)
						{

							Hist_ngamma->Fill(ngamma);
							Hist_ngam->Fill(ngam);

							// Fill Histograms
							HistEgam->Fill(egam);
							HistCMenergy->Fill(cmeg);
							HistcosTkpCM->Fill(cosTkpCM);
							//HistWVsCos->Fill(cosTkpCM,cmeg);

							HistMMsq->Fill((P4miss).M2());
							HistMM->Fill(xmmall);
							HistMMKp->Fill(xmmkp);
							HistM2pi->Fill(xm2pi);
							HistM2piK->Fill(xmall);
							HistMMvM2pi->Fill(xmm2pi,xm2pi);
							HistMMvM2piK->Fill(xmmall,xmall);

								// exclusive reaction (MM=0)
								if( xmmall > 0.9 && xmmall < 1.0 ) 
								{
									Hist1MMsq->Fill((P4miss).M2());
									Hist1M2pi->Fill(xm2pi);
									Hist1MM2pi->Fill(xmm2pi);
									Hist1MMvM2pi->Fill(xmm2pi,xm2pi);
									
									//remove K0 events
									if( fabs(xm2pi-0.495) > 0.015 )
									{
										HistWVsCos->Fill(cosTkpCM,cmeg);
										
										Hist1MMkpim->Fill(xmmkppim);
										Hist1MMkpip->Fill(xmmkppip);

										if(xmmkppim > minsigp && xmmkppim < maxsigp)
										{
											if (fabs(xmmkppim - sigp) < fabs(xmmkppip - sigp))
											{
												Hist12MMKp->Fill(xmmkp);
												HistGlb12MMKp->Fill(xmmkp);
												
												if (cmeg > WMIN && cmeg <= WMAX)
												{
													if (cosTkpCM > CMIN && cosTkpCM <= CMAX)
													{
														w = int(abs((cmeg - WMIN) * WBINS / (WMAX - WMIN)));
														c = int(abs((cosTkpCM - CMIN) / CBINSIZE));
														HistWCKPM[CBINS*w+c]->Fill(xmmkp); 
													}
												}

												if (evt == evt)
												{
													++npho1sigp;
												}
												if (evt == evt -1)
												{
													++npho12sigp;
												}
											}
										}

										if(xmmkppip > minsigm && xmmkppip < maxsigm)
										{
											if (fabs(xmmkppip - sigm) < fabs(xmmkppim - sigm))
											{
												Hist13MMKp->Fill(xmmkp);
												HistGlb13MMKp->Fill(xmmkp);
																
												if (cmeg > WMIN && cmeg <= WMAX)
												{
													if (cosTkpCM > CMIN && cosTkpCM <= CMAX)
													{
														w = int(abs((cmeg - WMIN) * WBINS / (WMAX - WMIN)));
														c = int(abs((cosTkpCM - CMIN) / CBINSIZE));
														HistWCKPP[CBINS*w+c]->Fill(xmmkp); 
													}
												}

												if (evt == evt)
												{
													++npho1sigm;
												}
												if (evt == evt -1)
												{
													++npho12sigm;
												}
											}
										}

										//************* fiducial and paddle Histos with all cuts **********************//
										// "nominal" fiducial cut
										Hist_fidcutkp->Fill(P4kp1.Theta(), P4kp1.Phi());
										Hist_fidcutpip->Fill(P4pip1.Theta(), P4pip1.Phi());
										Hist_fidcutpim->Fill(P4pim1.Theta(), P4pim1.Phi());
								
										Hist_fidcutkpmomphi->Fill(P3kp.Mag(), P4kp1.Phi());
										Hist_fidcutpipmomphi->Fill(P3pip.Mag(), P4pip1.Phi());
										Hist_fidcutpimmomphi->Fill(P3pim.Mag(), P4pim1.Phi());
								
										Hist_fidcutkpang->Fill(P4kp1.Theta()*sin(P4kp1.Phi()), P4kp1.Theta()* cos(P4kp1.Phi()));
										Hist_fidcutpipang->Fill(P4pip1.Theta()*  sin(P4pip1.Phi()), P4pip1.Theta()* cos(P4pip1.Phi()));
										Hist_fidcutpimang->Fill(P4pim1.Theta()*sin(P4pim1.Phi()), P4pim1.Theta()* cos(P4pim1.Phi()));
																
										for(sect = 1; sect < 7 ; ++sect)
										{
											if(sec[0] == sect)	Histfidcutpipmomthetasec[sect-1]->Fill(P3pip.Mag(), P4pip1.Theta());
											if(sec[1] == sect)	Histfidcutpimmomthetasec[sect-1]->Fill(P3pim.Mag(), P4pim1.Theta());
											if(sec[2] == sect)	Histfidcutkpmomthetasec[sect-1]->Fill(P3kp.Mag(), P4kp1.Theta());
										}
										
										// Paddle cut histos with cut12
										for(sect = 1; sect < 7 ; ++sect)
										{
											if(sec[0] == sect)	Histpadcutpipmomsec[sect-1]->Fill(npd[0], P3pip.Mag());
											if(sec[1] == sect)	Histpadcutpimmomsec[sect-1]->Fill(npd[1], P3pim.Mag());
											if(sec[2] == sect)	Histpadcutkpmomsec[sect-1]->Fill(npd[2], P3kp.Mag());
										}
										//************* fiducial and paddle Histos with all cuts **********************//

										//******************** Track Efficiency Histos without any cuts************//
										Histkp1Effallcut->Fill(Efficiency[0]);
										Histpip1Effallcut->Fill(Efficiency[1]);
										Histpim1Effallcut->Fill(Efficiency[2]);
										// Histpip2Effallcut->Fill(Efficiency[3]);
										// Histpim2Effallcut->Fill(Efficiency[4]);
										//******************** Track Efficiency Histos without any cuts************//
										++totsigmsigp;
									}
								}

								// exclusive (MM = Lambda)
								if( xmmall > 1.09 && xmmall < 1.15 )
								{
									Hist2MMsq->Fill((P4miss).M2());
									Hist2M2pi->Fill(xm2pi);
									Hist2MM2pi->Fill(xmm2pi);
									Hist2MMvM2pi->Fill(xmm2pi,xm2pi);
									if( fabs(xm2pi-0.495) > 0.015 )
									{
										Hist2MMkpim->Fill(xmmkppim);
										Hist2MMkpip->Fill(xmmkppip);
										Hist2MMvMM->Fill(xmmkppim,xmmkppip);
										if(xmmkppim > 1.34 && xmmkppim < 1.42 )
										{
											Hist22MMKp->Fill(xmmkp);
										}
										if(xmmkppip > 1.34 && xmmkppip < 1.42 )
										{
											Hist23MMKp->Fill(xmmkp);
										}
									}
								}

								// exclusive (MM = Sigma)
								if( xmmall > 1.16 && xmmall < 1.25 )
								{
									Hist3MMsq->Fill(P4miss.M2());
									Hist3M2pi->Fill(xm2pi);
									Hist3MM2pi->Fill(xmm2pi);
									Hist3MMvM2pi->Fill(xmm2pi,xm2pi);
									if( fabs(xm2pi-0.495) > 0.015 )
									{
										Hist3MMkpim->Fill(xmmkppim);
										Hist3MMkpip->Fill(xmmkppip);
										Hist3MMvMM->Fill(xmmkppim,xmmkppip);
										if(xmmkppim > 1.62 && xmmkppim < 1.7 )
										{
											Hist32MMKp->Fill(xmmkp);
										}
										if(xmmkppip > 1.62 && xmmkppip < 1.7 )
										{
											Hist33MMKp->Fill(xmmkp);
										}
									}
								}
						
						}	// allcuts
						
					}	// Loop over 1 and 2 piminus
					}	// Loop over 1 and 2 piplus
				}	// Loop over 1 and 2 photons
			//} while(j-1 == evt);
			}//for entries (FORloop)
			cout << "-----nentries = " << nentries << "-------" << endl;
			cout << "-----nphoton1 = " << npho1 << "-------" << endl;
			cout << "-----nphoton1&2 = " << npho12 << "-------" << endl;
			cout << "-----photon multiplicity = " << 1.*npho12/(npho1)*100.0 << " % -------" << endl;
			cout << "------total = " << total << endl;
			cout << "-----nphoton1sigp = " << npho1sigp << "-------" << endl;
			cout << "-----nphoton1&2sigp = " << npho12sigp << "-------" << endl;
			cout << "-----nphoton1sigm = " << npho1sigm << "-------" << endl;
			cout << "-----nphoton1&2sigm = " << npho12sigm << "-------" << endl;
			cout << "-----neventsigmsigp = " << totsigmsigp << "-------" << endl;
			//} while(evt == evt - 1);
			//cout << nentries << endl;//cout << evnt << endl;
		}//if TTree *p1 = (TTree*)

	    else
	    {
			cout << "File has no TTree or p1 does not exist!!" << endl;
	    }

	    cout << n_arg << "/" << __argc - 1 << endl;  
   	}//n_arg
//*************************       MAIN PROGRAM     *****************************//

   outFile.Write(); // write to the output file
   outFile.Close(); // close the output file

   cout << "Time taken: " << (double)(clock() - start_s)/CLOCKS_PER_SEC << " s" << endl;

}


