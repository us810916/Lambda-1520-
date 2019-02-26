/*******************************************************************************
	This program is for analysis of g12 data for the K+ pi+ pi- skim.
	Lambda(1520) analysis: 
	Utsav Shrestha
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
#include "/home/shrestha/project/others/g12/G12_Corrections-master/g12_corrections.hpp"

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

//const Float_t vlight=30.0;

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
	//char *outFileName = "g12Kpipitr1fd1pd12momcut.root";// output filename
	//char *outFileName = "g12Kpipinewskim.root";// output filename
	//char *outFileName = "g12Kpipioujlab3.root";// output filename
	
	char *outFileName = "g12Kpipioujlab3allcut.root";// output filename
 	extern int optind;   
 	TROOT troot();  
	TFile outFile(outFileName,"recreate");  

	Int_t q[3],id[3];
	Int_t sec[3],npd[3];
	Int_t EID, TID;
	int ngamma, ngam, ngam1, ngm[2];
	
	Float_t egam[2],px[3],py[3],pz[3],pe[3],vx[3],vy[3],vz[3],tof[3],path[3]; 
	Float_t xmmall,xmm2pi,xmmkp,xm2pi,xmall,xmmkppim,xmmkppip;
	Float_t cosTkpCM, cmeg, s_value;
	TLorentzVector P4pim1,P4pip1,P4kp1,P4miss,P4target,P4pho,P4inc,P4all,P4pipi,P4kpCM,P4phoCM; 
	TVector3 V3pim1, V3pip1, V3kp, Bcm;
	TVector3 P3pim, P3pip, P3kp;
	Float_t mompip, mompim, momkp;
	int tFlag = -9;	// 3 particles detected in 3 different sectors.

 	const Float_t xmp =0.938272;
	
	TH1F *HistNgamma = new TH1F("Ngamma", "Ngamma",20,0.,200.);
	TH1F *HistNgam = new TH1F("Ngam", "Ngam",100,-20.,150.);
	TH1F *HistEgam = new TH1F("Egam","E_{#gamma}",150,0.8,3.8);
	TH1F *Histngm0 = new TH1F("ngm0","ngm0",5,0.0,5.);
	TH1F *Histngm1 = new TH1F("ngm1","ngm1",5,0.0,5.);
	TH1F *HistCMenergy = new TH1F("CMenergy", "CM energy",100, 1.0,4.0);
	TH1F *HistcosTkpCM = new TH1F("cosTkpCM", "Boost K+ Pz/P", 200, -1.0, 1.0);
	TH2F *HistEgamVsEID = new TH2F("EgamVsEID", "E_{#gamma} Vs EID", 1000,1.0,1001.0,150,0.8,3.8);
	TH2F *HistWVsCos = new TH2F("WVsCos", "WVsCos", 300, -1.05, 1.05, 300, 1.5, 3.5);
	
	

	// histos with no cuts
	TH1F *HistMMsq = new TH1F("MMsq","MM^{2}(#pi #pi K+)",800, -1.,4.);
	TH1F *HistMM = new TH1F("MM","MM(#pi #pi K+)",200, 0.0,2.);
	TH1F *HistM2pi = new TH1F("M2pi","M(#pi^{+} #pi^{-})",200, 0.,2.);
	TH1F *HistM2piK = new TH1F("M2piK","M(#pi #pi K+)",200, 0.,3.);
	TH2F *HistMMvM2pi = new TH2F("MMvM2pi","MM(#pi #pi) v. M(#pi^{+}#pi^{-})",200,0.8,2.8,200,0.46,0.52);
	TH1F *HistMMKp = new TH1F("MMKp","MM(K+)",300,0.0,2.5);
	TH2F *HistMMvM2piK = new TH2F("MMvM2piK","MM(#pi #pi K+) v. M(#pi #pi K+)",200,0.0,2.0,200,0.5,2.5);

   	// histos with cut 1 (exclusive, MM = proton)
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
	
	//Momentum Plot
	//no cut
	TH1F *Histmomkp = new TH1F("momkp","momkp",100, -1.0,5.0);
	TH1F *Histmompip = new TH1F("mompip", "mompip",100, -1.0,5.0);
	TH1F *Histmompim = new TH1F("mompim", "mompim",100, -1.0,5.0);
	//n cut
	TH1F *Hist1momkp = new TH1F("1momkp","1momkp",100, -1.0,5.0);
	TH1F *Hist1mompip = new TH1F("1mompip", "1mompip",100, -1.0,5.0);
	TH1F *Hist1mompim = new TH1F("1mompim", "1mompim",100, -1.0,5.0);
	//sigma cut
	TH1F *Hist12momkp = new TH1F("12momkp","12momkp",100, -1.0,5.0);
	TH1F *Hist12mompip = new TH1F("12mompip", "12mompip",100, -1.0,5.0);
	TH1F *Hist12mompim = new TH1F("12mompim", "12mompim",100, -1.0,5.0);
	TH1F *Hist13momkp = new TH1F("13momkp","13momkp",100, -1.0,5.0);
	TH1F *Hist13mompip = new TH1F("13mompip", "13mompip",100, -1.0,5.0);
	TH1F *Hist13mompim = new TH1F("13mompim", "13mompim",100, -1.0,5.0);
	
	// histos with cut 2 (exclusive, MM = Lambda)
	TH1F *Hist2MMsq = new TH1F("2MMsq","MM^{2}(#pi #pi K+)",200,0.,4.);
	TH1F *Hist2M2pi = new TH1F("2M2pi","M(#pi^{+} #pi^{-})",300, 0.,1.5);
	TH1F *Hist2MM2pi = new TH1F("2MM2pi","MM(#pi #pi)",200, 1.4,2.4);
	TH2F *Hist2MMvM2pi = new TH2F("2MMvM2pi","MM(#pi #pi) v. M(#pi^{+}#pi^{-})",200,0.8,2.8,200,0.46,0.52);
	TH1F *Hist2MMkpim = new TH1F("2MMkpim","MM(#pi- K+)",200, 1.,2.);
	TH1F *Hist2MMkpip = new TH1F("2MMkpip","MM(#pi+ K+)",200, 1.,2.);
	TH2F *Hist2MMvMM = new TH2F("2MMvMM","MM(#pi- K+) v. MM(#pi+ K+)",200,1.2,2.0,200,1.2,2.0);
	TH1F *Hist22MMKp = new TH1F("22MMKp","MM(K+) (#Sigma^{*+} -K^{0})",300,1.4,2.9);
	TH1F *Hist23MMKp = new TH1F("23MMKp","MM(K+) (#Sigma^{*-} -K^{0})",300,1.4,2.9);

	// histos with cut 3 (exclusive, MM = Sigma)
	TH1F *Hist3MMsq = new TH1F("3MMsq","MM^{2}(#pi #pi K+) (3)",200,0.,4.);
	TH1F *Hist3M2pi = new TH1F("3M2pi","M(#pi^{+} #pi^{-})",300, 0.,1.5);
	TH1F *Hist3MM2pi = new TH1F("3MM2pi","MM(#pi #pi)",200, 1.4,2.4);
	TH2F *Hist3MMvM2pi = new TH2F("3MMvM2pi","MM(#pi #pi) v. M(#pi^{+}#pi^{-})",200,0.8,2.8,200,0.46,0.52);
	TH1F *Hist3MMkpim = new TH1F("3MMkpim","MM(#pi- K+) (-K^{0})",200, 1.,2.);
	TH1F *Hist3MMkpip = new TH1F("3MMkpip","MM(#pi+ K+) (-K^{0})",200, 1.,2.);
	TH2F *Hist3MMvMM = new TH2F("3MMvMM","MM(#pi- K+) v. MM(#pi+ K+)",200,1.2,2.0,200,1.2,2.0);
	TH1F *Hist32MMKp = new TH1F("32MMKp","MM(K+) (#Sigma^{**+} -K^{0})",300,1.4,2.9);
	TH1F *Hist33MMKp = new TH1F("33MMKp","MM(K+) (#Sigma^{**-} -K^{0})",300,1.4,2.9);
	
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
	
	// Histos for paddle cut check
	TH2F *Histpaddle10 = new TH2F("Histpaddle10", "Histpaddle10", 50,1.,50.,500, 0., 5.);
	TH2F *Histpaddle11 = new TH2F("Histpaddle11", "Histpaddle11",50,1.,50.,500, 0., 5.);
	
	// histos for W and Cos
	int w, c;
	char wc[20];
	TH1F *HistWC[WBINS*CBINS];
	TH1F *HistWCKPM[WBINS*CBINS];
	TH1F *HistWCKPP[WBINS*CBINS];
	for(w = 0; w < WBINS; ++w)
	{
		for(c = 0; c < CBINS; ++c)
		{
			sprintf(wc, "W%dC%02d", w, c);
			HistWC[w*CBINS+c] = new TH1F(wc, "Data: Cos#theta binned in W(CMenery)", WCBIN, WCMIN, WCMAX);
			
			sprintf(wc, "W%dC%02d(Kppim #Sigma^{+})", w, c);
			HistWCKPM[w*CBINS+c] = new TH1F(wc, "Data: Cos#theta binned in W(CMenery)", WCBIN, WCMIN, WCMAX);

			sprintf(wc, "W%dC%02d(Kppip #Sigma^{-})", w, c);
			HistWCKPP[w*CBINS+c] = new TH1F(wc, "Data: Cos#theta binned in W(CMenery)", WCBIN, WCMIN, WCMAX);
		}
	}
	
	
	

//*************************       MAIN PROGRAM     *****************************

	for(int n_arg = optind; n_arg < __argc; n_arg++)
	{ 
		//loop over data files
   
		TFile inFile(__argv[n_arg]);

		// open each input file        
		if( TTree *p1 = (TTree*)inFile.Get("p0"))
		{
			p1->SetBranchAddress("Id", id);
			p1->SetBranchAddress("Q", q);
			p1->SetBranchAddress("Sector", sec);
			p1->SetBranchAddress("Paddle", npd);
			p1->SetBranchAddress("ngamma", &ngamma);
			p1->SetBranchAddress("e_gamma",egam);
			p1->SetBranchAddress("ngm", ngm);
			p1->SetBranchAddress("P_x",px);	
			p1->SetBranchAddress("P_y", py);
			p1->SetBranchAddress("P_z", pz);
			p1->SetBranchAddress("P_e", pe);
			p1->SetBranchAddress("Vx", vx);
			p1->SetBranchAddress("Vy", vy);
			p1->SetBranchAddress("Vz", vz);
			p1->SetBranchAddress("Tof", tof);
			p1->SetBranchAddress("Path", path);
			p1->SetBranchAddress("EID", &EID);
			p1->SetBranchAddress("TID", &TID);
			p1->SetBranchAddress("tFlag", &tFlag);
			
			Int_t nentries = (Int_t)p1->GetEntries();
			for (Int_t j=0;j<nentries;j++)
			{ 
				p1->GetEntry(j);	//Loop over entries
				
				//cout << "\t" << ngamma << "\t" << ngm[0] << "\t" << ngm[1]  << "\t" << egam[0] << "\t" << egam[1] << endl;
				//if(ngamma == 1) cout << ngam << "\t" << ngamma << "\t" << egam[0] << "\t" <<  egam[1] << endl;
				//if(ngamma > 1) cout << ngam << "\t" << ngamma << "\t" << egam[0] << "\t" <<  egam[1] << endl;
				//if(ngamma != 1) ngamma = 2;
				//cout << ngam << "\t" << ngamma << "\t" << egam[0] << "\t" << egam[1] << endl;
				//Histngm0->Fill(ngm[0]);
				//Histngm0->Fill(ngm[1]);
				//HistNgamma->Fill(ngamma);
								
				if(ngamma != 1) ngamma = 2;
				//cout << ngamma << "\t" << egam[0] << "\t" << egam[1] << endl;
				//if(ngamma >= 0)
				//{ cout << ngamma << endl;
				
				for (Int_t ngam = 0; ngam < ngamma; ngam++) //Loop over 1 and 2 photon case
				{
					// Set the target and photon values
					P4target.SetPxPyPzE(0.0,0.0,0.0,xmp);		
					P4pho.SetPxPyPzE(0.0,0.0,egam[ngam],egam[ngam]);
					//P4pho1.SetPxPyPzE(0.0,0.0,egam[1],egam[1]);
					P4inc =  P4target + P4pho;

					// Four-vector of the particles
					P4pip1.SetPxPyPzE(px[0],py[0],pz[0],pe[0]); 
					P4pim1.SetPxPyPzE(px[1],py[1],pz[1],pe[1]);  
					P4kp1.SetPxPyPzE(px[2],py[2],pz[2],pe[2]);

					// Vertices 
					V3pip1.SetXYZ(vx[0],vy[0],vz[0]);  
					V3pim1.SetXYZ(vx[1],vy[1],vz[1]);  
					V3kp.SetXYZ(vx[2],vy[2],vz[2]);

					//momentum vectors
					P3pip = P4pip1.Vect();
					P3pim = P4pim1.Vect();
					P3kp = P4kp1.Vect();
					
					//momentum magnitude
					momkp = P3kp.Mag();
					mompip = P3pip.Mag();
					mompim = P3pim.Mag();
					
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
					Bcm.SetZ((-egam[ngam]/(egam[ngam]+xmp)));	//eq46.4 pdg
					P4kpCM = P4kp1;
					P4phoCM = P4pho;
					P4kpCM.Boost(Bcm);		//Boost the k+ to CM frame
					P4phoCM.Boost(Bcm);		//Boost the photon to CM frame
					cosTkpCM = P4kpCM.CosTheta();		//cosTkpCM.Boost(Bcm);
				
					s_value = (P4inc).M2();		//Mandelstam s-variable
					cmeg = sqrt(s_value);

					// Fill Histograms
					
					ngam1 = ngam*1.0;
					
					HistNgam->Fill(ngam1);
					HistEgam->Fill(egam[ngam]);
					HistCMenergy->Fill(cmeg);
					HistcosTkpCM->Fill(cosTkpCM);
					HistEgamVsEID->Fill(EID,egam[ngam]);
					HistWVsCos->Fill(cosTkpCM,cmeg);

					HistMMsq->Fill((P4miss).M2());
					HistMM->Fill(xmmall);
					HistMMKp->Fill(xmmkp);
					HistM2pi->Fill(xm2pi);
					HistM2piK->Fill(xmall);
					HistMMvM2pi->Fill(xmm2pi,xm2pi);
					HistMMvM2piK->Fill(xmmall,xmall);
					
					//mom no cut
					Histmomkp->Fill(momkp);
					Histmompip->Fill(mompip);
					Histmompim->Fill(mompim);
					
					
					bool fiducial = false;
					
					if(clas::g12::g12_PosParticle_fiducial_cuts(P3kp.Mag(), P4kp1.Theta()*180/PI, P4kp1.Phi()*180/PI, "nominal")== true &&
						clas::g12::g12_PosParticle_fiducial_cuts(P3pip.Mag(), P4pip1.Theta()*180/PI, P4pip1.Phi()*180/PI, "nominal")== true &&
							clas::g12::g12_NegParticle_fiducial_cuts(P3pim.Mag(), P4pim1.Theta()*180/PI, P4pim1.Phi()*180/PI, "nominal")== true)
							{
								fiducial = true;
							}
					
					
					//Bad paddle knock out
					
					if(sec[1]==2)	Histpaddle10->Fill(npd[1],P3pip.Mag());
					
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
					
					if(padcut12 == true && sec[1]==2)	Histpaddle11->Fill(npd[1],P3pip.Mag());
					
					
					// Fill WC histograms
												
					if (cmeg > WMIN && cmeg < WMAX)
					{
						if (cosTkpCM > CMIN && cosTkpCM < CMAX)
						{
							w = int(abs((cmeg - WMIN) * WBINS / (WMAX - WMIN)));
							c = int(abs((cosTkpCM - CMIN) / CBINSIZE));
							HistWC[CBINS*w+c]->Fill(xmmkp); 
						}
					}
					
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
					
					// exclusive reaction (MM=0)
					if( xmmall > 0.9 && xmmall < 1.0 )
					{
						Hist1MMsq->Fill((P4miss).M2());
						Hist1M2pi->Fill(xm2pi);
						Hist1MM2pi->Fill(xmm2pi);
						Hist1MMvM2pi->Fill(xmm2pi,xm2pi);
						
						Hist1momkp->Fill(momkp);
						Hist1mompip->Fill(mompip);
						Hist1mompim->Fill(mompim);
						
						//Fiducal Histos without fiducial cut
						Hist_fidkp->Fill(P4kp1.Theta(), P4kp1.Phi());
						Hist_fidpip->Fill(P4pip1.Theta(), P4pip1.Phi());
						Hist_fidpim->Fill(P4pim1.Theta(), P4pim1.Phi());

						Hist_fidkpmomphi->Fill(P3kp.Mag(), P4kp1.Phi());
						Hist_fidpipmomphi->Fill(P3pip.Mag(), P4pip1.Phi());
						Hist_fidpimmomphi->Fill(P3pim.Mag(), P4pim1.Phi());
						
						for(sect = 1; sect < 7 ; ++sect)
						{
							if(sec[0] == sect)	Histfidpipmomthetasec[sect-1]->Fill(P3pip.Mag(), P4pip1.Theta());
							if(sec[1] == sect)	Histfidpimmomthetasec[sect-1]->Fill(P3pim.Mag(), P4pim1.Theta());
							if(sec[2] == sect)	Histfidkpmomthetasec[sect-1]->Fill(P3kp.Mag(), P4kp1.Theta());
						}
														
						Hist_fidkpang->Fill(P4kp1.Theta()*sin(P4kp1.Phi()), P4kp1.Theta()*cos(P4kp1.Phi()));
						Hist_fidpipang->Fill(P4pip1.Theta()* sin(P4pip1.Phi()), P4pip1.Theta()*cos(P4pip1.Phi()));
						Hist_fidpimang->Fill(P4pim1.Theta()*sin(P4pim1.Phi()), P4pim1.Theta()* cos(P4pim1.Phi()));
						
						//Paddle histos without cut
						for(sect = 1; sect < 7 ; ++sect)
						{
							if(sec[0] == sect)	Histpadpipmomsec[sect-1]->Fill(npd[0], P3pip.Mag());
							if(sec[1] == sect)	Histpadpimmomsec[sect-1]->Fill(npd[1], P3pim.Mag());
							if(sec[2] == sect)	Histpadkpmomsec[sect-1]->Fill(npd[2], P3kp.Mag());
						}
							
						//remove K0 events
						//if( fabs(xm2pi-0.495) > 0.015 && tFlag == 1 && fiducial == true && padcut12 == true && cutM == true) //tFlag for trigger condition 3 particles trigger in separate sector
						if( fabs(xm2pi-0.495) > 0.015 && tFlag == 1 && fiducial == true && padcut12 == true) //tFlag for trigger condition 3 particles trigger in separate sector
						//if( fabs(xm2pi-0.495) > 0.015)
						{
							Hist1MMkpim->Fill(xmmkppim);
							Hist1MMkpip->Fill(xmmkppip);
							
							
							// Fiducial Histos (with "nominal" fiducial cut) 
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
							
							
							if(xmmkppim > minsigp && xmmkppim < maxsigp)
							{
								if (fabs(xmmkppim - sigp) < fabs(xmmkppip - sigp))
								{
									Hist12MMKp->Fill(xmmkp);
									HistGlb12MMKp->Fill(xmmkp);
									
									Hist12momkp->Fill(momkp);
									Hist12mompip->Fill(mompip);
									Hist12mompim->Fill(mompim);
									
									if (cmeg > WMIN && cmeg <= WMAX)
									{
										if (cosTkpCM > CMIN && cosTkpCM <= CMAX)
										{
											w = int(abs((cmeg - WMIN) * WBINS / (WMAX - WMIN)));
											c = int(abs((cosTkpCM - CMIN) / CBINSIZE));
											HistWCKPM[CBINS*w+c]->Fill(xmmkp); 
										}
									}
								}
							}
							if(xmmkppip > minsigm && xmmkppip < maxsigm)
							{
								if (fabs(xmmkppip - sigm) < fabs(xmmkppim - sigm))
								{
									Hist13MMKp->Fill(xmmkp);
									HistGlb13MMKp->Fill(xmmkp);
									
									Hist13momkp->Fill(momkp);
									Hist13mompip->Fill(mompip);
									Hist13mompim->Fill(mompim);
									
									if (cmeg > WMIN && cmeg <= WMAX)
									{
										if (cosTkpCM > CMIN && cosTkpCM <= CMAX)
										{
											w = int(abs((cmeg - WMIN) * WBINS / (WMAX - WMIN)));
											c = int(abs((cosTkpCM - CMIN) / CBINSIZE));
											HistWCKPP[CBINS*w+c]->Fill(xmmkp); 
										}
									}
								}
							}
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
				}
				//}

			}//for entries (FORloop)
		}//if TTree *p1 = (TTree*)

		else
		{
			cout << "File has no TTree or p1 does not exist!!" << endl;
		} 
		//cout << __argv[n_arg] << endl; 
		cout << n_arg <<"/"<< __argc << endl;
	}//n_arg

	outFile.Write(); // write to the output file
	outFile.Close(); // close the output file
}
			
/*				if(ngamma == 2 && (egam[0] - egam[1]) < 0.001)*/
/*				{*/
/*					//cout << ngamma << endl;*/
/*					ngamma = 1;*/
/*				}*/
/*				else*/
/*				{*/
/*					cout << ngamma << endl;*/
/*					ngamma = ngamma;*/
/*				}*/
				
/*				for (Int_t ngam = 0; ngam < ngamma; ngam++) //Loop over 1 and 2 photon case*/
/*				{*/
/*					if(ngam == 2 && (egam[0] - egam[1]) < 0.001)*/
/*					{*/
/*						//cout << ngamma << endl;*/
/*						ngam = 1;*/
/*					}*/
/*					else*/
/*					{*/
/*						//cout << ngamma << endl;*/
/*						ngam = ngam;*/
/*					}*/
/*				}*/
