#include <iostream>
extern "C" {
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <errno.h>
#include <ntypes.h>
#include <bostypes.h>
#include <clas_cern.h>
#include <particleType.h>
#include <kinematics.h>
#include <pdgutil.h>
#include <pid.h>
#include <dataIO.h>
#include <itape.h>
#include <vertex.h>
#include <trk.h>
#include <makebanks.h>
}

#include <Vec.h>
#include <lorentz.h>
#include <pputil.h>
#include <clasEvent.h>

//Import g12 corrections
#include "g12_corrections/g12_ecor.hpp"
#include "g12_corrections/g12_pcor.hpp"
#include "g12_corrections/g12_tofres_knockout.hpp"
//#include "g12_corrections/g12_TOF_knockout.hpp"
#include "g12_corrections/g12_NegParticle_fiducial_cuts.hpp"
#include "g12_corrections/g12_PosParticle_fiducial_cuts.hpp"
#include "g12_corrections/g12_trackEfficiency_proton.hpp"
#include "g12_corrections/g12_trackEfficiency_pip.hpp"
#include "g12_corrections/g12_trackEfficiency_pim.hpp"

//Import root stuff
#include "TNtuple.h" 
#include "TTree.h" 
#include "TFile.h"
#include "TBranch.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"



/*bitwise drop flags*/
#define DROP_RAW BIT(0)
#define DROP_DC0 BIT(1)
#define DROP_DC1 BIT(2)
#define DROP_HBLA BIT(3)
#define DROP_TBLA BIT(4)
#define DROP_HBTB BIT(5)
#define DROP_SC BIT(6)
#define DROP_EC BIT(7)
#define DROP_HBID BIT(8)
#define DROP_CL01 BIT(9)
#define DROP_SEB BIT(10)
#define DROP_TBID BIT(11)
#define DROP_HDPL BIT(12)
#define DROP_LAC BIT(13)
#define DROP_CC BIT(14)
#define DROP_ST BIT(15)
#define DROP_DHCL BIT(16)
#define DROP_TAGR BIT(17)

//Define BOS banks
#define CC_BANKS "CCRCCC01"
#define SEB_BANKS "HEVTEVNTDCPBSCPBCCPBUNUSEVHBTRKSSTPBTGPBLCPB"
//#define TAGGER_BANKS "TAGRTAGIPSO "
#define SC_BANKS "SC1 SCR SCRC"
#define EC_BANKS "EC01ECHBECPIECPCECPBEC  EC1 "
#define ST_BANKS "ST1 STR "
#define REGION1_BANKS "RGLK"



#define BUFSIZE 200000

static int requestedRawData = 0;

/* ----------- Function prototypes ---------------- */

static void signalINT(int isig);
static void signalSEGV(int isig);
int StartRun(int runNo);
int EndRun(int run);
int ProcessEvent(clasEvent &event,int verbose,int silent);
void installFaultHandlers();
void DropList(int DropFlag);
void pHeader(clasEvent &event);

double beta(const fourVec &p);
double beta(const fourVec &p,double m);
double beta(const threeVec &p,double m);
double deltaBeta(double l,double t);
double chisq(const threeVec &p,double mass,double len,double t);
double gamma(double beta);
double mass(double p,double beta);
double twoThreeVecAngle(threeVec, threeVec);
int npho( clasTAGR_t *TAGR, float tpho, int nb);  //jb npho2
tagr_t *my_get_photon_tagr(clasTAGR_t *TAGR,clasBID_t *TBID, threeVec myVert); //jb v13
//double pcor(const double p, const double phi, const int part_id);
//fourVec pcorFourVec( fourVec particle, int part_id );
fourVec pcorFourVec_new( fourVec particle, int part_id );
void printView(Float_t *array);
int GetHelicity(clasHEVT_t *HEVT);


//Global definitions
TFile *f;
TNtuple *kpipi;
TTree *tree;
const int nData = 80;
Float_t nt_value[nData];
//Float_t p1_mat[5][5];
//Float_t p2_mat[5][5];
//TBranch *bp1_mat;
//TBranch *bp2_mat;

clas::g12::MomentumCorrection pcor((char*) "");

extern particleDataTable PDGtable;


//declare the bos common 
BOSbank bcs_;
BOSbank wcs_;
#define LIGHT_SPEED 29.979 // speed of light cm/nsec
#define RAD_TO_DEGREE 57.295820614 

int CurrentRun = 0;
int CurrentEvent = 0;
int rootStatus;

void PrintUsage(char *processName){
    cerr << processName << " [options] file1 file2 ..." << endl;
    //cerr << "\t[-cm#]\t\tSet mini-torus current to #" << endl;
    //cerr << "\t[-ct#]\t\tSet torus current to #" << endl;
    cerr << "\t[-L]\t\tPrint labels" << endl;
    //cerr << "\t[-v]\t\tverbose mode" << endl;
    cerr << "\t[-ofile]\twrite selected events to bos file" << endl;
    cerr << "\t[-rfile]\twrite to rootfile" << endl;
    //cerr << "\t[-d]\t\tdebug mode" << endl;
    //cerr << "\t[-D#]\t\tSet drop list to #" << endl;
    cerr << "\t[-M#]\t\tProcess only # events" << endl;
    //cerr << "\t[-s]\t\tsilent mode" << endl;
    //cerr << "\t[-N]\t\tDo not build TDPL bank (unless necessary)" <
    cerr << "\t[-h]\t\tPrint this message" << endl;
}


main(int argc,char **argv){

    FILE *fp = NULL;
    int max = 0;
    
    int verbose = 0;
    int Silent = 0;
    
    int writeMode = 0;
    int forceWrite = 0;
    
    int reBuild = -1;
    int doBuild = 0;
    int rebuildTBID = 0;
    
    int ret = 1;
    int rejectFlag = 0;
    int NoFault = 0;
    int pEvent;
    
    int Nevents = 0;
    int nfile = 0;
    int Nwrite = 0;
    int nProc = 0;
    time_t time1; 
    float rate; 
    int i;
    char *argptr;
    char *word;
    int Dispatch = 0;
    unsigned int triggerMask = 0;
    
    int useGPID = 0;
    int useSEB = 0; 
    int partbank = 1;
    char *BOSout = (char*) "output.bos";
    int debug = 0;
    
    // Dispatcher
    char *server = NULL;
    int dispatcherPipelinePrefill = 5;
    
    int allFlag = 0;
    int dropTDPL = 0;
    
    


    // drop sector list
    
    vector<int> dropSector;
    
    // itape stuff
    FILE *finput = stdin;
    //gDebug=2;
    // bos stuff
    rootStatus = 0;
    int textStatus = 0;
    char* rootfile;

    clasOutput coutput;
    int status = 0;
    int OutputUnitNo = 9;
    extern char *optarg;
    int c; 
    extern int optind;
    int MaxBanks = 1000; 
    char  out[300];
    
    int DropFlag = 0x0; /*don't drop anything*/
    
    for (int i = 0; i < argc; ++i)
        cerr << argv[i] << " ";
    cerr << endl;

    while ( (c = getopt(argc,argv, "o:LhM:r:")) != -1 ) {
        switch (c) {
                
            case 'G':
                useGPID = 1;
                break;
                
            case 'S':
                useSEB = 1;
                break;
                
            case 'P':
                partbank = atoi(optarg);
                break;
                
            case 'F':
                forceWrite = 1;
                cout<<"FORCE WRITE"<<endl;
                break;
                
            case 'R':
                // write out rejected events
                rejectFlag = 1;
                break;
                
            case 'v':
                verbose = 1;
                break;
                
            case 'o':
                BOSout = optarg;
                status = coutput.open(BOSout, OutputUnitNo);
                cout <<"Status:["<<status<<"] " <<"Output file: " << BOSout << endl;
                break;

            case 'r':  
                rootfile = optarg;
                rootStatus = 1;
                break;

            case 'f':
                // do not install fault handlers
                NoFault = 1;
                break;
                
            case 'd':
                debug = 1;
                break;
                
            case 'D':
                DropFlag = strtoul(optarg,NULL,0);
                break;
                
            case 'M':
                max = atoi(optarg);
                break;
                
            case 't':
                triggerMask = strtoul(optarg,NULL,0);
                break;
                
            case 's':
                Silent = 1;
                break;
                
            case 'h':
                PrintUsage(argv[0]);
                exit(0);
                break;
                
            default:
                cerr << "Unrecognized argument: " << argptr << endl;;
                break;
        }
    }

    PDGtable.initialize();

    //Initialize root file
    if(rootStatus){
        f = new TFile((char *) rootfile,"recreate");
        kpipi = new TNtuple("kpipi","data","run:evt:pxkp1pcor:pykp1pcor:pzkp1pcor:pxpip1pcor:pypip1pcor:pzpip1pcor:pxpim1pcor:pypim1pcor:pzpim1pcor:nn:n0:np:ebeamcor:ebeam:vtime:stvtime:delvt:kp1toflen:kp1beta:kp1Beta:pip1toflen:pip1beta:pip1Beta:pim1toflen:pim1beta:pim1Beta:vxkp1:vykp1:vzkp1:vxpip1:vypip1:vzpip1:vxpim1:vypim1:vzpim1:seckp1:secpip1:secpim1:padkp1:padpip1:padpim1:tpho:helicity:nkp:npip:npim:pxpip2pcor:pypip2pcor:pzpip2pcor:pxpim2pcor:pypim2pcor:pzpim2pcor:pip2toflen:pip2beta:pip2Beta:vxpip2:vypip2:vzpip2:secpip2:padpip2:pim2toflen:pim2beta:pim2Beta:vxpim2:vypim2:vzpim2:secpim2:padpim2:kp1Efficiency:pip1Efficiency:pim1Efficiency:pip2Efficiency:pim2Efficiency:normalfidcut:tightfidcut:loosefidcut:tofko1:tofko2");
    }
    
    if (!NoFault)
        installFaultHandlers();
    
    // Initialize BOS
    bnames_(&MaxBanks);
    initbos();
    configure_banks(stderr,0);

    for (i = 1; i < argc; ++i) {
        argptr = argv[i];
        if (*argptr != '-') {
            clasEvent event(argptr,&bcs_,1,0);
            cerr << "initialize:\t" << argptr << endl;
            
            Dispatch = isDispatcher(argptr);
            if (event.status()) {
                ret = 1;
                
                while ((max ? Nevents < max : 1) && ret) {
                    clasHEAD_t *HEAD;	  
                    ret = event.read(1);  // still one event?
                    
                    if (ret == DISIO_DATA) {
                        if (event.type() == 1) {
                            if (HEAD = (clasHEAD_t *)getBank(&bcs_, "HEAD")) {
                                int runno=HEAD->head[0].nrun;
                                CurrentRun = HEAD->head[0].nrun;
                                CurrentEvent =  HEAD->head[0].nevent;
                                StartRun(runno); 

                                pEvent = ProcessEvent(event,verbose,Silent);
                                nProc++;
                                if (status && pEvent) {
                                    coutput.write(&bcs_);
                                    Nwrite++;
                                }
                                
                            }
                            event.clean();
                            Nevents++;
                        }
                        
                    }
                    else if (ret == DISIO_COMMAND) {
                        cerr << "Message from Giant Head: " << getBuffer() << endl;;   
                    }
                    
                }
                cerr << "\nTotal number of itape records:\t" << Nevents << endl;
                cerr << "\tTotal number of records written:\t" << Nwrite << endl;
                
            }
            else {
                cerr << "Unable to open " << argptr << endl;
            }
        }
    }
    if (Dispatch)
        disIO_command("FINISHED");
    cerr << "\nTotal number of itape records:\t" << Nevents << endl;
    if(status) {
        coutput.write(&bcs_,(char*)"0");
        coutput.close();
    }
    if(rootStatus){
       kpipi->Write();
        f->Close();
    }
}

int StartRun(int runNo){  
    int static CurrentRun = -1;
    if (CurrentRun != runNo) {
        //    dc_xvst_init_(&runNo);
        //    vertex_brun(runNo);
        //  sc_begin_run(runNo);
        //   dc_begin_run(runNo);
        CurrentRun = runNo; 
    }
    return 0;
}





int ProcessEvent(clasEvent &evt,int verbose,int silent){


	fourVec kp1, pip1, pim1, pip2, pim2;
    lorentzTransform L;
	clasParticle kp1x, pip1x, pim1x, pip2x, pim2x;
	fourVec beam,target,beamecor;
    threeVec Cross;
	//fourVec X;
	int ret;
	double eBeamLo = 4.4;
	double eBeamHi = 7.0;
    float tpho;
	//beam1 = evt.getBeam(eBeamLo, eBeamHi);
    target = evt.target().get4P();
    clasTAGR_t *TAGR = (clasTAGR_t *)getBank(&bcs_,"TAGR"); 
    clasBID_t *TBID  =  (clasBID_t *)getBank(&bcs_,"TBID"); 
    clasHEVT_t *HEVT = (clasHEVT_t *)getBank(&bcs_,"HEVT");
    threeVec evtVert=evt.V();
    tagr_t *tagr = my_get_photon_tagr(TAGR, TBID, evtVert);  //jb v13
    double beamp = -1000.0;
    if(tagr){
        //eid = tagr->e_id;
        beamp = tagr->erg;
        //jb v12    wo_beamp = tagr->erg;
        tpho=tagr->tpho; //jb npho
    }

    beam.set(beamp,threeVec(0.0,0.0,beamp));
    if( (evt.N(KPlus)==1) && (evt.N(PiPlus)>=1) && (evt.N(PiMinus)>=1) ){
    kp1x = evt.cp(KPlus, 1);
    pip1x = evt.cp(PiPlus, 1);
    pim1x = evt.cp(PiMinus, 1);
    if(evt.N(PiPlus) == 2)
    {
        pip2x = evt.cp(PiPlus, 2);
    }
    if(evt.N(PiMinus) == 2)
    {
        pim2x = evt.cp(PiMinus, 2);
    }
    //
     ret=0;
     if(TAGR){
            for(int k=0;k<TAGR->bank.nrow;k++){
                if (abs(tpho - TAGR->tagr[k].tpho) < (0.05+2.004)) {
                //if (fabs(tpho - TAGR->tagr[k].tpho) < 1.00) {
		    ret=1;
                beam.set(TAGR->tagr[k].erg,threeVec(0.0,0.0,TAGR->tagr[k].erg));
                    kp1 = kp1x.p();
                    pip1 = pip1x.p();
                    pim1 = pim1x.p();
                    pip2 = pip2x.p();
                    pim2 = pim2x.p();
              double ebeamcor = clas::g12::Beam_correction(evt.run(),beam.t());
                        beamecor.set(ebeamcor,threeVec(0.0,0.0,ebeamcor));

                	//Run info
                        nt_value[0] = evt.run();
                        nt_value[1] = evt.event() ;

                        //Apply eLoss
                        evt.eLoss((char*)"g11a",1);

                        //Apply Momentum Corrections
                        fourVec kp1_pcor, pip1_pcor, pim1_pcor, pip2_pcor, pim2_pcor;
                        kp1_pcor = pcorFourVec_new(kp1, 11);
                        pip1_pcor = pcorFourVec_new(pip1, 8);
                        pim1_pcor = pcorFourVec_new(pim1, 9);
                        pip2_pcor = pcorFourVec_new(pip2, 8);
                        pim2_pcor = pcorFourVec_new(pim2, 9);

                        //lab angles of momentum corrected protons
                        nt_value[2] = kp1_pcor.V().x();
                        nt_value[3] = kp1_pcor.V().y();
                        nt_value[4] = kp1_pcor.V().z();
                        nt_value[5] = pip1_pcor.V().x();
                        nt_value[6] = pip1_pcor.V().y();
                        nt_value[7] = pip1_pcor.V().z();
                        nt_value[8] = pim1_pcor.V().x();
                        nt_value[9] = pim1_pcor.V().y();
                        nt_value[10] = pim1_pcor.V().z();

                        
                        //Number of negative, neutral, and positive tracks
                        nt_value[11] = evt.N(-1);    
                        nt_value[12] = evt.N(0);    
                        nt_value[13] = evt.N(1);    

                        //Beam energy w/correction
                        nt_value[14] = ebeamcor;

                        //Beam Energy
                        //nt_value[18] = beam.t();
                        nt_value[15] = TAGR->tagr[k].erg;

                        //Event Timing
                        nt_value[16] = evt.vtime();
                        nt_value[17] = evt.stVtime();
                        nt_value[18] = tpho - TAGR->tagr[k].tpho;

                        //Path Length, beta, and Beta.
                        nt_value[19] = kp1x.tofLength();
                        nt_value[20] = kp1x.beta();
                        nt_value[21] = kp1x.Beta();

                        nt_value[22] = pip1x.tofLength();
                        nt_value[23] = pip1x.beta();
                        nt_value[24] = pip1x.Beta();
                        

                        nt_value[25] = pim1x.tofLength();
                        nt_value[26] = pim1x.beta();
                        nt_value[27] = pim1x.Beta();

                        //Event Vertex
                        nt_value[28] = kp1x.x();
                        nt_value[29] = kp1x.y();
                        nt_value[30] = kp1x.z();
                        nt_value[31] = pip1x.x();
                        nt_value[32] = pip1x.y();
                        nt_value[33] = pip1x.z();
                        nt_value[34] = pim1x.x();
                        nt_value[35] = pim1x.y();
                        nt_value[36] = pim1x.z();

			//Sector and Paddle
                        nt_value[37] = kp1x.scPaddleId_sector();
                        nt_value[38] = pip1x.scPaddleId_sector();
                        nt_value[39] = pim1x.scPaddleId_sector();
                        nt_value[40] = kp1x.scPaddleId();
                        nt_value[41] = pip1x.scPaddleId();
                        nt_value[42] = pim1x.scPaddleId();

                        nt_value[43] = tpho;
                        nt_value[44] = -1000.;
			if(HEVT) {
                        nt_value[44] = GetHelicity(HEVT);
			}
                         //Number of KPlus, PiPlus, PiMinus
                        nt_value[45] = evt.N(KPlus);
                        nt_value[46] = evt.N(PiPlus);
                        nt_value[47] = evt.N(PiMinus);

                        nt_value[48] = pip2_pcor.V().x();
                        nt_value[49] = pip2_pcor.V().y();
                        nt_value[50] = pip2_pcor.V().z();

                        nt_value[51] = pim2_pcor.V().x();
                        nt_value[52] = pim2_pcor.V().y();
                        nt_value[53] = pim2_pcor.V().z();
                        nt_value[54] = -1000;
                        nt_value[55] = -1000;
                        nt_value[56] = -1000;
                        nt_value[57] = -1000;
                        nt_value[58] = -1000;
                        nt_value[59] = -1000;
                        nt_value[60] = -1000;
                        nt_value[61] = -1000;
                        

                        if(evt.N(PiPlus) == 2)
                        {
                            nt_value[54] = pip2x.tofLength();
                            nt_value[55] = pip2x.beta();
                            nt_value[56] = pip2x.Beta();
                        
                            nt_value[57] = pip2x.x();
                            nt_value[58] = pip2x.y();
                            nt_value[59] = pip2x.z();

                            nt_value[60] = pip2x.scPaddleId_sector();
                            nt_value[61] = pip2x.scPaddleId();
                        }

                        nt_value[62] = -1000;
                        nt_value[63] = -1000;
                        nt_value[64] = -1000;
                        nt_value[65] = -1000;
                        nt_value[66] = -1000;
                        nt_value[67] = -1000;
                        nt_value[68] = -1000;
                        nt_value[69] = -1000;
                       

                        if(evt.N(PiMinus) == 2)
                        {
                            nt_value[62] = pim2x.tofLength();
                            nt_value[63] = pim2x.beta();
                            nt_value[64] = pim2x.Beta();
                        
                            nt_value[65] = pim2x.x();
                            nt_value[66] = pim2x.y();
                            nt_value[67] = pim2x.z();

                            nt_value[68] = pim2x.scPaddleId_sector();
                            nt_value[69] = pim2x.scPaddleId();
                        }

                        double P_Vz = kp1x.z();
                        double P_Ptot = kp1.V().r();
                        double P_Phi = kp1.V().phi();
                        double P_Theta = kp1.V().theta();
                        char* Units = (char*)"radian";
                        char* Effdir = (char*)"g12_corrections/EFFICIENCY_HISTS/";
                        nt_value[70] = clas::g12::g12_proton_efficiency(P_Vz, P_Ptot, P_Phi, P_Theta, Units, Effdir);
                        P_Vz = pip1x.z();
                        P_Ptot = pip1.V().r();
                        P_Phi = pip1.V().phi();
                        P_Theta = pip1.V().theta();
                        nt_value[71] = clas::g12::g12_pip_efficiency(P_Vz, P_Ptot, P_Phi, P_Theta, Units, Effdir);
                        P_Vz = pim1x.z();
                        P_Ptot = pim1.V().r();
                        P_Phi = pim1.V().phi();
                        P_Theta = pim1.V().theta();
                        nt_value[72] = clas::g12::g12_pim_efficiency(P_Vz, P_Ptot, P_Phi, P_Theta, Units, Effdir);

                        nt_value[73] = -1000;
                        if(evt.N(PiPlus) == 2)
                        {
                            P_Vz = pip2x.z();
                            P_Ptot = pip2.V().r();
                            P_Phi = pip2.V().phi();
                            P_Theta = pip2.V().theta();
                            nt_value[73] = clas::g12::g12_pip_efficiency(P_Vz, P_Ptot, P_Phi, P_Theta, Units, Effdir);
                        }
                        nt_value[74] = -1000;
                        if(evt.N(PiMinus) == 2)
                        {
                            P_Vz = pim2x.z();
                            P_Ptot = pim2.V().r();
                            P_Phi = pim2.V().phi();
                            P_Theta = pim2.V().theta();
                            nt_value[74] = clas::g12::g12_pim_efficiency(P_Vz, P_Ptot, P_Phi, P_Theta, Units, Effdir);
                        }

                        // // nominal, tight and loose fiducial cuts
                        if ((clas::g12::g12_PosParticle_fiducial_cuts(kp1_pcor.r(), kp1_pcor.theta()*180/PI, kp1_pcor.phi()*180/PI, "nominal") == true) && (clas::g12::g12_PosParticle_fiducial_cuts(pip1_pcor.r(), pip1_pcor.theta()*180/PI,  pip1_pcor.phi()*180/PI, "nominal") == true || clas::g12::g12_PosParticle_fiducial_cuts(pip2_pcor.r(), pip2_pcor.theta()*180/PI,  pip2_pcor.phi()*180/PI, "nominal") == true) && (clas::g12::g12_NegParticle_fiducial_cuts(pim1_pcor.r(), pim1_pcor.theta()*180/PI, pim1_pcor.phi()*180/PI, "nominal") == true || clas::g12::g12_NegParticle_fiducial_cuts(pim2_pcor.r(), pim2_pcor.theta()*180/PI, pim2_pcor.phi()*180/PI, "nominal") == true))
                        {
                            nt_value[75] = 1;
                        }
                        else
                        {
                            nt_value[75] = 0;
                        }
                        if ((clas::g12::g12_PosParticle_fiducial_cuts(kp1_pcor.r(), kp1_pcor.theta()*180/PI, kp1_pcor.phi()*180/PI, "tight") == true) && (clas::g12::g12_PosParticle_fiducial_cuts(pip1_pcor.r(), pip1_pcor.theta()*180/PI,  pip1_pcor.phi()*180/PI, "tight") == true || clas::g12::g12_PosParticle_fiducial_cuts(pip2_pcor.r(), pip2_pcor.theta()*180/PI,  pip2_pcor.phi()*180/PI, "tight") == true) && (clas::g12::g12_NegParticle_fiducial_cuts(pim1_pcor.r(), pim1_pcor.theta()*180/PI, pim1_pcor.phi()*180/PI, "tight") == true || clas::g12::g12_NegParticle_fiducial_cuts(pim2_pcor.r(), pim2_pcor.theta()*180/PI, pim2_pcor.phi()*180/PI, "tight") == true))
                        {
                            nt_value[76] = 1;
                        }
                        else
                        {
                            nt_value[76] = 0;
                        }
                        if ((clas::g12::g12_PosParticle_fiducial_cuts(kp1_pcor.r(), kp1_pcor.theta()*180/PI, kp1_pcor.phi()*180/PI, "loose") == true) && (clas::g12::g12_PosParticle_fiducial_cuts(pip1_pcor.r(), pip1_pcor.theta()*180/PI,  pip1_pcor.phi()*180/PI, "loose") == true || clas::g12::g12_PosParticle_fiducial_cuts(pip2_pcor.r(), pip2_pcor.theta()*180/PI,  pip2_pcor.phi()*180/PI, "loose") == true) && (clas::g12::g12_NegParticle_fiducial_cuts(pim1_pcor.r(), pim1_pcor.theta()*180/PI, pim1_pcor.phi()*180/PI, "loose") == true || clas::g12::g12_NegParticle_fiducial_cuts(pim2_pcor.r(), pim2_pcor.theta()*180/PI, pim2_pcor.phi()*180/PI, "loose") == true))
                        {
                            nt_value[77] = 1;
                        }
                        else
                        {
                            nt_value[77] = 0;
                        }

                        //TOF knockout bad occupancy and resolution 
                        // if((clas::g12::pass_TOFKO(kp1x) == 1) && (clas::g12::pass_TOFKO(pip1x) == 1 || clas::g12::pass_TOFKO(pip2x) == 1) && (clas::g12::pass_TOFKO(pim1x) == 1 || clas::g12::pass_TOFKO(pim2x) == 1))
                        if((clas::g12::pass_TOFKO(kp1x) == 1) && (clas::g12::pass_TOFKO(pip1x) == 1) && (clas::g12::pass_TOFKO(pim1x) == 1))
                        {
                            nt_value[78] = 1;
                        }
                        else
                        {
                            nt_value[78] = 0;
                        }



                        nt_value[79] = -1000;

                        if(evt.N(PiMinus) == 2 || evt.N(PiMinus) == 2)
                        {
                            if((clas::g12::pass_TOFKO(kp1x) == 1 && clas::g12::pass_TOFKO(pip2x) == 1 && clas::g12::pass_TOFKO(pim1x) == 1) || (clas::g12::pass_TOFKO(kp1x) == 1 && clas::g12::pass_TOFKO(pip1x) == 1 && clas::g12::pass_TOFKO(pim2x) == 1) || (clas::g12::pass_TOFKO(kp1x) == 1 && clas::g12::pass_TOFKO(pip2x) == 1 && clas::g12::pass_TOFKO(pim2x) == 1))
                            {
                                nt_value[79] = 1;
                            }
                            // else
                            // {
                            //     nt_value[79] = 0;
                            // }
                        }
                        

                        if(rootStatus){
                            f->cd();
                            kpipi->Fill(nt_value);      
                        }else{
                            printView(nt_value);
                        }
                    }
                }
            }
        }
	return ret;
}

//*******************jb npho********************************************
int  npho( clasTAGR_t *TAGR, float tpho, int nb){  //jb npho func
    int nph=0; //jb npho
    for(int k=0;k<TAGR->bank.nrow;k++){
        if (abs(tpho - TAGR->tagr[k].tpho) < (0.05+2.004*nb)) {  //jb npho  
            nph=nph+1;
        }
    }
    return(nph);
}

int GetHelicity(clasHEVT_t *HEVT){
    int helicity = 0;
    int readout = HEVT->hevt[0].trgprs;
    if(readout > 0) helicity = -1;
    if(readout < 0) helicity = 1;
    return helicity;
}

//tagr_t *my_get_photon_tagr(clasTAGR_t *TAGR,clasBID_t *TBID){ 
tagr_t *my_get_photon_tagr(clasTAGR_t *TAGR,clasBID_t *TBID, threeVec myVert){  //jb v13
    
    /* This routine only works for time-based tracks! */  
    float best_diff=ST_TAG_COINCIDENCE_WINDOW;
    float tprop=0.0;
    tagr_t *tagr = NULL;
    clasTBTR_t *TBTR=(clasTBTR_t *)getBank(&bcs_,"TBTR");
    float g11targz=-10; //this was for g11
    float g12targz=-90; //jb v11
    int i, j; 
    
    
    
    /* Exit from function if missing the requisite banks... */
    if(!TAGR || !TBTR || !TBID) return(NULL);
    
    for (i=0;i<TBID->bank.nrow;i++){
        int trk_ind=TBID->bid[i].track;
        if(trk_ind){     
            //jb v13      tprop=(TBTR->tbtr[trk_ind-1].vert.z - g12targz)/SPEED_OF_LIGHT*1e9;  //jb v11 i think i need to take this out
            tprop=(myVert.z() - g12targz)/LIGHT_SPEED;  //jb v13
            //cout<<"myVert.z()=" << myVert.z() <<" :  TBTR->tbtr[trk_ind-1].vert.z="<< TBTR->tbtr[trk_ind-1].vert.z<<endl; //jb v13 just a check
            if (TBID->bid[i].st.stat){       
                
                float mytime=-1000; //jb v14
                float myenergy=-1000; //jb v14
                
                for(j=0;j<TAGR->bank.nrow;j++){
                    float diff=fabs(TBID->bid[i].st.vtime-(TAGR->tagr[j].tpho+tprop));
                    //jb v14     if (diff<ST_TAG_COINCIDENCE_WINDOW&&diff<best_diff
                    //jb v14	&& (TAGR->tagr[j].stat==7 || TAGR->tagr[j].stat==15)){    	   
                    //cout <<abs(TAGR->tagr[j].tpho - tagr->tpho)<<endl;
                    if (( diff<ST_TAG_COINCIDENCE_WINDOW&&diff<best_diff )  &&  (TAGR->tagr[j].stat==7 || TAGR->tagr[j].stat==15)){  //jb v14
                        best_diff=diff;
                        tagr=&(TAGR->tagr[j]);
                        mytime=tagr->tpho; //jb v14
                        myenergy=tagr->erg; //jb v14
                        
                    }
                    
                } 
                
            } 
        } 
    }
    return(tagr);
}



double twoThreeVecAngle(threeVec v1, threeVec v2){
    double dot_product = v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
    double magnitude = v1.r() * v2.r();                                      
    return dot_product/magnitude;
}

int EndRun(int run){
	return 0;
}

int GetDispatcherCommand(){
    int ret;
    int maxSend = 2;
    int OkToSendData = 1;
    int WaitForALLFINISHED = 0;
    int runNo;
    
    /* wait for command from Dispatcher */
    
    fd_set rset;
    fd_set wset;
    struct timeval timeout;
    int maxfd = 0;
    
    FD_ZERO(&rset);
    FD_ZERO(&wset);
    
    FD_SET(disIO_socket,&rset); if (disIO_socket > maxfd) maxfd = disIO_socket;
    
    if (OkToSendData && (requestedRawData > 0)){
        FD_SET(disIO_socket,&wset); if (disIO_socket > maxfd) maxfd = disIO_socket;
    }
    
    timeout.tv_sec  = 1;
    timeout.tv_usec = 0;
    
    if (OkToSendData && (requestedRawData > 0)){
        timeout.tv_sec  = 0;
        timeout.tv_usec = 0;
    }
    
    ret = select(maxfd + 1, &rset, &wset, NULL, &timeout);
    if (ret < 0){
        cerr << "DisFeedDAQ: Error: select() returned " <<  ret << "errno: " << errno << " " << strerror(errno) << endl;
        //	  exitFlag = 1;
        exit(0);
    }    
    
    /* ok, we got a command. Now parse it */
    static char *msg = NULL;
    static int msglen = 0;
    char *cmd0;
    char *word;
    int maybeUnexpected = 0;
    
    if (msg)
        free(msg);
    
    msg = NULL;
    
    ret = disIO_readRAW_alloc((void **)&msg,&msglen,0);
    
    if (msg) {
        
        word = strtok(msg,":");
        cmd0 = word;
        if (word) {
            cerr << "COMMAND: " << cmd0 << endl;
            
            if (strcmp(word,"NOP") == 0){
                /* do nothing */
            }
            else if (strcmp(word,"PING") == 0){
                printf("DisFeedDAQ: Command from Dispatcher: %s\n",word);
                disIO_command("PING-REPLY");
            }
            else if (strcmp(word,"REQUEST_DATA") == 0){
                int imore;
                
                /* fprintf(stderr,"DisFeedDAQ: Command from Dispatcher: %s\n",word); */
                
                word = strtok(NULL,":");
                
                if (word)
                    imore = strtoul(word,NULL,0);
                else
                    imore = 1;
                
                /* printf("REQUEST_DATA: more: %d, requested events: %d, sent: %d\n",imore,requestedRawData,isent); */
                
                requestedRawData += imore;
            }
            else if (strcmp(word,"MAXSEND") == 0)
            {
                cerr << "DisFeedDAQ: Command from Dispatcher: " << word << endl;
                
                word = strtok(NULL,":");
                
                if (word)
                    maxSend = strtoul(word,NULL,0);
                
                printf("DisFeedDAQ: New maxSend is %d\n",maxSend);
            }
            else if (strcmp(word,"ALLFINISHED") == 0)
            {
                if (WaitForALLFINISHED)
                {
                    cerr << "DisFeedDAQ: Command ALLFINISHED from Dispatcher: " << word << endl;
                    //	    SendBeginTape(runNo);
                    OkToSendData = 1;
                    WaitForALLFINISHED = 0;
                }
                else
                {
                    fprintf(stderr,"DisFeedDAQ: Unexpected command ALLFINISHED from Dispatcher was ignored.\n");
                }
            }
            else if (strcmp(word,"QUIT") == 0)
            {
                fprintf(stderr,"DisFeedDAQ: Command QUIT from Dispatcher: %s\n",word);
                exit(0);
            }
            else if (strcmp(word,"EXIT") == 0)
            {
                fprintf(stderr,"DisFeedDAQ: Command EXIT from Dispatcher: %s\n",word);
                exit(0);
            }
            else
            {
                fprintf(stderr,"DisFeedDAQ: Unexpected command from Dispatcher: [%s] was ignored.\n",msg);
            }
        }
    }
    else {
        cerr << "Received empty message from the Dispatcher" << endl;
    }
    
    
    
    return(ret);
}
extern "C" {
    int ir_isnan_(float *x){
        return(0);
    }
    
    int ir_isinf_(float *x){ 
        return(0);
    }
    
    
}


void pHeader(clasEvent &event)
{
    int trig = event.trig() & (int) 0xf; 
    cout << event.run() << " " << event.event() << " " << trig << " " ;
    
}


void DropList(int DropFlag)
{
    
    /* Mask off banks according to DropFlag*/
    
    if (DropFlag & DROP_RAW) bankList(&bcs_, (char*)"E-", (char*) "R");  
    if (DropFlag & DROP_DC0) bankList(&bcs_, (char*)"E-", (char*) "DC0 ");
    if (DropFlag & DROP_DC1) bankList(&bcs_, (char*)"E-", (char*) "DC1 ");
    if (DropFlag & DROP_HBLA) bankList(&bcs_, (char*)"E-", (char*) "HBLA");
    if (DropFlag & DROP_TBLA) bankList(&bcs_, (char*)"E-", (char*) "TBLA");
    if (DropFlag & DROP_HBTB) bankList(&bcs_, (char*)"E-", (char*) "HBTB");
    if (DropFlag & DROP_SC) bankList(&bcs_, (char*)"E-", (char*)SC_BANKS);
    if (DropFlag & DROP_EC) bankList(&bcs_, (char*)"E-", (char*)EC_BANKS);
    if (DropFlag & DROP_HBID) bankList(&bcs_, (char*)"E-", (char*) "HBID");
    if (DropFlag & DROP_CL01) bankList(&bcs_, (char*)"E-", (char*) "CL01");
    if (DropFlag & DROP_SEB) bankList(&bcs_, (char*)"E-", (char*)SEB_BANKS);
    if (DropFlag & DROP_TBID) bankList(&bcs_, (char*)"E-", (char*) "TBIDPARTTBERTBTR");
    if (DropFlag & DROP_HDPL) bankList(&bcs_, (char*)"E-", (char*) "HDPL");
    if (DropFlag & DROP_LAC) bankList(&bcs_, (char*)"E-", (char*)"EC1R");
    if (DropFlag & DROP_CC) bankList(&bcs_, (char*)"E-", (char*)CC_BANKS);
    if (DropFlag & DROP_ST) bankList(&bcs_, (char*)"E-", (char*)ST_BANKS);
    if (DropFlag & DROP_DHCL) bankList(&bcs_, (char*)"E-", (char*) "DHCL");
    if (DropFlag & DROP_TAGR) bankList(&bcs_, (char*)"E-", (char*)TAGGER_BANKS);
    
    
}

/* Fault handlers */

void installFaultHandlers()
{
    signal(SIGINT,signalINT);  
    signal(SIGSEGV,signalSEGV);
    //  signal(SIGABRT,signalSEGV);
    //  signal(SIGBUS,signalSEGV); 
    //  signal(SIGPIPE,signalPIPE);
    cerr << "Fault handlers installed" << endl;
}

static void signalINT(int isig){
    clasEvent e;
    char mess[100];
    static int count = 0;
    
    cerr << "Run/Event: " << CurrentRun << " " << CurrentEvent << endl;
    
    count++;
    if (count < 4){
        signal(isig,signalINT);
    } 
    else {
        exit(0);
    }
}
static void signalSEGV(int isig){
    clasEvent e;
    static int icount = 0;
    
    cerr << "signalSEGV: Caught signal " << endl;
    cerr << "Run/Event: " << CurrentRun << " " << CurrentEvent << endl;
    exit(1);
}


void printView(Float_t *array){
    cout << "PP ";
    for(int i=0; i<nData; i++){
        cout << array[i] << " ";
    }
    cout << endl;
}

fourVec pcorFourVec_new( fourVec particle, int part_id ){
    float proton_mass = 0.938272;
    float pion_mass = 0.139570;
    float kaon_mass = 0.493667;
    //Assign mass based on pid.
    double PID_mass =-1000;
    PID_mass = (part_id == 14                )?proton_mass:PID_mass;
    PID_mass = (part_id == 11|| part_id == 12)?kaon_mass  :PID_mass;
    PID_mass = (part_id == 8 || part_id == 9 )?pion_mass  :PID_mass;
    
       float id = part_id; // proton (geant 3 ID codes)
        float p =  particle.r(); // total momentum
        float px =  particle.x();
        float py =  particle.y();

        float pz = sqrt(p*p - px*px - py*py);
        float phi = atan2(py,px);
        float theta = acos(pz/p);

        /// Momentum correction ////////////////////////////////////////
        float new_p = p + pcor.pcor(phi,id);
        /// ////////////////////////////////////////////////////////////

        float new_px = (new_p / p) * px;
        float new_py = (new_p / p) * py;
        float new_pz = (new_p / p) * pz;



    
    double pE = sqrt( PID_mass*PID_mass + new_p*new_p );
    fourVec correctedFourVec;
    correctedFourVec.set(pE, new_px, new_py, new_pz);
    return correctedFourVec; 
}
