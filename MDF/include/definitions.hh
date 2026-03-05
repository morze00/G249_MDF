#ifndef _DEFINITIONS_HH_
#define _DEFINITIONS_HH_

#include"libs.hh"
#include"R3BGladFieldMap.h"
#include"R3BTPropagator.h"
#include"R3BTrackingParticle.h"
#include <TMath.h>

//------------- Parameters of the training sample tracks -----------
const UInt_t    NEVENTS             = 2000; //Number of tracks for the training sample
const Double_t  AMU                 = 0.9314940038; //GeV 90 from AME 2016
const Double_t  BEAM_Z              = 9.; //primary Z
const Double_t  BEAM_A              = 25.;//Primary A
const Double_t  BEAM_EKIN_GEV_U     = 0.63;//per unit in GeV
const Double_t  BEAM_DTHETA_DEG     = 0.040*180./3.14;//max spread of theta angle in deg (+-)
const Double_t  BEAM_RADIUS_CM      = 1;//radius of the beamspot in cm 
const Double_t  SCALE_MIN           = 0.5;//scaling for subtraction: p_min = Ptot-Ptot*scale_min
const Double_t  SCALE_MAX           = 0.7; // scaling for adition:   p_max = Ptot+Ptot*scale_max
const Double_t  GLAD_CURRENT        = -2443.0;//23F @ 630 MeV/u to 18 deg
//Double_t BEAM_MASS_GEV_C2 = 0;//in GeV, initialized in the main function from FairIon

//Nominal detector positions in G249
//Positions in [cm] w.r.t. Turning Point (TP) of GLAD along 18 deg line
const Double_t Angle        = 18.0 * TMath::DegToRad();//18 deg line
const Double_t Angle_det    = TMath::Pi()/2. - Angle;//angle of the detector

const Double_t Target_to_FOOT5 = 24.5;//first FOOT after the target

const Double_t Target_to_GLAD_flange = 120.31;//cm
const Double_t Target_to_TP = Target_to_GLAD_flange + 165.4;//cm, TP = Turning point
const Double_t Target_to_RP = Target_to_GLAD_flange +  54.05 - 0.55580;//cm

const Double_t TP_to_ExitFlange  = 253.01834;//cm, along 18 deg line
const Double_t TP_to_F32  = TP_to_ExitFlange + 175.;//cm, first X
const Double_t TP_to_F30  = TP_to_F32 + 27.5;//cm, first Y
const Double_t TP_to_F31  = TP_to_F30 + 82.8;//cm second X WX
const Double_t TP_to_F33  = TP_to_F30 + 68.;//cm second X Messel
const Double_t TP_to_TOFD  = 560;
const Double_t F32_to_TOFD  = 277.6;

const UInt_t  FIT_VALUE = 5; //0=PoQ, 1=TX0, 2=TY0, 3=FlightPath, 4=TX1 after GLAD, 5=TY1 after glad, 6=Y after glad

//------------- Set Multi Dimensional Fit Class for Tracking -----------

const UInt_t   nVars   	= 9;
const UInt_t   nReduct 	= 0;//numbers of vars to ignore in PCA 
const UInt_t   Polynom_Type     = 0; //0=kChebyshev,1 =kLegendre(fails for MDFWrapper, don't use!),2-kMonomials 
const UInt_t   Power   	        = 8;
const UInt_t   MDFMaxFunctions  = 20000;
const UInt_t   MDFMaxStudy      = 400000;
const UInt_t   MDFMaxTerms      = 50; //PoQ
const UInt_t   MDFPowerLimit    = 1;
const Double_t MDFMinAngle      = 0.1;
const Double_t MDFMaxAngle      = 0.1;
const Double_t MDFMinError      = 1e-19;

//------------- Plotting ranges for the trajectories in [cm] -------------

const  Double_t Ymin = -100;
const  Double_t Ymax =  100;
const  Double_t Xmin = -1050;
const  Double_t Xmax =  1050;
const  Double_t Zmin = -300;
const  Double_t Zmax =  1800;
const  UInt_t NbinsX = 1000;
const  UInt_t NbinsY = 1000;
const  UInt_t NbinsZ = 1000;

//-------------- Data container for a given track --------------

struct Data
{
    Double_t edata[nVars];//data before PCA
    Double_t pdata[nVars - nReduct];//output data from PCA
    Double_t value;//main fit value
    R3BTrackingParticle primary;//original primary particle
};

//------------ Declaration of all fucntions ------------ 

R3BTrackingParticle Generate_Primary();
bool ExtractData(Data &d, R3BTPropagator* fPropagator);
void Print_MDF_params(TMultiDimFit* mdf);
void Save_MDF_params(TMultiDimFit* mdf, const char * filename);
void Save_PCA_params(TPrincipal* pca, const char * filename);
void run_MDF(Bool_t use_PCA, R3BTPropagator * fPropagator);

#endif




