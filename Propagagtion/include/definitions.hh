#ifndef _DEFINITIONS_HH_
#define _DEFINITIONS_HH_

#include"libs.hh"

const UInt_t    NEVENTS             = 50; //Number of particles to propagate
const Double_t  AMU                 = 0.9314940038; //GeV 90 from AME 2016
const Int_t     BEAM_Z              = 9; //primary Z
const Int_t     BEAM_A              = 25;//Primary A
const Double_t  BEAM_EKIN_GEV_U     = 0.63;//per unit in GeV
const Double_t  BEAM_DTHETA_DEG     = 0.040*180./3.14;//max spread of theta angle in deg (+-)
const Double_t  BEAM_RADIUS_CM      = 1;//radius of the beamspot in cm 
const Double_t  SCALE_MIN           = 0.5;//scaling for subtraction: p_min = Ptot-Ptot*scale_min
const Double_t  SCALE_MAX           = 0.7; // scaling for adition:   p_max = Ptot+Ptot*scale_max
//const Double_t  GLAD_CURRENT        = -2668.0;//25F @ 630 MeV/u to 18 deg
const Double_t  GLAD_CURRENT        = -2443.0;//23F @ 630 MeV/u to 18 deg


//Plotting ranges in cm
const  Double_t Ymin = -100;
const  Double_t Ymax =  100;
const  Double_t Xmin = -350;
const  Double_t Xmax =  350;
const  Double_t Zmin = -100;
const  Double_t Zmax =  1100;
//const  UInt_t NbinsX = 1050;
const  UInt_t NbinsX = 500;
//const  UInt_t NbinsY = 1000;
const  UInt_t NbinsY = 100;
const  UInt_t NbinsZ = 1000;
//const  UInt_t NbinsZ = 1050;

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

const Double_t TP_to_RPC  = 350;

const Double_t Angle = -18. * TMath::DegToRad();
const Double_t Angle_RPC = -37. * TMath::DegToRad();

#endif

