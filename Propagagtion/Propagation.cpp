#include"libs.hh"
#include"definitions.hh"
#include"R3BTrackingParticle.h"
#include"R3BGladFieldMap.h"
#include"R3BTPropagator.h"
#include"R3BMCTrack.h"
#include"FairIon.h"
using namespace std;

R3BTrackingParticle Generate_Primary()
{
    FairIon* tparticle = new FairIon(TString::Format("FairIon_%d_%d_%d", BEAM_Z, BEAM_A, BEAM_Z), BEAM_Z, BEAM_A, BEAM_Z);
    double BEAM_MASS_GEV_C2 = tparticle->GetMass();//in GeV
    cout << "\nIon Mass: " << BEAM_MASS_GEV_C2/AMU << " AMU";

    const Double_t  Ekin_Total = BEAM_EKIN_GEV_U * BEAM_MASS_GEV_C2/AMU;
    Double_t P, Theta, RCos, Phi, vx, vy, vz, Rand_R, beta;
    TVector3 PVector;
    P     = sqrt(Ekin_Total*(Ekin_Total+2*BEAM_MASS_GEV_C2)); //initial momentum
    auto Ptot  = P;
    P       = (Ptot+gRandom->Uniform((-1)*Ptot*SCALE_MIN, Ptot*SCALE_MAX));//min and max momentum spread

    while(1)
    {
        RCos    = gRandom->Uniform(TMath::Cos(BEAM_DTHETA_DEG*TMath::Pi()/180.),1.); 
        Theta   = TMath::ACos(RCos); 
        Phi     = gRandom->Uniform(0., 2.*TMath::Pi());
        Rand_R  = gRandom->Uniform(0, BEAM_RADIUS_CM);
        vx	    = Rand_R * cos(Phi);
        vy	    = Rand_R * sin(Phi);
        vz = 0; //somehow R3BTPropagator fails at z=0
        PVector.SetMagThetaPhi(P, Theta, Phi);
        cout << "\nX: " << vx + (Target_to_RP-vz) * PVector.X()/PVector.Z();
        cout << "\nY: " << vy + (Target_to_RP-vz) * PVector.Y()/PVector.Z();
        //Do not use events which may potentially hit the wall of GLAD chamber
        if(fabs(vx + (Target_to_RP-vz) * PVector.X()/PVector.Z() )<13. && 
                fabs(vy + (Target_to_RP-vz) * PVector.Y()/PVector.Z() )<13.)
            break;
    }
    beta = 1. / TMath::Sqrt(1 + TMath::Power(BEAM_MASS_GEV_C2 / PVector.Mag(), 2));

    R3BTrackingParticle ion(BEAM_Z,
            vx,
            vy,
            vz,
            PVector.X(),
            PVector.Y(),
            PVector.Z(),
            beta,
            BEAM_MASS_GEV_C2);
    return ion;
}
R3BTrackingParticle Generate_Primary(int Z, int A, double AMeV_min, double AMeV_max, double dTheta)
{
    Double_t vx, vy, vz, Rand_R, beta;
    TVector3 position_vec;
    TVector3 momentum_vec;
    
    //Get mass of the tparticle
    FairIon* tparticle = new FairIon(TString::Format("FairIon_%d_%d_%d", Z, A, Z), Z, A, Z);
    double mass_tparticle = tparticle->GetMass();//in GeV
    const Double_t  Ekin_Total_min = A * AMeV_min * 0.001; // MeV to GeV
    const Double_t  Ekin_Total_max = A * AMeV_max * 0.001; // MeV to GeV
    
    //Get energy and 3-momentum of the tparticle 
    Double_t P_min = sqrt(Ekin_Total_min*(Ekin_Total_min+2*mass_tparticle)); //minimum momentum in GeV/c
    Double_t P_max = sqrt(Ekin_Total_max*(Ekin_Total_max+2*mass_tparticle)); //maximum momentum in GeV/c
    Double_t P = gRandom->Uniform(P_min,P_max);

    //Randomize angles of the particle
    Double_t RCos    = gRandom->Uniform(TMath::Cos(dTheta),1.); 
    Double_t Theta   = TMath::ACos(RCos); 
    Double_t Phi     = gRandom->Uniform(0., 2.*TMath::Pi());
    momentum_vec.SetMagThetaPhi(P, Theta, Phi);

    //Get start position of the tparticle
    while(1)
    {
        Phi     = gRandom->Uniform(0., 2.*TMath::Pi());
        Rand_R  = gRandom->Uniform(0, BEAM_RADIUS_CM);
        vx	    = Rand_R * cos(Phi);
        vy	    = Rand_R * sin(Phi);
        vz = 0;
        cout << "\nX: " << vx + (Target_to_RP-vz) * momentum_vec.X()/momentum_vec.Z();
        cout << "\nY: " << vy + (Target_to_RP-vz) * momentum_vec.Y()/momentum_vec.Z();
        cout << "\nMass: " << mass_tparticle << endl;
        //Do not use events which may potentially hit the wall of GLAD chamber
        if(fabs(vx + (Target_to_RP-vz) * momentum_vec.X()/momentum_vec.Z() )<13. && 
                fabs(vy + (Target_to_RP-vz) * momentum_vec.Y()/momentum_vec.Z() )<13.)
            break;
    }
    beta = momentum_vec.Mag()/sqrt(momentum_vec.Mag2() + pow(mass_tparticle,2));
    R3BTrackingParticle ion(Z,
            vx,
            vy,
            vz,
            momentum_vec.X(),
            momentum_vec.Y(),
            momentum_vec.Z(),
            beta,
            mass_tparticle);
    return ion;
}
R3BTrackingParticle Generate_Primary(int Z, int A, double AMeV, double sigma_Px, double sigma_Py, double sigma_Pz)
{
    Double_t P, Theta, RCos, Phi, vx, vy, vz, Rand_R, beta;
    TVector3 position_vec;
    TVector3 momentum_vec;

    //Get mass of the tparticle
    FairIon* tparticle = new FairIon(TString::Format("FairIon_%d_%d_%d", Z, A, Z), Z, A, Z);
    double mass_tparticle = tparticle->GetMass();//in GeV
    const Double_t  Ekin_Total = A * AMeV * 0.001; // MeV to GeV

    //Get energy and 3-momentum of the tparticle 
    P = sqrt(Ekin_Total*(Ekin_Total+2*mass_tparticle)); //initial momentum in GeV/c
    momentum_vec.SetXYZ(0,0,P);
    momentum_vec.SetX(gRandom->Gaus(momentum_vec.X(),sigma_Px));//GeV units
    momentum_vec.SetY(gRandom->Gaus(momentum_vec.Y(),sigma_Py));//GeV units
    momentum_vec.SetZ(gRandom->Gaus(momentum_vec.Z(),sigma_Pz));//GeV units

    //Get start position of the tparticle
    while(1)
    {
        Phi     = gRandom->Uniform(0., 2.*TMath::Pi());
        Rand_R  = gRandom->Uniform(0, BEAM_RADIUS_CM);
        vx	    = Rand_R * cos(Phi);
        vy	    = Rand_R * sin(Phi);
        vz = 0;
        cout << "\nX: " << vx + (Target_to_RP-vz) * momentum_vec.X()/momentum_vec.Z();
        cout << "\nY: " << vy + (Target_to_RP-vz) * momentum_vec.Y()/momentum_vec.Z();
        cout << "\nMass: " << mass_tparticle << endl;
        //Do not use events which may potentially hit the wall of GLAD chamber
        if(fabs(vx + (Target_to_RP-vz) * momentum_vec.X()/momentum_vec.Z() )<13. && 
                fabs(vy + (Target_to_RP-vz) * momentum_vec.Y()/momentum_vec.Z() )<13.)
            break;
    }
    beta = momentum_vec.Mag()/sqrt(momentum_vec.Mag2() + pow(mass_tparticle,2));
    R3BTrackingParticle ion(Z,
            vx,
            vy,
            vz,
            momentum_vec.X(),
            momentum_vec.Y(),
            momentum_vec.Z(),
            beta,
            mass_tparticle);
    return ion;
}
void Propagation()
{
    TApplication* theApp = new TApplication("App", 0, 0);
    gStyle->SetOptStat(0);
    cout << "\n\t*************************************************" << endl;
    cout << "\t*                                               *" << endl;
    cout << "\t*  Propagating particles through GLAD field     *" << endl;
    cout << "\t*                                               *" << endl;
    cout << "\t*************************************************" << endl;

    cout << "\n-- Initializing GLAD field" << endl;
    //R3BGladFieldMap* magField = new R3BGladFieldMap("R3BGladMap","A");
    //R3BGladFieldMap* magField = new R3BGladFieldMap("R3BGladMap_Bxyz_X-3to3_Y-1to1_Z-3to5_step10mm","R");
    R3BGladFieldMap* magField = new R3BGladFieldMap("R3BGladMap_Bxyz_X-3to3_Y-1to1_Z-4to13_step10mm","R");
    magField->SetPosition(0.7871, 1.75-1.526, Target_to_GLAD_flange + 54.05 - 0.55580);//x,y,z in cm 
    magField->SetXAngle(-0.113); //deg
    magField->SetYAngle(-14.08); //deg
    magField->SetZAngle(-0.83); //deg
    magField->SetScale(GLAD_CURRENT/3583.81);
    magField->Init();
    magField->Print();

    cout << "\n-- Initializing field propagator";
    R3BTPropagator * fPropagator = new R3BTPropagator(magField, kFALSE);

    cout << "\n-- Plotting field map";
    double x,y,z,B;
    TH2D* h2 = new TH2D("h2","h2", NbinsZ, Zmin, Zmax, NbinsX, Xmin, Xmax);
    for (auto i = 0; i < NbinsZ; ++i)
    {
        z = Zmin + (i + 0.5) * (Zmax-Zmin)/NbinsZ;
        y=0;
        for (auto j = 0; j < NbinsX; ++j)
        {
            x = Xmin + (j + 0.5) * (Xmax-Xmin)/NbinsX;
            B = 0.1 * sqrt(pow(magField->GetBx(x, y, z),2) + pow(magField->GetBy(x, y, z),2)+ pow(magField->GetBz(x, y, z),2));//0.1 for kG->Tesla
            h2->Fill(z, x, B);
        }
    }

    //=========== Drawing setup configuration ======================
    TCanvas * c1 = new TCanvas("c1","c1",1200,700);
    TPaletteAxis *palette;
    c1->cd(1);
    h2->Draw("colz Z");
    h2->GetZaxis()->SetRangeUser(1e-6, 3);
    c1->Update();
    c1->SetLogz();
    palette=(TPaletteAxis*)h2->FindObject("palette");
    palette->Draw();

    //Beam axis (Z-axis)
    TLine* line_beamline = new TLine(Zmin, 0., Zmax, 0.);
    line_beamline->SetLineWidth(1);
    line_beamline->SetLineColor(kBlue);
    line_beamline->SetLineStyle(9);
    line_beamline->Draw();


    TVector3 p0, p1, p2;

    //Drawing 14 deg line
    //TLine* line_18deg = new TLine(Target_to_TP, 0., Zmax, (Zmax - Target_to_TP) * TMath::Tan(-14.*TMath::DegToRad()));
    //line_18deg->SetLineWidth(2);
    //line_18deg->SetLineColor(kBlue);
    //line_18deg->Draw();
    //c1->Update();

    //Drawing entrance flange
    TLine* line_entrance_flange = new TLine(Target_to_GLAD_flange, -50., Target_to_GLAD_flange, 50.);
    line_entrance_flange->Draw();
    line_entrance_flange->SetLineWidth(1);

    //Drawing EXIT flange
    p0.SetXYZ(TP_to_ExitFlange * TMath::Sin(-14*TMath::DegToRad()), 0, Target_to_TP + TP_to_ExitFlange * TMath::Cos(-14*TMath::DegToRad()));
    p1.SetXYZ(-200., 0., 0.);
    p2.SetXYZ(200., 0., 0.);
    p1.RotateY(-14.*TMath::DegToRad());
    p2.RotateY(-14.*TMath::DegToRad());
    p1 += p0;
    p2 += p0;
    TLine* line_exit_flange = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_exit_flange->SetLineWidth(1);
    line_exit_flange->SetLineColor(kBlack);
    line_exit_flange->Draw();

    //Drawing Rotation point
    TLine* line_rp = new TLine(Target_to_RP, 0., Target_to_RP, 20.);
    line_rp->SetLineWidth(2);
    line_rp->SetLineColor(kBlue);
    line_rp->Draw();

    //Drawing Turning point 
    TLine* line_tp = new TLine(Target_to_TP, 0., Target_to_TP, 20.);
    line_tp->SetLineWidth(2);
    line_tp->SetLineColor(kBlue);
    line_tp->Draw();

    //Drawing F32 (first X coord)
    p0.SetXYZ(TP_to_F32 * TMath::Sin(Angle), 0, Target_to_TP + TP_to_F32 * TMath::Cos(Angle));
    p1.SetXYZ(-25., 0., 0.);
    p2.SetXYZ(25., 0., 0.);
    p1.RotateY(Angle);
    p2.RotateY(Angle);
    p1 += p0;
    p2 += p0;
    TLine* line_F32 = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_F32->SetLineWidth(2);
    line_F32->Draw();

    //Drawing F30 (first Y coord)
    p0.SetXYZ(TP_to_F30 * TMath::Sin(Angle), 0, Target_to_TP + TP_to_F30 * TMath::Cos(Angle));
    p1.SetXYZ(-25., 0., 0.);
    p2.SetXYZ(25., 0., 0.);
    p1.RotateY(Angle);
    p2.RotateY(Angle);
    p1 += p0;
    p2 += p0;
    TLine* line_F30 = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_F30->SetLineWidth(2);
    line_F30->Draw();

    //Drawing F31 (second X coord, WX)
    p0.SetXYZ(TP_to_F31 * TMath::Sin(Angle) - 25, 0, Target_to_TP + TP_to_F31 * TMath::Cos(Angle));
    p1.SetXYZ(-25., 0., 0.);
    p2.SetXYZ(25., 0., 0.);
    p1.RotateY(Angle);
    p2.RotateY(Angle);
    p1 += p0;
    p2 += p0;
    TLine* line_f31 = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_f31->SetLineWidth(2);
    line_f31->Draw();

    //Drawing F33 (second X coord, Messel)
    p0.SetXYZ(TP_to_F33 * TMath::Sin(Angle) + 25, 0, Target_to_TP + TP_to_F33 * TMath::Cos(Angle));
    p1.SetXYZ(-25., 0., 0.);
    p2.SetXYZ(25., 0., 0.);
    p1.RotateY(Angle);
    p2.RotateY(Angle);
    p1 += p0;
    p2 += p0;
    TLine* line_f33 = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_f33->SetLineWidth(2);
    line_f33->Draw();

    //Drawing TOFD
    p0.SetXYZ((TP_to_F32 + F32_to_TOFD) * TMath::Sin(Angle), 0, Target_to_TP + (TP_to_F32 +  F32_to_TOFD) * TMath::Cos(Angle));
    p1.SetXYZ(-60., 0., 0.);
    p2.SetXYZ(60., 0., 0.);
    p1.RotateY(Angle);
    p2.RotateY(Angle);
    p1 += p0;
    p2 += p0;
    TLine* line_tofd = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_tofd->SetLineWidth(2);
    line_tofd->Draw();


    //Drawing RPC
    p0.SetXYZ(TP_to_RPC * TMath::Sin(Angle_RPC) - 50, 0, Target_to_TP + TP_to_RPC * TMath::Cos(Angle));
    p1.SetXYZ(-75., 0., 0.);
    p2.SetXYZ(75., 0., 0.);
    p1.RotateY(Angle_RPC);
    p2.RotateY(Angle_RPC);
    p1 += p0;
    p2 += p0;
    TLine* line_RPC = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_RPC->SetLineWidth(2);
    line_RPC->Draw();
    c1->Update();

    cout << "\n-- Creating particle tracks" << endl;
    TVector3 momentum, position;
    TVector3 track_pos;
    TVector3 v1, v2, v3;//points on the the plane for propagation

    //Trajectory object
    double Xtraj[NbinsZ];
    double Ztraj[NbinsZ];
    TGraph * trajectory;
    R3BTrackingParticle * particle;
    R3BTrackingParticle primary;

    particle = new R3BTrackingParticle(0,0,0,0,0,0,0,0,0);

    for(int ev=0; ev<NEVENTS; ev++)
    {
        for(auto k=0; k<NbinsZ; k++)
        {
            Xtraj[k]=0;
            Ztraj[k]=0;
        }

        //cout << "\r-- Working on entry " << ev << flush;
        cout << "\n-- Working on entry " << ev;

        //primary = Generate_Primary(9,25,620,660,0.01);
        primary = Generate_Primary();

        position = primary.GetStartPosition();
        momentum = primary.GetStartMomentum();
        Double_t position_array[3];
        position_array[0] = position.X();
        position_array[1] = position.Y();
        position_array[2] = position.Z();

        particle->SetCharge(primary.GetCharge());
        particle->SetMass(primary.GetMass());
        particle->SetStartPosition(position);
        particle->SetStartMomentum(momentum);
        particle->SetStartBeta(primary.GetStartBeta());
        particle->SetMomentum(momentum);
        particle->SetPosition(position_array);
        particle->SetBeta(primary.GetBeta());

        for (auto i = 0; i < NbinsZ; i++)
        {
            z = (i + 1) * Zmax / NbinsZ;
            if(z<position.Z()) continue;

            v1.SetXYZ(0,0,z);
            v2.SetXYZ(1,0,z);
            v3.SetXYZ(0,1,z);

            double tof=0;
            //if(!fPropagator->PropagateToPlaneRK(particle,v1,v2,v3))
            if(!fPropagator->PropagateToPlaneRK_eloss(particle,v1,v2,v3,tof,0.1,false)){
                cout << "\nFailed propagation...";
                return;
            }

            track_pos = particle->GetPosition();
            momentum = particle->GetMomentum();
            position_array[0] = track_pos.X();
            position_array[1] = track_pos.Y();
            position_array[2] = track_pos.Z();

            double beta = particle->GetBeta();

            if(isnan(track_pos.X()) || isnan(track_pos.Z()))
            {
                cout << "\n\nERROR! Propagated NaN values! @ Z=" << z << "cm\n";
                cout << "\nStart position X, Y, Z: " << position.X() << "\t" << position.Y() << "\t" << position.Z()<< " cm";
                cout << "\nStart angles TX, TY: " << momentum.X()/momentum.Z() << "\t" << momentum.Y()/momentum.Z() << " \n";
                return;
            }

            Xtraj[i]=track_pos.X();
            Ztraj[i]=track_pos.Z();

            //Update the particle 
            particle->SetStartPosition(track_pos);
            particle->SetStartMomentum(momentum);
            particle->SetStartBeta(beta);
        }
        //delete particle;

        trajectory  = new TGraph(NbinsZ,Ztraj,Xtraj);
        trajectory->SetMarkerColor(kRed);
        trajectory->SetMarkerStyle(kFullCircle);
        trajectory->SetMarkerSize(0.1);
        trajectory->Draw("same P");
        c1->Update();
    }

    cout << "\n-- FINISH --" << endl;
    c1->Update();
    theApp->Run();
    return;
}
int main(Int_t argc, Char_t* argv[])
{
    gRandom = new TRandom3();
    gRandom->SetSeed(0);
    gStyle->SetPalette(kRainBow);
    Propagation();
    return 0;
}
