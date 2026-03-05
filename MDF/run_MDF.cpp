#include"libs.hh"
#include"definitions.hh"
#include"R3BGladFieldMap.h"
#include"R3BTPropagator.h"
#include"R3BTrackingParticle.h"
#include"R3BMDFWrapper.h"
#include"FairIon.h"
using namespace std;

Double_t BEAM_MASS_GEV_C2;//in GeV, initialized in the main function from FairIon

R3BTrackingParticle Generate_Primary()
{
    const Double_t  Ekin_Total = BEAM_EKIN_GEV_U * BEAM_MASS_GEV_C2/AMU;
    Double_t P, Theta, RCos, Phi, vx, vy, vz, Angle, Rand_R, beta;
    TVector3 PVector;
    P     = sqrt(Ekin_Total*(Ekin_Total+2*BEAM_MASS_GEV_C2)); //initial momentum
    auto Ptot  = P;
    P       = (Ptot+gRandom->Uniform((-1)*Ptot*SCALE_MIN, Ptot*SCALE_MAX));//min and max momentum spread
    while(1)
    {
        RCos    = gRandom->Uniform(TMath::Cos(BEAM_DTHETA_DEG*TMath::Pi()/180.),1.); 
        Theta   = TMath::ACos(RCos); 
        Phi     = gRandom->Uniform(0., 2.*TMath::Pi());
        const double Angle_r   = gRandom->Uniform(0., 2.*TMath::Pi());
        Rand_R  = gRandom->Uniform(0, BEAM_RADIUS_CM);
        vx	    = Rand_R * cos(Angle_r);
        vy	    = Rand_R * sin(Angle_r);
        vz      = gRandom->Uniform(0);//somehow R3BTPropagator fails at z=0
        //vz      = gRandom->Uniform(-BEAM_THICK_CM/2.,BEAM_THICK_CM/2.);//somehow R3BTPropagator fails at z=0
        PVector.SetMagThetaPhi(P, Theta, Phi);
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

//----------- Obtain track points as measured in the experiment ------------
bool ExtractData(Data &d, R3BTPropagator* fPropagator)
{
    double x,x1,x2,z,z1,z2,y,y1,y2,det_to_TP;//coordinates of the detector plane 
    TVector3 position = d.primary.GetStartPosition();
    TVector3 momentum = d.primary.GetStartMomentum();
    TVector3 v1, v2, v3;//points on the the plane for propagation
    R3BTrackingParticle * particle;
    TVector3 points[5];
    TVector3 mom_track[5];
    Double_t length[5];
    Double_t tof[5];
    int i =0;
    particle = new R3BTrackingParticle(BEAM_Z,
            position.X(),
            position.Y(),
            position.Z(),
            momentum.X(),
            momentum.Y(),
            momentum.Z(),
            d.primary.GetBeta(),
            d.primary.GetMass());

    //Assuming 5 detectors where the trajectories are measured
    for(i=0; i<4; i++)
    {
        //Caluclate x, and z to define the detector plane
        switch(i)
        {
            case 0 ://FOOT5, first after the target
                z1 = Target_to_FOOT5 + gRandom->Uniform(-2.,2.);
                x1 = 0;
                z2 = z1;
                x2 = Xmax;
                break;

            case 1 ://First Fi32 (X) 
                det_to_TP = TP_to_F32 + gRandom->Uniform(-6,6);
                z1 = Target_to_TP + det_to_TP * TMath::Cos(Angle);
                x1 = 0. - (det_to_TP * TMath::Sin(Angle));
                z2 = z1 + 10 * TMath::Cos(Angle_det);
                x2 = x1 + 10 * TMath::Sin(Angle_det);
                break;

            case 2 ://Second Fi30 (Y)
                det_to_TP = TP_to_F30 + gRandom->Uniform(-6,6);
                z1 = Target_to_TP + det_to_TP * TMath::Cos(Angle);
                x1 = 0. - (det_to_TP * TMath::Sin(Angle));
                //2nd point is 10 cm from the central point along det surface
                z2 = z1 + 10 * TMath::Cos(Angle_det);
                x2 = x1 + 10 * TMath::Sin(Angle_det);
                break;

            case 3 ://Last two fibers
                det_to_TP = (TP_to_F31 + TP_to_F33)/2. +  gRandom->Uniform(-20,20);
                z1 = Target_to_TP + det_to_TP * TMath::Cos(Angle);
                x1 = 0. - (det_to_TP * TMath::Sin(Angle));
                //2nd point is 10 cm from the central point along det surface
                z2 = z1 + 10 * TMath::Cos(Angle_det);
                x2 = x1 + 10 * TMath::Sin(Angle_det);
                break;
        }

        double tof_value = 0;
        v1.SetXYZ(x1,0,z1); //central det point 
        v2.SetXYZ(x1,Ymax,z1);//vertical det point above central
        v3.SetXYZ(x2,0,z2);//3rd point to define the plane

        //if(i == 2){
        if(!fPropagator->PropagateToPlaneRK_eloss(particle,v1,v2,v3,tof_value,0.1,false)) 
        {
            cout << "\n\t-- WARNING! --> Failed propagation...";
            return false;
        }
        /*}
          else{
          if(!fPropagator->PropagateToPlaneRK_eloss(particle,v1,v2,v3,tof_value,0.05,false)) 
          {
          cout << "\n\t-- WARNING! --> Failed propagation...";
          return false;
          }
          }*/

        points[i] = particle->GetPosition();
        length[i] = particle->GetLength();
        mom_track[i] = particle->GetMomentum();
        tof[i] = tof_value;
        if(isnan(points[i].X()) || isnan(points[i].Z()) || isnan(points[i].Y()))
        {
            cout << "\n\t-- WARNING! --> Propagated NaN values!";
            return false;
        }
    }
    i=0;

    //FOOT5 XYZ point
    d.edata[i]   =  points[0].X(); 
    d.edata[++i] =  points[0].Y();
    d.edata[++i] =  points[0].Z();
    d.edata[++i] = mom_track[0].X()/mom_track[0].Z();
    d.edata[++i] = mom_track[0].Y()/mom_track[0].Z();

    d.edata[++i] =  points[1].X(); 
    d.edata[++i] =  points[1].Z();
    d.edata[++i] = (points[3].X() - points[1].X())/(points[3].Z() - points[1].Z());
    d.edata[++i] = (points[2].Y() - points[0].Y())/(points[2].Z() - points[0].Z());

    Double_t FlightPath = length[1] - length[0];

    for(i=0; i<nVars; ++i)
    {
        if(TMath::IsNaN(d.edata[i]))
        {
            cout << "\nNaN values for tracking!!";
            return false;
        }
    }

    switch(FIT_VALUE)
    {
        case 0 ://PoQ
            d.value = momentum.Mag()/BEAM_Z; 
            break;
        case 1 ://TX0
            d.value = mom_track[0].X()/mom_track[0].Z();
            break;
        case 2 ://TY0
            d.value = mom_track[0].Y()/mom_track[0].Z();
            break;
        case 3 ://Flight Path
            d.value = FlightPath;
            break;
        case 4 ://TX1 after glad
            d.value = mom_track[1].X()/mom_track[1].Z();
            break;
        case 5 ://TY1 after glad
            d.value = mom_track[1].Y()/mom_track[1].Z();
            break;
        case 6 ://Y after glad
            d.value = points[2].Y();
            break;
        default :
            cout << "\n-- ERROR! Choose correct fit value!";
            return false;
    }
    delete particle;
    return true;
}

void Print_MDF_params(TMultiDimFit* mdf)
{
    cout << "\n\n--Real MDF parameters:" << endl;
    Int_t i, j;

    Int_t     gNVariables               = mdf->GetNVariables();
    Int_t     gNCoefficients            = mdf->GetNCoefficients();
    Int_t     gNMaxFunctions            = mdf->GetMaxFunctions();
    Double_t  gDMean                    = mdf->GetMeanQuantity();
    Int_t     gPolyType                 = mdf->GetPolyType();

    const TVectorD * gXMean             = mdf->GetMeanVariables();
    const TVectorD * gXMin              = mdf->GetMinVariables();
    const TVectorD * gXMax              = mdf->GetMaxVariables();
    const TVectorD * gCoefficient       = mdf->GetCoefficients();
    const TVectorD * gCoefficientRMS    = mdf->GetCoefficientsRMS();
    const Int_t    * gPower             = mdf->GetPowers();
    Int_t          * gPowerIndex        = mdf->GetPowerIndex();

    //outfile << fixed << setprecision (9);
    cout << setprecision(20) << gNVariables << " " << gNCoefficients<< " " << gNMaxFunctions 
        << " " << gDMean << " " << gPolyType << " ";
    for(i=0; i<gNVariables;    i++) cout << setprecision(20) << (*gXMean)(i) << " ";
    for(i=0; i<gNVariables;    i++) cout << setprecision(20) << (*gXMin)(i) << " ";
    for(i=0; i<gNVariables;    i++) cout << setprecision(20) << (*gXMax)(i) << " ";
    for(i=0; i<gNCoefficients; i++) cout << setprecision(20) << (*gCoefficient)(i) << " ";
    for(i=0; i<gNCoefficients; i++) cout << setprecision(20) << (*gCoefficientRMS)(i) << " ";
    for (i = 0; i < gNCoefficients; i++)
    {
        for (j = 0; j < gNVariables; j++)
        {
            cout << gPower[gPowerIndex[i] * gNVariables + j] << " ";
        }
    }
    return;
}

//Do not change! This matches to R3BMDFWrapper
void Save_MDF_params(TMultiDimFit* mdf, const char * filename)
{
    Int_t i, j;
    ofstream outfile;
    outfile.open(filename);

    Int_t     gNVariables               = mdf->GetNVariables();
    Int_t     gNCoefficients            = mdf->GetNCoefficients();
    Int_t     gNMaxFunctions            = mdf->GetMaxFunctions();
    Double_t  gDMean                    = mdf->GetMeanQuantity();
    Int_t     gPolyType                 = mdf->GetPolyType();

    const TVectorD * gXMean             = mdf->GetMeanVariables();
    const TVectorD * gXMin              = mdf->GetMinVariables();
    const TVectorD * gXMax              = mdf->GetMaxVariables();
    const TVectorD * gCoefficient       = mdf->GetCoefficients();
    const TVectorD * gCoefficientRMS    = mdf->GetCoefficientsRMS();
    const Int_t    * gPower             = mdf->GetPowers();
    Int_t          * gPowerIndex        = mdf->GetPowerIndex();

    //outfile << fixed << setprecision (9);
    outfile << setprecision(20) << gNVariables << " " << gNCoefficients<< " " << gNMaxFunctions 
        << " " << gDMean << " " << gPolyType << " ";
    for(i=0; i<gNVariables;    i++) outfile << setprecision(20) << (*gXMean)(i) << " ";
    for(i=0; i<gNVariables;    i++) outfile << setprecision(20) << (*gXMin)(i) << " ";
    for(i=0; i<gNVariables;    i++) outfile << setprecision(20) << (*gXMax)(i) << " ";
    for(i=0; i<gNCoefficients; i++) outfile << setprecision(20) << (*gCoefficient)(i) << " ";
    for(i=0; i<gNCoefficients; i++) outfile << setprecision(20) << (*gCoefficientRMS)(i) << " ";
    for (i = 0; i < gNCoefficients; i++)
    {
        for (j = 0; j < gNVariables; j++)
        {
            outfile << setprecision(20) << gPower[gPowerIndex[i] * gNVariables + j] << " ";
        }
    }
    outfile.close();
    return;

}

//Do not change! This matches to R3BMDFWrapper
void Save_PCA_params(TPrincipal* pca, const char * filename)
{
    Int_t i, j;
    ofstream outfile;
    outfile.open(filename);
    Int_t gNVariables               = nVars;
    const TMatrixD* gEigenVectors   = pca->GetEigenVectors();
    const TVectorD* gEigenValues    = pca->GetEigenValues();
    const TVectorD* gMeanValues     = pca->GetMeanValues();
    const TVectorD* gSigmas         = pca->GetSigmas();

    outfile << gNVariables << " ";
    for (i = 0; i < gNVariables; i++) {
        for (j = 0; j < gNVariables; j++) {
            outfile << (*gEigenVectors)(i,j) << " ";
        }
    }

    for (i = 0; i < gNVariables; i++) outfile << (*gEigenValues)(i) << " ";
    for (i = 0; i < gNVariables; i++) outfile << (*gMeanValues)(i) << " ";
    for (i = 0; i < gNVariables; i++) outfile << (*gSigmas)(i) << " ";

    outfile.close();
    return;
}

//Main function for MDF fit called from main()
void run_MDF(Bool_t use_PCA, R3BTPropagator * fPropagator)
{

    FairIon* tparticle = new FairIon(TString::Format("FairIon_%d_%d_%d", BEAM_Z, BEAM_A, BEAM_Z), BEAM_Z, BEAM_A, BEAM_Z);
    //initialize ion mass
    BEAM_MASS_GEV_C2 = tparticle->GetMass();//GeV/c2, used inside GenerateParticle

    TApplication* theApp = new TApplication("App", 0, 0);
    TFile* output = new TFile("output/MDF.root", "RECREATE");

    TH2F * h_beam_x_vs_y   = new TH2F("h_beam_x_vs_y","h_beam_x_vs_y",1000,-20,20,1000,-20,20);
    TH2F * h_beam_tx_vs_ty = new TH2F("h_beam_tx_vs_ty","h_beam_tx_vs_ty",1000,0,10,1000,0,10);

    TH2F * h_track_vs_fit;
    TH2F * h_residual_vs_track;
    TH1F * h_residual;
    TH1F * h_track;

    Int_t    hist_Nbins = 1000;
    Double_t hist_res_min;
    Double_t hist_res_max;
    Double_t hist_value_min;
    Double_t hist_value_max;

    switch(FIT_VALUE)
    {
        case 0 ://PoQ
            hist_res_min = -0.1; hist_res_max = 0.1; hist_value_min = 0; hist_value_max = 10;
            break;
        case 1 ://TX0
            hist_res_min = -0.01; hist_res_max = 0.01; hist_value_min = -0.08; hist_value_max = 0.08;
            break;
        case 2 ://TY0
            hist_res_min = -0.01; hist_res_max = 0.01; hist_value_min = -0.08; hist_value_max = 0.08;
            break;
        case 3 ://FlightPath
            hist_res_min = -0.5; hist_res_max = 0.5; hist_value_min = 600; hist_value_max = 750;
            break;
        case 4 ://TX1 after GLAD
            hist_res_min = -0.01; hist_res_max = 0.01; hist_value_min = (-30) * 3.14/180. ; hist_value_max = (-5)*3.14/180.;
            break;
        case 5 ://TY1 after GLAD
            hist_res_min = -0.01; hist_res_max = 0.01; hist_value_min = -0.08; hist_value_max = 0.08;
            break;
        case 6 ://Y after GLAD
            hist_res_min = -5; hist_res_max = 5; hist_value_min = -50; hist_value_max = 50;
            break;
        default :
            cout << "\n-- ERROR! Choose correct fit value!";
            return;
    }

    h_track_vs_fit   = new TH2F("h_track_vs_fit","h_track_vs_fit",
            hist_Nbins, hist_value_min, hist_value_max,
            hist_Nbins, hist_value_min, hist_value_max);

    h_residual       = new TH1F("h_residual","h_residual",
            hist_Nbins, hist_res_min, hist_res_max);

    h_track          = new TH1F("h_track","h_track",
            hist_Nbins, hist_value_min, hist_value_max);

    h_residual_vs_track = new TH2F("h_residual_vs_track","h_residual_vs_track",
            hist_Nbins, hist_value_min, hist_value_max,
            hist_Nbins, hist_res_min,   hist_res_max);

    TH2F * h_lab_xz = new TH2F("h_lab_xz","h_lab_xz",1000,0,800,1000,-500,500);

    cout << "\n\n\t*************************************************" << endl;
    cout << "\t*                                               *" << endl;
    cout << "\t*  Multidimensional Fit for tracking in R3BRoot *" << endl;
    cout << "\t*                                               *" << endl;
    cout << "\t*************************************************" << endl;
    if(use_PCA) cout << "\n\n-- Using PCA variables for tracking";
    else cout << "\n\n-- No PCA, using real variables for tracking";

    //------------- Calculating correct number of variables for MDF ------------/

    Int_t nv;
    if(use_PCA){ nv = nVars - nReduct;  } 
    else{
        if(nReduct!=0) cout << "\n-- nReduct is not 0! --> will be ignored if --pca flag is not set...";
        nv = nVars;
    }
    const Int_t nVars_MDF = nv;

    //------------- Set Multi Dimensional Fit Class for Tracking -----------/

    TMultiDimFit* MDFf;
    switch(Polynom_Type)
    {
        case 0 : MDFf=new TMultiDimFit(nVars_MDF, TMultiDimFit::kChebyshev,"v");
                 cout << "\n-- Using kChebyshev polynoms for MDF";
                 break;
        case 1 : MDFf=new TMultiDimFit(nVars_MDF, TMultiDimFit::kLegendre,"v"); 
                 cout << "\n-- Using kLegendre polynoms for MDF";
                 break;
        case 2 : MDFf=new TMultiDimFit(nVars_MDF, TMultiDimFit::kMonomials,"v"); 
                 cout << "\n-- Using kMonomials polynoms for MDF";
                 break;
        default:
                 MDFf=new TMultiDimFit(nVars_MDF, TMultiDimFit::kChebyshev,"v");
                 cout << "\n-- Using kChebyshev polynoms for MDF";
                 break;
    }
    Int_t mPowers[nVars_MDF];
    MDFf->SetMaxFunctions(MDFMaxFunctions);
    MDFf->SetMaxStudy(MDFMaxStudy);
    MDFf->SetMaxTerms(MDFMaxTerms);
    MDFf->SetPowerLimit(MDFPowerLimit);
    MDFf->SetMinAngle(MDFMinAngle);
    MDFf->SetMaxAngle(MDFMaxAngle);
    MDFf->SetMinRelativeError(MDFMinError);
    for(UInt_t i=0; i<nVars_MDF; i++) mPowers[i]=Power;
    MDFf->SetMaxPowers(mPowers);

    TPrincipal * PCA = new TPrincipal(nVars,"ND");//main PCA object
    UInt_t counter, ev, i, j;
    Double_t track_value;
    std::vector<Data> data_array;
    Data data;

    //---------------- Collecting training sample data --------------/

    cout << "\n-- Collecting track data.....\n";
    counter = 0;
    R3BTrackingParticle Primary;
    while(counter<NEVENTS)
    {
        cout << "\r-- Track # : " << counter << flush;
        data.primary = Generate_Primary(); 
        if(!ExtractData(data, fPropagator))
        {
            cout << "\n\t-- WARNING! --> Failed to extract track data, trying different primary... \n";
            continue;
        }
        if(use_PCA) PCA->AddRow(data.edata);
        data_array.push_back(data);
        counter++;
    }
    cout<<"\n-- Collected: " << data_array.size() << " tracks";
    if(use_PCA)
    {
        std::cout<<"\n\n#######  Principal Components:\n ";
        PCA->MakePrincipals();
        PCA->MakeHistograms();
        PCA->Print();
        PCA->MakeCode("output/PCA");
    }

    //-------------------- Filling up MDF data -------------------/

    cout << "\n-- Getting data for MDF\n";
    for(ev=0; ev<data_array.size(); ev++)
    {
        cout << "\r-- Track # :" << ev << std::flush;
        if(use_PCA)
        {
            PCA->X2P(data_array[ev].edata,data_array[ev].pdata);
            MDFf->AddRow(data_array[ev].pdata, data_array[ev].value, 1e-19);
        }
        else
        {
            MDFf->AddRow(data_array[ev].edata, data_array[ev].value, 1e-19);
        }
    }
    MDFf->MakeHistograms();
    cout << "\n-- Number of events in the training sample: " << ev; 
    cout << "\n-- Running MDF parameterization:\n";
    MDFf->FindParameterization();
    MDFf->Print("rm");

    //--------------- Collecting data and perfroming MDF fitting ---------------/

    cout << "\n-- Performing MDF fitting.... ";
    Double_t *xMax = new Double_t[MDFf->GetNVariables()];
    Double_t *xMin = new Double_t[MDFf->GetNVariables()];
    for (i = 0; i < MDFf->GetNVariables(); i++) {
        xMax[i] = (*MDFf->GetMaxVariables())(i);
        xMin[i] = (*MDFf->GetMinVariables())(i);
    }
    counter=0;
    //while(counter<NEVENTS)
    while(1)
    {
        cout << "\r-- Track # : " << counter << flush;
        data.primary = Generate_Primary(); 
        if(!ExtractData(data, fPropagator))
        {
            cout << "\n\t-- WARNING! --> Failed to extract track data, trying different primary... \n";
            continue;
        }
        if(use_PCA)
        {
            PCA->X2P(data.edata,data.pdata);//converting to PCA variables
            for (i = 0; i < MDFf->GetNVariables(); i++)
                if (data.pdata[i] < xMin[i] || data.pdata[i] > xMax[i]) break; // Check if the variables are in the range
            if (i != MDFf->GetNVariables()){ cout << "\n-- Variables out of range, skip this event...\n"; continue; }
            MDFf->AddTestRow(data.pdata, data.value, 1e-19);
        }
        else
        {
            for (i = 0; i < MDFf->GetNVariables(); i++)
                if (data.edata[i] < xMin[i] || data.edata[i] > xMax[i]) break;
            if (i != MDFf->GetNVariables()){ cout << "\n-- Variables out of range, skip this event...\n"; continue; }
            MDFf->AddTestRow(data.edata, data.value, 1e-19);
        }

        counter++;
        if(counter>(MDFf->GetNCoefficients()*50))
        {
            cout << "\n-- Stopped after getting " << counter <<" tracks";
            break;
        }
    }

    //--------------- Performing  final MDF fitting ---------------/
    cout << "\n-- Performing final fit on " << counter << " collected samples" << "\n\n";
    MDFf->Fit("M");
    //MDFf->MakeCode("output/MDFcode");

    //======= Save MDF and PCA params to be read by R3Broot ========
    MDFf->Print("PSCRFKM v");
    //Print_MDF_params(MDFf);
    Save_MDF_params(MDFf, "output/MDF_params.txt");
    if(use_PCA) Save_PCA_params(PCA, "output/PCA_params.txt");

    R3BMDFWrapper * mdfw = new R3BMDFWrapper();
    mdfw->InitMDF("output/MDF_params.txt");
    if(use_PCA) mdfw->InitPCA("output/PCA_params.txt");
    mdfw->PrintMDF();

    //--------------- Performing MDF quality check ---------------/

    cout << "\n-- Performing MDF quality check after fitting.... " << endl;
    counter =0;
    while(counter<NEVENTS)
    {
        cout << "\r-- Track # : " << counter << flush;
        data.primary = Generate_Primary(); 

        if(!ExtractData(data, fPropagator))
        {
            cout << "\n\t-- WARNING! --> Failed to extract track data, trying different primary... \n";
            continue;
        }
        if(use_PCA)
        {
            mdfw->X2P(data.edata,data.pdata);
            track_value = mdfw->MDF(data.pdata);
        }
        else
        {
            track_value = mdfw->MDF(data.edata);
        }

        counter++;

        h_track_vs_fit->Fill(track_value,data.value);
        h_residual->Fill(track_value-data.value);
        h_residual_vs_track->Fill(data.value, data.value-track_value);
        h_track->Fill(data.value);
        //Fill MDF input data
    }

    cout << "\n-- FINISH --" << endl;
    TCanvas * c1 = new TCanvas("c1","c1",1000,1000);
    c1->Divide(2,2);
    c1->cd(1);
    h_track_vs_fit->Draw("colz");
    c1->cd(2);
    h_residual->Draw();
    h_residual->Fit("gaus");
    c1->cd(3);
    h_residual_vs_track->Draw("colz");
    c1->cd(4);
    h_track->Draw();
    c1->Write();
    //h_lab_xz->Write();

    output->Write();
    //output->Close();
    //delete output;
    theApp->Run();
    return;
}
