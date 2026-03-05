#include"libs.hh"
#include"definitions.hh"
#include"R3BMDFWrapper.h"
#include"R3BTrackingParticle.h"
#include"R3BGladFieldMap.h"
#include"R3BTPropagator.h"

using namespace std;

int main(Int_t argc, Char_t* argv[])
{
    TStopwatch t1; 
    t1.Start(); 
    Bool_t is_PCA = kFALSE;
    Bool_t NeedHelp = kTRUE;

    if (argc > 1){
        for (Int_t i = 0; i < argc; i++)
        {
            if (strncmp(argv[i],"--pca",5) == 0)
            {
                is_PCA = kTRUE;
                //NeedHelp = kFALSE;
            }
        }
    }
    //if (NeedHelp)
    //{
    //    cout << "\nUse the following flags: ";
    //    cout << "\n--pca  :  (optional) to use PCA variables. Default tracking is without PCA \n" << endl;
    //    return 0;
    //}
    
    gRandom = new TRandom3();
    gRandom->SetSeed(0);
    gROOT->Macro("~/rootlogon.C");
    gStyle->SetPalette(kRainBow);

    //R3BGladFieldMap* magField = new R3BGladFieldMap("R3BGladMap","A");
    //R3BGladFieldMap* magField = new R3BGladFieldMap("R3BGladMap_Bxyz_X-3to3_Y-1to1_Z-4to13_step5mm","R");
    R3BGladFieldMap* magField = new R3BGladFieldMap("R3BGladMap_Bxyz_X-3to3_Y-1to1_Z-4to13_step10mm","R");
    //R3BGladFieldMap* magField = new R3BGladFieldMap("test_field","R");
    //R3BGladFieldMap* magField = new R3BGladFieldMap("R3BGladMap_Bxyz_X-3to3_Y-1to1_Z-3to5_step10mm","R");
    magField->SetPosition(0.7871, 1.75-1.526, Target_to_GLAD_flange + 54.05 - 0.55580);//x,y,z in cm 
    magField->SetXAngle(-0.113); //deg
    magField->SetYAngle(-14.08); //deg
    magField->SetZAngle(-0.83); //deg
    magField->SetScale(GLAD_CURRENT/3583.81);
    magField->Init();

    cout << "\n-- Readout time of the field map file:\n";
    t1.Print();

    //magField->Print();
    R3BTPropagator * fPropagator = new R3BTPropagator(magField, kFALSE);

    run_MDF(is_PCA,fPropagator);

    return 0;
}
