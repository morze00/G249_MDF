#ifndef _TRACK_HH_
#define _TRACK_HH_

#include"libs.hh"
#include"definitions.hh"
#include"R3BMCTrack.h"
#include"R3BFibPoint.h"

class MDFTrack
{
    public:

        MDFTrack();
        void Clear();
        bool SetTrack(TClonesArray** FibArray, TClonesArray *MCTrack);
        bool GetData(Double_t *edata, bool exact);
        Double_t GetPoQ(){ return _TargetMom.Mag()/Qion; }
        //inline Double_t GetFitValue(){ return _TargetMom.X()/_TargetMom.Z(); }
        //inline Double_t GetFitValue(){ return _TargetPos.X(); }
        inline Double_t GetFitValue(){ return _TargetMom.Mag()/Qion; }

    private:

        TVector3 _TargetPos;// [cm, ns]
        TVector3 _TargetMom;// GeV

        TVector3 _Fib3aPos;
        TVector3 _Fib3aMom;

        TVector3 _Fib10Pos;
        TVector3 _Fib10Mom;

        TVector3 _Fib12Pos;
        TVector3 _Fib12Mom;

        bool isTarget;
        bool isFib3a;
        bool isFib10;
        bool isFib12;
};


//Track variable is container for an event
struct Data
{
    Double_t edata[nVars];
    Double_t pdata[nVars - nReduct];
    Double_t value;
};


void Print_MDF_params(TMultiDimFit* mdf);
void Save_MDF_params(TMultiDimFit* mdf, const char * outfile);
void Save_PCA_params(TPrincipal* pca, const char *outfile);
void run_MDF(Char_t*filename, Bool_t using_pca);

#endif
