/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Basic histogram implimentation via AliJHistogramInterface.
// author: O. Saarimaki, D.J. Kim (dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla, Finland
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#ifndef ALIJCDIJETHISTOS_H
#define ALIJCDIJETHISTOS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliJHistogramInterface.h"
#include "AliJCard.h"


using namespace std;

class AliJJetJtHistos : public AliJHistogramInterface
{

  public:

    AliJJetJtHistos(AliJCard *card); //constructor
    AliJJetJtHistos(const AliJJetJtHistos& obj); // copy constructor
    virtual ~AliJJetJtHistos();    //destructor
    AliJJetJtHistos& operator=(const AliJJetJtHistos& obj); // equal sign operator

    void CreateEventTrackHistos();
    void SetCard( AliJCard * card ){fCard=card;}

    static int GetCentralityClass(Double_t);
    void SetCentralityBinsHistos( vector<double> centralityClasses ) {
      CentBin=centralityClasses;
      fNCentBin = CentBin.size();
    }
    static int fNCentBin;
    static vector<double> CentBin;//[NCENT+1]; //8
    double fJetEtaCut;
    AliJCard              * fCard; ///< Pointer to the configuration card
    TVector  *fJetTriggPtBorders; ///< Jet pT bin borders
    TVector  *fJetConstPtLowLimits; ///<  Comment needed
    TVector  *fJetAssocPtBorders; ///< Jet constituent pT bin borders
    TVector *fJetLeadPtBorders; ///< Leading track pT bin borders
    TVector *fJetMultBorders; ///< Jet multiplicity bin borders
    TVector  *fDeltaRBorders;
    TVector *fCentralityBorders; ///< Jet multiplicity bin borders
    TVector *fXlongBorders; ///< Xlong bin borders

    AliJHistManager * fHMG; //! Histogram manager
    AliJBin fHistCentBin;   //! Centrality bin
    AliJBin fJetBin;        //! Jet bin
    //===================================================
    // Event/Track histograms
    //===================================================
    AliJTH1D fh_events;     //! // for counting events, jets, dijets and so on.
    AliJTH1D fh_info;       //! // General information about the run.
    AliJTH1D fh_centrality; //! // centrality histogram
    AliJTH1D fh_zvtx;       //! // z-vertex histogram

    AliJTH1D fh_pt;     //! // for pt dist of tracks
    AliJTH1D fh_eta;    //! // for eta dist of tracks
    AliJTH2D fh_etaPhi; //! // for (eta,phi) dist of tracks
    AliJTH1D fh_phi;    //! // for phi dist of tracks

    AliJTH1D fh_rho;        //! // for event pt density
    AliJTH1D fh_rhom;       //! // for event mt density

    AliJTH1D fh_jetPt;      //! // for pt dist of jets
    AliJTH1D fh_jetEta;     //! // for eta dist of jets
    AliJTH1D fh_jetPhi;     //! // for phi dist of jets
    AliJTH2D fh_jetEtaPhi;  //! // for (eta,phi) dist of jets
    AliJTH1D fh_jetArea;    //! // for jet area spectrum
    AliJTH1D fh_jetAreaRho; //! // for jet area*pt spectrum

    AliJTH1D fh_dijetInvM;                //! // for dijet invariant mass
    AliJTH1D fh_dijetPtPair;              //! // for dijet pt
    AliJTH1D fh_dijetDeltaPhi;            //! // for dijet deltaPhi
    AliJTH1D fh_dijetPtPairDeltaPhiCut;   //! // for dijet pt after deltaPhi cut
    AliJTH1D fh_dijetInvMDeltaPhiCut;     //! // for dijet invariant mass after deltaPhi cut

    AliJBin fJetFinderBin;
    AliJBin fJetTriggerBin;
    AliJBin fTrkPtBin;
    AliJBin fTrkLimPtBin;
    AliJBin fJetLeadPtBin;
    AliJBin fJetMultBin;
    AliJBin fdRBin;
    AliJBin fiHist;
    AliJBin fCentralityBin;
    AliJBin fXlongBin;
    AliJBin fktFinderBin;
    AliJBin fDeltaPhiCutBin;
    AliJBin fJetFinderBinMC;
    AliJBin fJetTriggerBinMC;
    AliJBin fTrkPtBinMC;
    AliJBin fTrkLimPtBinMC;
    AliJBin fJetLeadPtBinMC;
    AliJBin fJetMultBinMC;
    AliJBin fdRBinMC;
    AliJBin fiHistMC;
    AliJBin fCentralityBinMC;
    AliJBin fBgTypeBin;
    AliJBin fiHist2MC;

    //===================================================
    // Jt Histograms
    // ==================================================
    AliJTH1D fhBgTrkNumber; // ! // Comment
    AliJTH1D fhBgTrkNumberBin; // ! // Comment
    AliJTH1D fhJetPtBin; // ! // Comment
    AliJTH1D fhEventJtWithPtCutWeightBinBin; // ! // Comment
    AliJTH1D fhEventJtWeightBin; // ! // Comment
    AliJTH1D fhEventJtBin; // ! // Comment
    AliJTH1D fhJetConeTrkPt; // ! // Comment
    AliJTH1D fhJetConeTrkPtBin; // ! // Comment
    AliJTH1D fhJetConeTrkPtWeightBin; // ! // Comment
    AliJTH1D fhJetConeZ; // ! // Comment
    AliJTH1D fhJetConeZBin; // ! // Comment
    AliJTH1D fhJetConeJt; // ! // Comment
    AliJTH1D fhJetConeJtBin; // ! // Comment
    AliJTH1D fhJetConeJtWeightBin; // ! // Comment
    AliJTH2D fhJetConeJtWeight2D; // ! // Comment
    AliJTH1D fhJetConeJtWeightWithTrackCutBinBin; // ! // Comment
    AliJTH1D fhJetConeJtWeightWithMultiplicityCutBinBin; // ! // Comment
    AliJTH1D fhJetConeJtLeadingRefBin; // ! // Comment
    AliJTH1D fhJetConeJtWeightLeadingRefBin; // ! // Comment
    AliJTH1D fhJetConeJtWeightLeadingRefWithTrackCutBinBin; // ! // Comment
    AliJTH1D fhJetConeJtWithPtCutWeightBinBin; // ! // Comment
    AliJTH1D fhBgTrkPt; // ! // Comment
    AliJTH1D fhBgTrkPtBin; // ! // Comment
    AliJTH1D fhBgTrkPtWeightBin; // ! // Comment
    AliJTH1D fhBgRBin; // ! // Comment
    AliJTH1D fhBgZ; // ! // Comment
    AliJTH1D fhBgZBin; // ! // Comment
    AliJTH1D fhBgJt; // ! // Comment
    AliJTH1D fhBgJtBin; // ! // Comment
    AliJTH1D fhBgJtWeightBin; // ! // Comment
    AliJTH1D fhBgJtWithPtCutWeightBinBin; // ! //
};

#endif //ALIJCDIJETHISTOS_H
