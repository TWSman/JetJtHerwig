/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Container class for histograms needed in the analysis.

#include "AliJJetJtHistos.h"
#include <TGrid.h>
#include "AliJCard.h"
#include <TPRegexp.h>

//Double_t AliJJetJtHistos::pttJacek[74+16] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500};
//UInt_t AliJJetJtHistos::NpttJacek = sizeof(AliJJetJtHistos::pttJacek)/sizeof(AliJJetJtHistos::pttJacek[0])-1;
vector<double> AliJJetJtHistos::CentBin;
int AliJJetJtHistos::fNCentBin;

//______________________________________________________________________________
AliJJetJtHistos::AliJJetJtHistos(AliJCard *card) :
	fHMG(NULL),
	fHistCentBin(),
	fJetBin(),
  fCard(card),
	fh_events(),
	fh_centrality(),
	fh_zvtx(),
	fh_pt(),
	fh_eta(),
	fh_phi(),
	fh_rho(),
	fh_rhom(),
	fh_jetPt(),
	fh_jetEta(),
	fh_jetPhi(),
	fh_jetEtaPhi(),
	fh_jetArea(),
	fh_jetAreaRho(),
  fh_dijetInvM(),
  fh_dijetPtPair(),
  fh_dijetDeltaPhi(),
  fh_dijetPtPairDeltaPhiCut(),
  fh_dijetInvMDeltaPhiCut(),
  fJetFinderBin(),
  fJetTriggerBin(),
  fTrkPtBin(),
  fTrkLimPtBin(),
  fJetLeadPtBin(),
  fJetMultBin(),
  fdRBin(),
  fiHist(),
  fCentralityBin(),
  fXlongBin(),
  fDeltaPhiCutBin(),
  fBgTypeBin(),
  fhBgTrkNumber(),
  fhBgTrkNumberBin(),
  fhJetPtBin(), //Jet pt as a fn of jet pt
  fhEventJtWithPtCutWeightBinBin(),
  fhEventJtWeightBin(),
  fhEventJtBin(),
  fhJetConeTrkPt(),
  fhJetConeTrkPtBin(),
  fhJetConeTrkPtWeightBin(),
  fhJetConeZ(),
  fhJetConeZBin(),
  fhJetConeJt(),
  fhJetConeJtBin(),
  fhJetConeJtWeightBin(),
  fhJetConeJtWeight2D(),
  fhJetConeJtWeightWithTrackCutBinBin(),
  fhJetConeJtWeightWithMultiplicityCutBinBin(),
  fhJetConeJtLeadingRefBin(),
  fhJetConeJtWeightLeadingRefBin(),
  fhJetConeJtWeightLeadingRefWithTrackCutBinBin(),
  fhJetConeJtWithPtCutWeightBinBin(),
  fhBgTrkPt(),
  fhBgTrkPtBin(),
  fhBgTrkPtWeightBin(),
  fhBgRBin(),
  fhBgZ(),
  fhBgZBin(),
  fhBgJt(),
  fhBgJtBin(),
  fhBgJtWeightBin(),
  fhBgJtWithPtCutWeightBinBin()
{

}

//______________________________________________________________________________
AliJJetJtHistos::AliJJetJtHistos(const AliJJetJtHistos& obj) :
  fHMG(obj.fHMG),
  fHistCentBin(obj.fHistCentBin),
  fJetBin(obj.fJetBin),
  fCard(obj.fCard),
  fh_events(obj.fh_events),
  fh_centrality(obj.fh_centrality),
  fh_zvtx(obj.fh_zvtx),
  fh_pt(obj.fh_pt),
  fh_eta(obj.fh_eta),
  fh_phi(obj.fh_phi),
  fh_rho(obj.fh_rho),
  fh_rhom(obj.fh_rhom),
  fh_jetPt(obj.fh_jetPt),
  fh_jetEta(obj.fh_jetEta),
  fh_jetPhi(obj.fh_jetPhi),
  fh_jetEtaPhi(obj.fh_jetEtaPhi),
  fh_jetArea(obj.fh_jetArea),
  fh_jetAreaRho(obj.fh_jetAreaRho),
  fh_dijetInvM(obj.fh_dijetInvM),
  fh_dijetPtPair(obj.fh_dijetPtPair),
  fh_dijetDeltaPhi(obj.fh_dijetDeltaPhi),
  fh_dijetPtPairDeltaPhiCut(obj.fh_dijetPtPairDeltaPhiCut),
  fh_dijetInvMDeltaPhiCut(obj.fh_dijetInvMDeltaPhiCut),
  fJetFinderBin(obj.fJetFinderBin),
  fJetTriggerBin(obj.fJetTriggerBin),
  fTrkPtBin(obj.fTrkPtBin),
  fTrkLimPtBin(obj.fTrkLimPtBin),
  fJetLeadPtBin(obj.fJetLeadPtBin),
  fJetMultBin(obj.fJetMultBin),
  fdRBin(obj.fdRBin),
  fiHist(obj.fiHist),
  fCentralityBin(obj.fCentralityBin),
  fXlongBin(obj.fXlongBin),
  fDeltaPhiCutBin(obj.fDeltaPhiCutBin),
  fBgTypeBin(obj.fBgTypeBin),
  fhBgTrkNumber(obj.fhBgTrkNumber),
  fhBgTrkNumberBin(obj.fhBgTrkNumberBin),
  fhJetPtBin(obj.fhJetPtBin),
  fhEventJtWithPtCutWeightBinBin(obj.fhEventJtWithPtCutWeightBinBin),
  fhEventJtWeightBin(obj.fhEventJtWeightBin),
  fhEventJtBin(obj.fhEventJtBin),
  fhJetConeTrkPt(obj.fhJetConeTrkPt),
  fhJetConeTrkPtBin(obj.fhJetConeTrkPtBin),
  fhJetConeTrkPtWeightBin(obj.fhJetConeTrkPtWeightBin),
  fhJetConeZ(obj.fhJetConeZ),
  fhJetConeZBin(obj.fhJetConeZBin),
  fhJetConeJt(obj.fhJetConeJt),
  fhJetConeJtBin(obj.fhJetConeJtBin),
  fhJetConeJtWeightBin(obj.fhJetConeJtWeightBin),
  fhJetConeJtWeight2D(obj.fhJetConeJtWeight2D),
  fhJetConeJtWeightWithTrackCutBinBin(obj.fhJetConeJtWeightWithTrackCutBinBin),
  fhJetConeJtWeightWithMultiplicityCutBinBin(obj.fhJetConeJtWeightWithMultiplicityCutBinBin),
  fhJetConeJtLeadingRefBin(obj.fhJetConeJtLeadingRefBin),
  fhJetConeJtWeightLeadingRefBin(obj.fhJetConeJtWeightLeadingRefBin),
  fhJetConeJtWeightLeadingRefWithTrackCutBinBin(obj.fhJetConeJtWeightLeadingRefWithTrackCutBinBin),
  fhJetConeJtWithPtCutWeightBinBin(obj.fhJetConeJtWithPtCutWeightBinBin),
  fhBgTrkPt(obj.fhBgTrkPt),
  fhBgTrkPtBin(obj.fhBgTrkPtBin),
  fhBgTrkPtWeightBin(obj.fhBgTrkPtWeightBin),
  fhBgRBin(obj.fhBgRBin),
  fhBgZ(obj.fhBgZ),
  fhBgZBin(obj.fhBgZBin),
  fhBgJt(obj.fhBgJt),
  fhBgJtBin(obj.fhBgJtBin),
  fhBgJtWeightBin(obj.fhBgJtWeightBin),
  fhBgJtWithPtCutWeightBinBin(obj.fhBgJtWithPtCutWeightBinBin)
{
  // copy constructor
}

//______________________________________________________________________________
AliJJetJtHistos& AliJJetJtHistos::operator=(const AliJJetJtHistos& obj){
  // copy constructor
  return *this;
}

//______________________________________________________________________________
AliJJetJtHistos::~AliJJetJtHistos() {
  // destructor
  delete fHMG;
}


//______________________________________________________________________________
void AliJJetJtHistos::CreateEventTrackHistos(){
  fJetTriggPtBorders = fCard->GetVector("JetTriggPtBorders");
  fJetTriggPtBorders->Print();
  fJetConstPtLowLimits = fCard->GetVector("JetConstPtLowLimits");
  fJetAssocPtBorders = fCard->GetVector("JetAssocPtBorders");
  fJetLeadPtBorders = fCard->GetVector("JetLeadPtBorders");
  fJetMultBorders = fCard->GetVector("JetMultBorders");
  fDeltaRBorders = fCard->GetVector("DeltaRBorders");
  fJetEtaCut = fCard-> Get("JetEtaCut");
  fXlongBorders = fCard->GetVector("xEBorders");
  // Create basic event histograms
  fHMG = new AliJHistManager("AliJCDijetHistManager","jcdijet");
  // set AliJBin here //
  fHistCentBin.Set("CentBin","CentBin","Cent:",AliJBin::kSingle).SetBin(fNCentBin);
  fJetBin.Set("JetBin","JetBin","Jet bin:",AliJBin::kSingle).SetBin(5);
  fJetFinderBin.Set("JetFinderBin","JetFinderBin","Jet Finder Bin:",AliJBin::kSingle).SetBin(3);

  fJetTriggerBin.Set("JetTriggerBin","JetPt","p_{T,jet} : %.1f - %.1f").SetBin(fCard->GetVector("JetTriggPtBorders"));
  fTrkPtBin .Set("TrkPtBin","TrkPt","p_{T,constituent}:%.1f-%.1f").SetBin(fCard->GetVector("JetAssocPtBorders"));
  fTrkLimPtBin .Set("TrkLimitPtBin","TrkLimitPt","p_{T,Limit}<%.1f", AliJBin::kSingle).SetBin(fJetConstPtLowLimits->GetNoElements());
  fJetLeadPtBin.Set("LeadPtBin","LeadPt","p_{T,leading}:%.1f-%.1f").SetBin(fCard->GetVector("JetLeadPtBorders"));
  fJetMultBin.Set("MultBin","Multiplicity","N_{const.}:%.1f-%.1f").SetBin(fCard->GetVector("JetMultBorders"));
  fdRBin.Set("dRBin","dR","dR : %.1f - %.1f ").SetBin(fCard->GetVector("DeltaRBorders"));
  fiHist.Set("iHist","iHist","iHist : %d ", AliJBin::kSingle).SetBin(10);
  fCentralityBin.Set("EventCentariltyBin","NCent","Centrality: %d - %d ").SetBin(fCard->GetVector("CentBinBorders"));
  fDeltaPhiCutBin.Set("DeltaPhiCutBin","DPhi","DPhi:%d", AliJBin::kSingle).SetBin(2); // 0 for no cut, 1 for deltaPhi cut
  fXlongBin.Set("XlongBin","Xlong","Xlong: %.1f - %.1f ").SetBin(fCard->GetVector("xEBorders"));

  // fh_events counts several things:
  // 0:  Number of events
  // 1:  Number of ch. particles
  // 2:  Number of accepted ch. particles
  // 3:  Number of events with no rho calculations
  // 4:  Number of events with proper rho calculations
  // 5:  Number of jets
  // 6:  Number of accepted jets
  // 7:  Number of accepted jets after const. cut
  // 8:  Number of accepted bg subtracted jets
  // 9:  Number of accepted bg subtracted const. cut jets
  // 10: Number of kt-jets
  // 11: Number of accepted kt-jets
  // 12: Number of jets that drop under leading pt cut after bg subtraction
  // 13: Number of jets that drop under subleading pt cut after bg subtraction
  // 14: Number of raw dijets
  // 15: Number of raw dijets after leading pt cut
  // 16: Number of accepted raw dijets
  // 17: Number of accepted raw dijets with delta phi cut
  // 14: Number of bg subtr. dijets
  // 15: Number of bg subtr. dijets after leading pt cut
  // 16: Number of accepted bg subtr. dijets
  // 17: Number of accepted bg subtr. dijets with delta phi cut
  // 18: Number of bg subtr. const. cut dijets
  // 19: Number of bg subtr. const. cut dijets after leading pt cut
  // 20: Number of accepted bg subtr. const. cut dijets
  // 21: Number of accepted bg subtr. const. cut dijets with delta phi cut
  // 22: Number of const. cut dijets
  // 23: Number of const. cut dijets after leading pt cut
  // 24: Number of accepted const. cut dijets
  // 25: Number of accepted const. cut dijets with delta phi cut
  // 26: Number of kt-dijets
  // 27: Number of kt-dijets after leading pt cut
  // 28: Number of accepted kt-dijets
  // 29: Number of accepted kt-dijets with delta phi cut
  fh_events
    << TH1D("h_events", "h_events", 40, 0.0, 40.0 )
    << fHistCentBin
    << "END" ;

  fh_info
    << TH1D("h_info", "h_info", 40, 0.0, 40.0 )
    << "END" ;

  fh_centrality
    << TH1D("h_centrality", "h_centrality", 100, 0.0, 100.0 )
    << "END" ;

  fh_zvtx
    << TH1D("h_zvtx", "h_zvtx", 40, -20.0, 20.0 )
    << "END" ;

  int NBINSJet=150;
  double LogBinsXJet[NBINSJet+1], LimLJet=0.1, LimHJet=500;
  double logBWJet = (log(LimHJet)-log(LimLJet))/NBINSJet;
  for(int ijetBin=0;ijetBin<=NBINSJet;ijetBin++) LogBinsXJet[ijetBin]=LimLJet*exp(ijetBin*logBWJet);

  // ============= CHARGED PARTICLE HISTOS ============= 
  fh_pt
    //<< TH1D("h_pt", "h_pt", AliJJetJtHistos::NpttJacek, AliJJetJtHistos::pttJacek)
    << TH1D("h_pt","h_pt",NBINSJet, LogBinsXJet )
    << fHistCentBin
    << "END" ;

  fh_eta
    << TH1D("h_eta", "h_eta", 100, -1.0, 1.0 )
    << fHistCentBin
    << "END" ;

  fh_phi
    << TH1D("h_phi", "h_phi", 100, -TMath::Pi(), TMath::Pi())
    << fHistCentBin
    << "END" ;

  fh_etaPhi
    << TH2D("h_etaPhi", "h_etaPhi", 100, -1.0, 1.0, 100, -TMath::Pi(), TMath::Pi())
    << fHistCentBin
    << "END" ;

  // ============= JET HISTOS ============= 
  fh_rho
    << TH1D("h_rho", "h_rho", NBINSJet, LogBinsXJet)
    << fHistCentBin
    << "END" ;

  fh_rhom
    << TH1D("h_rhom", "h_rhom", NBINSJet, LogBinsXJet)
    << fHistCentBin
    << "END" ;

  fh_jetPt
    //<< TH1D("h_jetPt", "h_jetPt", AliJJetJtHistos::NpttJacek, AliJJetJtHistos::pttJacek)
    << TH1D("h_jetPt","h_jetPt",NBINSJet, LogBinsXJet )
    << fHistCentBin << fJetBin
    << "END" ;

  fh_jetEta
    << TH1D("h_jetEta", "h_jetEta", 100, -1.0, 1.0)
    << fHistCentBin << fJetBin
    << "END" ;

  fh_jetPhi
    << TH1D("h_jetPhi", "h_jetPhi", 100, -TMath::Pi(), TMath::Pi())
    << fHistCentBin << fJetBin
    << "END" ;

  fh_jetEtaPhi
    << TH2D("h_jetEtaPhi", "h_jetEtaPhi", 100, -1.0, 1.0, 100, -TMath::Pi(), TMath::Pi())
    << fHistCentBin << fJetBin
    << "END" ;

  fh_jetArea
    //<< TH1D("h_jetArea", "h_jetArea", AliJJetJtHistos::NpttJacek, AliJJetJtHistos::pttJacek)
    << TH1D("h_jetArea", "h_jetArea", NBINSJet, LogBinsXJet )
    << fHistCentBin << fJetBin
    << "END" ;

  fh_jetAreaRho
    //<< TH1D("h_jetAreaRho", "h_jetAreaRho", AliJJetJtHistos::NpttJacek, AliJJetJtHistos::pttJacek)
    << TH1D("h_jetAreaRho", "h_jetAreaRho", NBINSJet, LogBinsXJet )
    << fHistCentBin << fJetBin
    << "END" ;

  int NBINSDijet=170;
  double logBinsXDijet[NBINSDijet+1], LimLDijet=0.1, LimHDijet=1000;
  double logBWDijet = (log(LimHDijet)-log(LimLDijet))/NBINSDijet;
  for(int iDijet=0;iDijet<=NBINSDijet;iDijet++) logBinsXDijet[iDijet]=LimLDijet*exp(iDijet*logBWDijet);

  // ============= DIJET HISTOS ============= 
  fh_dijetInvM
    << TH1D("h_dijetInvM", "h_dijetInvM", NBINSDijet, logBinsXDijet)
    << fHistCentBin << fJetBin
    << "END" ;

  fh_dijetPtPair
    //<< TH1D("h_dijetPtPair", "h_dijetPtPair", AliJJetJtHistos::NpttJacek, AliJJetJtHistos::pttJacek )
    << TH1D("h_dijetPtPair", "h_dijetPtPair", NBINSDijet, logBinsXDijet )
    << fHistCentBin << fJetBin
    << "END" ;

  fh_dijetDeltaPhi
    << TH1D("h_dijetDeltaPhi", "h_dijetDeltaPhi", 100, 0, 10)
    << fHistCentBin << fJetBin
    << "END" ;

  fh_dijetPtPairDeltaPhiCut
    //<< TH1D("h_dijetPtPair", "h_dijetPtPair", AliJJetJtHistos::NpttJacek, AliJJetJtHistos::pttJacek )
    << TH1D("h_dijetPtPairDeltaPhiCut", "h_dijetPtPairDeltaPhiCut", NBINSDijet, logBinsXDijet )
    << fHistCentBin << fJetBin
    << "END" ;

  fh_dijetInvMDeltaPhiCut
    << TH1D("h_dijetInvMDeltaPhiCut", "h_dijetInvMDeltaPhiCut", NBINSDijet, logBinsXDijet)
    << fHistCentBin << fJetBin
    << "END" ;

  int NBINSJt=64;
  double LogBinsJt[NBINSJt+1], LimLJt=0.01, LimHJt=10;
  double logBWJt = (TMath::Log(LimHJt)-TMath::Log(LimLJt))/NBINSJt;
  for(int ij=0;ij<=NBINSJt;ij++) LogBinsJt[ij]=LimLJt*exp(ij*logBWJt);
  int NBINSJtW=64;
  double LimLJtW=TMath::Log(0.01), LimHJtW=TMath::Log(10);

  int NBINS=150;
  double LogBinsX[NBINS+1], LimL=0.1, LimH=500;
  double logBW = (log(LimH)-log(LimL))/NBINS;
  for(int ij=0;ij<=NBINS;ij++) LogBinsX[ij]=LimL*exp(ij*logBW);

  int NBINSNumber = 100;
  int LBinsNumber = 0;
  int HBinsNumber = 100;
  fhBgTrkNumber
    << TH1D("BgTrkNumber","",NBINSNumber,LBinsNumber,HBinsNumber)
    << fJetFinderBin
    <<"END";
  fhBgTrkNumberBin
    << TH1D("BgTrkNumberBin","",NBINSNumber,LBinsNumber,HBinsNumber)
    << fJetFinderBin << fJetTriggerBin
    <<"END";

  fhJetPtBin
    << TH1D("JetPtBin","",NBINS, LogBinsX )
    << fJetFinderBin << fJetTriggerBin
    <<"END";

  fhEventJtBin
    << TH1D("EventJtBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhEventJtWeightBin
    << TH1D("EventJtWeightBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin
    <<"END";

  fhEventJtWithPtCutWeightBinBin
    << TH1D("EventJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin << fTrkPtBin
    <<"END";

  int NBINSPt=64;
  double LogBinsPt[NBINSPt+1], LimLPt=0.01, LimHPt=50;
  double logBWPt = (TMath::Log(LimHPt)-TMath::Log(LimLPt))/NBINSPt;
  for(int ij=0;ij<=NBINSPt;ij++) LogBinsPt[ij]=LimLPt*exp(ij*logBWPt);

  fhJetConeTrkPt
    << TH1D("JetConeTrkPt","",NBINSPt,LogBinsPt)
    << fJetFinderBin
    <<"END";

  fhJetConeTrkPtBin
    << TH1D("JetConeTrkPtBin","",NBINSPt,LogBinsPt)
    << fJetFinderBin << fJetTriggerBin
    <<"END";

  fhJetConeTrkPtWeightBin
    << TH1D("JetConeTrkPtWeightBin","",NBINSPt,LogBinsPt)
    << fJetFinderBin << fJetTriggerBin
    <<"END";

  int NBINSZ=64;
  double LogBinsZ[NBINSZ+1], LimLZ=0.001, LimHZ=1.1;
  double logBWZ = (TMath::Log(LimHZ)-TMath::Log(LimLZ))/NBINSZ;
  for(int ij=0;ij<=NBINSZ;ij++) LogBinsZ[ij]=LimLZ*exp(ij*logBWZ);//

  fhJetConeZ
    << TH1D("JetConeZ","",NBINSZ, LogBinsZ )
    << fJetFinderBin
    <<"END";
  fhJetConeZBin
    << TH1D("JetConeZBin","",NBINSZ, LogBinsZ )
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhJetConeJt
    << TH1D("JetConeJt","",NBINSJt, LogBinsJt )
    << fJetFinderBin
    <<"END";
  fhJetConeJtBin
    << TH1D("JetConeJtBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhJetConeJtWeightBin
    << TH1D("JetConeJtWeightBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhJetConeJtWeightWithTrackCutBinBin
    << TH1D("JetConeJtWeightWithTrackCutBinBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin << fJetLeadPtBin
    <<"END";
  fhJetConeJtWeightWithMultiplicityCutBinBin
    << TH1D("JetConeJtWeightWithMultiplicityCutBinBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin << fJetMultBin
    <<"END";
  int NBINSC = fJetTriggPtBorders->GetNoElements();
  double BinsC[NBINSC];
  for(int i = 0 ; i < NBINSC ; i++){
    BinsC[i] = (*fJetTriggPtBorders)[i+1];
  }
  fhJetConeJtWeight2D
    << TH2D("JetConeJtWeight2D","",NBINSJt, LogBinsJt,NBINSC-1,BinsC )
    << fJetFinderBin
    <<"END";
  fhJetConeJtLeadingRefBin
    << TH1D("JetConeJtLeadingRefBin","",NBINSJt,LogBinsJt)
    << fJetFinderBin << fJetTriggerBin << fXlongBin
    <<"END";
  fhJetConeJtWeightLeadingRefBin
    << TH1D("JetConeJtWeightLeadingRefBin","",NBINSJt,LogBinsJt)
    << fJetFinderBin << fJetTriggerBin << fXlongBin
    <<"END";
  fhJetConeJtWeightLeadingRefWithTrackCutBinBin
    << TH1D("JetConeJtWeightLeadingRefWithTrackCutBinBin","",NBINSJt,LogBinsJt)
    << fJetFinderBin << fJetTriggerBin << fJetLeadPtBin << fXlongBin
    << "END";
  fhJetConeJtWithPtCutWeightBinBin
    << TH1D("JetConeJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin << fTrkPtBin
    <<"END";
  fhBgTrkPt
    << TH1D("BgTrkPt","",NBINSPt,LogBinsPt)
    << fJetFinderBin
    <<"END";
  fhBgTrkPtBin
    << TH1D("BgTrkPtBin","",NBINSPt,LogBinsPt)
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhBgTrkPtWeightBin
    << TH1D("BgTrkPtWeightBin","",NBINSPt,LogBinsPt)
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  int NBINSR = 100;
  double LimLR = 0;
  double LimHR = 2;
  fhBgRBin
    << TH1D("BgRBin","",NBINSR,LimLR,LimHR)
    << fJetFinderBin << fJetTriggerBin
    << "END";
  fhBgZ
    << TH1D("BgZ","",NBINSZ, LogBinsZ )
    << fJetFinderBin
    <<"END";
  fhBgZBin
    << TH1D("BgZBin","",NBINSZ, LogBinsZ )
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhBgJt
    << TH1D("BgJt","",NBINSJt, LogBinsJt )
    << fJetFinderBin
    <<"END";
  fhBgJtBin
    << TH1D("BgJtBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin
    <<"END";
  fhBgJtWeightBin
    << TH1D("BgJtWeightBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin
    <<"END";

  fhBgJtWithPtCutWeightBinBin
    << TH1D("BgJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt )
    << fJetFinderBin << fJetTriggerBin << fTrkPtBin
    <<"END";
}

int AliJJetJtHistos::GetCentralityClass(Double_t fCent){
  for(int iCbin = 0; iCbin < fNCentBin; iCbin++){
    if(fCent > CentBin[iCbin] && fCent < CentBin[iCbin+1])
      return iCbin;
  }
  return -1;
}

