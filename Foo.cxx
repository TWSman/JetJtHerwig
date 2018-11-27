// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Foo class.
//

#include "Foo.h"
#include "TFile.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "AliJBaseTrack.h"
#include "AliJJetJtHistos.h"
#include "AliJCard.h"
#include "AliJJet.h"



using namespace JtAnalysis;

Foo::Foo():
  fJTracks("AliJBaseTrack",1000),
  fJJets("AliJJet",100)
  {
  // Default constructor

}

Foo::~Foo() {}



#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

void Foo::analyze(tEventPtr event, long ieve, int loop, int state) {
  vector<fastjet::PseudoJet> chparticles;
  vector<fastjet::PseudoJet> jets;
  fastjet::PseudoJet fTrack;
  fJTracks.Clear();
  fJJets.Clear();
  int iContainer = 0;
  int itrack = 0;
  int doLeadingRef = 0;
  double z; double jt; double jtleading; double zleading;
  double phi, eta, pt;
  double lJetCone = 0.4;
  double lParticlePtCut = 0.15;
  double lParticleEtaCut = 0.8;
  double etaMaxCutForJet = 0.25;
  TLorentzVector  vOrtho;
  TLorentzVector summedJet;
  int moveJet = 1;
  int iBin = 0;
  int leadingTrackIndex;
  //AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  if ( loop > 0 || state != 0 || !event ) return;
  /** create local variable to store the multiplicity */
  int mult(0);
  /** get the final-state particles */
  tPVector particles=event->getFinalState();
  /** loop over all particles */
  for (tPVector::const_iterator pit = particles.begin();
      pit != particles.end(); ++pit){
    /** Select only the charged particles */
    if( ChargedSelector::Check(**pit) ){
      Particle *part = *pit;
      Lorentz5Momentum ph = part->momentum();
      pt = ph.perp()/GeV;
      phi = ph.phi();
      double px  = ph.x()/GeV;
      double py = ph.y()/GeV;
      double pz = ph.z()/GeV;
      double E = ph.t()/GeV;
      //cout << "Eta: " << part->eta() << " Phi: " << phi << " pT: " << pt << " px: " << px << " py: " << py << " pz " << pz << " E: " << E << endl;
      //cout << *part << endl;
      //AliJBaseTrack *particle = new AliJBaseTrack(px,py,pz,E,itrack,0,1);
      new (fJTracks[itrack]) AliJBaseTrack(px,py,pz,E, itrack,0,1);
      AliJBaseTrack * particle = static_cast<AliJBaseTrack*>(fJTracks[itrack]);
      particle->SetLabel(part->number());
      //fJTracks[itrack] = particle;
      fTrack = fastjet::PseudoJet(px,py,pz,E);
      fTrack.set_user_index(itrack);
      chparticles.push_back(fTrack);
      ++mult;
      itrack++;
    }
  }


  if(chparticles.size()==0) return; // We are not intereted in empty events.
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4, fastjet::pt_scheme); //Other option: fastjet::E_scheme

  fastjet::ClusterSequence cs(chparticles, jet_def);
  jets = sorted_by_pt(cs.inclusive_jets(5));
  //cout << "Clustering with " << jet_def.description() << endl;
  //cout <<   "        pt y phi" << endl;
  for (unsigned i = 0; i < jets.size(); i++) {
    //cout << "jet " << i << ": "<< jets[i].pt() << "  "  << jets[i].rap() << " " << jets[i].phi() << endl;
    AliJJet *jet = new AliJJet(jets[i].px(),jets[i].py(), jets[i].pz(), jets[i].E(),i,0,0);
    vector<fastjet::PseudoJet> constituents = jets[i].constituents();
    for (unsigned j = 0; j < constituents.size(); j++) {
      int itrack = jets[i].constituents()[j].user_index();
      jet->AddConstituent(fJTracks[itrack]);
    }

    double conPtMax =0;
    eta = jet->Eta();
    // anti-kt-jet eta cut
    int doBkg = 1;
    int counter = 0;
    vOrtho.SetVect(jet->Vect());
    vOrtho.SetE(jet->E());
    vOrtho.SetPhi(jet->Phi()+TMath::Pi()/2);
    if(TMath::Abs(eta) < etaMaxCutForJet) {
      pt = jet->Pt();
      iBin = GetBin(fhistos->fJetTriggPtBorders,pt); // fill jetPt histos
      phi = jet->Phi();
      for (int icon = 0; icon<jet->GetConstituents()->GetEntries(); icon++){
        AliJBaseTrack *con = jet->GetConstituent(icon);
        if (con->Pt()>conPtMax) conPtMax = con->Pt();
        leadingTrackIndex = icon;
      }
      AliJBaseTrack *leadingTrack = jet->GetConstituent(leadingTrackIndex);
      double leadingTrackPt = jet->LeadingParticlePt(); //FIXME? For MC tracks this is possibly a track with no charge
      int jetMult = jet->GetConstituents()->GetEntries();
      int iBin2 = GetBin(fhistos->fJetLeadPtBorders,leadingTrackPt);
      int iBin3 = GetBin(fhistos->fJetMultBorders,jetMult);
      fhistos->fhJetPtBin[iContainer][iBin]->Fill(pt);

      double maxconpt = 0;

      int noTracks = fJTracks.GetEntries();
      //cout << "number of tracks: " << noTracks << endl;
      for (int itrack = 0; itrack < noTracks; ++itrack) {//loop over all the particles in the event
        //cout << "itrack: " << itrack << endl;
        AliJBaseTrack *track = (AliJBaseTrack*)fJTracks[itrack];
        double pta = track->Pt();
        eta = track->Eta();
        if (!track){
          cout << "No track " << endl;
          continue;
        }
        /*if(track->GetCharge() == 0){
          cout << "No Charge" << endl;
          continue;
        }*/
        if (pta<lParticlePtCut || TMath::Abs(eta) > lParticleEtaCut) continue;
        phi = track->Phi();
        if (pta > maxconpt) maxconpt = pta;
        int iptaBin = GetBin(fhistos->fJetAssocPtBorders, pta);
        if( iptaBin < 0 ) continue;
        z = (track->Vect()*jet->Vect().Unit())/jet->P();
        jt = (track->Vect()-z*jet->Vect()).Mag();
        //jT for all tracks in the event
        double deltaR   = getDiffR(jet->Phi(),track->Phi(),jet->Eta(),track->Eta());
        //cout << "jet Phi: " << jet->Phi() << "Eta: " << jet->Eta() << " Track Phi: " << track->Phi() << " Eta: " << track->Eta() << endl;
        //cout << "deltaR: " << deltaR << endl;
        if(deltaR < TMath::Pi()/2){
          fhistos->fhEventJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]->Fill(jt, 1.0/jt);
          fhistos->fhEventJtWeightBin[iContainer][iBin]->Fill(jt, 1.0/jt);
          fhistos->fhEventJtBin[iContainer][iBin]->Fill(jt,1);
        }
        //Jet Cone Jt here
        if ( deltaR < lJetCone){
          fhistos->fhJetConeTrkPt[iContainer]->Fill(pta,1);
          fhistos->fhJetConeTrkPtBin[iContainer][iBin]->Fill(pta,1);
          fhistos->fhJetConeTrkPtWeightBin[iContainer][iBin]->Fill(pta,1/pta);
          fhistos->fhJetConeZ[iContainer]->Fill( z , 1.0);
          fhistos->fhJetConeZBin[iContainer][iBin]->Fill( z , 1);
          fhistos->fhJetConeJt[iContainer]->Fill( jt , 1);
          fhistos->fhJetConeJtBin[iContainer][iBin]->Fill( jt , 1);
          fhistos->fhJetConeJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * 1);
          fhistos->fhJetConeJtWeight2D[iContainer]->Fill(jt,jet->Pt(),1.0/jt * 1);
          if(iBin2 > -1){
            fhistos->fhJetConeJtWeightWithTrackCutBinBin[iContainer][iBin][iBin2]->Fill( jt, 1.0/jt * 1);
          }
          if(iBin3 > -1){
            fhistos->fhJetConeJtWeightWithMultiplicityCutBinBin[iContainer][iBin][iBin3]->Fill( jt, 1.0/jt * 1);
          }
          if(pta < 0.99*leadingTrackPt && doLeadingRef){
            int xlongBin = GetBin(fhistos->fXlongBorders, pta/leadingTrackPt);
            if( xlongBin < 0 ) {
              continue;
            }
            zleading = (track->Vect()*leadingTrack->Vect().Unit())/leadingTrack->P();
            jtleading =  (track->Vect()-zleading*leadingTrack->Vect()).Mag();
            fhistos->fhJetConeJtLeadingRefBin[iContainer][iBin][xlongBin]->Fill(jtleading,1);
            fhistos->fhJetConeJtWeightLeadingRefBin[iContainer][iBin][xlongBin]->Fill(jtleading,1.0/jtleading);
            if(iBin2 > -1){
              fhistos->fhJetConeJtWeightLeadingRefWithTrackCutBinBin[iContainer][iBin][iBin2][xlongBin]->Fill(jtleading,1.0/jtleading);
            }
          }

          if (iptaBin < 0) continue;
          fhistos->fhJetConeJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
            ->Fill( jt, 1.0/jt);
        }

        if ( doBkg ){
          //Background jt
          deltaR   = getDiffR(vOrtho.Phi(),track->Phi(),vOrtho.Eta(),track->Eta());
          if ( deltaR < lJetCone){
            counter++;
            fhistos->fhBgTrkPt[iContainer]->Fill(pta,1);
            fhistos->fhBgTrkPtBin[iContainer][iBin]->Fill(pta,1);
            fhistos->fhBgTrkPtWeightBin[iContainer][iBin]->Fill(pta,1.0/pta);
            if(moveJet){
              summedJet = track->GetLorentzVector() + vOrtho;
              fhistos->fhBgRBin[iContainer][iBin]->Fill(getDiffR(summedJet.Phi(), track->Phi(), summedJet.Eta(), track->Eta()),1);
              z = (track->Vect()*summedJet.Vect().Unit())/summedJet.P();
              jt = (track->Vect()-z*summedJet.Vect()).Mag();
            }else{
              fhistos->fhBgRBin[iContainer][iBin]->Fill(deltaR,1);
              z = (track->Vect()*vOrtho.Vect().Unit())/vOrtho.P();
              jt = (track->Vect()-z*vOrtho.Vect()).Mag();
            }
            fhistos->fhBgZ[iContainer]->Fill( z , 1);
            fhistos->fhBgZBin[iContainer][iBin]->Fill( z , 1);
            fhistos->fhBgJt[iContainer]->Fill( jt , 1);
            fhistos->fhBgJtBin[iContainer][iBin]->Fill( jt , 1);
            fhistos->fhBgJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt);

            if (iptaBin < 0) continue;
            fhistos->fhBgJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
              ->Fill( jt, 1.0/jt);
          } //End of If
        } //End of background
      } //End of track loop
      if(doBkg){
        fhistos->fhBgTrkNumber[iContainer]->Fill(counter);
        fhistos->fhBgTrkNumberBin[iContainer][iBin]->Fill(counter);
      }
    }
  }//end of the anti-kt-jet loop
}

LorentzRotation Foo::transform(tcEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void Foo::analyze(const tPVector & particles, double weight) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void Foo::analyze(tPPtr, double weight) {}

void Foo::dofinish() {
  fout->Write();
  fout->Close();
  //AnalysisHandler::dofinish();
  // *** ATTENTION *** Normalize and post-process histograms here.
}

void Foo::doinitrun() {
  fout = new TFile("output.root","RECREATE");
  fout->cd();//opening of the output file
  char* cardName = "cardAlice_pp.input";
  AliJCard *fCard = new AliJCard(cardName);
  fhistos = new AliJJetJtHistos(fCard);
  fhistos->CreateEventTrackHistos();
  //AnalysisHandler::doinitrun();
  // *** ATTENTION *** histogramFactory().registerClient(this); // Initialize histograms.
  // *** ATTENTION *** histogramFactory().mkdirs("/SomeDir"); // Put histograms in specal directory.
}


IBPtr Foo::clone() const {
  return new_ptr(*this);
}

IBPtr Foo::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).



// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeNoPIOClass<Foo,AnalysisHandler>
  describeJtAnalysisFoo("JtAnalysis::Foo", "Foo.so");

void Foo::Init() {
  static ClassDocumentation<Foo> documentation
    ("There is no documentation for the Foo class");

}


double Foo::getDiffR(double phi1, double phi2, double eta1, double eta2){
  double diffPhi = TMath::Abs(phi1-phi2);
  if(diffPhi > TMath::Pi()){
    diffPhi = 2*TMath::Pi() - diffPhi;
  }
  return TMath::Sqrt(TMath::Power(diffPhi,2)+TMath::Power(eta1-eta2,2));
}

int Foo::GetBin(TVector *array, double val){
  int iBin=-1;
  for(int i=1; i< array->GetNoElements(); i++){
    if((*array)[i] <= val && val<(*array)[i+1]){
      iBin=i-1;
      break;
    }
  }
  return iBin;
}
