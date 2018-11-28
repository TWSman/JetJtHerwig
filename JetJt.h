// -*- C++ -*-
#ifndef JtAnalysis_JetJt_H
#define JtAnalysis_JetJt_H
//
// This is the declaration of the JetJt class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <TClonesArray.h>
#include "TFile.h"
#include "AliJJetJtAnalysis.h"

class AliJJetJtHistos;

namespace JtAnalysis {

using namespace ThePEG;

/**
 * Here is the documentation of the JetJt class.
 *
 * @see \ref JetJtInterfaces "The interfaces"
 * defined for JetJt.
 */
class JetJt: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  JetJt();

  /**
   * The destructor.
   */
  virtual ~JetJt();
  //@}
  Double_t getDiffR(double phi1, double phi2, double eta1, double eta2);
  int GetBin(TVector *array, double val);
  void CalculateJetsJt(TClonesArray *inList,
      int    lDebug,
      int    lCBin,
      double lParticleEtaCut,
      double lParticlePtCut,
      double lJetCone,
      double lktJetCone,
      int    lktScheme,
      bool   lusePionMassInkt,
      bool   luseDeltaPhiBGSubtr,
      double lConstituentCut,
      double lDeltaPhiCut,
      int iContainer,
      int doLeadingRef);


  public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);

  /**
   * Return a LorentzTransform which would put the event in the
   * desired Lorentz frame.
   * @param event a pointer to the Event to be considered.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tcEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   * @param weight the weight of the current event.
   */
  virtual void analyze(const tPVector & particles, double weight);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   * @param weight the weight of the current event.
   */
  virtual void analyze(tPPtr particle, double weight);
  //@}

  protected:

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();


  public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

  TClonesArray*         GetJTracks()  {return &fJTracks;}

  protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


  // If needed, insert declarations of virtual function defined in the
  // InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


  private:
  TFile *fout;                            //! Output file
  TH1F *histo;
  TH1F *hJetPt;
  TClonesArray              fJTracks;     //! tracks array
  TClonesArray              fJJets;     //! tracks array
  AliJJetJtHistos *fhistos;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  JetJt & operator=(const JetJt &);

};

}

#endif /* JtAnalysis_JetJt_H */
