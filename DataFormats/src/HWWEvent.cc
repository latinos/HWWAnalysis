/*
 * HWWEvent.cc
 *
 *  Created on: Feb 18, 2011
 *      Author: ale
 */
#include "HWWAnalysis/DataFormats/interface/HWWEvent.h"

//______________________________________________________________________________
HWWEvent::HWWEvent() : NEles(0), NMus(0), PFNJets(0) {
	// TODO Auto-generated constructor stub

	Clear();
}

//______________________________________________________________________________
HWWEvent::~HWWEvent() {
	// TODO Auto-generated destructor stub

}

//______________________________________________________________________________
void HWWEvent::Clear( Option_t* option ){
	Run            = 0;
	Event          = 0;
	LumiSection    = 0;
    Weight         = 1;

	PrimVtxPosition = math::XYZPoint();
	PrimVtxx       = 0;
	PrimVtxy       = 0;
	PrimVtxz       = 0;
	NVrtx          = 0;
	NPileUp        = 0;

    TCMet          = math::XYZTLorentzVector(); //TLorentzVector();
    PFMet          = math::XYZTLorentzVector(); //TLorentzVector();
    ChargedMet     = math::XYZTLorentzVector(); //TLorentzVector();


    ReducedPFMomenta.clear();

	NSoftMus       = 0;
	NBTaggedJets   = 0;
	NEles          = 0;
	NMus           = 0;
	PFNJets        = 0;

	Els.clear();
	Mus.clear();
	PFJets.clear();
    JetBtags.clear();
    HltPaths.clear();

    BTaggers.CombSecVrtxMVA.clear();
    BTaggers.SimpleSecVrtxHighEff.clear();
    BTaggers.SimpleSecVrtxHighPur.clear();
    BTaggers.JetBProb.clear();
    BTaggers.JetProb.clear();
    BTaggers.TkCntHighEff.clear();
    BTaggers.TkCntHighPur.clear();
}

