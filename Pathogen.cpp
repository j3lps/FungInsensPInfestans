#include "Pathogen.h"

CPathogen::CPathogen(){
	
	// When the objects are created, set these vectors to zero of the appropriate size
	Latent.assign(TOTALGENOTYPES, 0.0);
	Infectious.assign(TOTALGENOTYPES, 0.0);
	sporesEnteringOWPool.assign(TOTALGENOTYPES, 0.0);
}

CPathogen::CPathogen(const CPathogen& objectToCopy){
	Latent = objectToCopy.Latent;
	Infectious = objectToCopy.Infectious;
	sporesEnteringOWPool = objectToCopy.sporesEnteringOWPool;
}

void CPathogen::zero(){
	Latent.assign(TOTALGENOTYPES, 0.0);
	Infectious.assign(TOTALGENOTYPES, 0.0);
	sporesEnteringOWPool.assign(TOTALGENOTYPES, 0.0);
}

CPathogen CPathogen::multiply(double scalar){
	CPathogen tempPathogen;
	for(unsigned int nGeno = 0; nGeno != TOTALGENOTYPES; ++nGeno){
		tempPathogen.Latent[nGeno] = this->Latent[nGeno] * scalar;
		tempPathogen.Infectious[nGeno] = this->Infectious[nGeno] * scalar;
		tempPathogen.sporesEnteringOWPool[nGeno] = this->sporesEnteringOWPool[nGeno] * scalar;
	}

	return tempPathogen;
}

CPathogen CPathogen::add(const CPathogen& obj){
	CPathogen tempPathogen;
	for(unsigned int nGeno = 0; nGeno != TOTALGENOTYPES; ++nGeno){
		tempPathogen.Latent[nGeno] = this->Latent[nGeno] + obj.Latent[nGeno];
		tempPathogen.Infectious[nGeno] = this->Infectious[nGeno] + obj.Infectious[nGeno];
		tempPathogen.sporesEnteringOWPool[nGeno] = this->sporesEnteringOWPool[nGeno] + obj.sporesEnteringOWPool[nGeno];
	}

	return tempPathogen;

}

bool CPathogen::anyNegative(){
	bool anyNeg = false;
	for(unsigned int nGeno = 0; nGeno != TOTALGENOTYPES; ++nGeno){
		if(Latent[nGeno] < 0.0) anyNeg = true;
		if(Infectious[nGeno] < 0.0) anyNeg = true;
		if (sporesEnteringOWPool[nGeno] < 0.0) anyNeg = true;
	}
	return anyNeg;
}

unsigned int CPathogen::TOTALGENOTYPES = 0;
