#pragma once

#include <vector>

class CPathogen{
	
public:

	static unsigned int TOTALGENOTYPES;

	// Constructor - this is how an object is created if called as 'CPathogen object()'
	CPathogen();
	// This is a copy constructor, and how to copy one object into another
	CPathogen(const CPathogen &objectToCopy);

	// Latent density for each genotype
	std::vector<double> Latent;

	// Infectious density for each genotype
	std::vector<double> Infectious;

	// A vector storing the density of spores in the overwinter pool for the next season
	std::vector<double> sporesEnteringOWPool;

	// A function to zero all 3 vectors
	void zero();
	// A function to multiply all elements by a scalar
	CPathogen multiply(double);
	// Add all elements of a pathogen to this pathogen
	CPathogen add(const CPathogen&);
	// Check if there are any negative values
	bool anyNegative();

};


