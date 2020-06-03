#pragma once

class CFungicide{

public:

	// Constructor
	CFungicide();

	double currentDose;

	// Zero dose between seasons
	void resetBWSeasons();
};


