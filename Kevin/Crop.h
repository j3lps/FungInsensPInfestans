#pragma once

class CCrop{

public:

	CCrop();

	// Total leaf area index
	double totalAreaIndex;
	// Healthy leaf area index
	double healthyAreaIndex;

	// Zero all elements
	void zero();
	// Multiply all elements by a scalar
	CCrop multiply(double);
	// Add all elements of a pathogen to this pathogen
	CCrop add(const CCrop&);
	// Are any of the values negative
	bool anyNegative();

};



