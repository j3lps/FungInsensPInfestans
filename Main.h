#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include "Fungicide.h"
#include "Pathogen.h"
#include "Crop.h"

// SIMULATION VARIABLES
unsigned int processNumber;

// TIME VARIABLES
// How many time units in a season?
unsigned int tEOS;
// Default time step
double defTimeStep;

// GENES
// The total number of genotypes
unsigned int TOTALGENOTYPES;
// The total number of genes
unsigned int TOTALGENES;
// The number of virulence genes
unsigned int nVirulenceGenes;
// The number of resistant genes
unsigned int nResistanceGenes;

// This vector stores the mutation probability from genotype to genotype for a single diploid gene <parent> <offspring>
// Fungicide resistance mutation probability
std::vector< std::vector< double > > singleFungResGeneMutationArray;
// Virulence gene mutation propbability
std::vector< std::vector< double > > singleVirGeneMutationArray;

// This vector stores the proportion of each genotype in the primary inoculum from the end of the previous year
std::vector<double> primaryInocProp;

// This vector converts from the genotype number (of all genotypes) to the genotype of a single (specified) gene
std::vector<unsigned int> genotypeToDiploid1D;

// Pathogen
CPathogen oPathogen;
// Crop
CCrop oCrop;
// Fungicide(s)
std::vector<CFungicide> vecFungicide;

// Adjusted infection efficiency and sporulation rate depending on resistance genes and fungicide dose
std::vector<double> infectionEfficiency;
std::vector<double> sporulationRate;
std::vector<double> latentPeriod;
// Base infection efficiency when virulence genes are taken into account - this stays constant throughout the sim
std::vector<double> baseIE;
std::vector<double> baseSR;
std::vector<double> baseLP;

// Iterate through a year
void iterate(double tStart, double tEnd, unsigned int year);
// Perform a 4th order runge-kutta
bool rungeKutta(double tNow, double timeStep);
// Calculate the derivatives -- takes inputs of temporary input and derivative
void deriv(double tNow, CCrop HA_IN, CCrop& HA_DV, const CPathogen &PA_IN, CPathogen& PA_DV, 
	const std::vector<CFungicide> &FU_IN, std::vector<CFungicide> &FU_DV);

// Reset the variables between years
void resetBwYears();
// Store the variables at the end of the year
void storeYearResults(unsigned int year);
// Store the results each day
void storeResults(double tNow, unsigned int year);
// Initialize densities
void initialize();

// Write results to file (2 files; 1 for each day, 1 for each year)
void writeResultsToFile();

// Print base infection efficiency and latent period
void printIE();
void printLP();
void printSR();

// Set the biological parameters
void setParameters();

// Work out the total number of spores of each genotype, taking into account mutation
void mutateOrDontIDontCare(const CPathogen&, std::vector<double>&);
// Calculate the current infection efficiency, dependent on the fungicide dose
void calculateInfectionEfficiency(const std::vector<CFungicide>&);
void calculateSporulationRate(const std::vector<CFungicide>&);
void calculateLatentPeriod(const std::vector<CFungicide>&);
// Deprecated function
void oldCalculateSporulationRate(const std::vector<CFungicide>&);
void oldCalculateInfectionEfficiency(const std::vector<CFungicide>&);
// Calculate the base infection efficiency for each genotype according to virulence
void calculateBaseIE();
// Calculate the base sporulation rate for each genotype dependent on virulence
void calculateBaseSR();
// Calculate the latent period for each genotype
void calculateBaseLP();
// Create the matrices for mutation
void createMutationMatrices();

// Calculate the current severity
double calcSeverity();

// Check whether a spray needs to be applied
void checkSprayTimes(double time, double timeStep);

// ************* PARAMETERS ***************

// CROP
// Initial crop leaf area index
double cropStartingArea;

// Potato tuber growth rate
double tuberGrowthRate;

// Parameters used in growth and senescense functions
double aCropParam, bCropParam, cCropParam, mCropParam, nCropParam;

// A vector specifying if the crop has receptors against each of the pathogen virulence genes
std::vector<bool> cropReceptor;

// PATHOGEN
// The exponent for how fast primary inoculum decays with time
double PIb;
double PIxi;

// The proprtion of spores that are produced and move to the overwinter pool
double Vtnu;
double Vtt0;

// Natural mortality rates for latent and infectoius lesions
double latentMortalityRate;
double infectiousMortalityRate;

// Lifespan of latent and infectious lesions
double latentLifespan;
double infectiousLifespan;

// Default infection efficiency
double defInfEff;
// Default sporulation rate
double defSporRate;

// Proportion of spores that leave the field
double propSporesLeavingField;

// Fitness cost to infection efficiency - expressed as a proportional reduction in the infection efficiency - i.e. 0.2 = 20% reduction; 1.0 = 100% reduction
double fitnessCostIE;
// The dominance of the fitness cost to infection efficiency: if Dom == 1, SR = RR fitnessCost; if Dom==0, SR = SS = no fitness cost
double fitnessCostIEDom;
// Fitness cost to sporulatoin rate - expressed as a proportional reduction in the sporulation rate - i.e. 0.2 = 20% reduction; 1.0 = 100% reduction
double fitnessCostSR;
// The dominance of the fitness cost to sporulation rate: if Dom == 1, SR = RR fitnessCost; if Dom==0, SR = SS = no fitness cost
double fitnessCostSRDom;
// Fitness cost to latent period - expressed as a proportional *rise* in the latent period
// i.e. fitness cost = 0.2 -> latent period = latent period * 1.2; fitness cost = 1.0 -> latent period = latent period * 2.0
double fitnessCostLP;
// The dominance fo the fitness cost to the latent period
double fitnessCostLPDom;

// Proportional reduction in infection efficiency as a result of being avirulent
double AVIRReductionIE;
// Dominance of the avirulence reduction
double AVIRDomIE;
// Proportional reduction in sporulation rate as a result of being avirulent
double AVIRReductionSR;
// Dominance of the avirulence reduction
double AVIRDomSR;
// Proportional *rise* in latent period as a result of being avirulentd
// i.e. rise = 0.2 -> latent period = latent period * 1.2; rise = 1.0 -> latent period = latent period * 2.0
double AVIRReductionLP;
// Doimnance
double AVIRDomLP;

// Fitness cost in infection efficiency as a result of being virulent
double fCReductionIE;
// Dominance of the virulent fitness cost reduction
double fCDomIE;
// Fitness cost in sporulation rate as a result of being virulent
double fCReductionSR;
// Dominance of the virulent fitness cost reduction
double fCDomSR;
// Proportional *rise* in latent period as a result of having a fitness cost
// i.e. fitness cost = 0.2 -> latent period = latent period * 1.2; fitness cost = 1.0 -> latent period = latent period * 2.0
double fCReductionLP;
// Doimnance
double fCDomLP;

// Mutation rate of one allele for a fungicide resistance gene
double mutationRateFungR;
// Mutation rate of one allele for a virulence gene
double mutationRateVir;
// Big mutation matrix - [parent genotype][offspring genotype] - will be 3^n x 3^n i.e. 3^2n large
std::vector< std::vector< double > > mutationMatrix;

// FUNGICIDE
// The number of fungicides
unsigned int nFungicides;

// Fungicide spray times and doses - each spray consists of a spray time (time units), and a dose
std::vector<std::pair<unsigned int, double> > sprayFung1;
std::vector<std::pair<unsigned int, double> > sprayFung2;

// Whether each resistance gene confers resistance to each fungicide
std::vector< std::vector< bool > > confersResistanceToFung;

// Decay rate of the fungicide
double fungDecayRate;

// Fungicide ~ infection efficiency parameters
// Dominance of the fungicide resistance gene
double fungResDom;
// Shape of the curve determinng how much each genotype contributes to resistance
double fungResPi;
// alphaMax is the maximum possible reduction in the infection efficiency
std::vector<double> alphaMax;
// kappa is the coefficient to the dose
std::vector<double> kappa;

// alpha is the maximum possible reduction of the infection efficiency - deprecated
std::vector<double> alpha;
// beta is the curvature of the dose-response curve - deprecated
std::vector<double> beta;

// OUTPUT
// Day results
struct SDayResults{
	std::vector<unsigned int> Year;
	std::vector<unsigned int> DegreeDay;
	std::vector<CPathogen> pathogenDensity;
	std::vector<double> healthyAreaIndex;
	std::vector<double> leafAreaIndex;
	std::vector<double> tubarWeight;
	std::vector< std::vector<double> > fungicideDose;
	std::vector< std::vector< double > > geneFreq;
	std::vector<double> totalDensity;
	std::vector<double> totalLatent;
	std::vector<double> totalInfectious;
	std::vector<double> severity;
};

SDayResults DayResults;

// Year results
struct SYearResults{
	std::vector<double> Year;
	std::vector<std::vector<double> > geneFreq;
	std::vector<double> tubarWeight;
};

SYearResults YearResults;

