#include "Main.h"

int main(int argc, char** argv) {

	if (argc == 2) {
		std::istringstream ss(argv[argc - 1]);
		if (!(ss >> processNumber)) std::cerr << "Invalid input " << argv[argc - 1] << "\n";
		std::cout << "Running simulation " << processNumber << std::endl;
	}
	else processNumber = 0;

	printDailyResults = false;
	printYearResults = true;

	// Set all the parameters for this simulation
	setParameters();

	// First run the simulation for a year without any pathogen
	runModelNoPathogen();

	// Run the model for a single year with each fungicide, and each resistant cultivar
	//runModelNoControl();

	// Run the model with each control singly for 100 years
	runModelEachControl();

	// Run the model specified initially for 100 years
	runModel(100);

	// Write the results to a file
	writeResultsToFile();

	return 0;

}

void runModelNoPathogen() {

	// Turn on printing daily results
	bool tempPDR = printDailyResults;
	printDailyResults = true;

	// Initialize densities
	initialize();

	// Turn off all primary inoculum
	double tempPI = primaryInoculum;
	primaryInoculum = 0.0;

	// Run the model for a single year
	runModel(1);

	// Write the results to a file
	writeResultsToFile("NPath");

	// Switch daily results back to default
	printDailyResults = tempPDR;
	primaryInoculum = tempPI;

	// Reset the recording objects
	reset();

}

void runModelNoControl() {

	// Turn on printing daily results
	bool tempPDR = printDailyResults;
	printDailyResults = true;

	// Initialize densities
	initialize();

	// Turn off all forms of control
	// 1. Set fungicides to dose zero
	std::vector<std::vector<std::pair<unsigned int, double>>> tempSF = sprayFung;
	for (std::vector<std::vector<std::pair<unsigned int, double>>>::iterator vecIt = sprayFung.begin(); vecIt != sprayFung.end(); ++vecIt) {
		for (std::vector<std::pair<unsigned int, double>>::iterator vecIt2 = vecIt->begin(); vecIt2 != vecIt->end(); ++vecIt2) {
			vecIt2->second = 0.0;
		}
	}
	// 2. Set the cultivar to be susceptible
	double tempDIE = defInfEff;
	setLifecycleParms(1);

	// Run the model for a single year
	runModel(1);

	// Write the results to a file
	writeResultsToFile("NoControl");

	// Switch daily results back to default
	sprayFung = tempSF;
	double posCult = tempDIE / (0.0165 * 0.01373 * 0.921);
	if (posCult < 0.9) {
		setLifecycleParms(1);
	}
	else if (posCult > 1.1) {
		setLifecycleParms(3);
	}
	else {
		setLifecycleParms(2);
	}

	// Reset the recording objects
	reset();

}

void runModelEachControl() {

	// Set all fungicide doses to zero
	std::vector<std::vector<std::pair<unsigned int, double>>> tempSF = sprayFung;
	for (std::vector<std::vector<std::pair<unsigned int, double>>>::iterator vecIt = sprayFung.begin(); vecIt != sprayFung.end(); ++vecIt) {
		for (std::vector<std::pair<unsigned int, double>>::iterator vecIt2 = vecIt->begin(); vecIt2 != vecIt->end(); ++vecIt2) {
			vecIt2->second = 0.0;
		}
	}
	// Store the initial default infection efficiency (so we can reset the cultivar)
	double tempDIE = defInfEff;

	// Run each cultivar by itself
	for (unsigned int iCult = 1; iCult <= 3; ++iCult) {

		setLifecycleParms(iCult);

		initialize();

		runModel(100);

		writeResultsToFile("Cult" + std::to_string(iCult));

		reset();

	}

	// Loop over each cultivar and fungicide combination
	for (unsigned int iCult = 1; iCult <= 3; ++iCult) {
		setLifecycleParms(iCult);
		for (unsigned int iFung = 0; iFung != nFungicides; ++iFung) {

			// Add the dose for a single fungicide
			for (std::vector<std::pair<unsigned int, double>>::iterator vecIt = sprayFung[iFung].begin(); vecIt != sprayFung[iFung].end(); ++vecIt) {
				vecIt->second = 1.6;
			}

			initialize();

			runModel(100);

			writeResultsToFile("Cult" + std::to_string(iCult) + "Fung" + std::to_string(iFung));

			reset();

		}
	}

	// Reset the fungicides and cultivars to how they were before
	sprayFung = tempSF;
	double posCult = tempDIE / (0.0165 * 0.01373 * 0.921);
	if (posCult < 0.9) {
		setLifecycleParms(1);
	}
	else if (posCult > 1.1) {
		setLifecycleParms(3);
	}
	else {
		setLifecycleParms(2);
	}

}

void reset() {

	// Reset between years
	resetBwYears(NULL);

	DayResults.reset();
	YearResults.reset();

}

void runModel(unsigned int NYEARS) {

	// Initialize densities
	initialize();

	// Iterate over the number of years required
	for (unsigned int year = 1; year != (NYEARS + 1); ++year) {

		// Store the year's initial conditions
		storeResults(0.0, year);

		std::cout << "Year " << year << " commencing." << std::endl;

		// Iterate between the start time and the end of the season
		iterate(0.0, tEOS, year);

		storeYearResults(year);

		resetBwYears(year);
	}

}

void initialize() {

	oPathogen.zero();
	primaryInocProp.assign(TOTALGENOTYPES, 0.0);

	// Put all density as the fully avirulent, fully susceptible individual.
	primaryInocProp[0] = 1.0;

	oCrop.totalAreaIndex = oCrop.healthyAreaIndex = cropStartingArea;
	vecFungicide.assign(nFungicides, CFungicide());

}

void setParameters() {

	defTimeStep = 1.0;

	// Specify the number of virulence genes the pathogen has. This is the same as the number of resistance genes in the crop.
	nVirulenceGenes = 0; // No virulence, because virulence is not evolving
	// Specify the number of fungicide resistant genes the pathogen has
	nResistanceGenes = 1; // Target-site resistance to fungicides, so one resistance gene
	// Total number of genes
	TOTALGENES = nVirulenceGenes + nResistanceGenes;
	// Calculate the total number of genotypes
	TOTALGENOTYPES = pow(3.0, int(nVirulenceGenes + nResistanceGenes));
	CPathogen::TOTALGENOTYPES = TOTALGENOTYPES;

	// Set up genotype to diploid array
	genotypeToDiploid1D.assign(TOTALGENOTYPES * TOTALGENES, 0);
	for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene) {
		for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {
			unsigned int genotype = iGeno;
			for (unsigned int jGene = iGene + 1; jGene != TOTALGENES; ++jGene) genotype = div(genotype, 3).quot;
			genotype = div(genotype, 3).rem;
			genotypeToDiploid1D[iGeno + iGene * TOTALGENOTYPES] = genotype;
		}
	}

	// Crop starting density - set so that maximum healthy area = 6.044 - same as Femke's. was 0.0000000663
	cropStartingArea = 5e-3; //5e-3 LLanilar; 1.15e-2 Auchincruive
	aCropParam = 6;//6 Llanilar; 6 Auchincruive
	bCropParam = 90; //90 Llanilar; 90 Auchincruive
	cCropParam = 633; //633 Llanilar; 560 Auchincruive
	mCropParam = 1900; //1900 Llanilar; 1650 Auchincruive
	nCropParam = 80; //80 Llanilar; 80 Auchincruive

	// Timed so that healthy area finishes at 1.3 leaf index. was originally 2750, we updated host growth
	tEOS = 1700; // 1708 for Llanilar; 1408 for Auchincruive

	// Dictate whether the crop has receptors against each of the virulence genes
	// If it does and the pathogen is avirulent then the IE, LP and SR are reduced
	// If the cultivar has receptor and the pathogen is virulent then the IE, LP and SR are not reduced
	// If the cultivar does not have receptor the virulent have a fitness cost
	cropReceptor.assign(nVirulenceGenes, false);

	// Pathogen parameters

	// Quantity of primary inoculum
	primaryInoculum = 0.1;
	// Shape of primary inoculum curve (k)
	PIxi = 3;
	// Another shape of primary inoculum curve (theta)
	PIb = 100;

	// Specify the number of fungicides
	nFungicides = 1;

	// Specify the spray times:
	std::vector<double> sprayTimes;
	// Modified Llanilar spray program
	sprayTimes.push_back(500); sprayTimes.push_back(600); sprayTimes.push_back(700); sprayTimes.push_back(800); sprayTimes.push_back(900); sprayTimes.push_back(1000);
	sprayTimes.push_back(1100); sprayTimes.push_back(1200); sprayTimes.push_back(1300); sprayTimes.push_back(1400); sprayTimes.push_back(1500); sprayTimes.push_back(1600);
	// Specify the dose of each fungicide - normally the dose applied is the same each spray:
	std::vector<double> dose; dose.push_back(1.6); // Infinito
	// Now create the spray program for each fungicide. You can adjust this manually later if you want.
	for (unsigned int iF = 0; iF != nFungicides; ++iF) {
		std::vector<std::pair<unsigned int, double> > tempVecPair;
		for (std::vector<double>::iterator vecIt = sprayTimes.begin(); vecIt != sprayTimes.end(); ++vecIt) {
			std::pair<unsigned int, double> tempPair(*vecIt, dose[iF]);
			// If doing a mixture then only spray half dose
			if (iF == 1) tempPair.second *= 0.5;
			tempVecPair.push_back(tempPair);
		}
		sprayFung.push_back(tempVecPair);
	}

	// Fungicide ~ infection efficiency parameters
	alphaMax.assign(nFungicides, 1.0);//0.25 * std::div(processNumber,5).rem);      // default is 1.0, numbers above 1 will not be biologically realistic.  
	kappa.assign(nFungicides, 0.1); //5.86 * 0.25 * std::div(processNumber, 5).rem
	// Make the 2nd fungicide less effective
	if (nFungicides == 2) alphaMax[1] = 0.5;

	// Decide, for each fungicide resistance gene, whether it confers resistance to each fungicide
	confersResistanceToFung.assign(nResistanceGenes, std::vector<bool>(nFungicides, true));
	if (nFungicides == 2) confersResistanceToFung[0][1] = false;

	// survival to next season
	Vtnu = 20;		//default is 20
	Vtt0 = 1600;   // 200gdd, 2 weeks, before burnoff of crop. 

	// Default infection efficiency
	defInfEff = 0.0165; // default is 0.0165

	// Natural mortality rates for latent and infectoius lesions
	latentMortalityRate = 0.0;     //default is 0
	infectiousMortalityRate = 0.0; //default is 0

	// Latent period
	latentLifespan = 50.0; //default is 50.66542  
	// Infectious period
	infectiousLifespan = 100.0;  //default is 102

	// Default sporulation rate, Rho0. 
	defSporRate = 200; // Default is 200

	// Site-specific transmission coefficient - 0.0364 is zeta, 0.901 is for Llanilar relative to Auchincruive
	transCoef = 0.0364 * 0.901;

	// The coefficients for each cultivar: 0.296 for Lady Balfour; 0.294 for Maris Piper; 0.251 for Shepody
	std::vector<double> cultRes; cultRes.push_back(0.296); cultRes.push_back(0.294); cultRes.push_back(0.251);

	// Choose the cultivar that you want to model
	unsigned int iCult = 0;

	// AVIRReduction = 0.0 implies no crop resistance; Sarpo Mira implies that AVIRReduction = 0.3 (1.0 - 0.8/1.15)
	AVIRReductionLP = cultRes[iCult];
	AVIRReductionIE = cultRes[iCult];
	AVIRReductionSR = cultRes[iCult];

	// Infection efficiency, latent period and sporulation rate vectors (for each genotype)
	baseIE.assign(TOTALGENOTYPES, defInfEff);
	baseLP.assign(TOTALGENOTYPES, latentLifespan);
	baseSR.assign(TOTALGENOTYPES, defSporRate);

	// Initialise infection efficiency and latent period vectors - these change during each season due to fungicide doses.
	infectionEfficiency.assign(TOTALGENOTYPES, defInfEff);
	latentPeriod.assign(TOTALGENOTYPES, latentLifespan);
	sporulationRate.assign(TOTALGENOTYPES, defSporRate);

	// Way to adjust fungicide insensitivity for partial insensitivity
	fungResPi.assign(nFungicides, 1.0);		// 1.0 = complete insensitivity; 0.0 = complete sensitivity of the insensitive strain; 0.6 for partial resistance.

	// Fitness costs to virulence, default is 0.002, dominance is 0.5
	// Fitness cost in infection efficiency as a result of being virulent (should be less than result of host resistance above)
	fCReductionIE = 0.002;// default is 0.002
	fCDomIE = 0.5;
	// Proportional *rise* in latent period as a result of having a fitness cost (should be less than result of host resistance above)
	// i.e. fitness cost = 0.2 -> latent period = latent period * 1.2; fitness cost = 1.0 -> latent period = latent period * 2.0
	fCReductionLP = 0.002;
	fCDomLP = 0.5;
	// Fitness cost and dominance of sporulation rate
	fCReductionSR = 0.002;
	fCDomSR = 0.5;

	// Proportion of spores leaving the field
	propSporesLeavingField = 0.0; //default is 0

	// Fitness costs to fungicide insensitivity, default is 0.002, dom is 0.5
	// Fitness costs of fungicide resistance (proportional reduction in infection efficiency for each RR gene)
	fitnessCostIE = 0.002;  // default is 0.002 in the QTL baed model, 0.001 in the allele based. 
	// Dominance of fitness cost to infection efficiency - i.e. defines fitness cost for SR
	fitnessCostIEDom = 0.5; //default is 0.5 under QTL based model, 1 under allele based model. See "descrbing partial resistance in a diploid system.doc". 
	// Fitness cost to latent period (multiply the latent period by this - make it longer)
	fitnessCostLP = 0.002;
	fitnessCostLPDom = 0.5;
	fitnessCostSR = 0.002;
	fitnessCostSRDom = 0.5;

	// theta dominance. default 0.5
	// Dominance of the avirulence reduction
	AVIRDomIE = 0.5;
	// Dominance of the avirulence reduction
	AVIRDomSR = 0.5;
	// Doimnance
	AVIRDomLP = 0.5;

	// Mutation rate of a single allele. Normally either 0.000001 or 0. 
	mutationRateFungR = 1e-10; //1e-10
	mutationRateVir = 1e-10; //1e-10

	fungResDom = 0.5; //in QTL model this is dominance, default is 0.5 (QTL model), based on literature. In allele model it is not dominance but expression level, which is 1.0 by default. 

	// Decay rate of the fungicide, 6 Jdays. 
	fungDecayRate = 0.007967; // default is 0.007967209, 6 julian days.

	// Functions that adjust the latent period, infection efficiency and sporulation rate for each genotype depending on the cultivar resistance
	calculateBaseLP();
	calculateBaseIE();
	calculateBaseSR();

	// Create a matrix specifying mutation between different genotypes
	createMutationMatrices();

}

void iterate(double tStart, double tEnd, unsigned int year) {

	double tNow = tStart;
	double tPrev = 0.0;

	// This is the default time step - if things go negative we might have to reduce the size of the time step
	double tempTimeStep = defTimeStep;

	while (floor(tNow + 0.5 * tempTimeStep) < floor(tEnd + 0.5 * tempTimeStep)) {

		// Check if a spray needs to be applied
		checkSprayTimes(tNow, tempTimeStep);

		while (!rungeKutta(tNow, tempTimeStep)) {
			tempTimeStep = tempTimeStep / 10.0;
			std::cout << "Trying a smaller timestep: " << tempTimeStep << "\n";
			if (tempTimeStep < 1E-20) {
				std::cerr << "Time step is too small: " << tempTimeStep << ". Please do something." << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		tPrev = tNow;
		tNow += tempTimeStep;

		// If a day has passed, record results
		if (floor(tNow + 0.5 * tempTimeStep) > floor(tPrev + 0.5 * tempTimeStep)) {
			storeResults(floor(tNow + 0.5 * tempTimeStep), year);
			tempTimeStep = defTimeStep;
		}

	}

}

void checkSprayTimes(double time, double timeStep) {

	// Loop over all the fungicides
	for (unsigned int iF = 0; iF != nFungicides; ++iF) {
		if (sprayFung[iF].size() > 0) {
			for (unsigned int iSpray = 0; iSpray != sprayFung[iF].size(); ++iSpray) {

				if (floor(time - 0.5 * timeStep) < floor(time + 0.5 * timeStep) && floor(time) == sprayFung[iF][iSpray].first) {
					if (sprayFung[iF][iSpray].second > 0.0) std::cout << "Spraying fungicide " << iF << " at dose " << sprayFung[iF][iSpray].second << "!" << std::endl;
					vecFungicide[iF].currentDose += sprayFung[iF][iSpray].second;
				}

			}
		}
	}

}

bool rungeKutta(double tNow, double timeStep) {

	// Make necessary vectors, and set everything to zero or the initial values.
	// For each variable, need to have: initial value (IV), temporary input (TI), k1, k2, k3, k4, derivative (DV)

	double inputTime = tNow;

	// *** Crop ***
	// One variable: the healthy area (HA)
	CCrop HA_IV, HA_TI, HA_k1, HA_k2, HA_k3, HA_k4, HA_DV;

	// *** Pathogen ***
	// Variables: Latent, Infectious of each genotype
	CPathogen PA_IV, PA_TI, PA_k1, PA_k2, PA_k3, PA_k4, PA_DV;
	PA_IV.zero(); PA_TI.zero(); PA_k1.zero(); PA_k2.zero(); PA_k3.zero(); PA_k4.zero(); PA_DV.zero();

	// *** Fungicide ***
	// Variables: Potentially the dose of two fungicides
	std::vector<CFungicide> FU_IV(nFungicides, CFungicide()), FU_TI(nFungicides, CFungicide()), FU_k1(nFungicides, CFungicide()),
		FU_k2(nFungicides, CFungicide()), FU_k3(nFungicides, CFungicide()), FU_k4(nFungicides, CFungicide()), FU_DV(nFungicides, CFungicide());

	// Set initial values
	// *** Crop ***
	HA_IV = oCrop;
	// *** Pathogen ***
	PA_IV = oPathogen;
	// *** Fungicide ***
	FU_IV = vecFungicide;

	// ~-~-~-~~~-~-~-~ Runge-Kutta 4th order ~-~-~-~~~-~-~-~

	// _-*-_-*-_-*-_ Calculate k1 _-*-_-*-_-*-_

	// Set inputs - for k1 this is simply the initial values
	inputTime = tNow;
	// *** Crop ***
	HA_TI = HA_IV;
	// *** Pathogen ***
	PA_TI = PA_IV;
	// *** Fungicide ***
	FU_TI = FU_IV;

	// Re-zero all derivative vectors
	// *** Crop ***
	HA_DV.zero();
	// *** Pathogen ***
	PA_DV.zero();
	// *** Fungicide ***
	FU_DV.assign(nFungicides, CFungicide());

	// Calculate derivatives
	deriv(inputTime, HA_TI, HA_DV, PA_TI, PA_DV, FU_TI, FU_DV);

	// Store k1
	// *** Crop ***
	HA_k1 = HA_DV;
	// *** Pathogen ***
	PA_k1 = PA_DV;
	// *** Fungicide ***
	FU_k1 = FU_DV;

	// _-*-_-*-_-*-_ Calculate k2 _-*-_-*-_-*-_

	// Set inputs - for k2 this = IV + 0.5 * h * k1
	inputTime = tNow + timeStep / 2.0;
	// *** Crop ***
	HA_TI = HA_k1.multiply(0.5);
	HA_TI = HA_TI.multiply(timeStep);
	HA_TI = HA_TI.add(HA_IV);
	// *** Pathogen ***
	PA_TI = PA_k1.multiply(0.5);
	PA_TI = PA_TI.multiply(timeStep);
	PA_TI = PA_TI.add(PA_IV);
	// *** Fungicide ***
	for (std::size_t iFun = 0; iFun != vecFungicide.size(); ++iFun) {
		FU_TI[iFun].currentDose = FU_IV[iFun].currentDose + 0.5 * timeStep * FU_k1[iFun].currentDose;
	}

	// Re-zero all derivative vectors
	// *** Crop ***
	HA_DV.zero();
	// *** Pathogen ***
	PA_DV.zero();
	// *** Fungicide ***
	FU_DV.assign(nFungicides, CFungicide());

	// Calculate derivatives
	deriv(inputTime, HA_TI, HA_DV, PA_TI, PA_DV, FU_TI, FU_DV);

	// Store k2
	// *** Crop ***
	HA_k2 = HA_DV;
	// *** Pathogen ***
	PA_k2 = PA_DV;
	// *** Fungicide ***
	FU_k2 = FU_DV;

	// _-*-_-*-_-*-_ Calculate k3 _-*-_-*-_-*-_

	// Set inputs -- = IV + 0.5 * timeStep * k2
	inputTime = tNow + timeStep / 2.0;
	// *** Crop ***
	HA_TI = HA_k2.multiply(0.5);
	HA_TI = HA_TI.multiply(timeStep);
	HA_TI = HA_TI.add(HA_IV);
	// *** Pathogen ***
	PA_TI = PA_k2.multiply(0.5);
	PA_TI = PA_TI.multiply(timeStep);
	PA_TI = PA_TI.add(PA_IV);
	// *** Fungicide ***
	for (std::size_t iFun = 0; iFun != vecFungicide.size(); ++iFun) {
		FU_TI[iFun].currentDose = FU_IV[iFun].currentDose + 0.5 * timeStep * FU_k2[iFun].currentDose;
	}

	// Re-zero all derivative vectors
	// *** Crop ***
	HA_DV.zero();
	// *** Pathogen ***
	PA_DV.zero();
	// *** Fungicide ***
	FU_DV.assign(nFungicides, CFungicide());

	// Calculate derivatives
	deriv(inputTime, HA_TI, HA_DV, PA_TI, PA_DV, FU_TI, FU_DV);

	// Store k3
	// *** Crop ***
	HA_k3 = HA_DV;
	// *** Pathogen ***
	PA_k3 = PA_DV;
	// *** Fungicide ***
	FU_k3 = FU_DV;

	// _-*-_-*-_-*-_ Calculate k4 _-*-_-*-_-*-_

	// Set inputs -- = IV + k3 * timeStep
	inputTime = tNow + timeStep;
	// *** Crop ***
	HA_TI = HA_k3.multiply(timeStep);
	HA_TI = HA_TI.add(HA_IV);
	// *** Pathogen ***
	PA_TI = PA_k3.multiply(timeStep);
	PA_TI = PA_TI.add(PA_IV);
	// *** Fungicide ***
	for (std::size_t iFun = 0; iFun != vecFungicide.size(); ++iFun) {
		FU_TI[iFun].currentDose = FU_IV[iFun].currentDose + timeStep * FU_k3[iFun].currentDose;
	}

	// Re-zero all derivative vectors
	// *** Crop ***
	HA_DV.zero();
	// *** Pathogen ***
	PA_DV.zero();
	// *** Fungicide ***
	FU_DV.assign(nFungicides, CFungicide());

	// Calculate derivatives
	deriv(inputTime, HA_TI, HA_DV, PA_TI, PA_DV, FU_TI, FU_DV);

	// Store k4
	// *** Crop ***
	HA_k4 = HA_DV;
	// *** Pathogen ***
	PA_k4 = PA_DV;
	// *** Fungicide ***
	FU_k4 = FU_DV;


	// Calculate final values (put in derivative (DV) variables, since these aren't needed anymore)
	// Final value = IV + (timeStep / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

	// Going to check that none of the values are negative
	bool noNeg = true;
	// *** Crop ***
	HA_k2 = HA_k2.multiply(2.0);
	HA_k3 = HA_k3.multiply(2.0);
	HA_k1 = HA_k1.add(HA_k2);
	HA_k1 = HA_k1.add(HA_k3);
	HA_k1 = HA_k1.add(HA_k4);
	HA_k1 = HA_k1.multiply(timeStep / 6.0);
	HA_DV = HA_IV.add(HA_k1);
	noNeg = !HA_DV.anyNegative();
	// *** Pathogen ***
	// First work out the multiple of ks
	PA_k2 = PA_k2.multiply(2.0);
	PA_k3 = PA_k3.multiply(2.0);
	PA_k1 = PA_k1.add(PA_k2);
	PA_k1 = PA_k1.add(PA_k3);
	PA_k1 = PA_k1.add(PA_k4);
	PA_k1 = PA_k1.multiply(timeStep / 6.0);
	PA_DV = PA_IV.add(PA_k1);
	noNeg = !PA_DV.anyNegative();
	// *** Fungicide ***
	for (size_t iFC = 0; iFC != vecFungicide.size(); ++iFC) {
		FU_DV[iFC].currentDose = FU_IV[iFC].currentDose + (timeStep / 6.0) * (FU_k1[iFC].currentDose + 2.0 * FU_k2[iFC].currentDose
			+ 2.0 * FU_k3[iFC].currentDose + FU_k4[iFC].currentDose);
		if (FU_DV[iFC].currentDose < 0.0) noNeg = false;
	}

	// If none of the final variables are negative then set them
	if (noNeg) {
		// *** Crop ***
		oCrop = HA_DV;
		// *** Pathogen ***
		oPathogen = PA_DV;
		// *** Fungicide ***
		for (unsigned int iFC = 0; iFC != vecFungicide.size(); ++iFC) vecFungicide[iFC].currentDose = FU_DV[iFC].currentDose;
	}
	else {
		std::cout << ".";
	}

	return noNeg;

}

void deriv(double tNow, CCrop HA_IN, CCrop& HA_DV, const CPathogen& PA_IN, CPathogen& PA_DV,
	const std::vector<CFungicide>& FU_IN, std::vector<CFungicide>& FU_DV) {

	// *** Crop stuff
	// Work out total crop area
	double totalCropArea = HA_IN.healthyAreaIndex;
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {
		totalCropArea += PA_IN.Latent[iGeno];
		totalCropArea += PA_IN.Infectious[iGeno];
	}

	double LAI = aCropParam / (1 + exp(-(tNow - cCropParam) / bCropParam));
	double SAI = aCropParam / (1 + exp(-(tNow - mCropParam) / nCropParam));
	double growthRate = ((aCropParam * exp(-(tNow - cCropParam) / bCropParam)) / (bCropParam * (1 + exp(-(tNow - cCropParam) / bCropParam)) * (1 + exp(-(tNow - cCropParam) / bCropParam)))) / (LAI - SAI);
	double senesRate = ((aCropParam * exp(-(tNow - mCropParam) / nCropParam)) / (nCropParam * (1 + exp(-(tNow - mCropParam) / nCropParam)) * (1 + exp(-(tNow - mCropParam) / nCropParam)))) / (LAI - SAI);

	HA_DV.totalAreaIndex += growthRate * HA_IN.totalAreaIndex;;
	HA_DV.healthyAreaIndex += growthRate * HA_IN.totalAreaIndex;;
	HA_DV.healthyAreaIndex -= senesRate * HA_IN.healthyAreaIndex;

	// Calculate the infection efficiency of each genotype at t = tNow
	newCalculateInfectionEfficiency(FU_IN);
	// Calculate the latent period
	newCalculateLatentPeriod(FU_IN);

	// Calculate the number of spores of each genotype - function of spores produced and mutation rate
	std::vector<double> secondarySpores(TOTALGENOTYPES, 0.0);
	mutateOrDontIDontCare(PA_IN, secondarySpores);

	// For each fungicide, make the dose decay
	for (size_t i = 0; i != FU_DV.size(); ++i) {
		FU_DV[i].currentDose -= fungDecayRate * FU_IN[i].currentDose;
	}

	// Calculate the total amount of primary inoculum at tNow - gamma distribution multiplied by pInoc
	double pI = primaryInoculum * pow(tNow, PIxi - 1) * exp(-tNow / PIb) / (tgamma(PIxi) * pow(PIb, PIxi));

	// Calculate the deposition probability
	double depoProb = (1 - propSporesLeavingField) * (1 - exp(-HA_IN.totalAreaIndex));

	// Loop through the genotypes and calculate change in the density of lesions
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {

		// Store a proportion of spores for the next season
		PA_DV.sporesEnteringOWPool[iGeno] += secondarySpores[iGeno] / (1.0 + exp((-(tNow - Vtt0)) / Vtnu));

		// New latent from primary inoculum
		double newLatentThisGeno = depoProb * (pI * primaryInocProp[iGeno] + secondarySpores[iGeno]) * HA_IN.healthyAreaIndex / HA_IN.totalAreaIndex * infectionEfficiency[iGeno];
		PA_DV.Latent[iGeno] += newLatentThisGeno;
		// Remove this area from the healthy area
		HA_DV.healthyAreaIndex -= newLatentThisGeno;

		// Infectious from latent - includes a fitness cost to latent period
		PA_DV.Latent[iGeno] -= PA_IN.Latent[iGeno] / latentPeriod[iGeno];
		PA_DV.Infectious[iGeno] += PA_IN.Latent[iGeno] / latentPeriod[iGeno];

		// Infectious to removed
		PA_DV.Infectious[iGeno] -= PA_IN.Infectious[iGeno] / infectiousLifespan;

		// Natural mortality of lesions
		PA_DV.Latent[iGeno] -= latentMortalityRate * PA_IN.Latent[iGeno];
		PA_DV.Infectious[iGeno] -= infectiousMortalityRate * PA_IN.Infectious[iGeno];

		// Senesce lesions if necessary
		PA_DV.Latent[iGeno] -= senesRate * PA_IN.Latent[iGeno];
		PA_DV.Infectious[iGeno] -= senesRate * PA_IN.Infectious[iGeno];

	}

}

void mutateOrDontIDontCare(const CPathogen& inPathogen, std::vector<double>& outSpores) {
	// NB: outSpores is a vector of zeros - probably ought to check this but meh.

	std::vector<double> infectiousLesions = inPathogen.Infectious;

	// Loop through each genotype for the parent
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {

		// Need to work out the proportion of spores produced by the parent genotype that goes to each offspring genotype
		// Create a vector that will store the proportion of spores of each genotype
		std::vector<double> propSpores(TOTALGENOTYPES, 0.0);

		// Loop through each potential offspring genotype
		for (unsigned int jGeno = 0; jGeno != TOTALGENOTYPES; ++jGeno) {

			// Store the proportion of iGeno that mutates (or doesn't) into jGeno
			double proportionjGeno = 1.0;

			if (mutationMatrix.size() > 0) {
				proportionjGeno = mutationMatrix[iGeno][jGeno];
			}
			else {
				// Loop through each gene
				for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene) {

					unsigned int parentGenotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];
					unsigned int offspringGenotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + jGeno];

					// If we're looking at a virulence gene, then use virulence gene mutation array
					if (iGene < nVirulenceGenes) proportionjGeno *= singleVirGeneMutationArray[parentGenotype][offspringGenotype];
					else proportionjGeno *= singleFungResGeneMutationArray[parentGenotype][offspringGenotype];
				}
			}

			outSpores[jGeno] += infectiousLesions[iGeno] * defSporRate * proportionjGeno;

		}

	}

}

void newCalculateInfectionEfficiency(const std::vector<CFungicide>& fungicides) { // this needs to contain confersResistanceToFung (16.08.2017)
																			   // Set sporulation rate array back to default for every genotype
	infectionEfficiency = baseIE;

	// Loop through each genotype, and for each fungicide resistance genotype reduce the sporulation rate accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {

		double proportionReduction = 1.0;

		// Loop through each fungicide
		for (size_t iFung = 0; iFung != fungicides.size(); ++iFung) {

			// This stores the consequence of resistance
			double coefficient = 1.0;

			// The following genes are fungicide resistance genes. N.B. This is not meant to go from zero, unless there are no virulence genes. It is correct to go from nVirulenceGenes.
			for (unsigned int iGene = nVirulenceGenes; iGene != TOTALGENES; ++iGene) {
				// Work out the diploid for this gene
				unsigned int genotype = 0;
				// Only work out the genotype if this gene confers resistance to this fungicide, otherwise it's susceptible
				if (confersResistanceToFung[iGene - nVirulenceGenes][iFung]) genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];
				switch (genotype) {
					// 0 = SS - if susceptible, there is no resistance - don't need to reduce alpha_max
				case 0:
					coefficient *= (1.0 - 0.0);
					break;
					// 1 = SR - partial reduction in infection efficiency by fungicide
				case 1:

					coefficient *= (1.0 - fungResPi[iFung] * fungResDom);

					break;
					// 2 = RR - no reduction in the infection efficiency
				case 2:

					coefficient *= (1.0 - fungResPi[iFung]);
					break;
				}
			}
			proportionReduction *= (1 - alphaMax[iFung] * coefficient * (1 - exp(-kappa[iFung] * fungicides[iFung].currentDose)));
		}
		infectionEfficiency[iGeno] *= proportionReduction;
	}
}

void newCalculateSporulationRate(const std::vector<CFungicide>& fungicides) { // this needs to contain confersResistanceToFung (16.08.2017)
																			   // Set sporulation rate array back to default for every genotype
	sporulationRate = baseSR;

	// Loop through each genotype, and for each fungicide resistance genotype reduce the sporulation rate accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {

		double proportionReduction = 1.0;

		// Loop through each fungicide
		for (size_t iFung = 0; iFung != fungicides.size(); ++iFung) {

			// This stores the consequence of resistance
			double coefficient = 1.0;

			// The following genes are fungicide resistance genes. N.B. This is not meant to go from zero, unless there are no virulence genes. It is correct to go from nVirulenceGenes.
			for (unsigned int iGene = nVirulenceGenes; iGene != TOTALGENES; ++iGene) {
				// Work out the diploid for this gene
				unsigned int genotype = 0;
				// Only work out the genotype if this gene confers resistance to this fungicide, otherwise it's susceptible
				if (confersResistanceToFung[iGene - nVirulenceGenes][iFung]) genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];
				switch (genotype) {
					// 0 = SS - if susceptible, there is no resistance - don't need to reduce alpha_max
				case 0:
					coefficient *= (1.0 - 0.0);
					break;
					// 1 = SR - partial reduction in sporulation rate by fungicide
				case 1:

					coefficient *= (1.0 - fungResPi[iFung] * fungResDom);

					break;
					// 2 = RR - no reduction in the sporulation rate
				case 2:

					coefficient *= (1.0 - fungResPi[iFung]);
					break;
				}
			}
			proportionReduction *= (1 - alphaMax[iFung] * coefficient * (1 - exp(-kappa[iFung] * fungicides[iFung].currentDose)));
		}
		sporulationRate[iGeno] *= proportionReduction;
	}
}

void newCalculateLatentPeriod(const std::vector<CFungicide>& fungicides) {

	// Set sporulation rate array back to default for every genotype
	latentPeriod = baseLP;

	// Loop through each genotype, and for each resistance genotype extend the latent period accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {

		double proportionExtension = 1.0;

		// Loop through each fungicide
		for (size_t iFung = 0; iFung != fungicides.size(); ++iFung) {
			// The mixing partner doesn't affect the latent period
			if (iFung == 1) continue;

			// This stores the consequence of resistance
			double coefficient = 1.0;
			// The following genes are fungicide resistance genes
			for (unsigned int iGene = nVirulenceGenes; iGene != TOTALGENES; ++iGene) {
				// Work out the diploid for this gene
				unsigned int genotype = 0;
				// Only work out the genotype if this gene confers resistance to this fungicide, otherwise it's susceptible
				if (confersResistanceToFung[iGene - nVirulenceGenes][iFung]) genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];
				switch (genotype) {
					// 0 = SS - if susceptible, there is no resistance - don't need to reduce alpha_max
				case 0:
					coefficient *= (1.0 - 0.0);
					break;
					// 1 = SR - partial reduction in sporulation rate by fungicide
				case 1:

					coefficient *= (1.0 - fungResPi[iFung] * fungResDom);

					break;
					// 2 = RR - no reduction in the sporulation rate
				case 2:

					coefficient *= (1.0 - fungResPi[iFung]);
					break;
				}
			}
			proportionExtension *= 1 - alphaMax[iFung] * coefficient * (1 - exp(-kappa[iFung] * fungicides[iFung].currentDose));
		}
		// Convert so that maximum dose results in a doubling of the latent period
		latentPeriod[iGeno] *= (2 - proportionExtension);
	}
}

void calculateBaseIE() {

	baseIE.assign(TOTALGENOTYPES, defInfEff);

	// Loop through each genotype, and for each resistance genotype reduce the infection efficiency accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {

		// Calculate the fitness cost for this genotype
		double fitnessCost = 1.0;

		for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene) {

			// Work out the diploid for this gene
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			// If a virulence gene, then account for reduction in infection efficiency for avirulent alleles
			if (iGene < nVirulenceGenes) {

				// Store the proportional reduction for this gene
				double proportionReduction = 1.0;

				switch (genotype) {
					// 0 = AA - avirulent; infection efficiency is reduced
				case 0:
					if (cropReceptor[iGene]) proportionReduction *= (1 - AVIRReductionIE);
					else proportionReduction *= 1.0;
					break;
					// 1 = Aa - partially virulent
				case 1:
					if (cropReceptor[iGene]) proportionReduction *= (1 - fCReductionIE * fCDomIE);
					else proportionReduction *= (1 - fCReductionIE * fCDomIE);
					break;
					// 2 = aa - virulent; use default infection efficiency
				case 2:
					if (cropReceptor[iGene]) proportionReduction *= (1 - fCReductionIE);
					else proportionReduction *= (1 - fCReductionIE);
					break;

				}

				baseIE[iGeno] *= proportionReduction;
			}
			// Otherwise account for fitness cost of having fungicide resistance
			else {

				switch (genotype) {
					// 0 = SS - no fitness cost of having fungicide resistance
				case 0:
					break;
					// 1 = SR - partial reduction in infection efficiency due to fitness cost 
				case 1:
					fitnessCost *= (1 - fitnessCostIE * fitnessCostIEDom);
					break;
					// 2 = RR - full fitness cost of having fungicide resistance
				case 2:
					fitnessCost *= (1 - fitnessCostIE);
					break;
				}

			}

		}

		baseIE[iGeno] *= fitnessCost;
	}

	//printIE();

}

void calculateBaseSR() {

	baseSR.assign(TOTALGENOTYPES, defSporRate);

	// Loop through each genotype, and for each resistance genotype reduce the sporulation rate accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {

		// Calculate the fitness cost for this genotype
		double fitnessCost = 1.0;

		for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene) {

			// Work out the diploid for this gene
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			// If a virulence gene, then account for reduction in sporulation rate for avirulent alleles
			if (iGene < nVirulenceGenes) {

				// Store the proportional reduction for this gene
				double proportionReduction = 1.0;

				switch (genotype) {
					// 0 = AA - avirulent; infection efficiency is reduced
				case 0:
					if (cropReceptor[iGene]) proportionReduction *= (1 - AVIRReductionSR);
					else proportionReduction *= 1.0;
					break;
					// 1 = Aa - partially virulent
				case 1:
					if (cropReceptor[iGene]) proportionReduction *= (1 - fCReductionSR * fCDomSR);
					else proportionReduction *= (1 - fCReductionSR * fCDomSR);
					break;
					// 2 = aa - virulent; use default sporulation rate
				case 2:
					if (cropReceptor[iGene]) proportionReduction *= (1 - fCReductionSR);
					else proportionReduction *= (1 - fCReductionSR);
					break;

				}

				baseSR[iGeno] *= proportionReduction;
			}
			// Otherwise account for fitness cost of having fungicide resistance
			else {

				switch (genotype) {
					// 0 = SS - no fitness cost of having fungicide resistance
				case 0:
					break;
					// 1 = SR - partial reduction in infection efficiency due to fitness cost 
				case 1:
					fitnessCost *= (1 - fitnessCostSR * fitnessCostSRDom);
					break;
					// 2 = RR - full fitness cost of having fungicide resistance
				case 2:
					fitnessCost *= (1 - fitnessCostSR);
					break;
				}

			}

		}

		baseSR[iGeno] *= fitnessCost;
	}

}

void calculateBaseLP() {

	baseLP.assign(TOTALGENOTYPES, latentLifespan);

	// Loop through each genotype, and for each resistance genotype reduce the infection efficiency accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {

		// Calculate the fitness cost for this genotype
		double fitnessCost = 1.0;

		for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene) {

			// Work out the diploid for this gene
			unsigned int genotype = genotypeToDiploid1D[(iGene * TOTALGENOTYPES) + iGeno];

			double proportionalExtension = 1.0;

			// If virulence gene, then extend the latent period if it's avirulent
			if (iGene < nVirulenceGenes) {

				switch (genotype) {
					// 0 = AA - avirulent; infection efficiency is reduced
				case 0:
					if (cropReceptor[iGene]) proportionalExtension *= (1 + AVIRReductionLP);
					else proportionalExtension *= 1.0;
					break;
					// 1 = Aa - partially virulent
				case 1:
					if (cropReceptor[iGene]) proportionalExtension *= (1 + AVIRDomLP * AVIRReductionLP) * (1 + fCReductionLP * fCDomLP);
					else proportionalExtension *= (1 + fCReductionLP * fCDomLP);
					break;
					// 2 = aa - fully virulent; use default infection efficiency
				case 2:
					if (cropReceptor[iGene]) proportionalExtension *= (1 + fCReductionLP);
					else proportionalExtension *= (1 + fCReductionLP);
					break;
				}

				baseLP[iGeno] *= proportionalExtension;

			}
			else {

				switch (genotype) {
					// 0 = SS - no fitness cost
				case 0:
					break;
					// 1 = SR - heterozygote fitness cost
				case 1:
					fitnessCost *= (1 + fitnessCostLPDom * fitnessCostLP);
					break;
					// 2 = RR - no reduction in the infection efficiency
				case 2:
					fitnessCost *= (1 + fitnessCostLP);
					break;
				}
			}

		}

		baseLP[iGeno] *= fitnessCost;
	}

}

void createMutationMatrices() {

	// Create a mutation matrix 3x3, for the probability of mutating from one genotype to another
	singleFungResGeneMutationArray.assign(3, std::vector<double>(3, 1.0));
	// I'm just going to write this out for now, can program it if necessary - it's not hard
	singleFungResGeneMutationArray[0][0] = 1.0 - 2.0 * mutationRateFungR + mutationRateFungR * mutationRateFungR;
	singleFungResGeneMutationArray[0][1] = 2.0 * mutationRateFungR - 2.0 * mutationRateFungR * mutationRateFungR;
	singleFungResGeneMutationArray[0][2] = mutationRateFungR * mutationRateFungR;
	singleFungResGeneMutationArray[1][0] = -mutationRateFungR * mutationRateFungR + mutationRateFungR;
	singleFungResGeneMutationArray[1][1] = 1.0 + 2.0 * mutationRateFungR * mutationRateFungR - 2.0 * mutationRateFungR;
	singleFungResGeneMutationArray[1][2] = -mutationRateFungR * mutationRateFungR + mutationRateFungR;
	singleFungResGeneMutationArray[2][0] = mutationRateFungR * mutationRateFungR;
	singleFungResGeneMutationArray[2][1] = 2.0 * mutationRateFungR - 2.0 * mutationRateFungR * mutationRateFungR;
	singleFungResGeneMutationArray[2][2] = 1.0 - 2.0 * mutationRateFungR + mutationRateFungR * mutationRateFungR;

	// Create a mutation matrix 3x3, for the probability of mutating from one genotype to another
	singleVirGeneMutationArray.assign(3, std::vector<double>(3, 1.0));
	// I'm just going to write this out for now, can program it if necessary - it's not hard
	singleVirGeneMutationArray[0][0] = 1.0 - 2.0 * mutationRateVir + mutationRateVir * mutationRateVir;
	singleVirGeneMutationArray[0][1] = 2.0 * mutationRateVir - 2.0 * mutationRateVir * mutationRateVir;
	singleVirGeneMutationArray[0][2] = mutationRateVir * mutationRateVir;
	singleVirGeneMutationArray[1][0] = -mutationRateVir * mutationRateVir + mutationRateVir;
	singleVirGeneMutationArray[1][1] = 1.0 + 2.0 * mutationRateVir * mutationRateVir - 2.0 * mutationRateVir;
	singleVirGeneMutationArray[1][2] = -mutationRateVir * mutationRateVir + mutationRateVir;
	singleVirGeneMutationArray[2][0] = mutationRateVir * mutationRateVir;
	singleVirGeneMutationArray[2][1] = 2.0 * mutationRateVir - 2.0 * mutationRateVir * mutationRateVir;
	singleVirGeneMutationArray[2][2] = 1.0 - 2.0 * mutationRateVir + mutationRateVir * mutationRateVir;

	// If there's enough space, create a mutation matrix for each genotype
	std::vector<double> test;
	std::cout << "Vector maximum size is: " << test.max_size() << std::endl;
	std::cout << "Asking for total size to be: " << pow(3.0, double(TOTALGENES) * 2) << std::endl;
	if (test.max_size() > pow(3.0, double(TOTALGENES) * 2)) {
		for (unsigned int parentGeno = 0; parentGeno != TOTALGENOTYPES; ++parentGeno) {
			mutationMatrix.push_back(std::vector<double>(TOTALGENOTYPES, 1.0));
			for (unsigned int offspringGeno = 0; offspringGeno != TOTALGENOTYPES; ++offspringGeno) {

				// Store the proportion of parentGeno that mutates (or doesn't) into offspringGeno
				double proportionOffspringGeno = 1.0;

				// Loop through each gene
				for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene) {

					unsigned int parentGenotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + parentGeno];
					unsigned int offspringGenotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + offspringGeno];

					// If we're looking at a virulence gene, then use virulence gene mutation array
					if (iGene < nVirulenceGenes) proportionOffspringGeno *= singleVirGeneMutationArray[parentGenotype][offspringGenotype];
					else proportionOffspringGeno *= singleFungResGeneMutationArray[parentGenotype][offspringGenotype];
				}

				mutationMatrix[parentGeno][offspringGeno] = proportionOffspringGeno;

			}
		}
	}

}

void storeYearResults(unsigned int year) {

	YearResults.Year.push_back(year);
	// Work out the frequency of each gene
	std::vector<double> geneFreq(TOTALGENES, 0.0);
	// Loop through each genotype, and add to each gene if heterozygote or resistant
	double totalDensity = 0.0;
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {
		totalDensity += oPathogen.Infectious[iGeno];

		for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene) {
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			switch (genotype) {
				// case 0 = SS
			case 0:
				break;
			case 1:
				geneFreq[iGene] += 0.5 * oPathogen.Infectious[iGeno];
				break;
			case 2:
				geneFreq[iGene] += oPathogen.Infectious[iGeno];
			}
		}
	}
	if (totalDensity > 0.0) for (unsigned int i = 0; i != TOTALGENES; ++i) geneFreq[i] /= totalDensity;
	YearResults.geneFreq.push_back(geneFreq);
	// Calculate the HAD for this year
	double sumHAI = 0.0;
	for (std::size_t iDD = 0; iDD != DayResults.healthyAreaIndex.size(); ++iDD) if (DayResults.Year[iDD] == year) sumHAI += DayResults.healthyAreaIndex[iDD];
	double HAD = sumHAI / 14; // Dividing by 14 to convert from degree days to real days
	YearResults.HAD.push_back(HAD);
	// Calculate the AUDPC for this year - sum of severity
	double sumSeverity = 0.0;
	for (std::size_t iDD = 0; iDD != DayResults.severity.size(); ++iDD) if (DayResults.Year[iDD] == year) sumSeverity += DayResults.severity[iDD];
	double AUDPC = 100.0 * sumSeverity / 14; // Dividing by 14 to convert from degree days to real days
	YearResults.AUDPC.push_back(AUDPC);
	// Get the end of year severity
	YearResults.Severity.push_back(*(DayResults.severity.end() - 1));
}

void storeResults(double tNow, unsigned int year) {

	DayResults.DegreeDay.push_back(floor(tNow));
	DayResults.Year.push_back(year);
	DayResults.pathogenDensity.push_back(oPathogen);
	DayResults.healthyAreaIndex.push_back(oCrop.healthyAreaIndex);
	DayResults.leafAreaIndex.push_back(oCrop.totalAreaIndex);
	DayResults.severity.push_back(calcSeverity());
	std::vector<double> tempDose;
	for (unsigned int iFC = 0; iFC != nFungicides; ++iFC) tempDose.push_back(vecFungicide[iFC].currentDose);
	DayResults.fungicideDose.push_back(tempDose);
	// Work out the frequency of each gene
	std::vector<double> geneFreq(TOTALGENES, 0.0);
	// Loop through each genotype, and add to each gene if heterozygote or resistant
	double totalDensity = 0.0, totalLatent = 0.0, totalInfectious = 0.0;
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno) {
		totalDensity += oPathogen.Infectious[iGeno];
		totalInfectious += oPathogen.Infectious[iGeno];
		totalLatent += oPathogen.Latent[iGeno];
		totalDensity += oPathogen.Latent[iGeno];

		for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene) {
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			switch (genotype) {
				// case 0 = SS
			case 0:
				break;
			case 1:
				geneFreq[iGene] += 0.5 * oPathogen.Infectious[iGeno];
				geneFreq[iGene] += 0.5 * oPathogen.Latent[iGeno];
				break;
			case 2:
				geneFreq[iGene] += oPathogen.Infectious[iGeno];
				geneFreq[iGene] += oPathogen.Latent[iGeno];
			}
		}
	}
	if (totalDensity > 0.0) for (unsigned int i = 0; i != TOTALGENES; ++i) geneFreq[i] /= totalDensity;
	DayResults.totalDensity.push_back(totalDensity);
	DayResults.totalLatent.push_back(totalLatent);
	DayResults.totalInfectious.push_back(totalInfectious);
	DayResults.geneFreq.push_back(geneFreq);

}

void resetBwYears(unsigned int year) {

	double totalInfectiousDensity = 0.0;

	// Calculate the proportion of each genotype, and store in primaryInocProp
	// Sum the total amount of spores in the overwinter pool
	double totalSpores = 0.0;
	for (std::vector<double>::const_iterator cVecIt = oPathogen.sporesEnteringOWPool.begin(); cVecIt != oPathogen.sporesEnteringOWPool.end(); ++cVecIt) {
		totalSpores += *cVecIt;
	}
	if (totalSpores > 0.0) {
		// Make the proportion = to the sporesEnteringOWPool
		primaryInocProp = oPathogen.sporesEnteringOWPool;
		// Divide each genotype by the total Spores to make a proportion
		for (std::vector<double>::iterator vecIt = primaryInocProp.begin(); vecIt != primaryInocProp.end(); ++vecIt) (*vecIt) /= totalSpores;
	}
	else {
		primaryInocProp.assign(TOTALGENOTYPES, 0.0);
	}

	// Now rezero the pathogen, and crop
	oPathogen.zero();
	oCrop.zero();
	oCrop.healthyAreaIndex = oCrop.totalAreaIndex = cropStartingArea;

	// Zero the fungicide
	for (std::vector<CFungicide>::iterator vecIt = vecFungicide.begin(); vecIt != vecFungicide.end(); ++vecIt) vecIt->resetBWSeasons();

	// If we've implemented a distinct strategy for each year, then change the resistance level and fungicide dose.
	if (stratProgram.size() > 0) {
		AVIRReductionLP = stratProgram[year].first;
		AVIRReductionIE = stratProgram[year].first;
		calculateBaseIE();
		calculateBaseLP();
		// For each fungicide, store the spray times and doses. <Fungicide><sprayIndex><time,dose>
		for (std::vector<std::vector<std::pair<unsigned int, double> > >::iterator it = sprayFung.begin(); it != sprayFung.end(); ++it) {
			for (std::vector<std::pair<unsigned int, double>>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
				it2->second = stratProgram[year].second;
			}
		}
	}

	//// Add a constant amount to default dose of fungicide for the next year (sprayFung is <fungicide><sprayTime><Pair>)
	//bool printed = false;
	//for (std::vector<std::vector<std::pair<unsigned int, double>>>::iterator vecIt = sprayFung.begin(); vecIt != sprayFung.end(); ++vecIt) {
	//	for (std::vector<std::pair<unsigned int, double>>::iterator vecIt2 = vecIt->begin(); vecIt2 != vecIt->end(); ++vecIt2) {
	//		if (vecIt2->second < 1.6) {
	//			vecIt2->second += 0.05;
	//			if (vecIt2->second > 1.6) vecIt2->second = 1.6;
	//		}
	//		if (!printed) {
	//			std::cout << "Dose is now " << vecIt2->second << std::endl;
	//			printed = true;
	//		}
	//	}
	//}

}

double calcSeverity() {

	double severity;

	double infected = 0.0, total = 0.0;
	for (std::vector<double>::const_iterator cVecIt = oPathogen.Infectious.begin(); cVecIt != oPathogen.Infectious.end(); ++cVecIt) {
		infected += *cVecIt;
		total += *cVecIt;
	}
	for (std::vector<double>::const_iterator cVecIt = oPathogen.Latent.begin(); cVecIt != oPathogen.Latent.end(); ++cVecIt) {
		total += *cVecIt;
	}
	total += oCrop.healthyAreaIndex;

	if (oCrop.totalAreaIndex <= 0.0) {
		std::cerr << "Trying to calculate severity when totalAreaIndex is zero. Returning zero." << std::endl;
		return(0.0);
	}
	else {
		severity = infected / total;
	}

	return severity;

}

void writeResultsToFile(const std::string& fileNameAddition) {

	std::ofstream myfile;
	if (printYearResults) {
		// Create and open a file
		myfile.open("Year" + std::to_string(processNumber) + fileNameAddition + ".csv");
		// Write a header
		myfile << "Year";
		for (unsigned int i = 0; i != TOTALGENES; ++i) myfile << ", Gene" << i;
		myfile << ", HAD, AUDPC, Severity\n";
		// Now write results
		for (unsigned int yyyy = 0; yyyy != YearResults.Year.size(); ++yyyy) {
			myfile << YearResults.Year[yyyy];
			for (unsigned int i = 0; i != TOTALGENES; ++i) myfile << ", " << YearResults.geneFreq[yyyy][i];
			myfile << ", " << YearResults.HAD[yyyy];
			myfile << ", " << YearResults.AUDPC[yyyy];
			myfile << ", " << YearResults.Severity[yyyy];
			myfile << "\n";
		}
		// Close file
		myfile.close();
	}

	if (printDailyResults) {
		// And also for the daily file
		myfile.open("Day" + std::to_string(int(processNumber + 0.5)) + fileNameAddition + ".csv");
		// Write a header
		myfile << "Year, DDay, TotalLeafArea, HealthyArea, Severity";
		for (unsigned int iFC = 0; iFC != nFungicides; ++iFC) myfile << ", FC" << iFC;
		for (unsigned int i = 0; i != TOTALGENOTYPES; ++i) myfile << ", L" << i;
		for (unsigned int i = 0; i != TOTALGENOTYPES; ++i) myfile << ", I" << i;
		for (unsigned int i = 0; i != TOTALGENES; ++i) myfile << ", GF" << i;
		myfile << ", TotalLatent, TotalInf, Total";
		myfile << "\n";
		// Now write results
		for (unsigned int counter = 0; counter != DayResults.Year.size(); ++counter) {
			myfile << DayResults.Year[counter] << ", " << DayResults.DegreeDay[counter] << ", ";
			myfile << DayResults.leafAreaIndex[counter] << ", " << DayResults.healthyAreaIndex[counter] << ", "
				<< DayResults.severity[counter];
			for (unsigned int iFC = 0; iFC != nFungicides; ++iFC) myfile << ", " << DayResults.fungicideDose[counter][iFC];
			for (unsigned int i = 0; i != TOTALGENOTYPES; ++i) myfile << ", " << DayResults.pathogenDensity[counter].Latent[i];
			for (unsigned int i = 0; i != TOTALGENOTYPES; ++i) myfile << ", " << DayResults.pathogenDensity[counter].Infectious[i];
			for (unsigned int i = 0; i != TOTALGENES; ++i) myfile << ", " << DayResults.geneFreq[counter][i];
			myfile << ", " << DayResults.totalLatent[counter] << ", " << DayResults.totalInfectious[counter] << ", " << DayResults.totalDensity[counter];
			myfile << "\n";

		}
		// Close file
		myfile.close();
	}

}

void SDayResults::reset()
{
	Year.clear();
	DegreeDay.clear();
	pathogenDensity.clear();
	healthyAreaIndex.clear();
	leafAreaIndex.clear();
	fungicideDose.clear();
	geneFreq.clear();
	totalDensity.clear();
	totalLatent.clear();
	totalInfectious.clear();
	severity.clear();
}
void SYearResults::reset()
{
	Year.clear();
	geneFreq.clear();
	AUDPC.clear();
	HAD.clear();
	Severity.clear();
}
