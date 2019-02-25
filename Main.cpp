#include "Main.h"

int main(int argc, char**argv){

	if (argc == 2){
		std::istringstream ss(argv[argc - 1]);
		if (!(ss >> processNumber)) std::cerr << "Invalid input " << argv[argc - 1] << "\n";
		std::cout << "Running simulation " << processNumber << std::endl;
	}
	else processNumber = 000.0;

	// Set all the parameters for this simulation
	setParameters();

	// Timed so that healthy area finishes at 1.3 leaf index. was originally 2750, we updated host growth
	tEOS = 3759;

	// Initialize densities
	initialize();

	// Iterate over the number of years required
	for (unsigned int year = 1; year != 36; ++year){

		// Store the year's initial conditions
		storeResults(0.0, year);

		std::cout << "Year " << year << " commencing." << std::endl;

		// Iterate between the start time and the end of the season
		iterate(0.0, tEOS, year);

		storeYearResults(year);

		resetBwYears();
	}

	writeResultsToFile();

	return 0;
}

void writeResultsToFile(){

	// Create and open a file
	std::ofstream myfile;
	myfile.open("Year.csv");
	// Write a header
	myfile << "Year, TubarWeight";
	for (unsigned int i = 0; i != TOTALGENES; ++i) myfile << ", Gene" << i;
	myfile << "\n";
	// Now write results
	for (unsigned int yyyy = 0; yyyy != YearResults.Year.size(); ++yyyy){
		myfile << YearResults.Year[yyyy] << ", ";
		myfile << YearResults.tubarWeight[yyyy];
		for (unsigned int i = 0; i != TOTALGENES; ++i) myfile << ", " << YearResults.geneFreq[yyyy][i];
		myfile << "\n";
	}
	// Close file
	myfile.close();

	// And also for the daily file
	myfile.open("Day.csv");
	// Write a header
	myfile << "Year, DDay, TotalLeafArea, HealthyArea, TubarWeight, Severity";
	for (unsigned int iFC = 0; iFC != nFungicides; ++iFC) myfile << ", FC" << iFC;
	for (unsigned int i = 0; i != TOTALGENOTYPES; ++i) myfile << ", L" << i;
	for (unsigned int i = 0; i != TOTALGENOTYPES; ++i) myfile << ", I" << i;
	for (unsigned int i = 0; i != TOTALGENES; ++i) myfile << ", GF" << i;
	myfile << ", TotalLatent, TotalInf, Total";
	myfile << "\n";
	// Now write results
	for (unsigned int counter = 0; counter != DayResults.Year.size(); ++counter){
		myfile << DayResults.Year[counter] << ", " << DayResults.DegreeDay[counter] << ", ";
		myfile << DayResults.leafAreaIndex[counter] << ", " << DayResults.healthyAreaIndex[counter] << ", "
			<< DayResults.tubarWeight[counter] << ", " << DayResults.severity[counter];
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

void initialize(){

	oPathogen.zero();
	primaryInocProp.assign(TOTALGENOTYPES, 0.0);
	
	// Put all density as the fully avirulent, fully susceptible individual
	primaryInocProp[0] = 1.0;
	
	oCrop.totalAreaIndex = oCrop.healthyAreaIndex = cropStartingArea;
	vecFungicide.assign(nFungicides, CFungicide());

}

void setParameters(){

	defTimeStep = 0.1;

	// Specify the number of virulence genes the pathogen has (equivalent to the number of resistance genes in the crop plant)
	nVirulenceGenes = 3;
	// Specify the number of fungicide resistant genes the pathogen has
	nResistanceGenes = 1;
	// Total number of genes
	TOTALGENES = nVirulenceGenes + nResistanceGenes;
	// Calculate the total number of genotypes
	TOTALGENOTYPES = pow(3.0, int(nVirulenceGenes + nResistanceGenes));
	CPathogen::TOTALGENOTYPES = TOTALGENOTYPES;



	// Crop starting density - set so that maximum healthy area = 6.044 - same as Femke's. was 0.0000000663
	cropStartingArea = 4.412073e-10;
	aCropParam = 6;//default is 6
	bCropParam = 90; //default is 90
	cCropParam = 2100; //default 2100
	mCropParam = 3850; //default is 3850
	nCropParam = 80; //default is 80

	//tuberGrowthRate = 0.00467448; // redundant, no longer used - replaced by insolation and HAA. 
	// Insolation parameters for tuber production. 
	// these three parameters describe a Gaussian distriubiton where y=a*exp(-.5*((x-x0)/b)^2) where x is the gdd day (tNOW). This describes average insolation in the UK in units of kWh.  
	solA=5; //default is 5
	solB=1115; //default is 1115
	solX=2735; // longest day of year, in GDD, default is 2735. Feel free to change if you are interested in potatoes on Mars.  
	K=0.4; //default is 0.4, this is a constant for potato crops, and represents the capture of light by a canopy (Beers law). See Waggoner and Breger 1987

	// Dictate whether the crop has receptors against each of the virulence genes
	// If it does and the pathogen is avirulent then the IE and LP are reduced
	// If the cultivar has receptor and the pathogen is virulent then the IE and LP are not reduced
	// If the cultivar does not have receptor the virulent have a fitness cost
	cropReceptor.assign(nVirulenceGenes, true);

	// Specify the number of fungicides
	nFungicides = 1;
	// Set up fungicide parameters - can add more sprays for each fungicide as required
	sprayFung1.assign(1, std::pair<unsigned int, double>());
	// N.B. .first is the time unit (degree day) on which to spray; .second is the dose that is sprayed.
	sprayFung1[0].first = 1900; sprayFung1[0].second = 1.0; //default is day is 1900, then every two weeks after. 
	sprayFung1.push_back(std::pair<unsigned int, double>()); sprayFung1[1].first = 2100; sprayFung1[1].second = sprayFung1[0].second;
	sprayFung1.push_back(std::pair<unsigned int, double>()); sprayFung1[2].first = 2300; sprayFung1[2].second = sprayFung1[0].second;
	sprayFung1.push_back(std::pair<unsigned int, double>()); sprayFung1[3].first = 2500; sprayFung1[3].second = sprayFung1[0].second;
	sprayFung1.push_back(std::pair<unsigned int, double>()); sprayFung1[4].first = 2700; sprayFung1[4].second = sprayFung1[0].second;
	sprayFung1.push_back(std::pair<unsigned int, double>()); sprayFung1[5].first = 2900; sprayFung1[5].second = sprayFung1[0].second;
	sprayFung1.push_back(std::pair<unsigned int, double>()); sprayFung1[6].first = 3100; sprayFung1[6].second = sprayFung1[0].second;
	sprayFung1.push_back(std::pair<unsigned int, double>()); sprayFung1[7].first = 3300; sprayFung1[7].second = sprayFung1[0].second;
	sprayFung1.push_back(std::pair<unsigned int, double>()); sprayFung1[8].first = 3500; sprayFung1[8].second = sprayFung1[0].second;
	

	// Decide, for each resistance gene, whether it confers resistance to each fungicide
	confersResistanceToFung.assign(nResistanceGenes, std::vector<bool>(nFungicides, true));

	// The exponent for how fast primary inoculum decays with time, for the beta distribution release - old. these are not used as we now have a truncated normal distribution release curve.  
	//PIb = 0.0000000001; 
	//PIxi = 0.0000000001;
	

	// Default sporulation rate, Rho0. 
	defSporRate = 8; // default is 8

	// survival to next season
	Vtnu = 20;		//default is 20
	Vtt0 =  3559; // 200gdd, 2 weeks, before burnoff of crop. 


	// Natural mortality rates for latent and infectoius lesions
	latentMortalityRate = 0.0; //default is 0
	infectiousMortalityRate = 0.0; //default is 0

	// Lifespans of latent and infectious lesions
	latentLifespan = 50.66542; //default is 50.66542
	infectiousLifespan = 102; //default is 102

	// Default infection efficiency
	defInfEff = 0.01655833; //default 0.01655833

	// Infection efficiency and sporulation rate vectors (for each genotype)
	baseIE.assign(TOTALGENOTYPES, defInfEff);
	baseSR.assign(TOTALGENOTYPES, defSporRate);
	baseLP.assign(TOTALGENOTYPES, latentLifespan);

	infectionEfficiency.assign(TOTALGENOTYPES, defInfEff);
	sporulationRate.assign(TOTALGENOTYPES, defSporRate);
	latentPeriod.assign(TOTALGENOTYPES, latentLifespan);

	// Proportion of spores leaving the field
	propSporesLeavingField = 0; //default is 0

	//fittness costs to fungicide insensitivity, default is 0.002, dom is 0.5
	// Fitness costs of fungicide resistance (proportional reduction in infection efficiency for each RR gene)
	fitnessCostIE = 0.002; 
	// Dominance of fitness cost to infection efficiency - i.e. defines fitness cost for SR
	fitnessCostIEDom = 0.5; //default is 0.5
	// Fitness costs of fungicide resistance (proportional reduction in sporulation rate for each RR gene)
	fitnessCostSR = 0.002;
	// Dominance of fitness cost to sporulation rate - i.e. defines fitness cost for SR
	fitnessCostSRDom = 0.5;
	// Fitness cost to latent period (multiply the latent period by this - make it longer)
	fitnessCostLP = 0.002;
	fitnessCostLPDom = 0.5;

	/// theta, default is 0.05 and 0.5 for dominance
	// Proportional reduction in infection efficiency as a result of being avirulent - i.e. due to host resistance
	AVIRReductionIE = 0.05; 
	// Dominance of the avirulence reduction
	AVIRDomIE = 0.5;
	// Proportional reduction in sporulation rate as a result of being avirulent
	AVIRReductionSR = 0.05;
	// Dominance of the avirulence reduction
	AVIRDomSR = 0.5;
	// Proportional extention of latent period as a result of being avirulent
	AVIRReductionLP = 0.05;
	// Doimnance
	AVIRDomLP = 0.5;

	//fittness costs to virulence, default is 0.002, dom is 0.5
	// Fitness cost in infection efficiency as a result of being virulent (should be less than result of host resistance above)
	fCReductionIE = 0.002;
	// Dominance of the virulent fitness cost reduction
	fCDomIE = 0.5;
	// Fitness cost in sporulation rate as a result of being virulent (should be less than result of host resistance above)
	fCReductionSR = 0.002;
	// Dominance of the virulent fitness cost reduction
	fCDomSR = 0.5;
	// Proportional *rise* in latent period as a result of having a fitness cost (should be less than result of host resistance above)
	// i.e. fitness cost = 0.2 -> latent period = latent period * 1.2; fitness cost = 1.0 -> latent period = latent period * 2.0
	fCReductionLP = 0.002;
	// Doimnance
	fCDomLP = 0.5;

	// Set up genotype to diploid array
	genotypeToDiploid1D.assign(TOTALGENOTYPES * TOTALGENES, 0);
	for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene){
		for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){
			unsigned int genotype = iGeno;
			for (unsigned int jGene = iGene + 1; jGene != TOTALGENES; ++jGene) genotype = div(genotype, 3).quot;
			genotype = div(genotype, 3).rem;
			genotypeToDiploid1D[iGeno + iGene * TOTALGENOTYPES] = genotype;
		}
	}

	// Mutation rate of a single allele. Normally either 0.000001 or 0. 
	mutationRateFungR = 0.000001;
	mutationRateVir =   0.0;//00001;

	// Fungicide ~ infection efficiency parameters
	alphaMax.assign(nFungicides,1.0); // default is 1.0 
	kappa.assign(nFungicides,5); // from LINK data, default is 5
	
	fungResPi = 0.2; // default is 0.2. means 1 dose reduces IE, sporulation or lp to 20% its default value.
	fungResDom = 0.5; //default is 0.5, based on literature. 

	// Decay rate of the fungicide, 6 Jdays. 
	fungDecayRate = 0.007967209; // default is 0.007967209, 6 julian days.



	createMutationMatrices();

	// Functions that adjust the latent period, infection efficiency and sporulation rate for each genotype depending on the cultivar resistance
	// N.B. Needs to happen after genoToDiploid1D has been created - we use genoToDiploid1D in the equation
	calculateBaseLP();
	calculateBaseIE();
	calculateBaseSR();

}

void iterate(double tStart, double tEnd, unsigned int year){

	double tNow = tStart;
	double tPrev = 0.0;

	// This is the default time step - if things go negative we might have to reduce the size of the time step
	double tempTimeStep = defTimeStep;

	while (floor(tNow + 0.5*tempTimeStep) < floor(tEnd + 0.5*tempTimeStep)){

		// Check if a spray needs to be applied
		checkSprayTimes(tNow, tempTimeStep);

		while (!rungeKutta(tNow, tempTimeStep)){
			tempTimeStep = tempTimeStep / 10.0;
			std::cout << "Trying a smaller timestep: " << tempTimeStep << "\n";
			if (tempTimeStep < 1E-20){
				std::cerr << "Time step is too small: " << tempTimeStep << ". Please do something." << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		tPrev = tNow;
		tNow += tempTimeStep;

		// If a day has passed, record results
		if (floor(tNow + 0.5*tempTimeStep) > floor(tPrev + 0.5*tempTimeStep)){
			storeResults(floor(tNow + 0.5*tempTimeStep), year);
			tempTimeStep = defTimeStep;
		}

	}

}

void checkSprayTimes(double time, double timeStep){

	// Loop over the first fungicide
	if (sprayFung1.size() > 0){
		for (unsigned int iSpray = 0; iSpray != sprayFung1.size(); ++iSpray){

			if (floor(time - 0.5*timeStep) < floor(time + 0.5*timeStep) && floor(time) == sprayFung1[iSpray].first){
				if(sprayFung1[iSpray].second > 0.0) std::cout << "Spraying fungicide 1 at dose " << sprayFung1[iSpray].second << "!" << std::endl;
				vecFungicide[0].currentDose += sprayFung1[iSpray].second;
			}

		}
	}

	// Loop over the second fungicide
	if (sprayFung2.size() > 0){
		for (unsigned int iSpray = 0; iSpray != sprayFung2.size(); ++iSpray){

			if (floor(time - 0.5*timeStep) < floor(time) && floor(time) == sprayFung2[iSpray].first){
				if (sprayFung2[iSpray].second > 0.0) std::cout << "Spraying fungicide 2!" << std::endl;
				vecFungicide[1].currentDose += sprayFung2[iSpray].second;
			}

		}
	}

}

bool rungeKutta(double tNow, double timeStep){

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
	for (std::size_t iFun = 0; iFun != vecFungicide.size(); ++iFun){
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
	for (std::size_t iFun = 0; iFun != vecFungicide.size(); ++iFun){
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
	for (std::size_t iFun = 0; iFun != vecFungicide.size(); ++iFun){
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
	for (size_t iFC = 0; iFC != vecFungicide.size(); ++iFC){
		FU_DV[iFC].currentDose = FU_IV[iFC].currentDose + (timeStep / 6.0) * (FU_k1[iFC].currentDose + 2.0 * FU_k2[iFC].currentDose
			+ 2.0 * FU_k3[iFC].currentDose + FU_k4[iFC].currentDose);
		if (FU_DV[iFC].currentDose < 0.0) noNeg = false;
	}

	// If none of the final variables are negative then set them
	if (noNeg){
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

void deriv(double tNow, CCrop HA_IN, CCrop& HA_DV, const CPathogen &PA_IN, CPathogen& PA_DV,
	const std::vector<CFungicide> &FU_IN, std::vector<CFungicide> &FU_DV){

	// *** Crop stuff
	// Work out total crop area
	double totalCropArea = HA_IN.healthyAreaIndex;
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){
		totalCropArea += PA_IN.Latent[iGeno];
		totalCropArea += PA_IN.Infectious[iGeno];
	}

		// Work out growth rate and senesence rate of leaves
	double LAI = (aCropParam / (1 + exp(-(tNow - cCropParam) / bCropParam)));
	double SAI = (aCropParam / (1 + exp(-(tNow - mCropParam)/nCropParam)));
	double growthRate = (aCropParam * exp(-(tNow - cCropParam) / bCropParam) / ((bCropParam * (1 + exp(-(tNow - cCropParam) / bCropParam)))* (1 + exp(-(tNow - cCropParam) / bCropParam)))) / (LAI - SAI);
	double senesRate = (aCropParam * exp(-(tNow - mCropParam) / nCropParam) / (nCropParam * (1 + exp(-(tNow - mCropParam) / nCropParam)) * (1 + exp(-(tNow - mCropParam) / nCropParam)))) / (LAI - SAI);
	HA_DV.totalAreaIndex += growthRate * HA_IN.totalAreaIndex;
	HA_DV.healthyAreaIndex += growthRate * HA_IN.totalAreaIndex;
	HA_DV.healthyAreaIndex -= senesRate * HA_IN.healthyAreaIndex;

	// Work out insolation at tNow (replicates gaussian distribution of annual average sunlight). 
	double insolation = 3.6*solA*exp(-0.5*pow(((tNow-solX)/solB),2));

	// Work out photosynthetic crop area (H+L) at tNow
	double photosyntheticArea = 0.0;
	photosyntheticArea += HA_IN.healthyAreaIndex;
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){
		photosyntheticArea += PA_IN.Latent[iGeno];
	}
	
	// work out Healty Area Absorbance(HAA) at tNow. This is the light absorbed by the healthy area 
	double HAA = 0.0;
	HAA = insolation*(photosyntheticArea/totalCropArea)*(1-exp(totalCropArea*K));

	// If tNow > 2247, then tubers are growing.   
	if (tNow > 2247){ 
		HA_DV.tuberWeight = HAA*-0.0001837431; // new model
		//HA_DV.tuberWeight += tuberGrowthRate * HA_IN.healthyAreaIndex; // old model
	}


	// Calculate the infection efficiency of each genotype at t = tNow
	calculateInfectionEfficiency(FU_IN);
	// Re-calculate sporulation rate based on current fungicide dose
	calculateSporulationRate(FU_IN);
	// Calculate the latent period
	calculateLatentPeriod(FU_IN);

	// Calculate the number of spores of each genotype - function of spores produced and mutation rate
	std::vector<double> secondarySpores(TOTALGENOTYPES, 0.0);
	mutateOrDontIDontCare(PA_IN, secondarySpores);

	// For each fungicide, make the dose decay
	for (size_t i = 0; i != FU_DV.size(); ++i){
		FU_DV[i].currentDose -= fungDecayRate * FU_IN[i].currentDose;
	}

	// Calculate the total amount of primary inoculum at tNow 
	//double primaryInoculum = PIxi * tNow * tNow * exp(-PIb * tNow); // beta distribution for release, old.
	double primaryInoculum= 0.01*exp(-0.5*pow(((tNow-2100)/50),2));// normal distribution for release, new. default is 0.01*exp(-0.5*pow(((tNow-2100)/50),2))
	
	double scale = 0.01;  // rounds primaryInoculum, i.e. round to nearest one-hundreth
	primaryInoculum = (int)(primaryInoculum / scale) * scale;

	// Calculate the deposition probability
	double depoProb = (1-propSporesLeavingField) * (1 - exp(-HA_IN.totalAreaIndex));
	
	// Loop through the genotypes and calculate change in the density of lesions
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){

		// Store a proportion of spores for the next season // TODO: Make these parameters into global parameters
		PA_DV.sporesEnteringOWPool[iGeno] += secondarySpores[iGeno] / (1.0 + exp((-(tNow - Vtt0))/Vtnu));
															
		// New latent from primary inoculum // TODO: What if totalCropArea = 0.0;
		double newPrimLatentThisGeno = depoProb * primaryInoculum * primaryInocProp[iGeno] * HA_IN.healthyAreaIndex / HA_IN.totalAreaIndex * infectionEfficiency[iGeno];
		PA_DV.Latent[iGeno] += newPrimLatentThisGeno;
		// Remove this area from the healthy area
		HA_DV.healthyAreaIndex -= newPrimLatentThisGeno;

		// New latent lesions from secondary infection
		double newSecondLatentThisGeno = depoProb * HA_IN.healthyAreaIndex / HA_IN.totalAreaIndex * infectionEfficiency[iGeno] * secondarySpores[iGeno];
		PA_DV.Latent[iGeno] += newSecondLatentThisGeno;
		HA_DV.healthyAreaIndex -= newSecondLatentThisGeno;

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

void mutateOrDontIDontCare(const CPathogen& inPathogen, std::vector<double>& outSpores){
	// NB: outSpores is a vector of zeros - probably ought to check this but meh.
	
	std::vector<double> infectiousLesions = inPathogen.Infectious;

	// Loop through each genotype for the parent
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){

		// Need to work out the proportion of spores produced by the parent genotype that goes to each offspring genotype
		// Create a vector that will store the proportion of spores of each genotype
		std::vector<double> propSpores(TOTALGENOTYPES, 0.0);

		// Loop through each potential offspring genotype
		for (unsigned int jGeno = 0; jGeno != TOTALGENOTYPES; ++jGeno){

			// Store the proportion of iGeno that mutates (or doesn't) into jGeno
			double proportionjGeno = 1.0;

			if (mutationMatrix.size() > 0){
				proportionjGeno = mutationMatrix[iGeno][jGeno];
			}
			else {
				// Loop through each gene
				for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene){

					unsigned int parentGenotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];
					unsigned int offspringGenotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + jGeno];

					// If we're looking at a virulence gene, then use virulence gene mutation array
					if (iGene < nVirulenceGenes) proportionjGeno *= singleVirGeneMutationArray[parentGenotype][offspringGenotype];
					else proportionjGeno *= singleFungResGeneMutationArray[parentGenotype][offspringGenotype];
				}
			}

			outSpores[jGeno] += infectiousLesions[iGeno] * sporulationRate[iGeno] * proportionjGeno;

		}

	}

}

void calculateInfectionEfficiency(const std::vector<CFungicide>& fungicides){
	// Set sporulation rate array back to default for every genotype
	infectionEfficiency = baseIE;

	// Loop through each genotype, and for each resistance genotype reduce the sporulation rate accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){

		double proportionReduction = 1.0;

		// Loop through each fungicide
		for (size_t iFung = 0; iFung != fungicides.size(); ++iFung){

	// This stores the part: (phi^nHomo * (1-(1-phi)*dom)^nHetero)
	double coefficient = 1.0;
	// The following genes are fungicide resistance genes
	for (unsigned int iGene = nVirulenceGenes; iGene != TOTALGENES; ++iGene){
	// Work out the diploid for this gene
	unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];
	switch (genotype){
	// 0 = SS - multiply by a full reduction in the sporulation rate
	case 0:
	coefficient *= fungResPi;
	break;
	// 1 = SR - partial reduction in sporulation rate by fungicide
	case 1:
	coefficient *= (1 - (1 - fungResPi) * fungResDom);
	break;
	// 2 = RR - no reduction in the sporulation rate
	case 2:
	break;
	}}
	proportionReduction *= 1 - alphaMax[iFung] * (1 - exp(-kappa[iFung] * fungicides[iFung].currentDose)) * (1 - coefficient);
	}infectionEfficiency[iGeno] *= proportionReduction;
}}



void calculateSporulationRate(const std::vector<CFungicide>& fungicides){

	// Set sporulation rate array back to default for every genotype
	sporulationRate = baseSR;

	// Loop through each genotype, and for each resistance genotype reduce the sporulation rate accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){

		double proportionReduction = 1.0;

		// Loop through each fungicide
		for (size_t iFung = 0; iFung != fungicides.size(); ++iFung){
			
// This stores the part: (phi^nHomo * (1-(1-phi)*dom)^nHetero)
double coefficient = 1.0;
// The following genes are fungicide resistance genes
for (unsigned int iGene = nVirulenceGenes; iGene != TOTALGENES; ++iGene){
// Work out the diploid for this gene
unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];
switch (genotype){
// 0 = SS - multiply by a full reduction in the sporulation rate
case 0:
coefficient *= fungResPi;
break;
// 1 = SR - partial reduction in sporulation rate by fungicide
case 1:
coefficient *= (1 - (1 - fungResPi) * fungResDom);
break;
// 2 = RR - no reduction in the sporulation rate
case 2:
break;
}}
proportionReduction *= 1 - alphaMax[iFung] * (1 - exp(-kappa[iFung] * fungicides[iFung].currentDose)) * (1 - coefficient);
}sporulationRate[iGeno] *= proportionReduction;
}}


		


void calculateLatentPeriod(const std::vector<CFungicide>& fungicides){

	// Set sporulation rate array back to default for every genotype
	latentPeriod = baseLP;

	// Loop through each genotype, and for each resistance genotype extend the latent period accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){

		double proportionExtension = 1.0;

		// Loop through each fungicide
		for (size_t iFung = 0; iFung != fungicides.size(); ++iFung){

// This stores the part: (phi^nHomo * (1-(1-phi)*dom)^nHetero)
double coefficient = 1.0;
// The following genes are fungicide resistance genes
for (unsigned int iGene = nVirulenceGenes; iGene != TOTALGENES; ++iGene){
// Work out the diploid for this gene
unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];
switch (genotype){
// 0 = SS - multiply by a full reduction in the sporulation rate
case 0:
coefficient *= fungResPi;
break;
// 1 = SR - partial reduction in sporulation rate by fungicide
case 1:
coefficient *= (1 - (1 - fungResPi) * fungResDom);
break;
// 2 = RR - no reduction in the sporulation rate
case 2:
break;


}}
proportionExtension *= 1 - alphaMax[iFung] * (1 - exp(-kappa[iFung] * fungicides[iFung].currentDose)) * (1 - coefficient);
}latentPeriod[iGeno] *= (1+proportionExtension);
}}


void calculateBaseIE(){

	baseIE.assign(TOTALGENOTYPES, defInfEff);

	// Loop through each genotype, and for each resistance genotype reduce the infection efficiency accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){

		// Calculate the fitness cost for this genotype
		double fitnessCost = 1.0;

		for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene){

			// Work out the diploid for this gene
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			// If a virulence gene, then account for reduction in infection efficiency for avirulent alleles
			if (iGene < nVirulenceGenes){

				// Store the proportional reduction for this gene
				double proportionReduction = 1.0;

				switch (genotype){
					// 0 = AA - avirulent; infection efficiency is reduced
				case 0:
				if (cropReceptor[iGene]) proportionReduction *= (1 - AVIRReductionIE);
					else proportionReduction *= 1.0;
					break;
					// 1 = Aa - partially virulent
				case 1:
				if (cropReceptor[iGene]) proportionReduction *= pow((1-fCReductionIE), fCDomSR)* (1-AVIRDomIE * AVIRReductionIE); // before 01/02/2016 this was (1 - fCReductionIE * fCDomSR);
					else proportionReduction *= pow((1-fCReductionIE), fCDomSR);
					break;
					// 2 = aa - virulent; use default infection efficiency
				case 2:
					if (cropReceptor[iGene]) proportionReduction *= (1-fCReductionIE);
					else proportionReduction *= (1 - fCReductionIE);
					break;

				}

				baseIE[iGeno] *= proportionReduction;
			}
			// Otherwise account for fitness cost of having fungicide resistance
			else {

				switch (genotype){
					// 0 = SS - no fitness cost of having fungicide resistance
				case 0:
					break;
					// 1 = SR - partial reduction in infection efficiency due to fitness cost
				case 1:
					fitnessCost *= pow((1 - fitnessCostIEDom), fitnessCostIE); //before 01/02/2016 this was (1 - fitnessCostIEDom * fitnessCostIE);
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

	printIE();

}

void calculateBaseSR(){

	baseSR.assign(TOTALGENOTYPES, defSporRate);

	// Loop through each genotype, and for each resistance genotype reduce the infection efficiency accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){

		// Calculate the fitness cost for this genotype
		double fitnessCost = 1.0;

		for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene){

			// Work out the diploid for this gene
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			// If a virulence gene, then account for reduction in infection efficiency for avirulent alleles
			if (iGene < nVirulenceGenes){

				// Store the proportional reduction for this gene
				double proportionReduction = 1.0;

				switch (genotype){
					// 0 = AA - avirulent; sporulation rate is reduced
				case 0:
					if (cropReceptor[iGene]) proportionReduction *= (1 - AVIRReductionSR);
					else proportionReduction *= 1.0;
					break;
					// 1 = Aa - partially virulent
				case 1:
					if (cropReceptor[iGene]) proportionReduction *= pow((1-fCReductionSR), fCDomSR)* (1 - AVIRDomSR * AVIRReductionSR); // before 01/02/2016 this was (1 - fCReductionSR * fCDomSR) * (1 - AVIRDomSR * AVIRReductionSR);
					else proportionReduction *= pow((1-fCReductionSR), fCDomSR); //before 01/02/2016 this was (1 - fCReductionSR * fCDomSR);
					break;
					// 2 = aa - fully virulent; default sporulation rate
				case 2:
					if (cropReceptor[iGene]) proportionReduction *= (1 - fCReductionSR);
					else proportionReduction *= (1 - fCReductionSR); 
					break;
				}

				baseSR[iGeno] *= proportionReduction;
			}
			// Otherwise account for fitness cost of having fungicide resistance
			else {

				switch (genotype){
					// 0 = SS - no fitness cost of having fungicide resistance
				case 0:
					break;
					// 1 = SR - partial reduction in infection efficiency due to fitness cost
				case 1:
					fitnessCost *= pow((1 - fitnessCostSRDom),fitnessCostSR); // before 01/02/2016 this was (1 - fitnessCostSRDom * fitnessCostSR)
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

	printSR();

}

void calculateBaseLP(){

	baseLP.assign(TOTALGENOTYPES, latentLifespan);

	// Loop through each genotype, and for each resistance genotype reduce the infection efficiency accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){

		// Calculate the fitness cost for this genotype
		double fitnessCost = 1.0;

		for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene){

			// Work out the diploid for this gene
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			double proportionalExtension = 1.0;

			// If virulence gene, then extend the latent period if it's avirulent
			if (iGene < nVirulenceGenes){

				switch (genotype){
					// 0 = AA - avirulent; infection efficiency is reduced
				case 0:
					if (cropReceptor[iGene]) proportionalExtension *= (1 + AVIRReductionLP);
					else proportionalExtension *= 1.0;
					break;
					// 1 = Aa - partially virulent 
				case 1:
					if (cropReceptor[iGene]) proportionalExtension *= (1 + AVIRDomLP * AVIRReductionLP) * pow((1 + fCReductionLP), fCDomLP); // before 01/02/2016 this was (1 + AVIRDomLP * AVIRReductionLP) * (1 + fCReductionLP * fCDomLP);
					else proportionalExtension *= pow((1 + fCReductionLP), fCDomLP); // before 01/02/2016 this was (1 + AVIRDomLP * AVIRReductionLP) * (1 + fCReductionLP * fCDomLP);
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

				switch (genotype){
					// 0 = SS - no fitness cost
				case 0:
					break;
					// 1 = SR - heterozygote fitness cost
				case 1:
					fitnessCost *= pow((1 + fitnessCostLPDom), fitnessCostLP); //before 01/02/2016 this was (1 + fitnessCostLPDom * fitnessCostLP);
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

	printLP();

}

void createMutationMatrices(){

	// Create a mutation matrix 3x3, for the probability of mutating from one genotype to another
	singleFungResGeneMutationArray.assign(3, std::vector<double>(3, 1.0));
	// I'm just going to write this out for now, can program it if necessary - it's not hard
	singleFungResGeneMutationArray[0][0] = 1.0 - 2.0*mutationRateFungR + mutationRateFungR * mutationRateFungR;
	singleFungResGeneMutationArray[0][1] = 2.0 * mutationRateFungR - 2.0 * mutationRateFungR * mutationRateFungR;
	singleFungResGeneMutationArray[0][2] = mutationRateFungR * mutationRateFungR;
	singleFungResGeneMutationArray[1][0] = -mutationRateFungR * mutationRateFungR + mutationRateFungR;
	singleFungResGeneMutationArray[1][1] = 1.0 + 2.0 * mutationRateFungR * mutationRateFungR - 2.0 * mutationRateFungR;
	singleFungResGeneMutationArray[1][2] = -mutationRateFungR * mutationRateFungR + mutationRateFungR;
	singleFungResGeneMutationArray[2][0] = mutationRateFungR * mutationRateFungR;
	singleFungResGeneMutationArray[2][1] = 2.0 *  mutationRateFungR - 2.0 * mutationRateFungR * mutationRateFungR;
	singleFungResGeneMutationArray[2][2] = 1.0 - 2.0 * mutationRateFungR + mutationRateFungR * mutationRateFungR;

	// Create a mutation matrix 3x3, for the probability of mutating from one genotype to another
	singleVirGeneMutationArray.assign(3, std::vector<double>(3, 1.0));
	// I'm just going to write this out for now, can program it if necessary - it's not hard
	singleVirGeneMutationArray[0][0] = 1.0 - 2.0*mutationRateVir + mutationRateVir * mutationRateVir;
	singleVirGeneMutationArray[0][1] = 2.0 * mutationRateVir - 2.0 * mutationRateVir * mutationRateVir;
	singleVirGeneMutationArray[0][2] = mutationRateVir * mutationRateVir;
	singleVirGeneMutationArray[1][0] = -mutationRateVir * mutationRateVir + mutationRateVir;
	singleVirGeneMutationArray[1][1] = 1.0 + 2.0 * mutationRateVir * mutationRateVir - 2.0 * mutationRateVir;
	singleVirGeneMutationArray[1][2] = -mutationRateVir * mutationRateVir + mutationRateVir;
	singleVirGeneMutationArray[2][0] = mutationRateVir * mutationRateVir;
	singleVirGeneMutationArray[2][1] = 2.0 *  mutationRateVir - 2.0 * mutationRateVir * mutationRateVir;
	singleVirGeneMutationArray[2][2] = 1.0 - 2.0 * mutationRateVir + mutationRateVir * mutationRateVir;

	// If there's enough space, create a mutation matrix for each genotype
	std::vector<double> test;
	std::cout << "Vector maximum size is: " << test.max_size() << std::endl;
	std::cout << "Asking for total size to be: " << pow(3.0, double(TOTALGENES) * 2) << std::endl;
	if (test.max_size() > pow(3.0, double(TOTALGENES) * 2)){
		for (unsigned int parentGeno = 0; parentGeno != TOTALGENOTYPES; ++parentGeno){
			mutationMatrix.push_back(std::vector<double>(TOTALGENOTYPES, 1.0));
			for (unsigned int offspringGeno = 0; offspringGeno != TOTALGENOTYPES; ++offspringGeno){

				// Store the proportion of parentGeno that mutates (or doesn't) into offspringGeno
				double proportionOffspringGeno = 1.0;

				// Loop through each gene
				for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene){

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

void storeYearResults(unsigned int year){

	YearResults.Year.push_back(year);
	YearResults.tubarWeight.push_back(oCrop.tuberWeight);
	// Work out the frequency of each gene
	std::vector<double> geneFreq(TOTALGENES, 0.0);
	// Loop through each genotype, and add to each gene if heterozygote or resistant
	double totalDensity = 0.0;
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){
		totalDensity += oPathogen.Infectious[iGeno];

		for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene){
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			switch (genotype){
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

}

void storeResults(double tNow, unsigned int year){

	DayResults.DegreeDay.push_back(floor(tNow));
	DayResults.Year.push_back(year);
	DayResults.pathogenDensity.push_back(oPathogen);
	DayResults.healthyAreaIndex.push_back(oCrop.healthyAreaIndex);
	DayResults.leafAreaIndex.push_back(oCrop.totalAreaIndex);
	DayResults.tubarWeight.push_back(oCrop.tuberWeight);
	DayResults.severity.push_back(calcSeverity());
	std::vector<double> tempDose;
	for (unsigned int iFC = 0; iFC != nFungicides; ++iFC) tempDose.push_back(vecFungicide[iFC].currentDose);
	DayResults.fungicideDose.push_back(tempDose);
	// Work out the frequency of each gene
	std::vector<double> geneFreq(TOTALGENES, 0.0);
	// Loop through each genotype, and add to each gene if heterozygote or resistant
	double totalDensity = 0.0, totalLatent = 0.0, totalInfectious = 0.0;
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){
		totalDensity += oPathogen.Infectious[iGeno];
		totalInfectious += oPathogen.Infectious[iGeno];
		totalLatent += oPathogen.Latent[iGeno];
		totalDensity += oPathogen.Latent[iGeno];

		for (unsigned int iGene = 0; iGene != TOTALGENES; ++iGene){
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			switch (genotype){
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

void resetBwYears(){

	double totalInfectiousDensity = 0.0;

	// Calculate the proportion of each genotype, and store in primaryInocProp
	// Sum the total amount of spores in the overwinter pool
	double totalSpores = 0.0;
	for (std::vector<double>::const_iterator cVecIt = oPathogen.sporesEnteringOWPool.begin(); cVecIt != oPathogen.sporesEnteringOWPool.end(); ++cVecIt){
		totalSpores += *cVecIt;
	}
	if (totalSpores > 0.0){
		// Make the proportion = to the sporesEnteringOWPool
		primaryInocProp = oPathogen.sporesEnteringOWPool;
		// Divide each genotype by the total Spores to make a proportion
		for (std::vector<double>::iterator vecIt = primaryInocProp.begin(); vecIt != primaryInocProp.end(); ++vecIt) (*vecIt) /= totalSpores;
	}
	else {
		primaryInocProp.assign(TOTALGENOTYPES, 0.0);
	}

	//primaryInocProp = oPathogen.Infectious;
	//double totalInoc = 0.0;
	//for (std::vector<double>::const_iterator cVecIt = oPathogen.Infectious.begin(); cVecIt != oPathogen.Infectious.end(); ++cVecIt){
	//	totalInoc += *cVecIt;
	//}
	//for (std::vector<double>::iterator vecIt = primaryInocProp.begin(); vecIt != primaryInocProp.end(); ++vecIt) *vecIt /= totalInoc;

	// Now rezero the pathogen, and crop
	oPathogen.zero();
	oCrop.zero();
	oCrop.healthyAreaIndex = oCrop.totalAreaIndex = cropStartingArea;

}

void printIE(){

	// Open a file
	std::ofstream myFile("IE.csv", std::ios::trunc);

	// Write a header
	myFile << "Genotype, IE" << std::endl;
	for (std::vector<double>::const_iterator cVecIt = baseIE.begin(); cVecIt != baseIE.end(); ++cVecIt){

		myFile << cVecIt - baseIE.begin() << ", " << *cVecIt << "\n";

	}
	myFile.close();

}

void printLP(){

	// Open a file
	std::ofstream myFile("LP.csv", std::ios::trunc);

	// Write a header
	myFile << "Genotype, LP" << std::endl;
	for (std::vector<double>::const_iterator cVecIt = baseLP.begin(); cVecIt != baseLP.end(); ++cVecIt){

		myFile << cVecIt - baseLP.begin() << ", " << *cVecIt << "\n";

	}
	myFile.close();

}

void printSR(){

	// Open a file
	std::ofstream myFile("SR.csv", std::ios::trunc);

	// Write a header
	myFile << "Genotype, SR" << std::endl;
	for (std::vector<double>::const_iterator cVecIt = baseSR.begin(); cVecIt != baseSR.end(); ++cVecIt){

		myFile << cVecIt - baseSR.begin() << ", " << *cVecIt << "\n";

	}
	myFile.close();

}

double calcSeverity(){

	double severity;

	double infected = 0.0;
	for (std::vector<double>::const_iterator cVecIt = oPathogen.Infectious.begin(); cVecIt != oPathogen.Infectious.end(); ++cVecIt){
		infected += *cVecIt;
	}
	for (std::vector<double>::const_iterator cVecIt = oPathogen.Latent.begin(); cVecIt != oPathogen.Latent.end(); ++cVecIt){
		infected += *cVecIt;
	}

	if (oCrop.totalAreaIndex <= 0.0){
		std::cerr << "Trying to calculate severity when totalAreaIndex is zero. Returning zero." << std::endl;
	}
	else { severity = infected / oCrop.totalAreaIndex; }

	return severity;

}

void oldCalculateSporulationRate(const std::vector<CFungicide>& fungicides){

	// This was how I modelled the change in the sporulation rate (and infection efficiency), but does lead to converging dose response curves

	// Set sporulation rate array back to default for every genotype
	sporulationRate = baseSR;

	// Loop through each genotype, and for each resistance genotype reduce the sporulation rate accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){

		// The following genes are fungicide resistance genes
		for (unsigned int iGene = nVirulenceGenes; iGene != TOTALGENES; ++iGene){

			// Work out the diploid for this gene
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			// Store the proportional reduction for this gene
			double proportionReduction = 1.0;

			switch (genotype){
				// 0 = SS - full reduction in sporulation rate by fungicide
			case 0:
				for (size_t iFung = 0; iFung != fungicides.size(); ++iFung){
					// This was the old way based on Peter's paper
					if (confersResistanceToFung[iGene - nVirulenceGenes][iFung]) proportionReduction *= (1 - alpha[iFung] * (1 - exp(-beta[iFung] * fungicides[iFung].currentDose)));
				}
				break;
				// 1 = SR - partial reduction in sporulation rate by fungicide - this doesn't seem to have dominance in...
			case 1:
				for (size_t iFung = 0; iFung != fungicides.size(); ++iFung){
					// This was the old way based on Peter's paper
					if (confersResistanceToFung[iGene - nVirulenceGenes][iFung]) proportionReduction *= (1 - alpha[iFung] * (1 - exp(-beta[iFung] * fungicides[iFung].currentDose)));
				}
				break;
				// 2 = RR - no reduction in the sporulation rate
			case 2:
				break;
			}

			sporulationRate[iGeno] *= proportionReduction;

		}

	}

}

void oldCalculateInfectionEfficiency(const std::vector<CFungicide>& fungicides){

	// Set infection efficiency array back to default for every genotype
	infectionEfficiency = baseIE;

	// Loop through each genotype, and for each resistance genotype reduce the infection efficiency accordingly
	for (unsigned int iGeno = 0; iGeno != TOTALGENOTYPES; ++iGeno){

		for (unsigned int iGene = nVirulenceGenes; iGene != TOTALGENES; ++iGene){

			// Work out the diploid for this gene
			unsigned int genotype = genotypeToDiploid1D[iGene * TOTALGENOTYPES + iGeno];

			// Store the proportional reduction for this gene
			double proportionReduction = 1.0;

			switch (genotype){
				// 0 = SS - full reduction in infection efficiency by fungicide
			case 0:
				for (size_t iFung = 0; iFung != fungicides.size(); ++iFung){
					if (confersResistanceToFung[iGene - nVirulenceGenes][iFung]) proportionReduction *= (1 - alpha[iFung] * (1 - exp(-beta[iFung] * fungicides[iFung].currentDose)));
				}
				break;
				// 1 = SR - partial reduction in infection efficiency by fungicide
			case 1:
				for (size_t iFung = 0; iFung != fungicides.size(); ++iFung){
					if (confersResistanceToFung[iGene - nVirulenceGenes][iFung]) proportionReduction *= (1 - alpha[iFung] * (1 - exp(-beta[iFung] * fungicides[iFung].currentDose)));
				}
				break;
				// 2 = RR - no reduction in the infection efficiency
			case 2:
				break;
			}

			infectionEfficiency[iGeno] *= proportionReduction;

		}

	}

}

