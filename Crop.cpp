#include "Crop.h"


CCrop::CCrop(){
	
	healthyAreaIndex = 0.0;
	totalAreaIndex = 0.0;
	tuberWeight = 0.0;

}

void CCrop::zero(){

	healthyAreaIndex = 0.0;
	totalAreaIndex = 0.0;
	tuberWeight = 0.0;

}

CCrop CCrop::multiply(double scalar){

	CCrop tempCrop;
	tempCrop.healthyAreaIndex = this->healthyAreaIndex * scalar;
	tempCrop.totalAreaIndex = this->totalAreaIndex * scalar;
	tempCrop.tuberWeight = this->tuberWeight * scalar;

	return tempCrop;

}

CCrop CCrop::add(const CCrop& obj){
	
	CCrop tempCrop;
	
	tempCrop.totalAreaIndex = this->totalAreaIndex + obj.totalAreaIndex;
	tempCrop.healthyAreaIndex = this->healthyAreaIndex + obj.healthyAreaIndex;
	tempCrop.tuberWeight = this->tuberWeight + obj.tuberWeight;

	return tempCrop;

}

bool CCrop::anyNegative(){
	
	bool anyNegative = false;
	if (totalAreaIndex < 0.0) anyNegative = true;
	if (healthyAreaIndex < 0.0) anyNegative = true;
	if (tuberWeight < 0.0) anyNegative = true;

	return anyNegative;
}

