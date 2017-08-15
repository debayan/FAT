
#pragma once

#include "Object.h"
#include "IBS.h"
#include "FuncRegion.h"

class SymGroup;

enum INTERACTION_TYPE{EMPTY, SINGLE, COMBINED_INIT, COMBINED_FEATURE, COMBINED_FULL};

class Interaction
{
public:
	Interaction(Scene* s, Object* object);
	Interaction(Scene* s, QVector<Object*> objects);
	Interaction(Scene* s, QVector<Interaction*> interactions, INTERACTION_TYPE type = COMBINED_FULL);
	Interaction(Scene* s, QVector<int> orgObjIdx);
	~Interaction();

public:
	void computeFeature();  // compute ibs and surface region features
	QVector<double> getFeature(int id); // 0-5: ibs_pfh, ibs->dirHist, ibs->distHist, region->pfh, region->dirHist, region->heightHist

private:
	void initialize(); 
	void construct();  // construct the interaction by compute IBS and region and their features
	void computeIBS();
	void findSurfaceRegion();

private:
	void computePairwiseIBS();
	void refinePairwiseIBS();
	void detectCollision();

	void computeIbsFeature();
	void computeRegionFeature();

	void storeOrigObj(QVector<Interaction*> interactions);  // store the index of combined original objects
	void generateCombinedObj(QVector<Interaction*> interactions);  // generate a new object which is the combination of all the objects in interactions
	void computeAverageFeatures(QVector<Interaction*> interactions); // compute the average feature among all the interactions
	void constructNewInteraction(QVector<Interaction*> interactions); // compute new combined object and its corresponding interaction

public:
	Scene* scene;

	Object* obj; // the interacting object, can be a combination of multiple objects
	IBS * ibs;
	FuncRegion* region;	
	
	// for multiple object, check whether there is some kind of pattern
	double mergeDist;
	bool isChecked;
	INTERACTION_TYPE type;
	int nodeIdx;
	int idx;

	// store the symmetry information: checked whether it covers any symmetry group and the corresponding interacting objects in
	QVector<SymGroup*> symGroups;
	QVector< QVector<int> > symObjs;

private:
	bool isCollisionFree;
	QVector<QVector<int>>  collisionFaceIdx; // the collision face idx for the interacting object and the central object
};