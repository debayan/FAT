#pragma once

#include "InterHierarchy.h"
#include "CommonSubtree.h"
#include "Eigen/Dense"
#include "SceneHierarchy.h"
#include "InterSet.h"
#include "Scene.h"
#include "SymGroup.h"

class DistMeasure
{
public:
	DistMeasure(){};
	DistMeasure(DistParameter para){distPara = para;}
	~DistMeasure(){};

public:
	double betweenPoses(Eigen::Matrix3Xd X, Eigen::Matrix3Xd Y);
	double betweenInteractions(Interaction* iter1, Interaction* iter2, bool insideScene = false);
	CommonSubtree betweenInterHierarchies(InterHierarchy* hier1, InterHierarchy* hier2, bool visualize, bool updateMatching);
	double betweenObjLevelFeature(ObjLevelFeature f1, ObjLevelFeature f2, int level);
	double betweenInterSets(InterSet set1, InterSet set2);
	double betweenGeoFeatures(ObjGeoFeature f1, ObjGeoFeature f2);
	Eigen::MatrixXd betweenGeoFeatures(QVector<ObjGeoFeature> features);
	double betweenSymGroups(SymGroup sg1, SymGroup sg2);

private:
	double computeIbsGeomDistance(IBS* ibs1, IBS* ibs2);
	int computeIbsTopoDistance(IBS* ibs1, IBS* ibs2);
	double computeRegionDistance(FuncRegion* r1, FuncRegion* r2);

	int symLevel(Interaction* iter1, Interaction* iter2); // detect the symmetry of interactions inside the same scene: 0: no symmetry, 1: outer symmetry, 2: inner symmetry, 3: pattern

	DistParameter distPara;
};