#pragma once
#include "Eigen/Dense"
#include "dataanalysis.h"

class HierarchicalCluster
{
public:
	HierarchicalCluster();
	~HierarchicalCluster();
   
	alglib::ahcreport* clustering(Eigen::MatrixXd distance);
};

