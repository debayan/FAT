#pragma once
#include "Object.h"

enum SYM_TYPE{TRANS, ROT, GRID};

class SymGroup
{
public:
	SymGroup(Object* obj);
	SymGroup(Object* obj, SYM_TYPE type);

public:
	Object* obj;

	SYM_TYPE type;

	// parameters: choose which to use based on the type, ignore others
	Eigen::Vector3d trans;

	Eigen::Vector3d rotCenter;
	double rotRadius;

	Eigen::Vector3d gridT1;
	Eigen::Vector3d gridT2;

	// store the order of the pattern inside each group
	QVector<int> pIdx1D; // for trans and rot
	QVector< QVector<int> > pIdx2D;  // for grid, inside vector corresponds to gridT1, outside vector corresponds to gridT2

	//Loading functions
	void LoadFromFile(QString file);

	// 
	void addOffsetToPatternIdx(int offset);
	int symLevel(int p1, int p2);
};