#pragma once

#include "Object.h"
#include "IBS.h"

class FuncRegion
{
public:
	FuncRegion(Object* obj);

public:
	void mapFromIbs(IBS* ibs);
	void computeFeature();

	void completeDetermingScore();

private:
	void computePFH();
	QVector<double> computePfhForSample(int sIdx);
	void computeDirHist();
	void computeHeightHist();

public:
	Object* object;
	IBS* ibs;
	QColor color;

	double maxScore;					// for illustration, normalize the weight within the region
	QVector<double> determiningScore;	// same size with that of samples in central object
	QVector<int> sampleIdxs;			// sample index that have score higher than zero
	QVector<bool> reverseNormal;		// same size with that of samples in object, but only those corresponding to samples in sampleIdxs are useful

	// Geometric features of the region
	QVector<double> pfh;			// point feature histogram
	QVector<double> dirHist;		// direction histogram
	QVector<double> heightHist;		// height histogram

	// store the symmetry information
	int patternIdx;
};