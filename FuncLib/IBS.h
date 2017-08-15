#pragma once

#include "SurfaceMeshModel.h"
#include "Sampler.h"
#include "RenderObjectExt.h"

class FuncRegion;
class Scene;
class Object;
using namespace SurfaceMesh;

class IBS
{

public:
	IBS();
	IBS(Scene *s);
	IBS(QVector<IBS*> ibsSet, QVector<bool> reverseNormal);
	~IBS();

public:
	void draw(bool drawIbsSample = false, bool drawIbsWeight = false, QColor color = Qt::red);
	void computeGeomFeatures();			// geometry features
	void computeTopoFeatures();			// topology features
	void computeSampleWeightForTri();

private:
	void sampling(int num);

	void computeDirHist();
	void computeDistHist();
	void computePFH();
	QVector<double> computePfhForSample(int sIdx, bool reverseNormal);

	void computeBettiNumbers();
	void computeBettiNumbers2();		// transfer original copy from Xi, extremely slow
	void ignoreSmallHoles();

	// implementation for community features
	QVector<double> combinedPFH(QVector<IBS*> ibsSet, QVector<bool> reverseNormal);
	QVector<double> computePfhForSample(int ibsIdx, int sIdx, QVector<IBS*> ibsSet, QVector<bool> reverseNormal);
	QVector<double> combinedDirHist(QVector<IBS*> ibsSet, QVector<bool> reverseNormal);
	QVector<double> combinedDistHist(QVector<IBS*> ibsSet);
	QVector<int> combinedBettiNumber(QVector<IBS*> ibsSet);

public:
	Scene *scene;
	FuncRegion* region;

	Object * obj1;
	Object * obj2; // normal points toward obj2 by default; for IBS beteween interacing object and central object, this is always the idx for the central object 

	bool pointToCentralObject; // true is objIdx2 is the central object (objIdx2 = ibsSetScene[i]->obj2-origIdx[0])

	SurfaceMeshModel* mesh;	
	QVector<QPair<int, int>> samplePairs;  // the pair of samples on the objects that corresponds the triangle

	// importance-based sampling
	double sampleRatio;  // for smaller IBS, the number of samples should be smaller
	QVector<SamplePoint> samples;
	starlab::PointSoup sampleRender;
	double maxWeight;
	double totalWeight;

	// Geometric features
	QVector<double> pfh;		// point feature histogram
	QVector<double> dirHist;	// direction histogram
	QVector<double> distHist;	// distance histogram

	// Topological features
	QVector<int> bettiNumbers;
};