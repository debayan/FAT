#pragma once

#include "IBS.h"
#include "Scene.h"
#include "Qhull.h"
#include "SubTreeSimilarity.h"

class IbsGenerator
{
public:
	IbsGenerator();
	~IbsGenerator();

public:
	QVector<IBS*> computeIBS(Scene * scene, QVector<Object*> objects);

private:
	void computeVoronoi();
	void findRidges();
	void buildIBS();
	
	QVector<Vec3d> getInputForVoronoi();	
	int findRidgesAroundVertex(vertexT *atvertex);  // find all unvisited Voronoi ridges for vertex (i.e., an input site)
	SurfaceMeshModel * buildIbsMesh(int i,  QVector<QPair<int, int>>& samplePairs);

private:
	Scene *scene;
	QVector<IBS*> ibsSet;
	//QVector<int> activeObjIdx;		// useless, can be commented
	QVector<Object*> objects;

private:
	orgQhull::Qhull* qhull;	
	QVector<Vec3d> voronoiVertices;

	// same size with ibsSet
	QMap<QPair<int, int>, int> objPair2IbsIdx;
	QVector< QVector<int> > ibsRidgeIdxs;	// using objPair2IbsIdx to index this vector

	QVector< QVector<int> > ridges;
	QVector< QVector<int> > ridgeSitePair;

	QVector<int> sampleObjIdx;				//the corresponding objectIdx for each sample points
	QVector<int> sampleLocalIdx;
};