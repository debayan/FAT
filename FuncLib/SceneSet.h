#pragma once
#include "Scene.h"

class SceneSet
{
public:
    SceneSet();
	~SceneSet();

public:
	void addScenes( QStringList filenames );
	void addScene( Scene* s );
	bool deleteScene(int i);	
	void clearAllScenes();

	QVector<int> getScenesWithNoCentralObject();
	void buildInteractionHierarchy();
	void loadInteractionHierarchy();

	QVector<QString> getModelNames(QString filename);
	void outputCentralObject(int idx);
	void outputCentralObjects();

	// compute all the models list inside the file with distPara and hierPara, note that each modelName should contain suffix
	void computeFeatures(QString filename, DistParameter distPara, HierParameter hierPara); 

	// For comparison to IBS paper:
	void computeIBSH();
	void computeIBSH(QString filename, bool useTopo);

	// For comparison to our hierarchical structure:
	// treat the interactions as a set and output the corresponding feature
	void outputInterSetFeature();

	// compute geometry features for central objects in each scene in the list
	void computeGeoFeature(QString filename);

public:
	QVector<Scene *> scenes;
};