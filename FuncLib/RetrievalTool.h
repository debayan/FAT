#pragma once

#include "InterHierarchy.h"
#include <QtCore>
#include "Eigen\Core"
#include "DistMeasure.h"
#include "UtilityGlobal.h"
#include "qfileinfo.h"
#include <fstream>
#include "LFD.h"
#include "SceneHierarchy.h"
#include "InterSet.h"

enum RETRIEVAL_MODE {NORMAL, SHAPE2POSE, PART};
enum FEATURE_TYPE {ICON, LFD, IBSH, LFDICON, POSE, ISET, GEO, GEOICON};

class RetrievalTool
{

public:
	RetrievalTool(QString filepath,int datatype = 0); //datatype = 0: our data. datatype = 1: shape2pose data
	~RetrievalTool(void);

	void computeDistance(FEATURE_TYPE featureType, QString paraString = "", double weight = 0.5);  // compute distance matrix for different type of feature; weight is only used for combining icon with lfd
	void computeFeature(FEATURE_TYPE featureType, QString featureFile);

	void computeICONFeatures(DistParameter distPara, HierParameter hierPara);
	void computeISETFeatures(DistParameter distPara);

	void analyzePR(FEATURE_TYPE featureType);  // evaluate the retrieval performance for each feature type
	bool prReady(); // check whether the pr-curve for the current featureType is ready
	void findFiles(QVector<QString> &modelPaths, const QString &path, const QString &fileext); //search and create .sl file

	QVector<int> returnRetrievalResult(int num);    // return the first num similar results for currHierarchyIdx	
	QStringList returnRetrievalResultName(int num); // return the name of first num similar results for currHierarchyIdx	

	Eigen::MatrixXd getPRShape(); // get the pr-curve for the retrieval result using current scene as query
	Eigen::MatrixXd getPRClass(); // get the pr-curve for the retrieval results of the category current scene in
	Eigen::MatrixXd getPRall(); // get the pr-curve for all the retrieval results

	Eigen::MatrixXd normalizePR(QVector<QVector<double>> &P,QVector<QVector<double>> &R,int idx,QString fn);

	void setCurrentRetrievalMode(RETRIEVAL_MODE mode);

private:
	void computeICONDistance();
	void computeLFDDistance();
	void computeLFDICONDistance(double weight);
	void computePOSEDistance();
	void computeIBSHDistance();
	void computeISETDistance();
	void computeGEODistance();
	void computeGEOICONDistance(double weight);

	int typeIdx(FEATURE_TYPE featureType);

public:
	FEATURE_TYPE currType;
	double combWeight;
	bool combUpdate;
	int currIdx;  // current selected scene index for retrieval
	int ibshDepth;
	double partCombWeight;
	bool partCombUpdate;

	RETRIEVAL_MODE currentRetrievalMode;				// 0 normal, 1 shape2pose

	QVector<QPixmap> sceneImages;			// the corresponding snapshot for each scene
	QVector<QString> scenePaths;			// the file path for each scene
	QVector<QString> sceneNames;			// the name for each scene
	QVector<QString> sceneLabels;			// the category label for each scene
	QVector<QPair<QString,int>*> labelNums; // number of scenes in each category

	QString dirName;  // the directory the output file in 
	QString filePath; // the file path to store output
	QString paraString; 

	QVector<InterHierarchy*> icons;  // the interaction hierarchies for all the scenes
	QVector<Eigen::Matrix3Xd> poses;
	QVector<ObjLevelFeature> ibshs;
	QVector<InterSet> sets;

	QVector< Eigen::MatrixXd > distance; // the distance between each pair of objects in the data set
	QVector< Eigen::MatrixXd > prShape;  // pr values for each shape
	QVector< Eigen::MatrixXd > prClass;  // pr value for each category of shapes
	QVector< Eigen::MatrixXd > prall;	 // pr value for each category of shapes
};