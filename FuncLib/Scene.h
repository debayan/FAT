#pragma once

#include "SceneHierarchy.h"
#include "InterHierarchy.h"
#include "dataanalysis.h"
#include <QDir>
#include <QString>

enum DIST_TYPE {L1, EMD};

enum CLUSTER_TYPE{AP, AHC};

enum IBS_DRAW_MODE{DRAW_PARTIAL, DRAW_ENTIRE};

typedef struct DrawParameter
{
	OBJECT_DRAW_MODE objMode;
	IBS_DRAW_MODE ibsMode;
	bool drawIbsSample;
	bool drawIbsWeight;
	bool drawClusterIbs;
	bool drawClusterIbsWeight;
	bool drawClusterIbsSample;
	bool drawInteraction;
	bool drawInteractionIBS;
	bool drawInteractionIR;
	bool drawInteractionObj;
} DrawParameter;

typedef struct DistParameter
{
	DIST_TYPE distType;
	bool useTopo;
	QVector<double> weight;
	QVector<double> ibsWeight;
	QVector<double> regionWeight;
	int symRatio;

	DistParameter()
	{
		useTopo = false;
		distType = L1;
		weight << 0.4 << 0.6;
		ibsWeight << 0.3 << 0.4 << 0.3;
		regionWeight << 0.3 << 0.4 << 0.3;
		symRatio = 1;
	}
	
	QString toString()
	{
		QString s = "(";

		if (distType == L1)
		{
			s += "L1,";
		}
		else
		{
			s += "EMD,";
		}

		if (useTopo)
		{
			s += "GT,";
		}
		else
		{
			s += "G,";
		}

		s += QString::number(symRatio) + ",";
		s += QString::number(weight[0]) + ",";
		s += QString::number(ibsWeight[0]) + "," + QString::number(ibsWeight[1]) + ",";
		s += QString::number(regionWeight[0]) + "," + QString::number(regionWeight[1]) + ")";

		return s;
	};

	void readFromString(QString para)
	{
		QStringList list = para.split(",");
		assert( list.size() ==  8);

		if (list[0] == "L1")
		{
			distType = L1;
		}
		else if (list[0] == "EMD")
		{
			distType = EMD;
		}

		if (list[1] == "GT")
		{
			useTopo = true;
		}
		else if (list[1] == "T")
		{
			useTopo = false;
		}

		symRatio = list[2].toInt();

		weight.clear();
		weight << list[3].toDouble();
		weight << 1.0 - weight[0];

		ibsWeight.clear();
		ibsWeight << list[4].toDouble() << list[5].toDouble();
		ibsWeight << 1.0 - ibsWeight[0] - ibsWeight[1];

		regionWeight.clear();
		regionWeight << list[6].toDouble() << list[7].toDouble();
		regionWeight << 1.0 - regionWeight[0] - regionWeight[1];
	}

}DistParameter;

typedef struct HierParameter
{
	bool methodtype;
	bool useMultipleTree;
	bool topDown;
	double similarThreshold;
	double lowerBound;
	double higherBound;

	HierParameter()
	{
		methodtype = true; // true if use ours; false if use recursive
		useMultipleTree = true;
		topDown = true;
		similarThreshold = 0.2;
		lowerBound = 0.5;
		higherBound = 0.6;
	}

	QString toString()
	{
		if(!methodtype)
		{
			QString s;

			return s;
		}

		QString s = "(";

		if (useMultipleTree)
		{
			s += "M,";
		}
		else
		{
			s += "S,";
		}

		if (topDown)
		{
			s += "T,";
		}
		else
		{
			s += "B,";
		}

		s += QString::number(similarThreshold) + ",";
		s += QString::number(lowerBound) + "," ;
		s += QString::number(higherBound) + ")";

		return s;
	};

	void readFromString(QString para)
	{
		QStringList list = para.split(",");
		if (list.size() != 5)
		{
			methodtype = false;
			return;
		}

		if (list[0] == "M")
		{
			useMultipleTree = true;
		}
		else if (list[0] == "S")
		{
			useMultipleTree = false;
		}

		if (list[1] == "T")
		{
			topDown = true;
		}
		else if (list[1] == "B")
		{
			topDown = false;
		}

		similarThreshold = list[2].toDouble();
		lowerBound = list[3].toDouble();
		higherBound = list[4].toDouble();
	}

}HierParameter;

class Scene
{
public:
	Scene();
	~Scene();

public:
	//void load(QDir dir);
	void load(QString filename);
	void loadSnapshot();

	void draw(DrawParameter para);

	bool hasCentralObj();
    int getSelectedObjectID(Eigen::Vector3d orig, Eigen::Vector3d dir);
	void setObjCentralState(QVector<bool> sel);
	void reverseObjCentralState(int objIdx);

	void setInteractionSelectState(QVector<bool> sel);
	void setInteractionHierSelectState(QVector<bool> sel);
	void resetInteractionSelectState();

	alglib::ahcreport* getAHC();

	QVector<int> getInteractions(int hierIdx);
	QVector<int> getInteractionsPreOrder(int hierIdx);
	QVector<int> getInterHierParentNodeIdx(int hierarchy);
	QVector<int> getInterHierParentNodeIdxPreOrder(int hierarchy);
	int getCandiHierNum();
	
	//Output
	void updateResultflodername();

	void outputCentricObject();

public:
	void computeInteractions();
	void buildInteractionHierarchy();
	void loadInteractionHierarchy();

	void analyzeInteractions();
	void generateCentralObject();

	void saveCentralIdx();

	void computeIBSH();
	void outputInterSetFeature();

	void computeGeometryFeature(int batnum);

private:
	
	void constructInteractions();
	void computeMaxHistForInteraction();

	void detectObjContact();

private:
	void initialize();
	void clearAllObjects();
	void loadFromTxt(QString filename);
	void loadFromObj(QString filename);
	void loadCentralIdx();

	void samplePoints();
	void updateInteractions();

	QVector<double> computeMaxHist(QVector<Interaction*> inters);

public:
	QString name;
	QString filename;
	QPixmap snapshot;
	QString resultFolder;

	QVector<QColor> colorMap;
	Eigen::AlignedBox3d bbox;
	Vec3d upright;				// upright direction of the scene

	QVector<Object*> objects;   // the input objects
	Object* centralObj;			// combine all the selected central objects into one object
	Object* allObj;

	QVector<IBS*> ibsSetScene;			// the ibs computed from the entire scene: for interaction identification
	Eigen::MatrixXi objPair2IbsSet;		// for IBS paper's scene hierarchy construction: if there is a IBS between two objects, store the IBS idx
	Eigen::MatrixXi objPairContact;		// for IBS paper's scene hierarchy construction: 1 if two objects are contacting each other, 0 otherwise

	QVector<int> obj2Inter;				// map each object to its corresponding interaction, -1 if no interaction exist
	QVector<int> interObjIdx;			// the interacting objects after central objects being selected
	QVector<Interaction*> interactions; // the interaction between each interacting object and the central object
	QVector<bool> isSeletedInteraction;	// for illustration; only show the seleted interactions
	QVector< QVector<int> > interIbsSet;// the scene ibs index corresponding to each interacting objects
	QVector<double> yMaxHistInteraction;
	QVector<QColor> interactionColors;
	QColor baseColor;

	//////////////////////////////////////////////////////////////////////////
	//for clustering
	QVector< Interaction* > clusteredInteractions; // Interactions after clustering
	int selectedClusterIdx;
	QVector<double> yMaxHistCluster;

	Eigen::MatrixXd affinity;
	Eigen::MatrixXd ibsAffinity;
	Eigen::MatrixXd regionAffinity;

	CLUSTER_TYPE clusterType;

	DistParameter distPara;
	HierParameter hierPara;
	bool visualizeHierarchy;

	//////////////////////////////////////////////////////////////////////////
	// for hierarchy
	alglib::ahcreport* ahcreport;
	InterHierarchy* hierarchy;
	QVector<Interaction*> allInteractions;		// all interactions including original interactions and merged interactions
	QVector<int> interactionsIdx;				// indexing interactions in allInteractions set
	QVector<int> intersParentIdx;				// keep track of parent node of each interactions

	QVector<bool> lastSelectedInters;
	QVector<bool> isSelectedInterGroup;

	//////////////////////////////////////////////////////////////////////////
	SceneHierarchy* sceneHierarchy;

	// statistics 
	int constructionTime;
	int interactionTime;
};