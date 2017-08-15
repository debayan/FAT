#pragma once

#include "Object.h"
#include "tree.h"
#include "IBS.h"
#include <QFile>

typedef struct ObjCluster
{
	int idx;
	QVector<int> origObjIdx;
	QSet<int> childClusterIdx;

	ObjCluster()
	{
		idx = -1;
		origObjIdx.clear();
		childClusterIdx.clear();
	}

	ObjCluster(int objIdx)
	{
		idx = objIdx;
		origObjIdx << objIdx;
	}
}ObjCluster;

typedef struct ObjLevelFeature
{
	QVector< QVector<IBS> > featureIBSs; // for each level, there are a set of IBSs

	void save(QString filename)
	{
		QFile file(filename + ".ibsh");
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream out(&file);

		// output features level by level
		for (int k=0; k<featureIBSs.size(); k++)
		{
			out << "level " << k << endl;
			out << "ibsNum " << featureIBSs[k].size() << endl;

			for (int j=0; j<featureIBSs[k].size(); j++)
			{
				out << "pfh";
				for (int i=0; i<featureIBSs[k][j].pfh.size(); i++)
				{
					out << " " << featureIBSs[k][j].pfh[i];
				}
				out << endl;

				out << "dirHist";
				for (int i=0; i<featureIBSs[k][j].dirHist.size(); i++)
				{
					out << " " << featureIBSs[k][j].dirHist[i];
				}
				out << endl;

				out << "distHist";
				for (int i=0; i<featureIBSs[k][j].distHist.size(); i++)
				{
					out << " " << featureIBSs[k][j].distHist[i];
				}
				out << endl;
			}
		}
		file.close();
	}

	void load(QString filename)
	{
		featureIBSs.clear();
		QFile file( filename + ".ibsh");
		if (file.open(QIODevice::ReadOnly | QIODevice::Text)) 
		{
			QTextStream in(&file);	
			while (!in.atEnd()) 
			{
				QString line = in.readLine();
				QStringList list = line.split(" ");
				if(list.isEmpty()) continue;

				if (list[0] == "level")
				{
					assert(list[1].toInt() == featureIBSs.size());

					QVector<IBS> levelIBSs;

					line = in.readLine();
					list = line.split(" ");
					assert(list[0] == "ibsNum");
					int n = list[1].toInt();

					for (int j=0; j<n; j++)
					{
						// load the features
						IBS ibs; 	

						line = in.readLine();
						list = line.split(" ");
						assert(list[0] == "pfh");
						for (int i=1; i<list.size(); i++)
						{
							ibs.pfh << list[i].toDouble();
						}
						assert(ibs.pfh.size() == 250);

						line = in.readLine();
						list = line.split(" ");
						assert(list[0] == "dirHist");
						for (int i=1; i<list.size(); i++)
						{
							ibs.dirHist << list[i].toDouble();
						}
						assert(ibs.dirHist.size() == 10);

						line = in.readLine();
						list = line.split(" ");
						assert(list[0] == "distHist");
						for (int i=1; i<list.size(); i++)
						{
							ibs.distHist << list[i].toDouble();
						}
						assert(ibs.distHist.size() == 10);

						levelIBSs << ibs;
					}

					featureIBSs << levelIBSs;
				}
			}
			file.close();
		}
	}

}ObjLevelFeature;

class SceneHierarchy
{
public:
	SceneHierarchy(Scene * s);
	~SceneHierarchy();

public:
	void construct();
	void computeFeature(int cObjIdx);

private:
	void computeobjPairWeight(); // compute the sampling weight on the IBS between each pair of objects, used to measure commitment
	double computeCommitment(QVector<int> objGroup1, QVector<int> objGroup2);
	double computeSimilarity(QVector<int> objGroup1, QVector<int> objGroup2);
	Eigen::MatrixXd computeSimilarityMatrix(QVector< ObjCluster > clusterSet);

	void findMergingClusters(QVector< ObjCluster > clusterSet, QVector< QSet<int> > &mergingClusters, QVector<bool> &isRemain);
	QVector< ObjCluster > getAllClusteringStep(); 

	void buildTree(QVector< ObjCluster > allClusterSet);
	void addChildren(tree<ObjCluster>::iterator parentNode, QVector< ObjCluster > allClusterSet);

	tree<ObjCluster>::leaf_iterator findLeafNode(int cIdx);
	void findProfiles(QVector< QVector< QVector<IBS*> > > &profiles, QVector< QVector< QVector<bool> > > &reverseIBS, tree<ObjCluster>::leaf_iterator leafNode);

public:
	Scene * scene;
	tree<ObjCluster> clusteringTree;

private:
	Eigen::MatrixXd objPairWeight;
	double threshold;
};
