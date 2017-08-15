#include "SceneHierarchy.h"
#include "Scene.h"
#include "Visualization.h"
#include "UtilityGlobal.h"

SceneHierarchy::SceneHierarchy(Scene * s)
{
	scene = s;
	threshold = 0.32;
}

SceneHierarchy::~SceneHierarchy()
{
}

void SceneHierarchy::construct()
{
	// compute the sampling weight on the IBS between each pair of objects
	computeobjPairWeight();

	// hierarchical agglomerative clustering
	QVector< ObjCluster > allClusterSet  = getAllClusteringStep();

	// build the tree structure
	buildTree(allClusterSet);
}

void SceneHierarchy::computeobjPairWeight()
{
	objPairWeight = Eigen::MatrixXd::Zero(scene->objects.size(), scene->objects.size());

	for (auto ibs : scene->ibsSetScene)
	{
		ibs->computeSampleWeightForTri();

		int objIdx1 = ibs->obj1->origIdx[0];
		int objIdx2 = ibs->obj2->origIdx[0];
		objPairWeight( objIdx1, objIdx2 ) = ibs->totalWeight;
		objPairWeight( objIdx2, objIdx1 ) = ibs->totalWeight;
	}
}

double SceneHierarchy::computeCommitment(QVector<int> objGroup1, QVector<int> objGroup2)
{
	double w1 = 0;			// object group 1 <-> object group 2
	QVector<bool> inGroup1(objPairWeight.cols(), false);
	for (auto i:objGroup1)
	{
		inGroup1[i] = true;
		for (auto j:objGroup2)
		{			
			w1 += objPairWeight(i, j);
		}
	}

	double w2 = 0;
	for (auto i:objGroup1)
	{
		for (int j=0; j<inGroup1.size(); j++)
		{
			if (!inGroup1[j])
			{
				w2 += objPairWeight(i, j);
			}
		}
	}

	if (w2 == 0)
	{
		assert(w1 == 0);
		return 1.0;
	}
	else
	{
		return w1/w2;						// equation 1 in Xi's paper
	}

}

double SceneHierarchy::computeSimilarity(QVector<int> objGroup1, QVector<int> objGroup2)
{
	double s = 0;

	s += computeCommitment(objGroup1, objGroup2);
	s += computeCommitment(objGroup2, objGroup1);

	return s;
}

Eigen::MatrixXd SceneHierarchy::computeSimilarityMatrix(QVector< ObjCluster > clusterSet)
{
	// contextual similarity
	int n = clusterSet.size();
	Eigen::MatrixXd M1 = Eigen::MatrixXd::Zero(n, n);
	for (int i=0; i<n; i++)
	{
		for (int j=i+1; j<n; j++)
		{
			M1(i, j) = computeSimilarity(clusterSet[i].origObjIdx, clusterSet[j].origObjIdx);
			M1(j, i) = M1(i, j);
		}
	}

	// contacting detection
	Eigen::MatrixXi M2 = Eigen::MatrixXi::Zero(n, n);
	for (int i=0; i<n; i++)
	{
		for (int j=i+1; j<n; j++)
		{
			bool contact = false;
			for (auto objIdx1 : clusterSet[i].origObjIdx)
			{
				for (auto objIdx2 : clusterSet[j].origObjIdx)
				{
					if (scene->objPairContact(objIdx1, objIdx2) == 1)
					{
						contact = true;
						M2(i, j) = 1;
						M2(j, i) = 1;
						break;
					}
				}
				if (contact)
				{
					break;
				}
			}
		}
	}

	// combined similarity
	Eigen::MatrixXd M3 = M1;
	if (M2.maxCoeff() > 0)
	{
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
			{
				M3(i, j) = M1(i, j) * M2(i, j);
			}
		}
	}

	return M3;
}

void SceneHierarchy::findMergingClusters(QVector< ObjCluster > clusterSet, QVector< QSet<int> > &mergingClusters, QVector<bool> &isRemain)
{
	// compute affinity matrix
	Eigen::MatrixXd S = computeSimilarityMatrix(clusterSet);	

	// find the pairs with similarity larger than the threshold
	int n = clusterSet.size();
	for (int i = 0; i<n; i++)
	{
		for (int j = i+1; j<n; j++)
		{
			if (S(i,j) > threshold)
			{
				bool found = false;
				QVector<int> clusterToMerge;
				for (int k=0; k<mergingClusters.size(); k++)
				{
					if (mergingClusters[k].contains(i) || mergingClusters[k].contains(j))
					{
						if (!found)
						{
							mergingClusters[k] << i << j;		// QSet does not allow duplicates
							found = true;
						}

						clusterToMerge << k;					// contain(i) or contain(j)
					}
				}

				if (!found)
				{
					QSet<int> cluster;
					cluster << i << j;
					mergingClusters << cluster;
				}

				isRemain[i] = false;
				isRemain[j] = false;

				// merge: for sets like (1 3) (2,4) connected  by (3,4)
				if (clusterToMerge.size() > 1)
				{
					assert(clusterToMerge.size() == 2);

					int idx1 = clusterToMerge[0];
					int idx2 = clusterToMerge[1];

					QVector< QSet<int> > mergingClusters_new;

					for (int k=0; k<mergingClusters.size(); k++)
					{
						if (k!=idx2)
						{
							mergingClusters_new << mergingClusters[k];
						}
					}

					for (auto cIdx : mergingClusters[idx2])
					{
						mergingClusters_new[idx1].insert(cIdx);
					}

					mergingClusters = mergingClusters_new;
				}
			}
		}
	}
}

QVector< ObjCluster > SceneHierarchy::getAllClusteringStep()
{
	QVector< ObjCluster > allClusterSet;
	for (int i=0; i<scene->objects.size(); i++)
	{
		allClusterSet << ObjCluster(i);
	}

	QVector< ObjCluster >  currClusterSet = allClusterSet;
	while (currClusterSet.size() > 1)
	{
		int n = currClusterSet.size();			// one single object forms a single cluster

		// find clusters need to be merged
		QVector< QSet<int> > mergingClusters;
		QVector<bool> isRemain(n, true);
		findMergingClusters(currClusterSet, mergingClusters, isRemain);

		// generate the new set of current clusters
		QVector< ObjCluster > remainingClusters;
		for (int i=0; i<n; i++)
		{
			if (isRemain[i])
			{
				remainingClusters << currClusterSet[i];
			}
		}

		QVector< ObjCluster > newClusterSet;
		for (int i=0; i<mergingClusters.size(); i++)
		{
			ObjCluster newCluster;

			newCluster.idx = allClusterSet.size();
			newCluster.origObjIdx.clear();
			newCluster.childClusterIdx.clear();
			for (auto j:mergingClusters[i])
			{
				newCluster.origObjIdx << currClusterSet[j].origObjIdx;
				newCluster.childClusterIdx << currClusterSet[j].idx;
			}

			newClusterSet << newCluster;
			allClusterSet << newCluster;
		}

		currClusterSet.clear();
		currClusterSet << remainingClusters << newClusterSet;
	}

	return allClusterSet;
}

void SceneHierarchy::buildTree(QVector< ObjCluster > allClusterSet)
{
	tree<ObjCluster>::iterator currNode = clusteringTree.begin();
	currNode = clusteringTree.insert(currNode, allClusterSet.last());
	addChildren(currNode, allClusterSet);

	if (scene->visualizeHierarchy)
	{
		Visualization::visualizeTree(&clusteringTree, scene->filename+"_ibsh");
	}
}

void SceneHierarchy::addChildren( tree<ObjCluster>::iterator parentNode, QVector< ObjCluster > allClusterSet)
{
	QSet<int>  childClusterIdx = allClusterSet[(*parentNode).idx].childClusterIdx;

	if ( childClusterIdx.isEmpty() )
	{
		return;
	}

	for (auto childIdx: childClusterIdx)
	{
		tree<ObjCluster>::iterator childNode = clusteringTree.append_child(parentNode, allClusterSet[childIdx]);
		addChildren(childNode, allClusterSet);		
	}
}

void SceneHierarchy::computeFeature(int cIdx)
{
	// find the central object in the tree
	tree<ObjCluster>::leaf_iterator leafNode = findLeafNode(cIdx);
	if (leafNode == clusteringTree.end_leaf())
	{
		return;
	}

	//debugBox("found leaf node");

	// find the IBSs in the community at each level from bottom up; number of levels = depth
	QVector< QVector< QVector<IBS*> > > profiles; // a set of IBS for each level
	QVector< QVector< QVector<bool> > > reverseIBS;
	findProfiles(profiles, reverseIBS, leafNode);

	//QStringList message;
	//message << "found profiles: level = " + QString::number(profiles.size()) + " ibsNum = " + QString::number(profiles[0].size());
	//for (int i=0; i<profiles[0].size(); i++)
	//	message << QString::number(profiles[0][i].size()); 
	//debugBoxList(message);

	// compute the features
	ObjLevelFeature feature;
	for (int i=0; i<profiles.size(); i++)
	{
		QVector<IBS> levelIBS;
		for (int j=0; j<profiles[i].size(); j++)
		{
			IBS newFeatureIBS(profiles[i][j], reverseIBS[i][j]);

			levelIBS << newFeatureIBS;
		}
		
		feature.featureIBSs << levelIBS;
	}
	feature.save(scene->filename);
}

tree<ObjCluster>::leaf_iterator SceneHierarchy::findLeafNode(int cIdx)
{
	tree<ObjCluster>::leaf_iterator leafNode = clusteringTree.begin_leaf();
	while (leafNode != clusteringTree.end_leaf() )
	{
		assert( (*leafNode).origObjIdx.size() == 1 );
		if ( (*leafNode).origObjIdx[0] == cIdx )
		{
			assert((*leafNode).idx == cIdx);
			break;
		}
		leafNode++;
	}
	if (leafNode == clusteringTree.end_leaf())
	{
		debugBox("Cannot find the central object! " + scene->name);
	}

	return leafNode;
}

void SceneHierarchy::findProfiles(QVector< QVector< QVector<IBS*> > > &profiles, QVector< QVector< QVector<bool> > > &reverseIBS, tree<ObjCluster>::leaf_iterator leafNode)
{
	int depth = clusteringTree.depth(leafNode);
	for (int i=0; i<depth; i++)
	{
		int cIdx = (*leafNode).idx;
		QVector<int>  objIdxSet1 = (*leafNode).origObjIdx;

		tree<ObjCluster>::leaf_iterator parent = clusteringTree.parent(leafNode);
		tree<ObjCluster>::sibling_iterator sibling = clusteringTree.begin(parent);

		QVector< QVector<IBS*> > profile;  // the IBSs between two communities 
		QVector< QVector<bool> > reverse;  // make sure all the normals on IBSs are pointing to the community leafnode in
		while (sibling != clusteringTree.end(parent))
		{
			if ((*sibling).idx != cIdx)
			{			
				QVector<IBS*> ibsSet;
				QVector<bool> needReverse;

				QVector<int>  objIdxSet2 = (*sibling).origObjIdx;
				for (auto i:objIdxSet1)
				{
					for (auto j:objIdxSet2)
					{
						if (scene->objPair2IbsSet(i,j) != -1)
						{
							ibsSet << scene->ibsSetScene[scene->objPair2IbsSet(i,j)];
							if (j > i)
							{
								needReverse << true;
							}
							else
							{
								needReverse << false;
							}
						}
					}
				}

				if (!ibsSet.isEmpty())
				{
					profile << ibsSet;
					reverse << needReverse;
				}
			}
			sibling++;
		}

		profiles << profile;
		reverseIBS << reverse;

		leafNode = parent;
	}
}