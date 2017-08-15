#include "InterHierarchy.h"
#include "Scene.h"
#include "Visualization.h"
#include "UtilityGlobal.h"
#include <fstream>

InterHierarchy::InterHierarchy(Scene* s)
{
	scene = s;
}

InterHierarchy::InterHierarchy(QString scene_name)
{
	scene = new Scene();
	scene->filename = scene_name;
	scene->name = QFileInfo(scene_name).baseName();
}

InterHierarchy::InterHierarchy( QString sceneName, QString paraString )
{
	scene = new Scene();
	scene->filename = sceneName;
	scene->name = QFileInfo(sceneName).baseName();

	QStringList list = paraString.split("(");
	scene->distPara.readFromString(list[1].left( list[1].size()-1 ));
	if (list.size() > 2)
	{
		scene->hierPara.readFromString(list[2].left( list[2].size()-1 ));
	}
	else
	{
		scene->hierPara.readFromString(QString());
	}
}

InterHierarchy::~InterHierarchy()
{
	deleteMergedInteractions();
}

void InterHierarchy::deleteMergedInteractions()
{
	for (int i=0; i<mergedInteractions.size(); i++)
	{
		delete mergedInteractions[i];
	}
	mergedInteractions.clear();
}

void InterHierarchy::construct()
{
	initializeBinaryTree();

	if (scene->hierPara.useMultipleTree)
	{
		generateMultipleTrees();
	}

	mergeBranches();

	addVirtualRoot();
}

void InterHierarchy::initializeBinaryTree()
{
	alglib::ahcreport* ahcreport = scene->getAHC();			// AHC clustering result is stored in ahcreport
	
	// generate all the interactions first including the individual interactions and also the clustered interactions
	QVector<Interaction*> interactions;
	interactions << scene->interactions;

	mergedInteractions.clear(); // store new generated merged interactions 
	for (int i=0; i<scene->interactions.size()-1; i++)
	{
		QVector<Interaction*> interactionToCluster;
		interactionToCluster << interactions[(int)ahcreport->z[i][0]];
		interactionToCluster << interactions[(int)ahcreport->z[i][1]];

		Interaction* clustedInteraction = new Interaction(scene, interactionToCluster, COMBINED_INIT);
		clustedInteraction->mergeDist = ahcreport->mergedist[i];
		interactions << clustedInteraction;
		mergedInteractions << clustedInteraction;
	}
	
	// build the tree
	tree<Interaction*>* binaryTree = new tree<Interaction*>();
	interTree << binaryTree;

	tree<Interaction*>::iterator currNode = binaryTree->begin();
	currNode = binaryTree->insert(currNode, interactions.last());
	addChildren(currNode, interactions, scene->interactions.size()-2);		// the last one is root node, so adding 
																			// children from the one before last one

	if (scene->visualizeHierarchy)
	{
		Visualization::visualizeTree(binaryTree, scene->filename + "_" + scene->distPara.toString() + "_tree_initial");
	}
}

inline bool readFromCSVfile(QString name, Eigen::MatrixXd &matrix)
{
	std::ifstream file(name.toStdString().c_str());
	if(!file)
		return false;
	std::string value,in,line;
	QVector<double> data;
	int rows = 0;
	while(file.good())
	{
		while(std::getline(file,line))
		{
			std::istringstream stream(line);
			while(std::getline(stream,value,','))
			{
				QString qstr = QString::fromStdString(value);
				data.push_back(qstr.toDouble());
			}
			rows++;
		}
	}
	int columns = data.size()/rows;
	matrix = Eigen::MatrixXd::Zero(rows,columns);
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
		{
			matrix(i,j) = data[i*columns + j];
		}
		file.close();
		return true;
}

void InterHierarchy::addChildren( tree<Interaction*>::iterator parentNode, QVector<Interaction*>& interactions, int idx )
{
	if ( idx < 0 || idx >= scene->interactions.size()-1)
	{
		return;
	}

	alglib::ahcreport* ahcreport = scene->getAHC();

	for (int i=0; i<2; i++)
	{
		int childIdx = ahcreport->z[idx][i];
		tree<Interaction*>::iterator childNode = interTree[0]->append_child(parentNode, interactions[childIdx]);

		if (childIdx >= scene->interactions.size())
		{
			addChildren(childNode, interactions, childIdx - scene->interactions.size());
		}
	}
}

QVector< QVector<int> > combination(int n, int m)
{
	QVector< QVector<int> > results;

	if( 0==n || 0==m || n<m )   
	{
		return results;
	}

	int size=n-1;       
	int position;
	QVector<int> vecOrder;
	for(int i=0; i<m; i++)
	{
		vecOrder.push_back(i);
	}

	position=m-1;
	while( 1 )
	{
		results << vecOrder;

		if ( vecOrder[0]==(size-m+1) )
		{
			break;
		}

		if( (position==m-1 && vecOrder[position]<size )|| (position<m-1 && vecOrder[position]+1<vecOrder[position+1]) ) //要求最右边的数字不等于size才加1&&不是最右边的数字加1后小于其右边的数
		{
			vecOrder[position]++;
			continue;
		}
		else
		{
			position--;
			if( vecOrder[position]+1 == size )
			{
				position--;
			}
			vecOrder[position]++;

			for(int i=position; i<m-1; i++)
			{
				vecOrder[i+1]=vecOrder[i]+1;
			}

			for(int i=m-1; i>position; i--)
			{
				if((i==m-1 && vecOrder[i]<size)||(i<m-1 && vecOrder[i]+1<vecOrder[i+1]))
				{
					position=i;
					break;
				}
			}
		}
	}

	int count = 1;
	for (int i=0; i<m; i++)
	{
		count *= (n-i);
	}
	for (int i=1; i<=m; i++)
	{
		count /= i;
	}

	assert(results.size() == count );

	return results;
}

QVector< QVector<int> > getAllCombination(int n)
{
	QVector< QVector<int> > results;
	for (int m=1; m<=n; m++)
	{
		results << combination(n, m);
	}
	return results;
}

void InterHierarchy::findFlexibleNodes( QVector< tree<Interaction*>::iterator >& flexibleNode, QVector< QVector<int> >& delCombs )
{
	unlabelAllInteractions();

	// 1. find all the nodes whose mergeDist is between threshold1 and threshold2
	tree<Interaction*>::iterator iter = interTree[0]->begin();			// interTree[0] represents the initial binary tree
	double maxMergeDist = (*iter)->mergeDist;
	double threshold1 = scene->hierPara.lowerBound * maxMergeDist;		// the merge dist of root node must be the max merge dist
	double threshold2 = scene->hierPara.higherBound * maxMergeDist;
	flexibleNode.clear();
	while (iter != interTree[0]->end())
	{
		if ( (*iter)->mergeDist > threshold1 && (*iter)->mergeDist < threshold2)
		{
			(*iter)->nodeIdx = flexibleNode.size();
			flexibleNode << iter;
		}
		iter++;
	}

	// 2. keep record of the posterity of each flexibleNode
	QVector< QVector<int> > posterity(flexibleNode.size());
	for (int i = 0; i<flexibleNode.size(); i++)
	{
		tree<Interaction*>::iterator node = flexibleNode[i];
		tree<Interaction*>::iterator parent = interTree[0]->parent(node); 
		while ( interTree[0]->is_valid(parent) )
		{
			if ( (*parent)->nodeIdx != -1 )				// if (*parent)->nodeIdx != -1, then the parent node is also flexible
			{
				posterity[(*parent)->nodeIdx] << i;
			}
			parent = interTree[0]->parent(parent); 
		}		
	}

	// 3. generated all possible trees by keeping or deleting the corresponding nodes
	QVector< QVector<int> > nodeToDelete = getAllCombination(flexibleNode.size());

	// 4. delete invalid combination: if it's ancestry has not been deleted, then the node should not be deleted
	delCombs.clear();
	for (int i=0; i<nodeToDelete.size(); i++)
	{
		// convert to binary vector
		QVector<bool> deleteFlag(flexibleNode.size(), false);
		for (int j=0; j<nodeToDelete[i].size(); j++)
		{
			deleteFlag[nodeToDelete[i][j]] = true;
		}

		// check whether if the node is kept then all it's posterity is kept
		bool valid = true;
		for (int j=0; j<flexibleNode.size(); j++)
		{
			if (!deleteFlag[j])
			{
				for (int k=0; k<posterity[j].size(); k++)
				{
					if (deleteFlag[posterity[j][k]])
					{
						valid = false;
						break;
					}
				}

				if (!valid)
				{
					break;
				}
			}
		}

		if (valid)
		{
			delCombs << nodeToDelete[i];
		}
	}
}

void InterHierarchy::generateMultipleTrees()
{
	QVector< tree<Interaction*>::iterator > flexibleNode;
	QVector< QVector<int> > delCombs;
	findFlexibleNodes(flexibleNode, delCombs);
	
	for (int i=0; i<delCombs.size(); i++)
	{
		unlabelAllInteractions();
		tree<Interaction*>* newTree = new tree<Interaction*>( *(interTree[0]) );

		tree<Interaction*>::post_order_iterator iter = findFirstUncheckedNodeBottomUp(newTree);
		while (iter != newTree->end_post())
		{
			(*iter)->isChecked = true;

			for (int j=0; j<delCombs[i].size(); j++)
			{
				if ( (*iter) == ( *(flexibleNode[delCombs[i][j]]) ))
				{
					iter = newTree->flatten( iter);
					newTree->erase( iter );
					break;
				}
			}

			iter = findFirstUncheckedNodeBottomUp(newTree);
		}

		interTree << newTree;

		if (scene->visualizeHierarchy)
		{
			Visualization::visualizeTree(newTree, scene->filename + "_" + scene->distPara.toString() + "_tree_possible_"+QString::number(i+1));
		}
	}
}

void InterHierarchy::mergeBranches()
{
	for (int i=0; i<interTree.size(); i++)
	{
		unlabelAllInteractions();

		double maxMergeDist = (*(interTree[i]->begin()))->mergeDist;
		double threshold1 = scene->hierPara.similarThreshold * maxMergeDist;
		double threshold2 = scene->hierPara.higherBound * maxMergeDist;

		if (scene->hierPara.topDown)
		{
			tree<Interaction*>::breadth_first_iterator iter = findFirstUncheckedNodeTopDown(interTree[i]);
			while (iter != interTree[i]->end_breadth_first())
			{
				(*iter)->isChecked = true;

				tree<Interaction*>::breadth_first_iterator parent = interTree[i]->parent(iter);
				if ( interTree[i]->is_valid(parent) )
				{
					bool similarNode = abs( (*iter)->mergeDist - (*parent)->mergeDist ) < threshold1;
					bool meaninglessNode = (*iter)->mergeDist > threshold2;
					if ( similarNode || meaninglessNode )
					{
						iter = interTree[i]->flatten(iter);
						interTree[i]->erase(iter);
					}
				}
				iter = findFirstUncheckedNodeTopDown(interTree[i]);
			}
		}
		else
		{
			tree<Interaction*>::post_order_iterator iter = findFirstUncheckedNodeBottomUp(interTree[i]);
			while (iter != interTree[i]->end_post())
			{
				(*iter)->isChecked = true;

				tree<Interaction*>::post_order_iterator parent = interTree[i]->parent(iter);
				if ( interTree[i]->is_valid(parent) )
				{
					bool similarNode = abs( (*iter)->mergeDist - (*parent)->mergeDist ) < threshold1;
					bool meaninglessNode = (*iter)->mergeDist > threshold2;
					if ( similarNode || meaninglessNode )
					{
						iter = interTree[i]->flatten(iter);
						interTree[i]->erase(iter);
					}
				}
				iter = findFirstUncheckedNodeBottomUp(interTree[i]);
			}
		}

		if (scene->visualizeHierarchy)
		{
			Visualization::visualizeTree( interTree[i], scene->filename + "_" + scene->distPara.toString() + "_tree_merged_" + QString::number(i) );
		}
	}	
}

void InterHierarchy::computeFeatures()
{
	for (int i=0; i<interTree.size(); i++)
	{
		tree<Interaction*>::iterator iter = interTree[i]->begin();
		while (iter != interTree[i]->end())
		{
			if ((*iter)->ibs == NULL)
			{
				(*iter)->computeFeature();
			}
			++iter;
		}
	}
}

void InterHierarchy::addVirtualRoot()
{
	if (interTree.size() == 1 && interTree[0]->size() == 1)				// virtual node
	{
		tree<Interaction*>::iterator newNode = interTree[0]->begin();
		interTree[0]->append_child(newNode, (*newNode));

		if (scene->visualizeHierarchy)
		{
			Visualization::visualizeTree( interTree[0], scene->filename + "_" + scene->distPara.toString() + "_tree");
		}
	}
}


void InterHierarchy::save()
{
	QString resultFolder = QFileInfo(scene->filename).dir().absolutePath() + "/../" + "(ICON)" + scene->distPara.toString() + scene->hierPara.toString();
	QDir dir(resultFolder);
	if (!dir.exists())
	{
		dir.mkpath(resultFolder);
	}

	QFile file(resultFolder + '/' + scene->name + ".icon");
	file.open(QIODevice::WriteOnly | QIODevice::Text);
	QTextStream out(&file);
	
	// find all useful interactions
	unlabelAllInteractions();
	QVector< Interaction* > allInteractions = scene->interactions;
	for (int i=0; i<mergedInteractions.size(); i++)
	{
		if (mergedInteractions[i]->ibs)
		{
			allInteractions.push_back(mergedInteractions[i]);
		}
	}
	for (int i=0; i<allInteractions.size(); i++)
	{
		allInteractions[i]->idx = i;
	}

	// output all the interactions first
	out << "interactionCount " << allInteractions.size() << endl;
	for (int k=0; k<allInteractions.size(); k++)
	{
		out << "interaction " << k << endl;

		out << "origObjIdx";
		for (int i=0; i<allInteractions[k]->obj->origIdx.size(); i++)
		{
			out << " " << allInteractions[k]->obj->origIdx[i];
		}
		out << endl;

		out << "ibs_pfh";
		for (int i=0; i<allInteractions[k]->ibs->pfh.size(); i++)
		{
			out << " " << allInteractions[k]->ibs->pfh[i];
		}
		out << endl;

		out << "ibs_dirHist";
		for (int i=0; i<allInteractions[k]->ibs->dirHist.size(); i++)
		{
			out << " " << allInteractions[k]->ibs->dirHist[i];
		}
		out << endl;

		out << "ibs_distHist";
		for (int i=0; i<allInteractions[k]->ibs->distHist.size(); i++)
		{
			out << " " << allInteractions[k]->ibs->distHist[i];
		}
		out << endl;

		out << "region_pfh";
		for (int i=0; i<allInteractions[k]->region->pfh.size(); i++)
		{
			out << " " << allInteractions[k]->region->pfh[i];
		}
		out << endl;

		out << "region_dirHist";
		for (int i=0; i<allInteractions[k]->region->dirHist.size(); i++)
		{
			out << " " << allInteractions[k]->region->dirHist[i];
		}
		out << endl;

		out << "region_heightHist";
		for (int i=0; i<allInteractions[k]->region->heightHist.size(); i++)
		{
			out << " " << allInteractions[k]->region->heightHist[i];
		}
		out << endl;
	}

	// output all the trees
	out << "treeCount " << interTree.size() << endl;
	for (int i=0; i<interTree.size(); i++)
	{
		out << "tree " << i << endl;
		out << "nodeCount " << interTree[i]->size() << endl;

		int idx = 0;
		tree<Interaction*>::breadth_first_iterator iter = interTree[i]->begin_breadth_first();
		while( iter != interTree[i]->end_breadth_first() )
		{
			out << "node " << idx << endl;

			int parentIdx = -1; 
			if ( idx > 0 )
			{
				parentIdx = ( *(interTree[i]->parent(iter)) )->nodeIdx;
				assert( parentIdx != -1 );
			}
			out << "parent " << parentIdx << endl;	

			out << "inter " << (*iter)->idx << endl;

			(*iter)->nodeIdx = idx;

			iter++;
			idx++;
		}
	}

	file.close();
}

void InterHierarchy::load()
{
	QString resultFolder = QFileInfo(scene->filename).dir().absolutePath() + "/../" + "(ICON)" + scene->distPara.toString() + scene->hierPara.toString();
	QString filename = resultFolder + '/' + scene->name;

	QFile file( filename + ".icon");
	if (file.open(QIODevice::ReadOnly | QIODevice::Text)) 
	{
		int interactionNum = 0;
		int treeNum = 0;
		int nodeNum = 0;
		tree<Interaction*>* newTree;
		QVector<tree<Interaction*>::iterator> nodes;

		QTextStream in(&file);	
		while (!in.atEnd()) {
			QString line = in.readLine();
			QStringList list = line.split(" ");
			if(list.isEmpty()) continue;

			if (list[0] == "interactionCount")
			{
				interactionNum = list[1].toInt();
			}
			else if (list[0] == "interaction")
			{
				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "origObjIdx");
				QVector<int> origObjIdx;
				for (int i=1; i<list.size(); i++)
				{
					origObjIdx << list[i].toInt();
				}

				Interaction* inter = new Interaction(scene, origObjIdx);
				mergedInteractions << inter;

				// load the features
				inter->ibs = new IBS(scene);
				inter->region = new FuncRegion(scene->centralObj);
				inter->region->ibs = inter->ibs;

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "ibs_pfh");
				for (int i=1; i<list.size(); i++)
				{
					inter->ibs->pfh << list[i].toDouble();
				}
				assert(inter->ibs->pfh.size() == 250);

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "ibs_dirHist");
				for (int i=1; i<list.size(); i++)
				{
					inter->ibs->dirHist << list[i].toDouble();
				}
				assert(inter->ibs->dirHist.size() == 10);

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "ibs_distHist");
				for (int i=1; i<list.size(); i++)
				{
					inter->ibs->distHist << list[i].toDouble();
				}
				assert(inter->ibs->distHist.size() == 10);

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "region_pfh");
				for (int i=1; i<list.size(); i++)
				{
					inter->region->pfh << list[i].toDouble();
				}
				assert(inter->region->pfh.size() == 250);

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "region_dirHist");
				for (int i=1; i<list.size(); i++)
				{
					inter->region->dirHist << list[i].toDouble();
				}
				assert(inter->region->dirHist.size() == 10);

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "region_heightHist");
				for (int i=1; i<list.size(); i++)
				{
					inter->region->heightHist << list[i].toDouble();
				}
				assert(inter->region->heightHist.size() == 10);
			}
			else if (list[0] == "treeCount")
			{
				treeNum = list[1].toInt();	
			}
			else if (list[0] == "tree")
			{
				assert(nodeNum == nodes.size());
				newTree = new tree<Interaction*>();
				interTree << newTree;
			}
			else if (list[0] == "nodeCount")
			{
				nodeNum = list[1].toInt();	
				nodes.clear();	
			}
			else if (list[0] == "node")
			{
				// add the new node
				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "parent");
				int parentIdx = list[1].toInt();	

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "inter");
				int interIdx = list[1].toInt();	

				tree<Interaction*>::iterator newNode;
				if (parentIdx == -1)
				{
					newNode = newTree->begin();
					newNode = newTree->insert(newNode, mergedInteractions[interIdx]);
				}
				else
				{
					newNode = newTree->append_child(nodes[parentIdx], mergedInteractions[interIdx]);
				}

				nodes << newNode;								
			}
		}

		assert(interactionNum == mergedInteractions.size());
		assert(treeNum == interTree.size());
		file.close();
	}
}

void InterHierarchy::unlabelAllInteractions()
{
	for (int i=0; i<mergedInteractions.size(); i++)
	{
		mergedInteractions[i]->isChecked = false;
		mergedInteractions[i]->nodeIdx = -1;
		mergedInteractions[i]->idx = -1;
	}
}

tree<Interaction*>::post_order_iterator InterHierarchy::findFirstUncheckedNodeBottomUp(tree<Interaction*>* t)
{
	tree<Interaction*>::post_order_iterator iter = t->begin_post();

	while (iter != t->end_post())
	{
		if (!(*iter)->isChecked && iter.number_of_children() > 0)
		{
			return iter;
		}

		iter++;
	}

	return iter;
}

tree<Interaction*>::breadth_first_iterator InterHierarchy::findFirstUncheckedNodeTopDown( tree<Interaction*>* t )
{
	tree<Interaction*>::breadth_first_iterator iter = t->begin_breadth_first();

	while (iter != t->end_breadth_first())
	{
		if (!(*iter)->isChecked && iter.number_of_children() > 0)
		{
			return iter;
		}

		iter++;
	}

	return iter;
}

//////////////////////////////////////////////////////////////////////////
// for visualizing merged interactions
void InterHierarchy::combineAllInteractionsForCurrentScene()
{
	// combine all useful interactions
	unlabelAllInteractions();
	if (!scene->allInteractions.isEmpty())
		scene->allInteractions.clear();

	scene->allInteractions = scene->interactions;
	for (int i = 0; i < mergedInteractions.size(); ++i)
	{
		if (mergedInteractions[i]->ibs)
		{
			scene->allInteractions.push_back(mergedInteractions[i]);
		}
	}
	for (int i = 0; i < scene->allInteractions.size(); ++i)
	{
		scene->allInteractions[i]->idx = i;
	}
}

// get current hierarchy node
QVector<int> InterHierarchy::getInteractions(int hierIdx)
{
	if (hierIdx < 0 || hierIdx > interTree.size())
	{
		return QVector<int>();
	}

	QVector<int> interactionsIdx;
	tree<Interaction*>::breadth_first_iterator iter = interTree[hierIdx]->begin_breadth_first();
	while( iter != interTree[hierIdx]->end_breadth_first() ) {
		interactionsIdx.push_back((*iter)->idx);
		++iter;
	}

	return interactionsIdx;
}

QVector<int> InterHierarchy::getInteractionsPreOrder(int hierIdx)
{
	if (hierIdx < 0 || hierIdx > interTree.size())
	{
		return QVector<int>();
	}

	QVector<int> interactionsIdx;
	tree<Interaction*>::pre_order_iterator iter = interTree[hierIdx]->begin();
	while (iter != interTree[hierIdx]->end())
	{
		interactionsIdx.push_back((*iter)->idx);
		++iter;
	}
	return interactionsIdx;
}

QVector<int> InterHierarchy::getInterHierParentIdx(int hierIdx)
{
	if (hierIdx < 0 || hierIdx > interTree.size())
	{
		return QVector<int>();
	}

	QVector<int> parentIdx;

	int idx = 0;
	tree<Interaction*>::breadth_first_iterator iter = interTree[hierIdx]->begin_breadth_first();
	while( iter != interTree[hierIdx]->end_breadth_first() )
	{
		int parentNode = -1; 
		if ( idx > 0 )
		{
			parentNode = ( *(interTree[hierIdx]->parent(iter)) )->nodeIdx;
			assert( parentNode != -1 );
		}

		parentIdx.push_back(parentNode);		
		(*iter)->nodeIdx = idx;

		iter++;
		idx++;
	}

	return parentIdx;
}

QVector<int> InterHierarchy::getInterHierParentIdxPreOrder(int hierIdx)
{
	if (hierIdx < 0 || hierIdx > interTree.size())
	{
		return QVector<int>();
	}

	QVector<int> parentIdx;

	int idx = 0;
	tree<Interaction*>::pre_order_iterator iter = interTree[hierIdx]->begin();
	while (iter != interTree[hierIdx]->end())
	{
		int parentNode = -1;
		if (idx > 0)
		{
			parentNode = (*(interTree[hierIdx]->parent(iter)))->nodeIdx;
			//assert( parentNode != -1 );
		}

		parentIdx.push_back(parentNode);
		(*iter)->nodeIdx = idx;

		iter++;
		idx++;
	}

	return parentIdx;
}

int InterHierarchy::getCandiHierNum()
{
	return interTree.size();
}