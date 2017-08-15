#include "CommonSubtree.h"
#include "DistMeasure.h"
#include "hungarian/hungarian.h"
#include "UtilityGlobal.h"

CommonSubtree::CommonSubtree( tree_inter* t1, tree_inter* t2 )
{
	tree1 = t1;
	tree2 = t2;

	mapping.clear();
	similarity = 0;

	totalWeight = 0;

	node n1 = tree1->begin();
	while( n1 != tree1->end()) 
	{
		totalWeight += pow(0.5, tree1->depth(n1));
		++n1;
	}

	node n2 = tree2->begin();
	while( n2 != tree2->end()) 
	{
		totalWeight += pow(0.5, tree2->depth(n2));
		++n2;
	}
}

CommonSubtree::~CommonSubtree()
{
}

void CommonSubtree::findMaxMatch()
{
	node root1 = tree1->begin();
	node root2 = tree2->begin();

	node n1 = root1;
	while( n1 != tree1->end()) 
	{
		CommonSubtree tmpSubtree = findAnchoredMaxMatch(n1, root2);
		if(tmpSubtree.similarity > similarity)
		{
			similarity = tmpSubtree.similarity;
			mapping = tmpSubtree.mapping;
			pairSimilarity = tmpSubtree.pairSimilarity;
		}
		++n1;
	}
	
	node n2 = root2;
	while( n2 != tree2->end()) 
	{
		CommonSubtree tmpSubtree = findAnchoredMaxMatch(root1, n2);
		if(tmpSubtree.similarity > similarity)
		{
			similarity = tmpSubtree.similarity;
			mapping = tmpSubtree.mapping;
			pairSimilarity = tmpSubtree.pairSimilarity;
		}
		++n2;
	}
}

CommonSubtree CommonSubtree::findAnchoredMaxMatch( node node1, node node2 )
{
	//// show the anchored nodes
	//QString nodeInformation = "find AnchoredMaxMatch for (";
	//for (int i=0; i<(*node1)->obj->origIdx.size()-1; i++)
	//{
	//	nodeInformation += QString::number((*node1)->obj->origIdx[i]) + ", ";
	//}
	//nodeInformation += QString::number((*node1)->obj->origIdx.last()) + ") and (";
	//for (int i=0; i<(*node2)->obj->origIdx.size()-1; i++)
	//{
	//	nodeInformation += QString::number((*node2)->obj->origIdx[i]) + ", ";
	//}
	//nodeInformation += QString::number((*node2)->obj->origIdx.last()) + ").";

	//qDebug() << nodeInformation;

	// start searching
	CommonSubtree subtree(tree1, tree2);
	DistParameter distPara = (*node1)->scene->distPara;
	DistMeasure dist(distPara);
	double weight = (pow(0.5, tree1->depth(node1)) + pow(0.5, tree2->depth(node2))) / 2.0;
	subtree.pairSimilarity.push_back(1 - dist.betweenInteractions(*node1, *node2));
	subtree.similarity = subtree.pairSimilarity.last()*weight;
	subtree.mapping.push_back(QPair<node, node>(node1, node2));

	int m = node1.number_of_children();
	int n = node2.number_of_children();
	if ( m != 0  && n != 0)
	{
		// find best matching for the children
		Eigen::MatrixXd S = Eigen::MatrixXd::Zero(m, n);
		QVector< QVector< QVector<QPair<node, node>> > >  M; // store all the mapping
		QVector< QVector< QVector<double> > > PS;			 // store all the pair similarity
 
		tree_inter::sibling_iterator child1 = node1.begin();
		while ( child1 != node1.end() )
		{
			int i = tree1->index(child1);
			QVector< QVector<QPair<node, node>> > M1;
			QVector< QVector<double> > PS1;

			tree_inter::sibling_iterator child2 = node2.begin();
			while ( child2 != node2.end() )
			{
				int j = tree2->index(child2);
				CommonSubtree temp = findAnchoredMaxMatch(child1, child2);		// recursively

				S(i,j) = temp.similarity;
				M1.push_back(temp.mapping);
				PS1.push_back(temp.pairSimilarity);
				child2++;
			}

			M.push_back(M1);
			PS.push_back(PS1);
			child1++;
		}

		QVector<int> assignment = findBestAssignment(S);

		// add the matching
		for (int i=0; i<assignment.size(); i++)
		{
			if (assignment[i] != -1)
			{
				subtree.mapping << M[i][assignment[i]];
				subtree.pairSimilarity << PS[i][assignment[i]];
				subtree.similarity += S(i,assignment[i]);
			}
		}
	}

	return subtree;
}

QVector<int> CommonSubtree::findBestAssignment( Eigen::MatrixXd S )
{
	int m = S.rows();
	int n = S.cols();

	Eigen::MatrixXi similarity = Eigen::MatrixXi::Zero(m, n);
	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
		{
			similarity(i,j) = S(i,j) * 10000; 
		}
	}

	QVector<int> assignment(m, -1);

	hungarian_t prob;
	if (m <= n)
	{
		Eigen::MatrixXi similarityT = similarity.transpose();

		hungarian_init(&prob, similarityT.data(), m, n, HUNGARIAN_MAX); 
		hungarian_solve(&prob);

		for (int i=0; i<m; i++)
		{
			assignment[i] = prob.a[i];
		}
	}
	else
	{
		hungarian_init(&prob, similarity.data(), n, m, HUNGARIAN_MAX);
		hungarian_solve(&prob);

		for (int i=0; i<n; i++)
		{
			assignment[prob.a[i]] = i;
		}
	}

	return assignment;
}

double CommonSubtree::distance()
{
	double d = 1.0 - similarity / (totalWeight - similarity);

	return d;
}
