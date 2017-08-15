#pragma once

#include "interHierarchy.h"
#include <QPair>

class CommonSubtree 
{

typedef tree<Interaction*> tree_inter;
typedef tree_inter::iterator node;

public:
	CommonSubtree(tree_inter* tree1, tree_inter* tree2);
	~CommonSubtree();

public:
	void findMaxMatch();
	double distance();

private:
	CommonSubtree findAnchoredMaxMatch(node n1, node n2);
	QVector<int> findBestAssignment(Eigen::MatrixXd S);

public:
	tree_inter* tree1;
	tree_inter* tree2;
	QVector< QPair<node, node> > mapping;
	QVector< double > pairSimilarity;
	double similarity;
	double totalWeight;
};
