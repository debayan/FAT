#pragma once

#include "tree.h"
#include <QVector>
#include "Eigen/Dense"
#include "QtCore/qfile.h"
#include <QtCore/QTextStream>
#include "hungarian/hungarian.h"


template<class TreeNode>
class SubTreeSimilarity
{
protected:
	typedef tree<TreeNode> tree_node;
public:
	SubTreeSimilarity();
	SubTreeSimilarity(tree_node tA,tree_node tB);
	~SubTreeSimilarity(void);

	double NodeSimilarity(TreeNode nA,TreeNode nB);    //node similarity measure
	double MaxSimilarityCommonSubtree();
	double AnchoredSimilarity(typename tree_node::iterator nA, typename tree_node::iterator nB);
	double Assignment(QVector<QVector<double>> DistAB);
	tree_node GetSubtreeA();
	tree_node GetSubtreeB();
private:
	tree_node TreeA;
	tree_node TreeB;
	tree_node CommonTreeA;
	tree_node CommonTreeB;
};

template<class TreeNode>
tree<TreeNode> SubTreeSimilarity<TreeNode>::GetSubtreeA()
{
	return CommonTreeA;
}

template<class TreeNode>
tree<TreeNode> SubTreeSimilarity<TreeNode>::GetSubtreeB()
{
	return CommonTreeB;
}

template<class TreeNode>
SubTreeSimilarity<TreeNode>::SubTreeSimilarity()
{
	//This function is just for testing
	TreeNode tmpN,tmpM,tmpP,tmpQ;
	tmpN = 1.2;
	tmpM = 2.1;
	tmpP = 3.3;
	tmpQ = 4.7;
	///////////////////////////////////////////////////////
	tree_node::iterator topA,secondA,thirdA;
	topA = TreeA.begin();
	secondA = TreeA.insert(topA, tmpN);
	thirdA = TreeA.append_child(secondA,tmpM);
	TreeA.append_child(secondA,tmpP);
	TreeA.append_child(secondA,tmpQ);
	TreeA.append_child(thirdA,tmpM);
	TreeA.append_child(thirdA,tmpQ);

	topA = TreeB.begin();
	secondA = TreeB.insert(topA, tmpM);
	thirdA = TreeB.append_child(secondA,tmpN);
	TreeB.append_child(secondA,tmpN);
	TreeB.append_child(secondA,tmpP);
	TreeB.append_child(thirdA,tmpN);
	TreeB.append_child(thirdA,tmpP);

	/*
	QFile file("a.dot");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
	out<<"digraph G{"<<endl;
	tree_node::sibling_iterator children;
	tree_node::iterator iterator;
	iterator = TreeA.begin();
	while(iterator!=TreeA.end())
	{
		children = TreeA.begin(iterator);
		while(children!=TreeA.end(iterator))
		{
			QString f = QString::number((*iterator));
			f = f + QString::number(tree_node::depth(iterator));
			QString c = QString::number((*children));
			c = c + QString::number(tree_node::depth(children));
			out<<f<<"->"<<c<<";"<<endl;
			++children;
		}
		++iterator;
	}
	out<<"}"<<endl;
	file.close();
	///////////////////////////////////////////////
	QString commandq = "dot -Tjpg a.dot -o a.jpg";
	QByteArray ba = commandq.toLocal8Bit();
	char* command = ba.data();
	std::system(command);*/
}

template<class TreeNode>
SubTreeSimilarity<TreeNode>::SubTreeSimilarity(tree_node tA,tree_node tB)
{
	TreeA = tA;
	TreeB = tB;
}

template<class TreeNode>
SubTreeSimilarity<TreeNode>::~SubTreeSimilarity(void)
{
}

template<class TreeNode>
double SubTreeSimilarity<TreeNode>::NodeSimilarity(TreeNode nA,TreeNode nB)
{
	return 1.0f/(abs(nA-nB)+1);
}

template<class TreeNode>
double SubTreeSimilarity<TreeNode>::MaxSimilarityCommonSubtree()
{
	double maxsim = 0;
	tree_node::iterator rootA = TreeA.begin();
	tree_node::iterator rootB = TreeB.begin();
	tree_node::sibling_iterator simA,simB,simAend,simBend;
	while(rootA!=TreeA.end()) 
	{
		double tmpsim = AnchoredSimilarity(rootA,rootB);
		if(tmpsim > maxsim)
		{
			maxsim = tmpsim;
			simA = TreeA.begin(rootA);
			simAend = TreeA.end(rootA);
			simB = TreeB.begin(rootB);
			simBend = TreeB.end(rootB);
		}
		++rootA;
	}
	rootA = TreeA.begin();
	while(rootB!=TreeB.end()) 
	{
		double tmpsim = AnchoredSimilarity(rootA,rootB);
		if(tmpsim > maxsim)
		{
			maxsim = tmpsim;
			simA = TreeA.begin(rootA);
			simAend = TreeA.end(rootA);
			simB = TreeB.begin(rootB);
			simBend = TreeB.end(rootB);
		}
		++rootB;
	}
	CommonTreeA = TreeA.subtree(simA,simAend);
	CommonTreeB = TreeA.subtree(simB,simBend);
	return maxsim;
}

template<class TreeNode>
double SubTreeSimilarity<TreeNode>::AnchoredSimilarity(typename tree_node::iterator iteratorA, typename tree_node::iterator iteratorB)
{	
	if(tree_node::number_of_children(iteratorA)==0||tree_node::number_of_children(iteratorB)==0)
		return NodeSimilarity(*iteratorA,*iteratorB);

	tree_node::iterator nodeAbegin = TreeA.begin(iteratorA); // return the children of iterationA or it's siblings?
	tree_node::iterator nodeBbegin = TreeB.begin(iteratorB);
	tree_node::iterator nodeAend = TreeA.end(iteratorA);
	tree_node::iterator nodeBend = TreeB.end(iteratorB);

	int dA = TreeA.depth(nodeAbegin);
	int dB = TreeB.depth(nodeBbegin);

	QVector<tree_node::iterator> ChildrenB;
	while(nodeBbegin!=nodeBend) {
     if(TreeB.depth(nodeBbegin)>dB) // it's sibling + children? why the depth will be different?
	 {
		 ++nodeBbegin;
		 continue;
	 }
	 ChildrenB.push_back(nodeBbegin);
     ++nodeBbegin;
     }

	QVector<QVector<double>> DistAB;
	DistAB.resize(tree_node::number_of_children(iteratorA));
	int count = 0;
	while(nodeAbegin!=nodeAend) {
     if(TreeA.depth(nodeAbegin)>dA)
	 {
		 ++nodeAbegin;
		 continue;
	 }
	 for each(tree_node::iterator tmpn in ChildrenB)
	 {
		 double disttmp = AnchoredSimilarity(nodeAbegin,tmpn);
		 DistAB[count].push_back(disttmp);
	 }
     ++nodeAbegin;
	 count++;
     }
	return NodeSimilarity(*iteratorA,*iteratorB) + Assignment(DistAB);
}

template<class TreeNode>
double SubTreeSimilarity<TreeNode>::Assignment(QVector<QVector<double>> DistAB)
{
	int n = DistAB.size();
	int m = DistAB[0].size();
	QVector<int> index;

	int *r = new int[m*n];
	if(n <= m)
	{
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < m; j++)
			{
				r[i*n + j] = DistAB[i][j]*100;
			}
		}
	}
	else
	{
		for(int i = 0; i < m; i++)
		{
			for(int j = 0; j < n; j++)
			{
				r[i*m + j] = DistAB[j][i]*100;
			}
		}
	}
	hungarian_t prob;
	if(n <= m)
		hungarian_init(&prob,(int*)r,n,m,HUNGARIAN_MAX);
	else
		hungarian_init(&prob,(int*)r,m,n,HUNGARIAN_MAX);
	hungarian_solve(&prob);
	return hungarian_benefit(&prob);
}

