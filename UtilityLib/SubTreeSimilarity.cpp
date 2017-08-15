#include "SubTreeSimilarity.h"
#include "hungarian/hungarian.h"

template<class TreeNode>
SubTreeSimilarity<TreeNode>::SubTreeSimilarity()
{
	TreeNode tmpN,tmpM,tmpP,tmpQ;
	tmpN = 1;
	tmpN = 2;
	tmpN = 3;
	tmpN = 4;
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
	std::system(command);
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
	return 1.0f;
}

template<class TreeNode>
double SubTreeSimilarity<TreeNode>::MaxSimilarityCommonSubtree()
{
	double maxsim = 0;
	tree_node::iterator rootA = TreeA.begin();
	tree_node::iterator rootB = TreeB.begin();
	tree_node::sibling iterator simA,simB;
	while(rootA!=TreeA.end()) 
	{
		double tmpsim = AnchoredSimilarity(rootA,rootB);
		if(tmpsim > maxsim)
		{
			maxsim = tmpsim;
			simA = TreeA.begin(rootA);
			simB = TreeB.begin(rootB);
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
			simB = TreeB.begin(rootB);
		}
		++rootB;
	}
	CommonTreeA = tree_node::subtree(simA,simA);
	CommonTreeB = tree_node::subtree(simB,simB);
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
				r[i*n + j] = DistAB[i][j];
			}
		}
	}
	else
	{
		for(int i = 0; i < m; i++)
		{
			for(int j = 0; j < n; j++)
			{
				r[i*m + j] = DistAB[j][i];
			}
		}
	}
	hungarian_t prob;
	hungarian_init(&prob,(int*)r,n,m,HUNGARIAN_MAX);
	hungarian_solve(&prob);
	return hungarian_benefit(&prob);
}