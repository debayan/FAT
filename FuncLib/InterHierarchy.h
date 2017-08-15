#pragma once

#include "Interaction.h"
#include "tree.h"

class InterHierarchy
{
public:
	InterHierarchy(Scene *s);
	InterHierarchy(QString scene_name);
	InterHierarchy(QString sceneName, QString paraString);
	~InterHierarchy();

public:
	void construct();
	void save();
	void load();

	void computeFeatures();      // compute the features for the intermediate node (clustered interactions) 
	void combineAllInteractionsForCurrentScene();
	QVector<int> getInteractions(int hierIdx);
	QVector<int> getInteractionsPreOrder(int hierIdx);
	QVector<int> getInterHierParentIdx(int hierIdx);
	QVector<int> getInterHierParentIdxPreOrder(int hierIdx);
	int getCandiHierNum();

private:
	void initializeBinaryTree();	// build the binary tree based on the AHC result
	void generateMultipleTrees();	// generate multiple possible tree if 
	void mergeBranches();			// merge branches to convert the binary tree into our general tree	
	void addVirtualRoot();			// if the scene has only one interaction, add a virtual  root

	void deleteMergedInteractions();

private:
	void addChildren(tree<Interaction*>::iterator parentNode,  QVector<Interaction*>& interactions, int idx); 
	void findFlexibleNodes(QVector< tree<Interaction*>::iterator >& nodes,  QVector< QVector<int> >& delCombs);

	void unlabelAllInteractions();
	tree<Interaction*>::post_order_iterator findFirstUncheckedNodeBottomUp(tree<Interaction*>* t);
	tree<Interaction*>::breadth_first_iterator findFirstUncheckedNodeTopDown(tree<Interaction*>* t);

public:
	Scene * scene;
	QVector< tree<Interaction*>* > interTree; 

private:
	QVector< Interaction* > mergedInteractions; // just for storing the merged interaction
};