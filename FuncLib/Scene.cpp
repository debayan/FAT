#include "Scene.h"
#include "Interaction.h"
#include "IbsGenerator.h"
#include "primitives.h"
#include <QFileInfo>
#include "SegMeshLoader.h"
#include "Colormap.h"
#include <QTime>
#include "DetectCollision.h"
#include "HierarchicalCluster.h"
#include "DistMeasure.h"
#include "InterSet.h"
#include "writeOBJ.h"

#define PQP_SUPPORT_SURFACEMESH
#include "PQPLib.h"

Scene::Scene()
{
	upright = Vec3d(0,0,1);					// initialized when created
	centralObj = NULL;
	allObj = NULL;

	clusterType = AP;
	selectedClusterIdx = -1;

	ahcreport = NULL;
	hierarchy = NULL;
	visualizeHierarchy = false;

	sceneHierarchy = NULL;
	baseColor = QColor(235,235,235,255);
}

Scene::~Scene()
{
	clearAllObjects();
}

void Scene::clearAllObjects()
{
	for(auto& inter:interactions)
	{
		if (inter)
		{
			delete inter;
			inter = NULL;
		}
	}
	interactions.clear();

	for (auto& cluster:clusteredInteractions)
	{
		if (cluster)
		{
			delete cluster;
			cluster = NULL;
		}
	}
	clusteredInteractions.clear();

	if (ahcreport)
	{
		delete ahcreport;
	}

	if (hierarchy)
	{
		delete hierarchy;
	}

	foreach(Object * obj, objects)
	{
		if (obj)
		{
			delete obj;
			obj = NULL;
		}
	}
	objects.clear();

	if (centralObj)
	{
		delete centralObj;
	}

	foreach(IBS * ibs, ibsSetScene)
	{
		if (ibs)
		{
			delete ibs;
			ibs = NULL;
		}
	}
	ibsSetScene.clear();

	if (sceneHierarchy)
	{
		delete sceneHierarchy;
	}
}

void Scene::initialize()
{
	// load snapshot
	loadSnapshot();

	// load selected object idx
	loadCentralIdx();

	// update bbox for rendering
	bbox = objects[0]->bbox;
	for (int i=1; i<objects.size(); i++)
	{
		bbox.extend(objects[i]->bbox);
	}

	// all the objects don't have corresponding interactions yet
	obj2Inter.clear();
	for (int i=0; i<objects.size(); i++)
	{
		obj2Inter.push_back(-1);
	}

	colorMap = rndColors2(objects.size()*10);
	for (int i=0; i<objects.size(); i++)
	{
		objects[i]->color = colorMap[i];
	}
}

void Scene::loadCentralIdx()
{
	QFile file(filename + ".ci");
	if (file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		QTextStream in(&file);	
		QString line = in.readLine();
		if (line == "x")
		{
			upright = Vec3d(1,0,0);
		}
		else if (line == "y")
		{
			upright = Vec3d(0,1,0);
		}
		else if (line == "z")
		{
			upright = Vec3d(0,0,1);
		}

		while (!in.atEnd()) {
			QString line = in.readLine();
			int idx = line.toInt();
			if (idx < 0 || idx >= objects.size())
			{
				continue;
			}
			objects[idx]->setCentral(true);
			objects[idx]->loadSymGroups();
		}	
		file.close();
	}
}

void Scene::saveCentralIdx()
{
	QFile file(filename + ".ci");

	if (file.open(QIODevice::WriteOnly | QIODevice::Text))
	{
		QTextStream out(&file);	
		if ( abs(upright[0] - 1) < 1.0e-6)
		{
			out << "x\n";
		}
		else if ( abs(upright[1] - 1) < 1.0e-6)
		{
			out << "y\n";
		}
		else if (abs(upright[2] - 1) < 1.0e-6)
		{
			out << "z\n";
		}


		for (int i=0; i<objects.size(); i++)
		{
			if (objects[i]->isCentral)
			{
				out << i << "\n";
			}
		}
	}
}

void Scene::load( QString filename )
{
	QFileInfo file(filename);
	QString extension = file.suffix();
	
	if (extension.toLower() == "txt")
		loadFromTxt(filename);				// from *.txt file
	else if (extension.toLower() == "obj")
		loadFromObj(filename);				// from SegMeshLoader

	this->name = file.baseName();
	this->filename = file.absolutePath() + "/" + name;

	initialize();
}

void Scene::loadSnapshot()
{
	snapshot = QPixmap(filename + ".png");
}

void Scene::loadFromTxt( QString filename )
{
	// open file to read
	QFile file(filename);
	if (file.open(QIODevice::ReadOnly | QIODevice::Text)) 
	{
		QFileInfo fileInfo(filename);
		QString dirName = fileInfo.dir().absolutePath();
		QString modelDir = dirName + "/../../models/";

		clearAllObjects();	
		QTextStream in(&file);	
		while (!in.atEnd()) {
			QString line = in.readLine();
			QStringList list = line.split(" ");
			if(list.isEmpty()) continue;

			if (list[0] == "modelCount")
			{
				//objects.resize(list[1].toInt()-1);
				continue;
			}
			else if (list[0] == "newModel")
			{
				int objectIdx = list[1].toInt() - 1;
				if (objectIdx == -1)
				{
					continue;
				}
				else
				{		
					QString modelName = modelDir + list[2] + ".obj";
					if (!QFile(modelName).exists())
					{
						modelName = modelDir + list[2] + ".off";
					}
					if (!QFile(modelName).exists())
					{
						debugBox("There is no existing model named " + list[2] + "!");
						continue;
					}

					Object * obj = new Object(this);
					obj->load(modelName);
					obj->origIdx.push_back(objects.size());		// object index ordered by .txt file
					objects.push_back(obj);
				}
			}
		}	
		file.close();
	}
}

void Scene::loadFromObj( QString filename )
{
	SurfaceMeshModel * entireMesh = new SurfaceMeshModel(filename);
	SegMeshLoader sml(entireMesh);
	QVector<SurfaceMeshModel*> models = sml.getSegMeshes();

	objects.clear();
	for (auto m:models)
	{
		Object * obj = new Object(this);
		obj->setMesh(m);
		obj->origIdx.push_back(objects.size());

		objects.push_back(obj);
	}

	delete entireMesh;
}

void Scene::draw(DrawParameter para)
{
	if (centralObj == NULL)
	{
		for (auto obj:objects)
		{
			obj->draw(para.objMode);
		}
	}
	else
	{
		if (para.drawInteraction)
		{
			centralObj->drawMesh(baseColor);

			for (int i = 0; i < interactionsIdx.size(); ++i)
			{
				if (isSelectedInterGroup[i])
				{
					int iterIdx = interactionsIdx[i];
					QVector<int> interObjGroup = allInteractions[iterIdx]->obj->origIdx;
					QColor color = colorMap[i];

					if (para.drawInteractionObj) {					// drawing interacting object
						for (int j = 0; j < interObjGroup.size(); ++j)
						{
							objects[interObjGroup[j]]->drawMesh(color);
						}
					}

					if (para.drawInteractionIBS) {					// draw IBS
						allInteractions[iterIdx]->ibs->draw(para.drawIbsSample, para.drawIbsWeight, color);
					}

					if (para.drawInteractionIR) {					// draw IR
						centralObj->drawIR(allInteractions[iterIdx]->region, color);
					}
				}
			}
		}
		else
		{
			for (auto obj:objects)
			{
				if (!obj->isCentral)
				{
					obj->draw(para.objMode);
				}
			}
			centralObj->draw(para.objMode);

			for (int i=0; i<interObjIdx.size(); i++)
			{
				if (isSeletedInteraction[i])
				{
					for (auto ibsIdx : interIbsSet[i])
					{
						if ( ibsIdx < ibsSetScene.size() )
						{
							ibsSetScene[ibsIdx]->draw();
						}	
					}
				}		
			}
		}
	}
}


// Ray tracking based object selection
int Scene::getSelectedObjectID( Eigen::Vector3d orig, Eigen::Vector3d dir )
{
    int selectedID = -1;

	Ray ray(orig, dir);
	double dist = DBL_MAX;

	for(int i=0; i<objects.size(); i++)
	{
		Eigen::AlignedBox3d aabb = objects[i]->bbox;
		Eigen::Vector3d extent = (aabb.max() - aabb.min()) / 2;
        BoundingBox box(aabb.center(), extent[0], extent[1], extent[2]);

		double t;
		if (box.intersectsWorking(ray, t) && t < dist)
		{
			dist = t;
			selectedID = i;
		}
	}
	return selectedID;
}

void Scene::setObjCentralState( QVector<bool> sel )
{
	for (int i=0; i<objects.size(); i++)
	{
		objects[i]->setCentral(sel[i]);
	}

	updateInteractions();
}

void Scene::reverseObjCentralState( int objIdx )
{
	objects[objIdx]->reverseCentralState();

	updateInteractions();
}

bool Scene::hasCentralObj()
{
	for (auto obj:objects)
	{
		if (obj->isCentral)
		{
			return true;
		}
	}

	return false;
}

void Scene::setInteractionSelectState( QVector<bool> sel )
{
	isSeletedInteraction = QVector<bool>(interObjIdx.size(), false);
	QVector<int> selectedInters;
	for (int i = 0; i < sel.size(); ++i)
	{
		if (sel[i]) {
			int interIdx = interactionsIdx[i];
			QVector<int> interObjsIdx = allInteractions[interIdx]->obj->origIdx;
			selectedInters << interObjsIdx;
		}
	}

	for (int i = 0; i < interObjIdx.size(); ++i) {
		if (selectedInters.contains(interObjIdx[i]))
			isSeletedInteraction[i] = true;
	}	
}

void Scene::setInteractionHierSelectState(QVector<bool> sel)
{
	assert(sel.size() == lastSelectedInters.size());
	
	int currSelectedIdx = -1;
	if (lastSelectedInters != sel)
	{
		for (int i = 0; i < sel.size(); ++i) {
			if (lastSelectedInters[i] != sel[i])
			{
				currSelectedIdx = i;
				break;
			}
		}
	}

	if (currSelectedIdx != -1)
	{
		// un-check children
		for (int i = 0; i < interactionsIdx.size(); ++i)
		{
			int parentIdx = intersParentIdx[i];
			while (parentIdx != -1) {				
				if (parentIdx == currSelectedIdx)
				{
					sel[i] = false;
					break;
				}
				parentIdx = intersParentIdx[parentIdx];
			}
		}

		// un-check parents
		int parentIdx = intersParentIdx[currSelectedIdx];
		while (parentIdx != -1) 
		{
			if (sel[parentIdx])
			{
				sel[parentIdx] = false;
			}
			parentIdx = intersParentIdx[parentIdx];
		}
	}

	lastSelectedInters = sel;
	isSelectedInterGroup = sel;
}

void Scene::resetInteractionSelectState()
{
	if (!lastSelectedInters.isEmpty())
	{
		for (int i = 0; i < lastSelectedInters.size(); ++i)
			lastSelectedInters[i] = false;
	}

	if (!isSelectedInterGroup.isEmpty())
	{
		for (int i = 0; i < isSelectedInterGroup.size(); ++i)
			isSelectedInterGroup[i] = false;
	}

	if (!isSeletedInteraction.isEmpty())
	{
		for (int i = 0; i < isSeletedInteraction.size(); ++i)
			isSeletedInteraction[i] = false;
	}
}

//////////////////////////////////////////////////////////////////////////
// find the interaction to the central object, 
// and represent each using IBS and interacting surface region
void Scene::computeInteractions()
{
	// compute IBS for the entire scene to get the interactions related to the central object
	analyzeInteractions();				// computing IBS

	// combine all the selected object into one central object
	generateCentralObject();			// including combing symmetry information

	// generate interaction between the central object and each interacting object
	constructInteractions();

	computeMaxHistForInteraction();		// for illustrating the features
	centralObj->prepareRegionDraw();	// for illustrating the contacting region

	// detect the symmetry pattern each region in
	//centralObj->identifyPatterEachRegionIn();			// TBD
}

void Scene::analyzeInteractions()
{
	for (auto& ibs : ibsSetScene) {
		if (ibs) {
			delete ibs;
			ibs = NULL;
		}
	}
	ibsSetScene.clear();

	//if (ibsSetScene.isEmpty())
	{
		samplePoints();

		IbsGenerator generator;
		ibsSetScene = generator.computeIBS(this, objects);

		// map each object pair to ibs
		objPair2IbsSet = Eigen::MatrixXi::Ones(objects.size(), objects.size()) * (-1);
		for (int i=0; i<ibsSetScene.size(); i++)
		{
			int objIdx1 = ibsSetScene[i]->obj1->origIdx[0];
			int objIdx2 = ibsSetScene[i]->obj2->origIdx[0];
			objPair2IbsSet( objIdx1, objIdx2 ) = i;
			objPair2IbsSet( objIdx2, objIdx1 ) = i;
		}
		
		updateInteractions();			// update interaction index only for interacting objects
	}
}

void Scene::samplePoints()
{
	if (objects.isEmpty())
	{
		return;
	}

	// sample based on surface area
	for (auto obj:objects)
	{
		if (obj->surfaceArea == 0)
		{
			obj->computeArea();
		}
	}

	double minArea = objects[0]->surfaceArea;
	double totalArea = 0;
	for (auto obj:objects)
	{
		if (obj->surfaceArea < minArea)
		{
			minArea = obj->surfaceArea;
		}
		totalArea += obj->surfaceArea;
	}

	int sampleNum = 2000 * objects.size();
	for (auto obj:objects)
	{
		int num = sampleNum * obj->surfaceArea / totalArea;
		if (num < 1000)
		{
			num = 1000;
		}
		if (obj->isCentral)
		{
			num = num * 2;
		}
		obj->sampling(num);
	}
}

void Scene::updateInteractions()
{
	for (int i=0; i<objects.size(); i++)
	{
		objects[i]->isInteracting = false;
		obj2Inter[i] = -1;
	}

	QMap<int, QVector<int>> interObj2IbsSet;
	for (int i=0; i<ibsSetScene.size(); i++)
	{
		if (ibsSetScene[i]->obj1->isCentral != ibsSetScene[i]->obj2->isCentral)
		{	
			int iIdx = ibsSetScene[i]->obj2->origIdx[0];
			ibsSetScene[i]->pointToCentralObject = false;

			if (ibsSetScene[i]->obj2->isCentral)
			{
				iIdx = ibsSetScene[i]->obj1->origIdx[0];
				ibsSetScene[i]->pointToCentralObject = true;
			}

			QMap<int, QVector<int>>::iterator iter = interObj2IbsSet.find(iIdx);
			if (iter != interObj2IbsSet.end())
			{
				iter.value().push_back(i);
			}
			else
			{
				QVector<int> newIbsSet;
				newIbsSet.push_back(i);
				interObj2IbsSet.insert(iIdx, newIbsSet);
				objects[iIdx]->isInteracting = true;
			}
		}
	}

	interObjIdx.clear();
	interIbsSet.clear();
	isSeletedInteraction.clear();
	QMap<int, QVector<int>>::const_iterator i = interObj2IbsSet.constBegin();
	while (i != interObj2IbsSet.constEnd()) {
		obj2Inter[i.key()] = interObjIdx.size();		// indexing interObjIdx (obj2Inter -> interIbsSet)
		interObjIdx.push_back(i.key());
		interIbsSet.push_back(i.value());
		isSeletedInteraction.push_back(false);
		++i;
	}
}

void Scene::generateCentralObject()
{
	saveCentralIdx();					// save *.ci

	if (centralObj != NULL)
	{
		delete centralObj;
	}

	// generate a new object to represent the central shape
	QVector<Object*> combinedObj;
	for (auto obj:objects)
	{
		if (obj->isCentral)
		{
			combinedObj.push_back(obj);
		}
	}

	centralObj = new Object(this);
	centralObj->combineObjects(combinedObj);
	centralObj->isCentral = true;
}

void Scene::constructInteractions()
{
	for ( auto& inter:interactions)
	{
		if (inter)
		{
			delete inter;
			inter = NULL;
		}
	}
	interactions.clear();
	interactionColors.clear();

	for (int i=0; i<interObjIdx.size(); i++)
	{
		Interaction *  newInteraction = new Interaction(this, objects[interObjIdx[i]]);
		interactions.push_back(newInteraction);
		interactionColors.push_back(colorMap[i]);
	}
}

void Scene::computeMaxHistForInteraction()
{
	yMaxHistInteraction = computeMaxHist(interactions);
}

QVector<double> Scene::computeMaxHist( QVector<Interaction*> inters )
{
	QVector<double> yMaxHist = QVector<double>(5,0);
	for (auto inter:inters)
	{
		for (auto v:inter->ibs->pfh)
		{
			yMaxHist[0] = (yMaxHist[0]>v) ? yMaxHist[0]:v;
		}

		for (auto v:inter->ibs->dirHist)
		{
			yMaxHist[1] = (yMaxHist[1]>v) ? yMaxHist[1]:v;
		}

		for (auto v:inter->ibs->distHist)
		{
			yMaxHist[2] = (yMaxHist[2]>v) ? yMaxHist[2]:v;
		}

		for (auto v:inter->region->pfh)
		{
			yMaxHist[3] = (yMaxHist[3]>v) ? yMaxHist[3]:v;
		}

		for (auto v:inter->region->dirHist)
		{
			yMaxHist[4] = (yMaxHist[4]>v) ? yMaxHist[4]:v;
		}
	}

	return yMaxHist;
}

//////////////////////////////////////////////////////////////////////////
// compute the functional descriptor: interaction hierarchy

void Scene::buildInteractionHierarchy()
{
	loadInteractionHierarchy();
}

alglib::ahcreport* Scene::getAHC()
{
	if (!ahcreport && !interactions.isEmpty())
	{
		DistMeasure dist(distPara);
		Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(interactions.size(), interactions.size());
		for (int i=0; i<interactions.size()-1; i++)
		{
			for (int j=i+1; j<interactions.size(); j++)
			{
				distance(i, j) = dist.betweenInteractions(interactions[i], interactions[j], true);
				distance(j, i) = distance(i, j);
			}
		}

		// build the binary hierarchy based on the geometry distance measure
		HierarchicalCluster ahc;
		ahcreport = ahc.clustering(distance); // I cannot delete centralObj->ahcreport, why?
	}

	return ahcreport;
}

void Scene::loadInteractionHierarchy()
{
	if (hierarchy)
	{
		delete hierarchy;
	}
	hierarchy = new InterHierarchy(this);

	resultFolder = QFileInfo(filename).dir().absolutePath() + "/../" + "(ICON)" + distPara.toString() + hierPara.toString();
	QString iconFile = resultFolder + '/' + name;

	QDir dir(resultFolder);
	if (!dir.exists())
	{
		dir.mkpath(resultFolder);
	}

	constructionTime = 0;
	interactionTime = 0;

	QTime t0, t1;
	t0.start();
	

	QFile file(iconFile+".icon");

	if (interactions.isEmpty())
	{
		t1.start();
		computeInteractions();
		interactionTime += t1.elapsed();
	}
	hierarchy->construct();

	t1.start();
	hierarchy->computeFeatures();
	interactionTime += t1.elapsed();

	hierarchy->save();

	constructionTime = t0.elapsed();

	// merge all interactions into scene->allInteractions
	hierarchy->combineAllInteractionsForCurrentScene();
}

void Scene::updateResultflodername()
{
	resultFolder = QFileInfo(filename).dir().absolutePath() + "/../" + distPara.toString() + hierPara.toString();
}

//////////////////////////////////////////////////////////////////////////
void Scene::computeIBSH()
{
	QVector<int> centralObjIdx;
	for (int i=0; i<objects.size(); i++)
	{
		if (objects[i]->isCentral)
		{
			centralObjIdx << i;
		}
	}
	if (centralObjIdx.isEmpty())
	{
		debugBox("Please select central object!");
		return;
	}
	if (centralObjIdx.size() > 1)
	{
		debugBox("Please select one central object only!");
		return;
	}
	
	detectObjContact();
	analyzeInteractions();				// compute ibs

	if (sceneHierarchy)
	{
		delete sceneHierarchy;
	}
	sceneHierarchy = new SceneHierarchy(this);
	sceneHierarchy->construct();
	sceneHierarchy->computeFeature(centralObjIdx[0]);		// centralObjIdx[0] refers to the central object
}

void Scene::outputInterSetFeature()
{
	// compute IBS for the entire scene to get the interactions related to the central object
	analyzeInteractions();

	// combine all the selected object into one central object
	generateCentralObject();

	// generate interaction between the central object and each interacting object
	constructInteractions();

	// store the set features
	InterSet interSet;
	interSet.interactions = interactions;

	resultFolder = QFileInfo(filename).dir().absolutePath() + "/../" + "(ISET)" + distPara.toString();
	QDir dir(resultFolder);
	if (!dir.exists())
	{
		dir.mkpath(resultFolder);
	}

	interSet.save(resultFolder + "/" + name);	
}

void Scene::detectObjContact()
{
	objPairContact = Eigen::MatrixXi::Zero(objects.size(), objects.size());

	// add all the involving meshes
	PQP::Manager manager(objects.size());
	for (auto obj : objects)
	{
		manager.addModel( makeModelPQP( obj->origMesh ) );
	}

	// detect interaction between objects
	double threshold = bbox.diagonal().norm() * 0.02;
	for (int i=0; i<objects.size(); i++)
	{
		for (int j=i+1; j<objects.size(); j++)
		{
			double dist = manager.getDistance(i, j); 
			if (dist < threshold)
			{
				objPairContact(i,j) = 1;
				objPairContact(j,i) = 1;
			}
		}		
	}
}

//////////////////////////////////////////////////////////////////////////
// compute geometry feature for central object
void Scene::computeGeometryFeature(int barNum)
{
	if (!centralObj)
	{
		generateCentralObject();
	}
	if (!allObj)
	{
		QVector<Object*> combinedObj;
		for (auto obj:objects)
			combinedObj.push_back(obj);

		allObj = new Object(this);
		allObj->combineObjects(combinedObj);
	}

	allObj->computeBS();
	centralObj->computeGeoFeature(barNum);
}

//////////////////////////////////////////////////////////////////////////
// get current hierarchy
QVector<int> Scene::getInteractions(int hierIdx)
{
	if (!hierarchy || hierIdx < 0 || hierIdx > getCandiHierNum())
		return QVector<int>();

	return hierarchy->getInteractions(hierIdx);
}

QVector<int> Scene::getInteractionsPreOrder(int hierIdx)
{
	if (!hierarchy || hierIdx < 0 || hierIdx > getCandiHierNum())
		return QVector<int>();

	return hierarchy->getInteractionsPreOrder(hierIdx);
}

QVector<int> Scene::getInterHierParentNodeIdx(int hierIdx)
{
	if (!hierarchy || hierIdx < 0 || hierIdx > getCandiHierNum())
		return QVector<int>();

	return hierarchy->getInterHierParentIdx(hierIdx);
}

QVector<int> Scene::getInterHierParentNodeIdxPreOrder(int hierIdx)
{
	if (!hierarchy || hierIdx < 0 || hierIdx > getCandiHierNum())
		return QVector<int>();

	return hierarchy->getInterHierParentIdxPreOrder(hierIdx);
}

int Scene::getCandiHierNum()
{
	if (hierarchy)
		return hierarchy->getCandiHierNum();
	else
		return 0;
}

void Scene::outputCentricObject()
{
	generateCentralObject();
	centralObj->origMesh->update_vertex_normals();
	std::string tmp = filename.toLocal8Bit().constData();
	
	QString centricObjFile = filename + "_centric.obj";
	writeOBJ::wirte(centralObj->origMesh, centricObjFile);
}