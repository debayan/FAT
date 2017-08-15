#include "Interaction.h"
#include "Scene.h"
#include "IbsGenerator.h"
#include "UtilityGlobal.h"
#include "DetectCollision.h"

Interaction::Interaction(Scene* s, Object* object)
{
	initialize();
	scene = s;
	type = SINGLE;
	obj = object;

	construct();
}

Interaction::Interaction( Scene* s, QVector<Object*> objects)
{
	initialize();
	scene = s;
	type = SINGLE;

	if (objects.size() == 1)
	{
		obj = objects[0];
	}
	else
	{
		obj = new Object(s);
		obj->combineObjects(objects);
	}

	construct();
}

Interaction::Interaction( Scene* s, QVector<Interaction*> interactions, INTERACTION_TYPE t)
{
	initialize();

	scene = s;
	type =  t;

	switch (type)
	{
	case COMBINED_INIT:
		storeOrigObj(interactions);
		break;
	case COMBINED_FULL:
		storeOrigObj(interactions);
		computeFeature();
		break;
	case SINGLE:
		constructNewInteraction(interactions);
		break;
	default:
		break;
	}
}

Interaction::Interaction( Scene* s, QVector<int> orgObjIdx )
{
	initialize();
	scene = s;	
	type = COMBINED_INIT;

	obj = new Object(s);
	obj->origIdx = orgObjIdx;
}

void Interaction::initialize()
{
	scene = NULL;
	obj = NULL;
	ibs = NULL;
	region = NULL;
	type = EMPTY;

	mergeDist = 0;
	isChecked = false;
	nodeIdx = -1;
	idx = -1;
}

Interaction::~Interaction()
{
	if (obj && obj->origIdx.size() > 1)
	{
		delete obj;
	}

	if (ibs)
	{
		delete ibs;
	}

	if (region)
	{
		delete region;
	}
}

void Interaction::construct()
{
	isCollisionFree = false;

	if (type != COMBINED_INIT)
	{
		computeIBS();					// compute IBS features
		findSurfaceRegion();			// construct IR and compute IR features
	}
}

void Interaction::computeIBS()
{
	computePairwiseIBS();

	if (type == SINGLE)					// refine the ibs by sampling more points only for the initial individual interaction
	{
		refinePairwiseIBS();
	}

	computeIbsFeature();
}

void Interaction::computePairwiseIBS()
{
	QVector<Object*> activeObj;
	activeObj.push_back(obj);
	activeObj.push_back(scene->centralObj);

	IbsGenerator generator;
	QVector<IBS*> newIbsSet = generator.computeIBS(scene, activeObj);

	if (newIbsSet.size() > 1)
	{
		debugBox("ERROR: find more than one IBS between two objects!") ;
	}
	else if ( newIbsSet[0]->obj1 != activeObj[0] || newIbsSet[0]->obj2 != activeObj[1])
	{
		debugBox( "ERROR: The object is wrong!");
	}
	newIbsSet[0]->pointToCentralObject = true;		// newIbsSet[0]->obj2 is central object

	if (ibs)
	{
		delete ibs;
	}
	ibs = newIbsSet[0];	
}

void Interaction::refinePairwiseIBS()
{
	detectCollision();			// refinement based on collision detection

	if ( !isCollisionFree )		// only recompute the IBS when there is intersection/collision
	{
		// sample more points on the intersecting objects
		// collisionFaceIdx[0] represents collision between ibs->mesh and ibs->obj1->mesh
		// collisionFaceIdx[1] represents collision between ibs->mesh and ibs->obj2->mesh
		if ( !collisionFaceIdx[0].isEmpty() )
		{
			obj->sampleMoreTri(collisionFaceIdx[0]);
		}	
		if ( !collisionFaceIdx[1].isEmpty() )
		{
			scene->centralObj->sampleMoreTri(collisionFaceIdx[1]);
		}

		// recompute all the pairwise IBS
		computePairwiseIBS();
	}
}

void Interaction::detectCollision()
{
	QVector<SurfaceMeshModel*> involvedMeshes;
	involvedMeshes.push_back(ibs->mesh);
	involvedMeshes.push_back(ibs->obj1->subdividedMesh);	
	involvedMeshes.push_back(ibs->obj2->subdividedMesh);

	// detect the intersection between IBS and it's corresponding 
	DetectCollision d(involvedMeshes);
	d.DoDetect(0, 1);
	d.DoDetect(0, 2);
	
	// get all the detecting faces 
	isCollisionFree = true;
	collisionFaceIdx.clear();
	for (int i=1; i<=2; i++)
	{
		QVector<int> faces = QVector<int>::fromStdVector(d.Get_contact_faces(i));
		collisionFaceIdx.push_back(faces);
		if (!faces.isEmpty())
		{
			isCollisionFree = false;
		}
	}
}

void Interaction::findSurfaceRegion()
{
	region = new FuncRegion(scene->centralObj);
	region->mapFromIbs(ibs);
	ibs->region = region;

	region->color = scene->colorMap[scene->centralObj->regions.size()].toHsv();
	scene->centralObj->regions.push_back(region);

	computeRegionFeature();
}

void Interaction::computeIbsFeature()
{
	ibs->computeGeomFeatures();
	ibs->computeTopoFeatures();
}

void Interaction::computeRegionFeature()
{
	region->computeFeature();
}

void Interaction::storeOrigObj( QVector<Interaction*> interactions )
{
	obj = new Object(scene);
	for (auto inter : interactions)
	{
		obj->origIdx << inter->obj->origIdx;
	}
}

void Interaction::generateCombinedObj( QVector<Interaction*> interactions )
{
	QVector<Object*> objects;
	for (auto inter : interactions)
	{
		objects << inter->obj;
	}

	obj = new Object(scene);
	obj->combineObjects(objects);
}

void Interaction::constructNewInteraction( QVector<Interaction*> interactions )
{
	generateCombinedObj(interactions);
	construct();
}

void Interaction::computeFeature()
{
	type = COMBINED_FULL;

	obj->combineSamples();			// some interactions are merged version, so combining samples is needed to re-compute ibs and ir
	//obj->generateNewMesh();	
	
	construct();
}

QVector<double> Interaction::getFeature( int id )
{

	QVector<double> distribution;

	switch (id)
	{
	case 0:
		distribution = ibs->pfh;
		break;
	case 1:
		distribution = ibs->dirHist;
		break;
	case 2:
		distribution = ibs->distHist;
		break;
	case 3:
		distribution = region->pfh;
		break;
	case 4:
		distribution = region->dirHist;
		break;
	case 5:
		distribution = region->heightHist;
		break;
	}	

	return distribution;
}

inline void addVectors(QVector<double> & vec1, QVector<double> vec2)
{
	if (vec1.isEmpty())
	{
		vec1 = vec2;
	}
	else
	{
		assert(vec1.size() == vec2.size());
		for (int i=0; i<vec1.size(); i++)
		{
			vec1[i] += vec2[i];
		}
	}
}

inline void normalizeVector(QVector<double> & vec)
{
	double sum = 0;
	for (auto v:vec)
	{
		sum += v;
	}

	for (auto &v:vec)
	{
		v /= sum;
	}
}

void Interaction::computeAverageFeatures( QVector<Interaction*> interactions )
{
	ibs = new IBS(scene);
	region = new FuncRegion(scene->centralObj);

	for (auto inter : interactions)
	{
		addVectors(ibs->pfh, inter->ibs->pfh);
		addVectors(ibs->dirHist, inter->ibs->dirHist);
		addVectors(ibs->distHist, inter->ibs->distHist);
		addVectors(region->pfh, inter->region->pfh);
		addVectors(region->dirHist, inter->region->dirHist);
		addVectors(region->heightHist, inter->region->heightHist);
	}

	normalizeVector(ibs->pfh);
	normalizeVector(ibs->dirHist);
	normalizeVector(ibs->distHist);
	normalizeVector(region->pfh);
	normalizeVector(region->dirHist);
	normalizeVector(region->heightHist);
}

