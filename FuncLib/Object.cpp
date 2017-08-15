#include "Object.h"
#include "Scene.h"
#include "FuncRegion.h"
#include "QuickMeshDraw.h"
#include "UtilityGlobal.h"
#include "writeOBJ.h"
#include "SymGroup.h"

Object::Object()
{
	scene = NULL;
	origMesh = NULL;
	subdividedMesh = NULL;
	surfaceArea = 0;

	isCentral = false;

	isInteracting = false;
	clusterIdx = -1;

	color = QColor(168, 168, 168);
	colors.resize(0);
}

Object::Object( Scene *s)
{
	scene = s;
	origMesh = NULL;
	subdividedMesh = NULL;
	surfaceArea = 0;

	isCentral = false;

	isInteracting = false;
	clusterIdx = -1;

	color = QColor(168, 168, 168);
	colors.resize(0);
}

Object::~Object()
{
	if (origMesh)
	{
		delete origMesh;

		if (subdividedMesh && subdividedMesh != origMesh)
		{
			delete subdividedMesh;
		}
	}
	else
	{
		for (auto comp : components)
		{
			if (comp)
			{
				delete comp;
			}
		}
	}

	for (auto sg:symGroups)
	{
		if (sg)
		{
			delete sg;
		}
	}
	symGroups.clear();
}

void Object::load( QString filename )
{
	clearAllComponents();

	if (origMesh)
	{
		delete origMesh;
	}

	origMesh = new SurfaceMeshModel();
	origMesh->read(filename.toStdString());	
	
	origMesh->updateBoundingBox();
	bbox = origMesh->bbox();	

	//subdividedMesh = origMesh;	// add this step to subdivide(), the next step
	subdivide();
	computeArea();

	determineLabelColor();
}

void Object::computeHeightRange()
{
	heightRange.clear();
	Surface_mesh::Vertex_property<Vector3> points = subdividedMesh->vertex_property<Vector3>("v:point");
	for (auto v:subdividedMesh->vertices())
	{
		double d = points[v].dot(scene->upright);

		if (heightRange.isEmpty())
		{
			heightRange << d << d;
		}
		else
		{
			heightRange[0] = (heightRange[0] < d)? heightRange[0] : d;
			heightRange[1] = (heightRange[1] > d)? heightRange[1] : d;
		}
	}
}

void Object::setMesh(SurfaceMeshModel * m)
{
	clearAllComponents();
	if (origMesh)
	{
		delete origMesh;
	}

	origMesh = m;

	origMesh->updateBoundingBox();
	bbox = origMesh->bbox();

	subdividedMesh = origMesh;
	subdivide();
	computeArea();
}

void Object::combineObjects( QVector<Object*> objects )
{
	origIdx.clear();
	int sampleNum = 0;
	QVector<SurfaceMeshModel*> meshes;
	for (auto obj:objects)
	{
		origIdx << obj->origIdx;
		meshes.push_back(obj->subdividedMesh);
		sampleNum += obj->samples.size();	
	}
	combineMeshes(meshes);
	sampling(sampleNum);

	// copy the symmetry information
	symGroups.clear();
	repeatedPatterns.clear();
	patternIdx.clear();
	patternToSym.clear();
	for (auto obj:objects)
	{
		if (!obj->symGroups.isEmpty())
		{
			for (int i=0; i<obj->patternToSym.size(); i++)
			{
				QVector<int> newPatterMap;
				newPatterMap << obj->patternToSym[i][0] + symGroups.size();
				newPatterMap << obj->patternToSym[i][1];

				patternToSym << newPatterMap;
			}			

			for (int i=0; i<obj->symGroups.size(); i++)
			{
				SymGroup* sg = new SymGroup(*(obj->symGroups[i]));
				sg->obj = this;
				sg->addOffsetToPatternIdx(repeatedPatterns.size());
				symGroups << sg;
			}

			for (int i=0; i<obj->patternIdx.size(); i++)
			{
				if (obj->patternIdx[i] == -1)
				{
					patternIdx << -1;
				}
				else
				{
					patternIdx << obj->patternIdx[i] + repeatedPatterns.size();
				}
			}

			repeatedPatterns << obj->repeatedPatterns;	
		}
	}
}

void Object::drawPatterns()
{
	if (isCentral && !repeatedPatterns.isEmpty())
	{
		if(!subdividedMesh->property("hasNormals").toBool())
		{
			subdividedMesh->update_face_normals();
			subdividedMesh->update_vertex_normals();
			subdividedMesh->setProperty("hasNormals",true);
		}
		glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		
		if(pattern_colors.size()==0)
		{
			pattern_colors.resize(repeatedPatterns.size());
			for(int i = 0; i < repeatedPatterns.size(); i++)
			{
				Eigen::Vector3d tmpc;
				tmpc(0) = (rand()%1000)/1000.0f;
				tmpc(1) = (rand()%1000)/1000.0f;
				tmpc(2) = (rand()%1000)/1000.0f;
				pattern_colors[i] = tmpc;
			}
		}
		Surface_mesh::Vertex_property<Vector3> points = subdividedMesh->vertex_property<Vector3>("v:point");
		Surface_mesh::Face_property<Vector3> fnormals = subdividedMesh->face_property<Vector3>("f:normal");

		Surface_mesh::Face_iterator fit, fend = subdividedMesh->faces_end();
		Surface_mesh::Vertex_around_face_circulator fvit, fvend;

		glBegin(GL_TRIANGLES);
		int count = 0;
		for (fit=subdividedMesh->faces_begin(); fit!=fend; ++fit){
			Eigen::Vector3d c;
			if(patternIdx[count]==-1)
				c = Eigen::Vector3d(0.5f,0.5f,0.5f);
			else
				c = pattern_colors[patternIdx[count]];
			glColor3d(c(0),c(1),c(2));

			glNormal3( fnormals[fit] );
			fvit = fvend = subdividedMesh->vertices(fit);
			do{ glVector3( points[fvit] ); } while (++fvit != fvend);
			count++;
		}
		glEnd();

		glEnable(GL_LIGHTING);
		glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	}
}

void Object::drawSegmentation()
{
	if (subdividedMesh&&isCentral)
	{
		if(!subdividedMesh->property("hasNormals").toBool())
		{
			subdividedMesh->update_face_normals();
			subdividedMesh->update_vertex_normals();
			subdividedMesh->setProperty("hasNormals",true);
		}
		glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		QVector<int> seg_labels;
//		QVector<Eigen::Vector3d> colors;
		int maxL = -1;
		QFile file(scene->filename+"_labels.seg");
		if (file.open(QIODevice::ReadOnly | QIODevice::Text))
		{
			QTextStream in(&file);
			while(!in.atEnd())
			{
				int l;
				in >> l;
				seg_labels.push_back(l);
				if(l > maxL)
					maxL = l;
			}
		}
		else
			return;
		if(colors.size()==0)
		{
			colors.resize(maxL+1);
			for(int i = 0; i < maxL+1; i++)
			{
				Eigen::Vector3d tmpc;
				tmpc(0) = (rand()%1000)/1000.0f;
				tmpc(1) = (rand()%1000)/1000.0f;
				tmpc(2) = (rand()%1000)/1000.0f;
				colors[i] = tmpc;
			}
		}
		Surface_mesh::Vertex_property<Vector3> points = subdividedMesh->vertex_property<Vector3>("v:point");
		Surface_mesh::Face_property<Vector3> fnormals = subdividedMesh->face_property<Vector3>("f:normal");

		Surface_mesh::Face_iterator fit, fend = subdividedMesh->faces_end();
		Surface_mesh::Vertex_around_face_circulator fvit, fvend;

		glBegin(GL_TRIANGLES);
		int count = 0;
		for (fit=subdividedMesh->faces_begin(); fit!=fend; ++fit){

			Eigen::Vector3d c = colors[seg_labels[count]];
			glColor3d(c(0),c(1),c(2));

			glNormal3( fnormals[fit] );
			fvit = fvend = subdividedMesh->vertices(fit);
			do{ glVector3( points[fvit] ); } while (++fvit != fvend);
			count++;
		}
		glEnd();

		glEnable(GL_LIGHTING);
		glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	}
}

void Object::combineMeshes( QVector<SurfaceMeshModel *> meshes )
{
	clearAllComponents();
	if (origMesh)
	{
		delete origMesh;
	}
	
	origMesh = new SurfaceMeshModel();
	int n = 0;
	Surface_mesh::Vertex_property<Vector3> points;
	Surface_mesh::Vertex_around_face_circulator fvit, fvend;
	for (auto m: meshes)
	{
		points = m->vertex_property<Vector3>("v:point");

		// add vertices
		for (auto v:m->vertices())
		{
			origMesh->add_vertex(points[v]);
		}

		// add faces
		for (auto f:m->faces())
		{
			std::vector<Vertex> face;
			fvit = fvend = m->vertices(f);
			do
			{ 
				face.push_back( Vertex(Vertex(fvit).idx() + n )); 
			} while (++fvit != fvend);

			origMesh->add_face(face);
		}

		n += m->vertices_size();
	}
	
	origMesh->updateBoundingBox();
	bbox = origMesh->bbox();	

	subdividedMesh = origMesh;			// subdividedMesh equals to origMesh
	computeHeightRange();
	computeArea();
}

void Object::draw(OBJECT_DRAW_MODE mode)
{
	// colors are controlled here
	QColor c(168, 168, 168);
	if (isCentral)
	{
		c = QColor(237, 125, 49);
	}
	else if (isInteracting)
	{
		c = QColor(153, 217, 234);
	}

	switch (mode)
	{
	case DRAW_MESH:
		drawMesh(c);
		break;

	case DRAW_SAMPLE:
		drawSamples(c);
		break;

	case DRAW_BBOX_OBJ:
		drawMeshBBox(c);
		break;

	case DRAW_BBOX_COMP:
		drawCompBBox(c);
		break;

	case DRAW_REGION_MESH:
		drawMeshRegion();
		break;

	case DRAW_REGION_SAMPLE:
		if (!regions.isEmpty())
		{
			drawBaseMesh();
			drawSampleRegion();
		}		
		break;

	case DRAW_CLUSTER_MESH:
		drawMeshCluster();
		break;

	case DRAW_CLUSTER_SAMPLE:		
		break;

	case DRAW_WIRE_ORIG:
		if (origMesh)
		{
			QuickMeshDraw::drawMeshWireFrame(origMesh);
		}
		break;

	case DRAW_WIRE_SUBD:
		if (subdividedMesh)
		{
//			QuickMeshDraw::drawMeshWireFrame(subdividedMesh);
			drawSegmentation();
		}
		break;

	case DRAW_MESH_COLOR:
		drawMesh(color);
		break;

	case DRAW_NONE:
		drawPatterns();
		break;
	}

}

void Object::drawBaseMesh()
{
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	drawMesh(QColor(235,235,235,255));		
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);	
}

void Object::drawMesh(QColor color)
{
	if (origMesh)
	{		
		glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		QuickMeshDraw::drawMeshSolid(origMesh, color);		
		glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	}
	else
	{
		for (auto comp : components)
		{
			glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
			QuickMeshDraw::drawMeshSolid(comp, color);
			glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
		}
	}
}

void Object::drawMeshRegion()
{
	if (subdividedMesh && faceDeterminingScore.is_valid())
	{
		if(!subdividedMesh->property("hasNormals").toBool())
		{
			subdividedMesh->update_face_normals();
			subdividedMesh->update_vertex_normals();
			subdividedMesh->setProperty("hasNormals",true);
		}

		glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		glEnable (GL_BLEND);
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		Surface_mesh::Vertex_property<Vector3> points = subdividedMesh->vertex_property<Vector3>("v:point");
		Surface_mesh::Face_property<Vector3> fnormals = subdividedMesh->face_property<Vector3>("f:normal");

		Surface_mesh::Face_iterator fit, fend = subdividedMesh->faces_end();
		Surface_mesh::Vertex_around_face_circulator fvit, fvend;

		glBegin(GL_TRIANGLES);
		for (fit=subdividedMesh->faces_begin(); fit!=fend; ++fit){

			QVector<Scalar> ibsWeight = faceDeterminingScore[fit];
			double maxWeight = 0;
			int idx = -1;
			for (int i=0; i<ibsWeight.size(); i++)
			{
				if (ibsWeight[i] > maxWeight)
				{
					maxWeight = ibsWeight[i];
					idx = i;
				}
			}

			QColor c(255,255,255,255);
			if (idx != -1)
			{
				c = regions[idx]->color;
			}
			c = c.toHsv();
			c.setHsvF(c.hsvHueF(), pow(maxWeight, 1.0/4.0)*0.8 + 0.2, 1.0);
			c = c.toRgb();

			glColor4d(c.redF(), c.greenF(), c.blueF(), c.alphaF());

			glNormal3( fnormals[fit] );
			fvit = fvend = subdividedMesh->vertices(fit);
			do{ glVector3( points[fvit] ); } while (++fvit != fvend);
		}
		glEnd();

		glEnable(GL_LIGHTING);
		glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	}
}

void Object::drawMeshCluster()
{
	QColor c(168, 168, 168); 
	if (isCentral)
	{
		c = QColor(237, 125, 49);
	}
	else if (isInteracting)
	{
		if (clusterIdx != -1)
		{
			c = scene->clusteredInteractions[clusterIdx]->region->color;
		}	
	}
	drawMesh(c.toRgb());
}

void Object::drawMeshBBox(QColor c)
{
	starlab::BoxSoup objBoxRender;
	objBoxRender.addBox(bbox, c);
	objBoxRender.draw();
}

void Object::drawCompBBox(QColor c)
{	
	if (origMesh)
	{
		drawMeshBBox(c);
	}
	else
	{
		starlab::BoxSoup compBoxRender;
		for (auto comp:components)
		{		
			compBoxRender.addBox(comp->bbox(), c);
		}
		compBoxRender.draw();
	}	
}

void Object::drawSamples(QColor c)
{
	starlab::PointSoup sampleRender;
	for (auto sample : samples)
	{
		//sampleRender.addPointNormal(sample.pos, sample.n, color);
		sampleRender.addPoint(sample.pos, c);
	}
	sampleRender.draw();
}

void Object::drawSampleRegion()
{
	if (regions.isEmpty())
	{
		return;
	}

	starlab::PointSoup  sampleRegionRender;
	for (int i=0; i<samples.size(); i++)
	{			
		int idx = dominateRegion[i];
		if (idx != -1)
		{
			QColor c = regions[idx]->color;
			double score = regions[idx]->determiningScore[i] / regions[idx]->maxScore;
			c.setHsvF(c.hsvHueF(), pow(score, 1.0/4.0)*0.8 + 0.2, 1.0);
			c = c.toRgb();	
			sampleRegionRender.addPoint(samples[i].pos, c);
		}				
	}

	sampleRegionRender.draw();
}

void Object::drawIR(FuncRegion* regions, QColor c)
{
	if (regions) {
		starlab::PointSoup sampleRegionRender;
		for (int i = 0; i < samples.size(); ++i)
		{
			double score = regions->determiningScore[i] / regions->maxScore;
			if (score > 0) {
				sampleRegionRender.addPoint(samples[i].pos, c);
			}
		}
		sampleRegionRender.draw();
	}
}

void Object::drawIR(int idx, QColor c)
{
	if (regions.isEmpty())
	{
		return;
	}

	starlab::PointSoup  sampleRegionRender;
	for (int i=0; i<samples.size(); i++)
	{		
		double score = regions[idx]->determiningScore[i] / regions[idx]->maxScore;
		if (score > 0)
		{
			//c.setHsvF(c.hsvHueF(), pow(score, 1.0/4.0)*0.8 + 0.2, 1.0);
			//c = c.toRgb();	
			sampleRegionRender.addPoint(samples[i].pos, c);		
		}				
	}

	sampleRegionRender.draw();
}

void Object::drawSampleCluster()
{
}

void Object::clearAllComponents()
{
	for(auto comp : components)
	{
		if (comp)
		{
			delete comp;
			comp = NULL;
		}
	}

	components.clear();
}

void Object::computeArea()
{
	surfaceArea = 0;

	if (subdividedMesh)				// when no components, subdivideMesh is necessary for computing area
	{		
		SurfaceMeshHelper h(subdividedMesh);
		ScalarFaceProperty farea = h.computeFaceAreas();

		for (auto f:subdividedMesh->faces())
		{
			surfaceArea += farea[f];
		}
	}
	else
	{
		for (auto comp:components)
		{
			SurfaceMeshHelper h(comp);
			ScalarFaceProperty farea = h.computeFaceAreas();

			for (auto f:comp->faces())
			{
				surfaceArea += farea[f];
			}
		}
	}
}

void Object::sampling(int num)
{
	Sampler s(subdividedMesh);						// subdividedMesh is used for sampling
	samples = s.getSamples(num, 0);

	sampleDensity = double(num) / surfaceArea;
}

void Object::sampleMoreTri( QVector<int> tIdx )
{
	ScalarFaceProperty farea = subdividedMesh->get_face_property<Scalar>(FAREA);

	ScalarFaceProperty fweight = subdividedMesh->get_face_property<Scalar>("f:weight");
	if (fweight.is_valid())
	{
		subdividedMesh->remove_face_property<Scalar>(fweight);
	}
	fweight = subdividedMesh->face_property<Scalar>("f:weight", 0);

	double totalArea = 0;
	for (auto i:tIdx)
	{
		fweight[Face(i)] = farea[Face(i)];
		totalArea += farea[Face(i)];
	}

	int sampleNum = sampleDensity * totalArea * 5;

	Sampler s(subdividedMesh, RANDOM_BARYCENTRIC_WEIGHTED);
	samples << s.getSamples(sampleNum, 0);
}

void Object::setCentral( bool s )
{
	isCentral = s;
}

void Object::reverseCentralState()
{
	isCentral = !isCentral;
}

void Object::subdivide()
{
	if (!subdividedMesh)
		subdividedMesh = origMesh;

	subFaceToOrig.clear();
	int i = 0;
	for (auto f : subdividedMesh->faces()) {
		subFaceToOrig.push_back(i);
		++i;
	}

	subdividedMesh->updateBoundingBox();

	//if (!subdividedMesh)
	//{
	//	return;
	//}
	//else
	//{
	//	subdividedMesh = origMesh;
	//	subFaceToOrig.clear();
	//	int i=0;
	//	for (auto f:subdividedMesh->faces())
	//	{
	//		subFaceToOrig << i;
	//		i++;
	//	}

	//	subdividedMesh->updateBoundingBox();
	//	return;
	//}

	//double edgeThreshold = subdividedMesh->bbox().diagonal().norm() / 9.0;

	//QVector<Vec3d> subdividedV;
	//QVector< QVector<int> > subdividedF;
	//subFaceToOrig.clear();

	//SurfaceMeshHelper h(subdividedMesh);
	//Vector3VertexProperty points = h.getVector3VertexProperty(VPOINT);
	//for (auto v:subdividedMesh->vertices())
	//{
	//	subdividedV.push_back(points[v]);
	//}

	//for (auto f:subdividedMesh->faces())
	//{
	//	Surface_mesh::Vertex_around_face_circulator fvit, fvend;
	//	fvit = fvend = subdividedMesh->vertices(f);

	//	QVector< QVector<int> > fStack;   //stack for split edge 

	//	QVector< int > firstFace;
	//	do{ firstFace.push_back( Vertex(fvit).idx() ); } while (++fvit != fvend);
	//	fStack.push_back(firstFace);

	//	while(fStack.size() != 0)
	//	{
	//		QVector<int> currFace = fStack.back();
	//		fStack.pop_back();

	//		int id1 = -1;     // id of first vertex of the longest edge
	//		int id2 = -1;
	//		double length = 0;

	//		//find longest edge
	//		for(int i = 0; i < currFace.size(); i++)
	//		{
	//			int j = (i+1) % currFace.size();
	//			double d = (subdividedV[currFace[i]] - subdividedV[currFace[j]]).norm();   //current edge length 

	//			if(length < d)
	//			{
	//				id1 = i;
	//				id2 = j;
	//				length = d;
	//			}
	//		}

	//		// if the face is small enough, add to the subdivided mesh
	//		if(length <= edgeThreshold)
	//		{				
	//			subdividedF.push_back(currFace);
	//			subFaceToOrig.push_back(f.idx());
	//			continue;
	//		}
	//		else // otherwise, split the edge
	//		{
	//			//add a point in the middle of the longest edge
	//			Vec3d mid_point = (subdividedV[currFace[id1]] + subdividedV[currFace[id2]]) / 2;
	//			subdividedV << mid_point;				

	//			//push all the result triangles
	//			int idx = subdividedV.size() - 1;
	//			int id3 = (id2+1) % currFace.size(); 

	//			QVector<int> t1;
	//			t1 << idx << currFace[id2] << currFace[id3];

	//			QVector<int> t2;
	//			t2 << idx << currFace[id3] << currFace[id1];

	//			fStack << t1 << t2;
	//		}
	//	}//while stack not empty
	//}	

	//// build the subdivided mesh
	//if (subdividedMesh != origMesh)
	//{
	//	delete subdividedMesh;
	//}
	//subdividedMesh = new SurfaceMeshModel();

	//// add vertices
	//for (auto v:subdividedV) 
	//{
	//	subdividedMesh->add_vertex(v);
	//}

	//// add faces
	//for(auto f:subdividedF)
	//{
	//	subdividedMesh->add_triangle(Vertex(f[0]), Vertex(f[1]), Vertex(f[2]));		
	//}

	//subdividedMesh->updateBoundingBox();
}

// this function is valid only for central object
void Object::prepareRegionDraw()
{
	completeRegionDeterminingScore();
	normalizeScore();
	propogateScore();
}

void Object::completeRegionDeterminingScore()
{
	for (auto r:regions)
	{
		r->completeDetermingScore();
	}
}

void Object::normalizeScore()			// normalization for all regions
{
	double maxScore = 0;
	for (auto r:regions)
	{
		for (auto sIdx:r->sampleIdxs)
		{
			maxScore = (maxScore > r->determiningScore[sIdx])? maxScore : r->determiningScore[sIdx];
		}
	}	

	for (auto r:regions)
	{
		for (auto sIdx:r->sampleIdxs)
		{
			r->determiningScore[sIdx] /= maxScore;
		}
	}
}

void Object::propogateScore()
{
	// map the score to the mesh
	if (regions.isEmpty())
	{
		return;
	}
		
	// get the dominate region idx for each sample, regions can be overlapping
	findDominateRegionForSample();

	// compute face score
	computeFaceScore();
}

void Object::findDominateRegionForSample()
{
	dominateRegion.clear();
	for (int i=0; i<samples.size(); i++)
	{
		int idx = -1;
		double maxWeight = 0;

		for (int j=0; j<regions.size(); j++)
		{
			if (regions[j]->determiningScore[i] > maxWeight)
			{
				maxWeight = regions[j]->determiningScore[i];
				idx = j;
			}
		}

		dominateRegion.push_back(idx);		// the dominant region for each sample
	}
}

void Object::computeFaceScore()
{
	int rNum = regions.size();
	faceDeterminingScore = subdividedMesh->face_property<QVector<Scalar>>("f:score", QVector<Scalar>(rNum, 0));
	for (int i=0; i<samples.size(); i++)
	{
		Face f(samples[i].findex);
		for (int j=0; j<rNum; j++)
		{
			faceDeterminingScore[f][j] += regions[j]->determiningScore[i];	// summation of contribution of all regions samples
		}
	}

	// normalize
	double highest = 0;
	for (auto f:subdividedMesh->faces())
	{
		for(int j=0; j<rNum; j++)
		{
			if (faceDeterminingScore[f][j] > highest)
			{
				highest = faceDeterminingScore[f][j];
			}
		}
	}

	for (auto f:subdividedMesh->faces())
	{
		for(int j=0; j<rNum; j++)
		{
			faceDeterminingScore[f][j] /= highest;			
		}
	}
}

void Object::generateNewMesh()
{
	if (origIdx.isEmpty())
	{
		return;
	}

	QVector<Object*> objects;
	for (int i=0; i<origIdx.size(); i++)
	{
		objects.push_back(scene->objects[origIdx[i]]);
	}

	combineObjects(objects);
}

void Object::combineSamples()
{
	if (origIdx.isEmpty())
	{
		return;
	}

	// bbox
	bbox = scene->objects[origIdx[0]]->bbox;
	for (int i=1; i<origIdx.size(); i++)
	{
		bbox.extend(scene->objects[origIdx[i]]->bbox);
	}

	// samples
	samples.clear();
	for (int i=0; i<origIdx.size(); i++)
	{
		//find the corresponding IBS
		int interIdx = scene->obj2Inter[origIdx[i]];
		IBS* origIbs = scene->interactions[interIdx]->ibs;
		Object* origObj = scene->interactions[interIdx]->obj;
		assert(origObj->origIdx.size() == 1);

		QSet<int> sampleIdx;
		for (int j=0; j<origIbs->samplePairs.size(); j++)
		{
			QPair<int,int> pair = origIbs->samplePairs[j];
			if (origIbs->pointToCentralObject)
			{
				sampleIdx.insert(pair.first);
			}
			else
			{
				sampleIdx.insert(pair.second);
			}
		}

		for (QSet<int>::iterator iter = sampleIdx.begin(); iter!=sampleIdx.end(); iter++)
		{
			SamplePoint s = origObj->samples[*iter];
			s.findex = -1;
			samples.push_back(s);
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// compute geometric features

void Object::computeGeoFeature(int barNum)
{
	if (samples.isEmpty())
	{
		sampling(1024);
	}

	computeBS();
//	geoFeature.partSize = bbox.diagonal().norm() / scene->bbox.diagonal().norm();
	geoFeature.partSize = ball.radious / scene->allObj->ball.radious;
//	computeHD(barNum);
	computeDD(barNum);
	computeLPS();

	geoFeature.save(scene->filename);
}

void Object::computeBS()
{
	Surface_mesh::Vertex_property<Vector3> points = origMesh->vertex_property<Vector3>("v:point");
	Surface_mesh::Vertex_iterator it_begin = origMesh->vertices_begin();
	Surface_mesh::Vertex_iterator it_end = origMesh->vertices_end();

	double maxX = -999.0f , maxY = -999.0f, maxZ = -999.0f , minX = 999.0f, minY = 999.0f, minZ = 999.0f ;  
    Surface_mesh::Vertex_iterator maxXi, maxYi, maxZi, minXi, minYi, minZi;

    //Find the max and min along the x-axie, y-axie, z-axie
	int count = 0;
    while(it_begin!=it_end)  
    {  
		if(points[it_begin].x() > maxX)
		{
			maxX = points[it_begin].x() ;
			maxXi = it_begin;
		}
		if(points[it_begin].y() > maxY)
		{
			maxY = points[it_begin].y() ;
			maxYi = it_begin;
		}
		if(points[it_begin].z() > maxZ)
		{
			maxZ = points[it_begin].z() ;
			maxZi = it_begin;
		}
		if(points[it_begin].x() < minX)
		{
			minX = points[it_begin].x() ;
			minXi = it_begin;
		}
		if(points[it_begin].y() < minY)
		{
			minY = points[it_begin].y() ;
			minYi = it_begin;
		}
		if(points[it_begin].z() < minZ)
		{
			minZ = points[it_begin].z() ;
			minZi = it_begin;
		}
        
		++it_begin;
		count++;
    }
  
    float x = 0;  
    Vector3 sub1 , sub2 ;  
    sub1[0] = points[maxXi].x() ; sub1[1] = points[maxXi].y() ; sub1[2] = points[maxXi].z() ;  
    sub2[0] = points[minXi].x() ; sub2[1] = points[minXi].y() ; sub2[2] = points[minXi].z() ;  
    sub1 = sub1 - sub2;  
	x = sub1.dot(sub1); 
  
    float y = 0 ;  
    sub1[0] = points[maxYi].x() ; sub1[1] = points[maxYi].y() ; sub1[2] = points[maxYi].z() ;  
    sub2[0] = points[minYi].x() ; sub2[1] = points[minYi].y() ; sub2[2] = points[minYi].z() ;  
    sub1 = sub1 - sub2;  
	y = sub1.dot(sub1);  
  
    float z = 0 ;  
    sub1[0] = points[maxZi].x() ; sub1[1] = points[maxZi].y() ; sub1[2] = points[maxZi].z() ;  
    sub2[0] = points[minZi].x() ; sub2[1] = points[minZi].y() ; sub2[2] = points[minZi].z() ;  
    sub1 = sub1 - sub2;  
	z = sub1.dot(sub1);  
  
    float dia = x ;  
    Surface_mesh::Vertex_iterator max = maxXi , min = minXi ;  
    if( z > x && z > y)  
    {  
        max = maxZi ;  
        min = minZi ;  
        dia = z ;  
    }else if(y > x && y > z)  
    {  
        max = maxYi ;  
        min = minYi ;  
        dia = y ;  
    }  
  
    //Compute the center point  
    ball.center[0] = 0.5 * (points[max].x() + points[min].x()) ;  
    ball.center[1] = 0.5 * (points[max].y() + points[min].y()) ;  
    ball.center[2] = 0.5 * (points[max].z() + points[min].z()) ;  
  
    //Compute the radious  
    ball.radious = 0.5 * sqrt(dia);  
  
    //Fix it
	it_begin = origMesh->vertices_begin();
    while(it_begin!=it_end)   
    {  
        Vector3 d;  
        d = points[it_begin] - ball.center;  
		float dist2 = d.dot(d);  
  
        if(dist2 > ball.radious * ball.radious)  
        {  
            float dist = sqrt(dist2);  
            float newRadious = (dist + ball.radious) * 0.5 ;  
            float k = (newRadious - ball.radious) / dist ;  
            ball.radious = newRadious ;  
            Vector3 temp  = d*k;
			ball.center = ball.center + temp; 
        }
		++it_begin;
    }
}

void Object::computeHD(int barNum)
{
	int uprightIdx = -1;	
	for (int i=0; i<3; i++)
	{
		if (abs(scene->upright[i] - 1) < 1.0e-6)
		{
			uprightIdx = i;
			break;
		}
	}
	double bottom = scene->bbox.min()[uprightIdx];
	double top = scene->bbox.max()[uprightIdx];

	geoFeature.heightDistribution = QVector<double>(barNum, 0);
	for (auto v:samples)
	{
		double d = v.pos.dot(scene->upright);
		d = (d - bottom)/(top - bottom);

		if ( abs(d) < 1.0e-8)
		{
			d = 0;
		}

		int idx = floor(d*barNum);
		if(idx == barNum)
			idx--;
		geoFeature.heightDistribution[idx] += 1;
	}
	
	for (int i=0; i<barNum; i++)
	{
		geoFeature.heightDistribution[i] /= samples.size();
	}
}

void Object::computeDD(int barNum)
{
	QVector<double> tmpD;
	double mind = 1e6;
	double maxd = 0;
	for (auto p:samples)
	{
		for (auto q:samples)
		{
			double d = (p.pos - q.pos).norm();
			tmpD.push_back(d);
			if(mind > d)
				mind = d;
			if(maxd < d)
				maxd = d;
		}
	}

	geoFeature.distDistribution = QVector<double>(barNum, 0);
	for (auto dist:tmpD)
	{
		double d = (dist - mind)/(maxd - mind);
		int idx = floor(d*barNum);
		if(idx == barNum)
			idx--;
		geoFeature.distDistribution[idx] += 1;
	}

	for (int i=0; i<barNum; i++)
	{
		geoFeature.distDistribution[i] /= samples.size()*samples.size();
	}
}

void Object::computeLPS()
{
	Eigen::Matrix3Xd points_data = Eigen::Matrix3Xd::Zero(3, samples.size());
	for (int i=0; i<samples.size(); i++)
	{
		SamplePoint p = samples[i];
		points_data(0,i) = p.pos.x();
		points_data(1,i) = p.pos.y();
		points_data(2,i) = p.pos.z();
	}

	Eigen::Vector3d points_mean;
	for(unsigned int i=0; i<3; ++i)
		points_mean(i) = (points_data.row(i).array()).sum()/points_data.cols();
	points_data.colwise() -= points_mean;
	/////////////////////////////////////////////////////
	Eigen::MatrixXd data = points_data*points_data.transpose();
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> s(data);
	Eigen::VectorXd val = s.eigenvalues();
	Eigen::MatrixXd vec = s.eigenvectors();

	int n = data.cols();
	for (int i = 0; i < n - 1; ++i) {
		int k;
		val.segment(i, n - i).maxCoeff(&k);
		if (k > 0) {
			std::swap(val[i], val[k + i]);
			vec.col(i).swap(vec.col(k + i));
		}
	}

	//for (int i = 0; i < n; ++i)
	//	geoFeature.LPS.push_back(vec.col(i));

	double sum = val[0] + val[1] + val[2];
	geoFeature.LPS[0] = (val[0] - val[1]) / sum;
	geoFeature.LPS[1] = 2 * (val[1] - val[2]) / sum;
	geoFeature.LPS[2] = 3 * val[2] / sum;
}

// assign color to objects based on the name
void Object::determineLabelColor()
{
	//QVector< QColor > colors;
	//colors << QColor(241, 95, 116);
	//colors << QColor(112, 173, 71);
	//colors << QColor(84, 129, 230);
	//colors << QColor(247, 216, 66);
	//colors << QColor(145, 60, 205);
	//colors << QColor(131, 144, 152);

	//QString filename = QFileInfo(origMesh->name).baseName();


}

void savePointsToMesh(QVector<Eigen::Vector3d> points, double radius, QString filename)
{
	SurfaceMeshModel* mesh = new SurfaceMeshModel();

	for (int i=0; i<points.size(); i++)
	{
		Eigen::Vector3d center = points[i];

		// add the points at the conner
		mesh->add_vertex(points[i] + radius * Eigen::Vector3d(-1, -1, -1));
		mesh->add_vertex(points[i] + radius * Eigen::Vector3d(1, -1, -1));
		mesh->add_vertex(points[i] + radius * Eigen::Vector3d(1, 1, -1));
		mesh->add_vertex(points[i] + radius * Eigen::Vector3d(-1, 1, -1));
		mesh->add_vertex(points[i] + radius * Eigen::Vector3d(1, -1, 1));
		mesh->add_vertex(points[i] + radius * Eigen::Vector3d(1, 1, 1));
		mesh->add_vertex(points[i] + radius * Eigen::Vector3d(-1, 1, 1));
		mesh->add_vertex(points[i] + radius * Eigen::Vector3d(-1, -1, 1));

		// add the faces
		int j = i*8;
		mesh->add_triangle(Vertex(j), Vertex(j+2), Vertex(j+1));
		mesh->add_triangle(Vertex(j), Vertex(j+3), Vertex(j+2));
		mesh->add_triangle(Vertex(j+4), Vertex(j+5), Vertex(j+7));
		mesh->add_triangle(Vertex(j+5), Vertex(j+6), Vertex(j+7));
		mesh->add_triangle(Vertex(j), Vertex(j+1), Vertex(j+7));
		mesh->add_triangle(Vertex(j+1), Vertex(j+4), Vertex(j+7));
		mesh->add_triangle(Vertex(j+1), Vertex(j+2), Vertex(j+5));
		mesh->add_triangle(Vertex(j+1), Vertex(j+5), Vertex(j+4));
		mesh->add_triangle(Vertex(j+2), Vertex(j+3), Vertex(j+5));
		mesh->add_triangle(Vertex(j+3), Vertex(j+6), Vertex(j+5));
		mesh->add_triangle(Vertex(j), Vertex(j+7), Vertex(j+3));
		mesh->add_triangle(Vertex(j+7), Vertex(j+6), Vertex(j+3));
	}

	mesh->update_face_normals();
	mesh->update_vertex_normals();

	writeOBJ::wirte(mesh, filename);	
	delete mesh;
}

void Object::saveIR(int idx)
{
	QVector<Eigen::Vector3d> points;
	for (int i=0; i<samples.size(); i++)
	{		
		double score = regions[idx]->determiningScore[i] / regions[idx]->maxScore;
		if (score > 0)
		{
			points << samples[i].pos;		
		}				
	}

	double radius = origMesh->bbox().diagonal().norm() * 0.001;

	//savePointsToMesh(points, radius, "IR_" + QString::number(idx) + "_0.obj");
	//savePointsToMesh(points, radius*2, "IR_" + QString::number(idx) + "_1.obj");

	savePointsToMesh(points, radius, scene->filename + "_IR_" + QString::number(idx) + "_0.obj");
	savePointsToMesh(points, radius*2, scene->filename + "_IR_" + QString::number(idx) + "_1.obj");
}

void Object::saveIR(FuncRegion* regions)
{
	QVector<Eigen::Vector3d> points;
	for (int i=0; i<samples.size(); i++)
	{		
		double score = regions->determiningScore[i] / regions->maxScore;
		if (score > 0)
		{
			points << samples[i].pos;		
		}				
	}

	double radius = origMesh->bbox().diagonal().norm() * 0.001;

	//savePointsToMesh(points, radius, "IR_" + QString::number(idx) + "_0.obj");
	//savePointsToMesh(points, radius*2, "IR_" + QString::number(idx) + "_1.obj");

	savePointsToMesh(points, radius, scene->filename + "_IR_" + "_0.obj");
	savePointsToMesh(points, radius*2, scene->filename + "_IR_" + "_1.obj");
}

void Object::loadSymGroups()
{
	// assign pattern index to each triangle (one triangle can only be in one pattern, right?)*/
	QString name = scene->filename + "_centric_labels.lb";
	patternIdx.resize(origMesh->faces_size());
	for(int i = 0; i < patternIdx.size(); i++)
		patternIdx[i] = -1;

	if(!QFile::exists(name))
		return;

	loadSymGroupsfromfile(name);	
	// update the patternIdx and repeatedPatterns for subdivided mesh
	if (!repeatedPatterns.isEmpty() && !subFaceToOrig.isEmpty())
	{
		QVector<int> temp = patternIdx;
		patternIdx.clear();
		for (int i=0; i<repeatedPatterns.size(); i++)
		{
			repeatedPatterns[i].clear();
		}
		for (auto i:subFaceToOrig)
		{
			if (temp[i] != -1)
			{
				repeatedPatterns[temp[i]] << patternIdx.size();
			}
			patternIdx << temp[i];
		}
	}

	// store map form pattern to symmetry group
	patternToSym.resize(repeatedPatterns.size());
	for (int i=0; i<patternToSym.size(); i++)
	{
		patternToSym[i] = QVector<int>(2, -1);
	}

	for (int k=0; k<symGroups.size(); k++)
	{
		SymGroup * sg = symGroups[k];
		for (int i=0; i<sg->pIdx1D.size(); i++)
		{
			patternToSym[ sg->pIdx1D[i] ][0] =  k;
			patternToSym[ sg->pIdx1D[i] ][1] =  0;
		}

		for (int i=0; i<sg->pIdx2D.size(); i++)
		{
			for (int j=0; j<sg->pIdx2D[i].size(); j++)
			{
				patternToSym[ sg->pIdx2D[i][j] ][0] = k;
				patternToSym[ sg->pIdx2D[i][j] ][1] = i;
			}
		}
	}

	// test symLevel
	Eigen::MatrixXd levels = Eigen::MatrixXd::Zero(repeatedPatterns.size(), repeatedPatterns.size());
	for (int i=0; i<repeatedPatterns.size(); i++)
	{
		for (int j=i; j<repeatedPatterns.size(); j++)
		{
			levels(i, j) = symLevel(i, j);
			levels(j, i) = levels(i, j);
		}
	}
	matrixToFile(levels, scene->filename + "_symLevel.csv");
}

void Object::loadSymGroupsfromfile(QString file)
{
	symGroups.clear();
	repeatedPatterns.clear();	QFile File(file);
	QVector<int> total_f;
	if(File.open(QFile::ReadOnly | QFile::Text))
	{
		QTextStream in(&File);
		QString line = in.readLine().simplified();
		int count = 0;
		SymGroup* sg = new SymGroup(this);
		while(!line.isNull())
		{			
			QStringList elements = line.split(" ");
			if(elements[0] == "end")
			{
				total_f.clear();
				symGroups.push_back(sg);
				line = in.readLine().simplified();
				continue;
			}
			if(elements[0] == "g")
			{
				sg = new SymGroup(this);
				count = 0;
				int typei = elements[1].toInt();
				switch (typei)
				{
				case 0:
					{
						sg->type = GRID;
						sg->gridT1[0] = elements[2].toDouble();
						sg->gridT1[1] = elements[3].toDouble();
						sg->gridT1[2] = elements[4].toDouble();
						sg->gridT2[0] = elements[5].toDouble();
						sg->gridT2[1] = elements[6].toDouble();
						sg->gridT2[2] = elements[7].toDouble();
						break;
					}
				case 1:
					{
						sg->type = TRANS;
						sg->trans[0] = elements[2].toDouble();
						sg->trans[1] = elements[3].toDouble();
						sg->trans[2] = elements[4].toDouble();
						break;
					}
				case 2:
					{
						sg->type = ROT;
						sg->rotCenter[0] = elements[2].toDouble();
						sg->rotCenter[1] = elements[3].toDouble();
						sg->rotCenter[2] = elements[4].toDouble();
						sg->rotRadius = elements[5].toDouble();
						break;
					}
				default:
					{
						QMessageBox::warning(NULL,"No this symmetry type!","symmetry type error!");
						break;
					}
				}
			}
			else if(elements[0] == "p")
			{
				if(elements[1].toInt() == 1)
				{
					for(int i = 0; i < elements[2].toInt(); i++)
						sg->pIdx1D.push_back(i+repeatedPatterns.size());
					//repeatedPatterns.resize(elements[2].toInt());
				}
				else
				{
					sg->pIdx2D.resize(elements[1].toInt());
					int pre_count = 0;
					for(int i =0; i < elements[1].toInt(); i++)
					{
						for(int j = 0; j < elements[i+2].toInt(); j++)
							sg->pIdx2D[i].push_back(j+pre_count+repeatedPatterns.size());
						pre_count += sg->pIdx2D[i].size();
					}
					//repeatedPatterns.resize(pre_count);
				}
			}
			else
			{
				QVector<int> tmp_pattern;
				for(int i = 0; i < elements.size(); i++)
				{
					int idx = elements[i].toInt();
					if(qFind(total_f.begin(),total_f.end(),idx)==total_f.end()||total_f.size()==0)
					{
						total_f.push_back(idx);
						tmp_pattern.push_back(idx);
						patternIdx[idx] = repeatedPatterns.size();
					}
					else
					{
						//const QString w = "Repeated faces!" + QString::number(idx);
						//QMessageBox::warning(NULL,"Warning",w);
					}				
				}
				repeatedPatterns.push_back(tmp_pattern);
			}
			line = in.readLine().simplified();
		}
	}
}

void Object::identifyPatterEachRegionIn()
{
	if (repeatedPatterns.isEmpty())
	{
		return;
	}

	QVector< QVector<double> > regionPatternWeight(regions.size());
	for (int i=0; i<regions.size(); i++)
	{
		regionPatternWeight[i] = QVector<double>(repeatedPatterns.size()+1, 0);
	}

	Surface_mesh::Face_iterator fit, fend = subdividedMesh->faces_end();
	for (fit=subdividedMesh->faces_begin(); fit!=fend; ++fit)
	{
		QVector<Scalar> ibsWeight = faceDeterminingScore[fit];
		int pIdx = patternIdx[Face(fit).idx()];

		for (int i=0; i<regions.size(); i++)
		{
			regionPatternWeight[i][pIdx+1] += ibsWeight[i];
		}			
	}

	QVector<int> rp;
	for (int i=0; i<regions.size(); i++)
	{
		int idx = -1;
		double maxWeight = 0;
		for (int j=0; j<regionPatternWeight[i].size(); j++)
		{
			if (regionPatternWeight[i][j] > maxWeight)
			{
				maxWeight = regionPatternWeight[i][j];
				idx = j;
			}
		}

		assert(idx >= 0);

		regions[i]->patternIdx = idx - 1;
		rp << idx - 1;
	}

	// test symLevel
	Eigen::MatrixXd levels = Eigen::MatrixXd::Zero(rp.size(), rp.size());
	for (int i=0; i<rp.size(); i++)
	{
		for (int j=i; j<rp.size(); j++)
		{
			levels(i, j) = symLevel(rp[i], rp[j]);
			levels(j, i) = levels(i, j);
		}
	}
	matrixToFile(levels, scene->filename + "_symLevel_obj.csv");
}

int Object::symLevel(int pIdx1, int pIdx2)
{
	int level = 0;

	if (pIdx1 < 0 || pIdx1 >= repeatedPatterns.size() || pIdx2 < 0 || pIdx2 >= repeatedPatterns.size())
	{
		return level;
	}
		
	if (pIdx1 == pIdx2) // same pattern 
	{
		level = 3;
	}
	else if(patternToSym[pIdx1][0] == patternToSym[pIdx2][0]) // same symmetry group: check whether they are in the inner or outer group
	{
		if (patternToSym[pIdx1][1] == patternToSym[pIdx2][1])
		{
			level = 2;
		}
		else
		{
			level = 1;
		}
	}

	return level;
}