#include "IbsGenerator.h"
#include "QhullFacetList.h"
#include "QhullVertexSet.h"
#include "OrientHelper.h"
#include "UtilityGlobal.h"

#define PI 3.14159265359

IbsGenerator::IbsGenerator()
{
	scene = NULL;
	qhull = NULL;
}

IbsGenerator::~IbsGenerator()
{
	if (qhull)
	{
		delete qhull;
	}	
}

QVector<IBS*> IbsGenerator::computeIBS(Scene * s, QVector<Object*> obj)
{
	scene = s;
	objects = obj;

	computeVoronoi();
	findRidges();
	buildIBS();

	return ibsSet;
}

void IbsGenerator::computeVoronoi()
{
	QVector<Vec3d> points = getInputForVoronoi();

	// use Qhull to create Voronoi diagram
	//qhull = new orgQhull::Qhull("", 3, (int)points.size(), points[0].data(), "v");
	qhull = new orgQhull::Qhull();
	qhull->runQhull("", 3, (int)points.size(), points[0].data(), "v");
	
	voronoiVertices.push_back(Vec3d(0, 0, 0)); // the index of vertices start from 1, 0 is for the infinite one
	for ( orgQhull::QhullFacet facet = qhull->firstFacet(); facet != qhull->endFacet(); facet=facet.next() )
	{
		orgQhull::QhullPoint p = facet.voronoiVertex(qhull->runId());
		voronoiVertices.push_back(Vec3d(p[0], p[1], p[2]));
	}

	//foreach(orgQhull::QhullFacet facet, qhull->facetList()) {
	//	orgQhull::QhullPoint p = facet.voronoiVertex(qhull->runId());
	//	voronoiVertices.push_back(Vec3d(p[0], p[1], p[2]));
	//}
}

QVector<Vec3d> IbsGenerator::getInputForVoronoi()
{
	sampleObjIdx.clear();
	sampleLocalIdx.clear();
	QVector<Vec3d> points;

	// 1. get the sample points on the objects
	for (int i=0; i<objects.size(); i++)
	{
		Object * obj = objects[i];

		int idx = 0;
		for (auto sample : obj->samples)
		{
			points.push_back(sample.pos);
			sampleObjIdx.push_back(i);			// index the corresponding object index
			sampleLocalIdx.push_back(idx++);	// index in each object
		}
	}

	// 2. add sample points on the bounding ball
	Eigen::AlignedBox3d bbox = objects[0]->bbox;
	for (auto obj:objects)
	{
		bbox.extend(obj->bbox);
	}
	double radius = bbox.diagonal().norm() * 2 / 3;
	Vec3d center = bbox.center();

	int M = 20; 
	int N = 40;
	double step_z = PI / M;
	double step_xy = 2*PI / N;

	double angle_z = 0.0;
	double angle_xy = 0.0;

	for(int i=0; i<M; i++)
	{
		angle_z = i * step_z;

		for(int j=0; j<N; j++)
		{
			angle_xy = j * step_xy;

			double x = radius * sin(angle_z) * cos(angle_xy);
			double y = radius * sin(angle_z) * sin(angle_xy);
			double z = radius * cos(angle_z);

			points.push_back(Vec3d(x + center[0], y + center[1], z + center[2]));
			sampleObjIdx.push_back(-1);				// not on the object, but the bounding sphere, so no local index
		}
	}

	//// 2. add samples on the bounding box
	//Vec3d v_min = bbox.center() - bbox.diagonal() * 1.1;
	//Vec3d v_max = bbox.center() + bbox.diagonal() * 1.1;
	//int stepNum = 10;
	//Vec3d step = (v_max - v_min) / stepNum;
	//
	//for (int i=0; i<=stepNum; i++)
	//{
	//	if (i==0 || i==stepNum)
	//	{
	//		for (int j=0; j<=stepNum; j++)
	//		{
	//			for (int k=0; k<=stepNum; k++)
	//			{
	//				points.push_back(Vec3d(v_min[0]+i*step[0], v_min[1]+j*step[1], v_min[2]+k*step[2]));
	//				objIdx.push_back(-1);
	//			}
	//		}
	//	}
	//	else
	//	{
	//		for (int j=0; j<=stepNum; j+=stepNum)
	//		{
	//			for (int k=0; k<=stepNum; k+=stepNum)
	//			{
	//				points.push_back(Vec3d(v_min[0]+i*step[0], v_min[1]+j*step[1], v_min[2]+k*step[2]));
	//				objIdx.push_back(-1);
	//			}
	//		}
	//	}
	//}

	return points;	
}

// refer the function, qh_printvdiagram in io.c in QHull library.

/*
 * visit all pairs of input sites (vertices) for selected Voronoi vertices
*/

void IbsGenerator::findRidges()
{	
	int numcenters;
	boolT isLower;
	facetT *facetlist = qhull->firstFacet().getFacetT();
	setT * vertices = qh_markvoronoi(facetlist, 0, true, &isLower, &numcenters);  // vertices are the input point sites, indexed by pointid

	vertexT *vertex;
	int vertex_i, vertex_n;

	FORALLvertices
		vertex->seen= False;

	int totcount = 0;
	objPair2IbsIdx.clear();
	FOREACHvertex_i_(vertices) 
	{
		if (vertex) 
		{
			if (qh GOODvertex > 0 && qh_pointid(vertex->point)+1 != qh GOODvertex)
				continue;
			totcount += findRidgesAroundVertex(vertex);
		}
	}
}

/// libqhull\io.c -> qh_eachvoronoi
int IbsGenerator::findRidgesAroundVertex(vertexT *atvertex)
{
	// parameters
	qh_RIDGE innerouter = qh_RIDGEall;
	printvridgeT printvridge = qh_printvridge;
	boolT inorder = true;
	boolT visitall = !qh_ALL;

	boolT unbounded;
	int count;
	facetT *neighbor, **neighborp, *neighborA, **neighborAp;
	setT *centers;
	setT *tricenters= qh_settemp(qh TEMPsize);

	vertexT *vertex, **vertexp;
	boolT firstinf;
	unsigned int numfacets= (unsigned int)qh num_facets;
	int totridges= 0;

	qh vertex_visit++;
	atvertex->seen= True;
	if (visitall) {
		FORALLvertices
			vertex->seen= False;
	}
	FOREACHneighbor_(atvertex) {
		if (neighbor->visitid < numfacets)
			neighbor->seen= True;
	}
	FOREACHneighbor_(atvertex) {
		if (neighbor->seen) {
			FOREACHvertex_(neighbor->vertices) {
				if (vertex->visitid != qh vertex_visit && !vertex->seen) {
					vertex->visitid= qh vertex_visit;
					count= 0;
					firstinf= True;
					qh_settruncate(tricenters, 0);
					FOREACHneighborA_(vertex) {
						if (neighborA->seen) {
							if (neighborA->visitid) {
								if (!neighborA->tricoplanar || qh_setunique(&tricenters, neighborA->center))
									count++;
							}else if (firstinf) {
								count++;
								firstinf= False;
							}
						}
					}

					if (count >= qh hull_dim - 1) { // Each ridge has to have at least hull_dim - 1 vertices
						if (firstinf) {
							if (innerouter == qh_RIDGEouter)
								continue;
							unbounded= False;
						}else {
							if (innerouter == qh_RIDGEinner)
								continue;
							unbounded= True;
						}
						totridges++;
						trace4((qh ferr, 4017, "qh_eachvoronoi: Voronoi ridge of %d vertices between sites %d and %d\n",
							count, qh_pointid(atvertex->point), qh_pointid(vertex->point)));
						if (printvridge) {
							if (inorder && qh hull_dim == 3+1) /* 3-d Voronoi diagram */
								centers= qh_detvridge3 (atvertex, vertex); // determine 3-d Voronoi ridge from 'seen' neighbors of atvertex and vertex, listed in adjacency order (not oriented)
							else
								centers= qh_detvridge(vertex);
							// centers : set of facets (i.e., Voronoi vertices)

							// save the new ridge to our own data structure
							{
								QVector< int > sites; // pair of input sites separated by this ridge
								sites.push_back( qh_pointid(atvertex->point) );
								sites.push_back( qh_pointid(vertex->point) );
								
								bool upperdelaunay = false;
								facetT *facet, **facetp;
								QVector<int> ridge;
								FOREACHfacet_(centers)			// the bisecting planer has more than one edges, which corresponds one ridge
								{
									ridge.push_back(facet->visitid);
									if ( facet->upperdelaunay )
									{
										upperdelaunay = true;
										break;
									}
								}

								// only add bounded
								if ( !upperdelaunay )
								{					
									int obj_id_1 = sampleObjIdx[sites[0]];
									int obj_id_2 = sampleObjIdx[sites[1]];	

									if ( obj_id_2 < obj_id_1 ) 
									{
										std::swap(obj_id_1, obj_id_2);
										std::swap(sites[0], sites[1]);
									}

									ridgeSitePair.push_back(sites);
									ridges.push_back(ridge);

									if ( obj_id_1 != -1 && obj_id_2 != -1 && ( obj_id_1 != obj_id_2 ) )
									{
										QPair<int, int> obj_id_pair(obj_id_1, obj_id_2);
										if ( objPair2IbsIdx.find(obj_id_pair) == objPair2IbsIdx.end() ) // check whether there is already ridge found between those two objects
										{
											// create a ridge list for a new ibs separating those two objects
											QVector<int> ridge_id_list;
											ibsRidgeIdxs.push_back(ridge_id_list);

											// map the new object pair to the ridge list
											objPair2IbsIdx[obj_id_pair] = (int)ibsRidgeIdxs.size()-1;
										}

										//update ridge_id_list in IBS
										int ibsIdx = objPair2IbsIdx[obj_id_pair];
										ibsRidgeIdxs[ibsIdx].push_back(ridges.size()-1);										
									}						
								}
								
							}

							qh_settempfree(&centers);
						}
					}
				}
			}
		}
	}
	FOREACHneighbor_(atvertex)
		neighbor->seen= False;
	qh_settempfree(&tricenters);
	return totridges;
} 

void IbsGenerator::buildIBS()
{
	ibsSet.clear();
	for (int i=0; i<ibsRidgeIdxs.size(); i++)
	{
		IBS * ibs = new IBS(scene);

		// 1. set the corresponding pair of objects
		QPair<int, int> objPair = objPair2IbsIdx.key(i);
		ibs->obj1 = objects[objPair.first];
		ibs->obj2 = objects[objPair.second];

		// 2. build the mesh & get the corresponding sample pairs
		ibs->mesh = buildIbsMesh(i, ibs->samplePairs);	

		ibsSet.push_back(ibs);
	}
}

SurfaceMeshModel* IbsGenerator::buildIbsMesh( int i,  QVector<QPair<int, int>>& samplePairs )
{
	//////////////////////////////////////////////////////////////////////////
	// 1. re-index the vertices and update the ridges
	QVector<int> vIdx;							// vertices
	QVector<int> vIdxNew(voronoiVertices.size(), -1);
	QVector< QVector<int> > ridgesNew;			// faces
	QVector<int> remainRidgeIdx;

	for (int j=0; j < ibsRidgeIdxs[i].size(); j++)
	{
		int ridgeIdx = ibsRidgeIdxs[i][j];
		QVector<int> ridge = ridges[ridgeIdx];
		QVector<int> ridgeWithNewId; 

		if(ridge.size() < 3)						// a face must contain more than 2 edges(ridges)
		{
			continue;
		}

		//find all the vertices on this ridge
		for(int k = 0; k < ridge.size(); k++)
		{
			int v_id = ridge[k];					// vertice_id(visited_id)
		
			if(vIdxNew[v_id] == -1)					// duplicated edges are not counted
			{
				vIdxNew[v_id] = vIdx.size();		// counting from 0, locally
				vIdx.push_back(v_id);				// preserve the original index of ridges[k]
			}

			// renew the vertex index
			ridgeWithNewId.push_back(vIdxNew[v_id]);
		}
		
		ridgesNew.push_back(ridgeWithNewId);
		remainRidgeIdx.push_back(ridgeIdx);			// remainRidgeIdx = ibsRidgeIbs[i]
	}
	
	//////////////////////////////////////////////////////////////////////////
	// 2. re-orient all faces coherently
	OrientHelper help;
	ridgesNew = help.reorient(ridgesNew, vIdx.size());	// ridgesNew stores new index

	//////////////////////////////////////////////////////////////////////////
	// 3. make the normal point from obj1 to obj2
	// get the first ridge triangle to see if its normal points to obj2
	QVector<Vec3d> fv;
	fv << voronoiVertices[vIdx[ridgesNew[0][0]]] << voronoiVertices[vIdx[ridgesNew[0][1]]] << voronoiVertices[vIdx[ridgesNew[0][2]]];
	Vec3d n = (fv[2] - fv[1]).cross(fv[0] - fv[1]).normalized();
	Vec3d center = (fv[0] + fv[1] + fv[2]) / 3.0;

	int ridge_id = ibsRidgeIdxs[i][0];
	QVector<int> pair = ridgeSitePair[ridge_id];
	int sIdx = pair[1];
	Vec3d s2 = objects[sampleObjIdx[sIdx]]->samples[sampleLocalIdx[sIdx]].pos;
	Vec3d d = (s2 - center).normalized();

	// n is the normal of Voronoi faces, and d is the direction pointing from center to s2
	// 因为上面已经翻转了所有的面，并保证了方向一致性，所以只需要计算一次flip就行
	// 然后所有的面的方向就一致了。
	bool flip = n.dot(d) < 0;	// if the normal is not pointing to objs, flip the mesh

	//////////////////////////////////////////////////////////////////////////
	// 4. build the mesh
	SurfaceMeshModel * mesh = new SurfaceMeshModel();

	// add vertices
	for (int j=0; j<vIdx.size(); j++) 
	{
		mesh->add_vertex(voronoiVertices[vIdx[j]]);
	}

	// add faces
	Face f;
	for(int j=0; j<ridgesNew.size(); j++)
	{
		QVector<int> pair = ridgeSitePair[remainRidgeIdx[j]];
		// add the triangle fan for each polygon
		for (int k=1; k<ridgesNew[j].size()-1; k++)
		{
			if (flip)
			{
				f = mesh->add_triangle(Vertex(ridgesNew[j][k+1]), Vertex(ridgesNew[j][k]), Vertex(ridgesNew[j][0]));
			}
			else
			{
				f = mesh->add_triangle(Vertex(ridgesNew[j][0]), Vertex(ridgesNew[j][k]), Vertex(ridgesNew[j][k+1]));
			}	

			if (f.is_valid())
			{	
				samplePairs.push_back(QPair<int, int>(sampleLocalIdx[pair[0]], sampleLocalIdx[pair[1]]));								
			}
		}
	}

	if (samplePairs.size() != mesh->faces_size())
	{
		debugBox( "ERROR: the number of sample pairs does not equal to that of triangles!!!  " + QString::number(samplePairs.size()) + " vs. " + QString::number(mesh->faces_size()));
	}

	return mesh;
}
