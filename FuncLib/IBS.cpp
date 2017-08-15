#include "IBS.h"
#include "Scene.h"
#include "QuickMeshDraw.h"
#include "RenderObjectExt.h"
#include "SurfaceMeshHelper.h"
#include <QTime>
#include "UtilityGlobal.h"

#define PI 3.1415926

IBS::IBS()
{
	scene = NULL;
	mesh = NULL;
	sampleRatio = 1.0;
	maxWeight = 0;
	totalWeight = 0;
	bettiNumbers << 0 << 0 << 0;
}


IBS::IBS(Scene *s)
{
	scene = s;
	mesh = NULL;
	sampleRatio = 1.0;
	maxWeight = 0;
	totalWeight = 0;
	bettiNumbers << 0 << 0 << 0;
}

IBS::~IBS()
{
	if (mesh)
	{
		delete mesh;
	}
}

void IBS::draw(	bool drawIbsSample, bool drawIbsWeight, QColor color)
{
	if (mesh)
	{
		ScalarFaceProperty fweight = mesh->get_face_property<Scalar>("f:weight");

		if (!drawIbsWeight || !fweight.is_valid())
		{
			color.setAlphaF(0.2);

			glDepthMask(GL_FALSE);
			QuickMeshDraw::drawMeshSolid(mesh, color, false);	
			glDepthMask(GL_TRUE);
		}
		else
		{
			if(!mesh->property("hasNormals").toBool())
			{
				mesh->update_face_normals();
				mesh->update_vertex_normals();
				mesh->setProperty("hasNormals",true);
			}

			glDepthMask(GL_FALSE);
			glEnable (GL_BLEND);
			glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glDisable(GL_LIGHTING);		
			
			Vector3VertexProperty points = mesh->vertex_property<Vector3>("v:point");
			Vector3FaceProperty fnormals = mesh->face_property<Vector3>("f:normal");

			Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
			Surface_mesh::Vertex_around_face_circulator fvit, fvend;

			glBegin(GL_TRIANGLES);
			for (fit=mesh->faces_begin(); fit!=fend; ++fit){
				QColor c = starlab::qtJetColor(fweight[fit], 0 , maxWeight);
				c.setAlphaF(0.2);
				glColorQt(c);

				glNormal3( fnormals[fit] );
				fvit = fvend = mesh->vertices(fit);
				do{ glVector3( points[fvit] ); } while (++fvit != fvend);
			}
			glEnd();

			glEnable(GL_LIGHTING);
			glDepthMask(GL_TRUE);
		}		

		if (drawIbsSample)
		{
			sampleRender.draw();
		}		
	}
}

void IBS::computeSampleWeightForTri()			// according to Xi's IBS paper
{
	ScalarFaceProperty fweight = mesh->get_face_property<Scalar>("f:weight");
	if (fweight.is_valid()) 
	{
		return;
	}	

	SurfaceMeshHelper h(mesh);
	Vector3VertexProperty points = h.getVector3VertexProperty(VPOINT);
	ScalarFaceProperty farea = h.computeFaceAreas();

	mesh->update_face_normals();
	Vector3FaceProperty fnormal = mesh->get_face_property<Vector3>(FNORMAL);
	ScalarFaceProperty fdist = mesh->add_face_property<Scalar>("f:dist");
	
	maxWeight = 0;
	totalWeight = 0;
	fweight = mesh->face_property<Scalar>("f:weight", 0);
	for (auto f : mesh->faces())
	{
		int idx = samplePairs[f.idx()].second;
		Vec3d s = obj2->samples[idx].pos;

		Vec3d center(0,0,0);
		Surface_mesh::Vertex_around_face_circulator fvit, fvend;
		fvit = fvend = mesh->vertices(f);
		do{ center += points[fvit]; } while (++fvit != fvend);
		center /= 3.0;

		Vec3d d = s - center;
		Vec3d n = fnormal[f];

		double dist = d.norm();
		double angle = acos(d.dot(n) / dist);	

		if ( angle >= PI / 2 )
		{
			fdist[f] = -dist;
			dist = 0;
			angle = PI - angle;
		}
		else
		{
			fdist[f] = dist;
		}
		
		//Eigen::AlignedBox3d bbox = scene->objects[objIdx1]->bbox;
		//bbox.extend(scene->objects[objIdx2]->bbox);
		//double distThreshold = bbox.diagonal().norm() * 0.5;

		double distThreshold = scene->bbox.diagonal().norm() * 0.5;
		double disWeight = (dist <= distThreshold)? pow(1 - (dist / distThreshold), 20) : 0;
			
		double angleThreshold = PI / 3;
		double angleWeight = (angle <= angleThreshold)? (1 - angle / angleThreshold) : 0 ;

		double areaWeight = farea[f];

		fweight[f] = disWeight * angleWeight * areaWeight;

		totalWeight += fweight[f];
		if (fweight[f] > maxWeight)
		{
			maxWeight = fweight[f];
		}
	}

	//for (auto f:mesh->faces())
	//{
	//	fweight[f] /= totalWeight;
	//}
	//maxWeight /= totalWeight;
}

void IBS::sampling( int num )
{
	if (samples.size() == num)
	{
		return;
	}

	computeSampleWeightForTri();

	Sampler s(mesh, RANDOM_BARYCENTRIC_WEIGHTED);
	samples = s.getSamples(num);
	
	sampleRender.clear();
	for (auto sample : samples)
	{
		//sampleRender.addPointNormal(sample.pos, sample.n, QColor(255, 255, 0));
		sampleRender.addPoint(sample.pos, QColor(255, 255, 0));
	}
}

void IBS::computeGeomFeatures()
{
	computePFH();
	computeDirHist();
	computeDistHist();
}

void IBS::computeTopoFeatures()
{
	if (scene->distPara.useTopo)
	{
		computeBettiNumbers();
	}
	
}

void IBS::computePFH()
{
	// sample points
	int sampleNum = 1000 * sampleRatio;
	sampling(sampleNum);

	// compute PFH, make sure normal points to interacting object
	QVector<QVector<double>> hist;
	for (int i = 0; i < samples.size(); ++i)
		hist << computePfhForSample(i, pointToCentralObject);

	QVector<double> mean(hist[0].size(), 0);
	for(int i=0; i<mean.size(); i++)
	{
		for (auto h:hist)
		{
			mean[i] += h[i];
		}

		mean[i] /= hist.size();
	}

	QVector<double> deviation(hist[0].size(), 0);		// standard deviation of the histogram
	for(int i=0; i<deviation.size(); i++)
	{
		for (auto h:hist)
		{
			deviation[i] += pow(h[i] - mean[i], 2);
		}

		deviation[i] /= hist.size();
		deviation[i] = sqrt(deviation[i]);
	}

	pfh.clear();
	pfh << mean << deviation;		
}

void IBS::computeDirHist()
{
	int sampleNum = 5000 * sampleRatio;
	sampling(sampleNum);

	int bNum = 10;
	QVector<double> h(bNum, 0);
	for (auto s:samples)
	{
		Vec3d n = s.n;
		if (pointToCentralObject)
		{
			n = -n;
		}

		double angle = acos(scene->upright.dot(n));

		int bIdx = angle / PI * 10;
		bIdx = (bIdx > 9)? 9:bIdx;

		h[bIdx]++;
	}

	for (auto& v:h)
	{
		v /= samples.size();
	}	

	dirHist = h;
}

void IBS::computeDistHist()
{
	int sampleNum = 5000 * sampleRatio;
	sampling(sampleNum);

	int bNum = 10;
	distHist.resize(bNum);
	for (int i=0; i<bNum; i++)
	{
		distHist[i] = 0;
	}

	Eigen::AlignedBox3d bbox = obj1->bbox;
	bbox.extend(obj2->bbox);
	double th = bbox.diagonal().norm() / 8.0;

	for (auto s:samples)
	{
		Vec3d site = obj1->samples[samplePairs[s.findex].first].pos;
		double d = (s.pos - site).norm();

		int bIdx = (int) (d / th * 10);
		bIdx = bIdx > 9? 9 : bIdx;

		distHist[bIdx]++;
	}

	for (auto& v:distHist)
	{
		v /= samples.size();
	}

}

QVector<double> IBS::computePfhForSample( int sIdx, bool reverseNormal )
{
	int bNum = 125;
	QVector<double> samplePFH(bNum, 0);
	Vec3d n1 = samples[sIdx].n;
	if (reverseNormal)
	{
		n1 = -n1;
	}
	for (int i=0; i<samples.size(); i++)
	{
		if (i != sIdx)
		{
			Vec3d n2 =  samples[i].n;
			if (reverseNormal)
			{
				n2 = -n2;
			}
			Vec3d p1p2 = samples[i].pos - samples[sIdx].pos;
			Vec3d v = p1p2.cross(n1).normalized();
			Vec3d w = n1.cross(v).normalized();
			Vec3d projected_n2 = Vec3d(w.dot(n2), n1.dot(n2), 0);

			double phi = acos(n1.dot(p1p2) / p1p2.norm());
			double alpha = acos(n2.dot(v));

			Vec3d local_e2(0,1,0);
			double theta = acos(local_e2.dot(projected_n2)/ projected_n2.norm());
			double cross = local_e2[0]*projected_n2[1] - local_e2[1]*projected_n2[0];
			if (cross < 0)
			{
				theta = 2*PI - theta;
			}

			double bWidth1 = PI / 5.0;
			double bWidth2 = 2 * bWidth1; 

			int phiIdx = phi / bWidth1;
			int alphaIdx = alpha / bWidth1;			// alpha \in [0, \pi]
			int thetaIdx = theta / bWidth2;			// theta \in [0, 2*\pi]

			int bIdx = alphaIdx * 25 + thetaIdx * 5 + phiIdx;

			samplePFH[bIdx]++;
		}
	}

	for (auto& h:samplePFH)			// normalized to [0, 1]
	{
		h /= (samples.size() - 1);
	}

	return samplePFH;
}

void IBS::computeBettiNumbers()
{
	QTime t;
	t.start();

	int positive, negative, positiveTri, negativeTri;
	positive = negative = positiveTri = negativeTri = 0;

	QVector< QVector<Edge> > complexEdgeIdx;
	QVector<int> complexOpenEdgeNumber;
	Surface_mesh::Edge_property<Integer> edgeComplexIdx = mesh->add_edge_property("e:complexIdx", -1);

	Surface_mesh::Halfedge_around_face_circulator feit, feend;
	for(auto f:mesh->faces())
	{
		QVector<Edge> consistingEdgeIdx;

		// find the complex idx for each edge
		feit = feend = mesh->halfedges(f);
		do 
		{
			Edge edge = mesh->edge(Halfedge(feit));
			consistingEdgeIdx << edge;			

			if (edgeComplexIdx[edge] != -1) 
			{
				complexOpenEdgeNumber[edgeComplexIdx[edge]]--;
			}
			else// if the edge haven't been checked, find the complex it belongs to
			{	
				int complexIdx = -1;

				// if any of the end point has been used in previous checked edges, copy the corresponding complexID
				QVector<int> vComplexIdx(2, -1);
				for (int vIdx=0; vIdx<2; vIdx++)
				{
					Vertex v = mesh->vertex(edge, vIdx);

					Surface_mesh::Halfedge_around_vertex_circulator vhit, vhend;
					vhit = vhend = mesh->halfedges(v);
					do 
					{
						if (edgeComplexIdx[mesh->edge(vhit)] != -1)
						{
							vComplexIdx[vIdx] = edgeComplexIdx[mesh->edge(vhit)];
							complexIdx = vComplexIdx[vIdx];
							break;
						}
					} while (++vhit != vhend);
				}

				
				if (complexIdx != -1)  // find connected complex
				{
					// add to one valid complex
					edgeComplexIdx[edge] = complexIdx;
					complexEdgeIdx[complexIdx].push_back(edge);
					complexOpenEdgeNumber[complexIdx]++;					

					if (vComplexIdx[0] == vComplexIdx[1]) // if two end points correspond to same complex 
					{
						positive++;
					}
					else// if two end points correspond to two different complexes
					{
						// if both two complexes are valid, merge them
						if (vComplexIdx[0]!=-1 && vComplexIdx[1]!=-1)
						{
							for (auto e:complexEdgeIdx[vComplexIdx[1]])
							{
								edgeComplexIdx[e] = vComplexIdx[0];
							}
							complexEdgeIdx[vComplexIdx[0]] << complexEdgeIdx[vComplexIdx[1]];
							complexOpenEdgeNumber[vComplexIdx[0]] += complexOpenEdgeNumber[vComplexIdx[1]];

							complexEdgeIdx[vComplexIdx[1]].clear();
							complexOpenEdgeNumber[vComplexIdx[1]] = -1;
						}

						negative++;
					}
				}
				else // create a new complex
				{
					QVector<Edge> newComplex;
					newComplex.push_back(edge);
					edgeComplexIdx[edge] = complexEdgeIdx.size();
					complexEdgeIdx.push_back(newComplex);
					complexOpenEdgeNumber.push_back(1);

					negative++;
				}
			}		

		} while (++feit != feend);	

		// If all three edges belong to the same complex and the complex has no open edge
		QVector<int> triComplexIdx;
		triComplexIdx << edgeComplexIdx[consistingEdgeIdx[0]] << edgeComplexIdx[consistingEdgeIdx[1]] << edgeComplexIdx[consistingEdgeIdx[2]];
		if (triComplexIdx[0] == triComplexIdx[1] && triComplexIdx[0] == triComplexIdx[2] && complexOpenEdgeNumber[triComplexIdx[0]] == 0)
		{
			positiveTri ++;
		}			
		else 
		{
			negativeTri++;
		}	
	}	

	bettiNumbers.clear();
	bettiNumbers.push_back(mesh->vertices_size() - 1 - negative);
	bettiNumbers.push_back(positive-negativeTri);
	bettiNumbers.push_back(positiveTri);

	// there might be small holes, set a threshold for computing b2
	int eThreshold = 10;

	int complexNum = 1;
	for( int i=0; i<complexOpenEdgeNumber.size(); i++ )
	{
		if( complexOpenEdgeNumber[i] >= 0 ) 
		{
			qDebug() << "Complex " << complexNum++ << ": have " << complexOpenEdgeNumber[i] << "open edges";

			if (complexOpenEdgeNumber[i]>0 && complexOpenEdgeNumber[i] <= eThreshold)
			{
				bettiNumbers[2]++;
			}
		}
	}

	qDebug("New code --- Time elapsed: %d ms", t.elapsed());

	//ignoreSmallHoles();
}

void IBS::ignoreSmallHoles()
{
	Surface_mesh::Halfedge_property<bool> hvisisted = mesh->halfedge_property<bool>("h:visisted", false);

	QVector< QVector<Halfedge> > holes;
	for (auto h:mesh->halfedges())
	{
		if( !mesh->is_boundary(h) || hvisisted[h] ) continue;

		QVector<Halfedge> hole;

		Halfedge hinit = h;
		h = mesh->next_halfedge(h);

		while( h != hinit ) {
			hole.push_back( h );
			hvisisted[h] = true;
			h = mesh->next_halfedge(h);
		}
		hole.push_back( h );

		holes.push_back(hole);
	}

	if (!(holes.size()==0 || holes.size() == bettiNumbers[1] + 1))
	{
		debugBox("Betti number computation error! holeNum = " + QString::number(holes.size()) + " vs Betti[1] = " + QString::number(bettiNumbers[1]));
	}
	else
	{
		int eThreshold = 10;
		for ( auto h: holes)
		{
			if ( h.size() < eThreshold)
			{
				bettiNumbers[1]--;
			}
		}
	}
}

// transfer original copy from Xi, extremely slow
void IBS::computeBettiNumbers2() 
{
	QTime t;
	t.start();

	int positive, negative, positiveTri, negativeTri;
	positive = negative = positiveTri = negativeTri = 0;

	QVector< QVector<int> > Edges; // Unique edges (v1, v2, complex id, next edge in complex, number of triangles using this edge, edge idx)
	QVector<int> Complexes; // Complexes created by edges, the index of latest edge added into the complex
	QVector<int> ComplexOpenEdges; // Open edges for each complex, number of open edges in that complex
	int TrigIDs[6] = {0,1,0,2,1,2}; // Vertex ID combinations, the i-th edge has two end points with idx TrigIDs[2*i] and TrigIDs[2*i+1]

	QVector<int> tmpEdge; // Temporal edge vector
	tmpEdge << -1 << -1 << -1 << -1 << 1 << 0;
	int TriangleEdges[3]; // the 3 edge idx for this triangle
	int tmpComplexID[2]; // Complex IDs for 2 vertices building up an edge

	Surface_mesh::Vertex_around_face_circulator fvit, fvend;
	for(auto f:mesh->faces())
	{
		QVector<int> vIdx;
		fvit = fvend = mesh->vertices(f);
		do{ vIdx.push_back(Vertex(fvit).idx()); } while (++fvit != fvend);

		for(int j=0; j<3; j++)
		{
			bool UniqueEdge = true;
			tmpComplexID[0] = tmpComplexID[1] = -1;
			for(int k=0; k < Edges.size(); k++)
			{
				// if any of the end point has been used in previous edges, copy the corresponding complexID
				if(Edges[k][0]==vIdx[TrigIDs[j*2]] || Edges[k][1]==vIdx[TrigIDs[j*2]]) 
				{
					tmpComplexID[0] = Edges[k][2];
				}
				if(Edges[k][0]==vIdx[TrigIDs[j*2+1]] || Edges[k][1]==vIdx[TrigIDs[j*2+1]])
				{
					tmpComplexID[1] = Edges[k][2];
				}

				// if the edge already exists, update the edge information
				if((Edges[k][0]==vIdx[TrigIDs[j*2]] && Edges[k][1]==vIdx[TrigIDs[j*2+1]]) ||
					(Edges[k][1]==vIdx[TrigIDs[j*2]] && Edges[k][0]==vIdx[TrigIDs[j*2+1]]))
				{
					UniqueEdge = false;
					TriangleEdges[j] = k;				

					// if this edge is shared by two triangles, then it will not be an open edge. The number of open edges in that complex --
					Edges[k][4]++;
					if(Edges[k][4]==2) 
					{
						ComplexOpenEdges[Edges[k][2]]--;
					}					
					break;
				}
			}
			if(UniqueEdge) // if this edge is new, add it to the Edges
			{
				// correspond this edge to existing complex if found any
				if(tmpComplexID[0]==-1)
				{
					tmpEdge[2] = tmpComplexID[1];
				}
				if(tmpComplexID[1]==-1)
				{
					tmpEdge[2] = tmpComplexID[0];
				}

				//  Add edge to a complex
				if(tmpEdge[2]>=0)   // if the complex does exist
				{
					// if two end points correspond to same complex 
					if(tmpComplexID[0] == tmpComplexID[1])
					{
						tmpEdge[2] = tmpComplexID[0];
						tmpEdge[3] = Complexes[tmpComplexID[0]];
						Complexes[tmpComplexID[0]] = Edges.size();
						ComplexOpenEdges[tmpEdge[2]]++;
						positive++;
					}
					else  // if two end points correspond to two different complexes
					{						
						if(tmpComplexID[0]==-1)  // Assign to complex 1
						{
							tmpEdge[2] = tmpComplexID[1];
							tmpEdge[3] = Complexes[tmpComplexID[1]];
							Complexes[tmpComplexID[1]] = Edges.size();
							ComplexOpenEdges[tmpComplexID[1]]++;
						}
						else if(tmpComplexID[1]==-1)  // Assign to complex 0
						{
							tmpEdge[2] = tmpComplexID[0];
							tmpEdge[3] = Complexes[tmpComplexID[0]];
							Complexes[tmpComplexID[0]] = Edges.size();
							ComplexOpenEdges[tmpComplexID[0]]++;
						}
						else  // merge complex 0 and complex 1: move edges from complex 1 to complex 0
						{
							// chase edges using Edges[k][3] , update the complex id
							int FirstEdgeOfComplex = Complexes[tmpComplexID[1]];
							int k = FirstEdgeOfComplex;
							while(k>=0)
							{
								Edges[k][2] = tmpComplexID[0];
								FirstEdgeOfComplex=k;
								k = Edges[k][3];
							}

							// connect those two complexes and transfer informations
							Edges[FirstEdgeOfComplex][3] = Complexes[tmpComplexID[0]];
							Complexes[tmpComplexID[0]] = Complexes[tmpComplexID[1]];
							Complexes[tmpComplexID[1]] = -1;
							ComplexOpenEdges[tmpComplexID[0]] += ComplexOpenEdges[tmpComplexID[1]];
							ComplexOpenEdges[tmpComplexID[1]] = -1;

							// add the new edge
							tmpEdge[2] = tmpComplexID[0];
							tmpEdge[3] = Complexes[tmpComplexID[0]];
							Complexes[tmpComplexID[0]] = Edges.size();
							ComplexOpenEdges[tmpComplexID[0]]++;
						}
						negative++;
					}
				}
				else  //if there is no such complex
				{
					// Add a new complex
					tmpEdge[2] = Complexes.size();
					tmpEdge[3] = -1;
					Complexes.push_back(Edges.size());
					ComplexOpenEdges.push_back(1);
					negative++;
				}
				TriangleEdges[j] = Edges.size();

				tmpEdge[0] = vIdx[TrigIDs[j*2]];
				tmpEdge[1] = vIdx[TrigIDs[j*2+1]];
				tmpEdge[5] = Edges.size();
				Edges.push_back(tmpEdge);
			}
		}

		// If all three edges belong to the same complex
		if(Edges[TriangleEdges[0]][2]==Edges[TriangleEdges[1]][2] && Edges[TriangleEdges[1]][2]==Edges[TriangleEdges[2]][2])
		{
			if(ComplexOpenEdges[Edges[TriangleEdges[2]][2]]==0) // if the complex has no open edges
			{
				positiveTri++;
			}
			else // check whether there is a same triangle has been visited 
			{
				bool SameTriangle = false;
				for(int j=0; j<f.idx()-1; j++)
				{
					QVector<int> vIdx2;
					fvit = fvend = mesh->vertices(Face(j));
					do{ vIdx2.push_back(Vertex(fvit).idx()); } while (++fvit != fvend);

					if((vIdx[0]==vIdx2[0] && vIdx[1]==vIdx2[1] && vIdx[2]==vIdx2[2]) ||
						(vIdx[0]==vIdx2[0] && vIdx[1]==vIdx2[2] && vIdx[2]==vIdx2[1]) ||
						(vIdx[0]==vIdx2[1] && vIdx[1]==vIdx2[0] && vIdx[2]==vIdx2[2]) ||
						(vIdx[0]==vIdx2[1] && vIdx[1]==vIdx2[2] && vIdx[2]==vIdx2[0]) ||
						(vIdx[0]==vIdx2[2] && vIdx[1]==vIdx2[1] && vIdx[2]==vIdx2[0]) ||
						(vIdx[0]==vIdx2[2] && vIdx[1]==vIdx2[0] && vIdx[2]==vIdx2[1]))
					{
						SameTriangle=true;
						break;
					}
				}
				if(SameTriangle)
				{
					positiveTri++;
				}
				else
				{
					negativeTri++;
				}
			}
		}
		else
		{
			negativeTri++;
		}
	}	

	bettiNumbers.clear();
	bettiNumbers.push_back(mesh->vertices_size() - 1 - negative);
	bettiNumbers.push_back(positive-negativeTri);
	bettiNumbers.push_back(positiveTri);

	qDebug() << "Triangles: " << mesh->faces_size() << "    Edges: " << mesh->edges_size() << "    Vertices: " << mesh->vertices_size();    
	qDebug() << "Open edges:";
	int obj=1;
	for( int i=0; i<(int)ComplexOpenEdges.size(); i++ )
	{
		if( ComplexOpenEdges[i]>=0 ) 
		{
			qDebug() << "Complex " << obj++ << ": " << ComplexOpenEdges[i];
		}
	}

	qDebug("Old code --- Time elapsed: %d ms", t.elapsed());
}

//////////////////////////////////////////////////////////////////////////
// implementation for community features

IBS::IBS(QVector<IBS*> ibsSet, QVector<bool> reverseNormal)
{
	scene = ibsSet[0]->scene;
	mesh = NULL;
	sampleRatio = 1.0;
	maxWeight = 0;
	totalWeight = 0;

	// make sure that the features of each consisting IBS have been computed
	for (auto ibs : ibsSet)
	{
		if (ibs->dirHist.isEmpty())
		{
			ibs->computeGeomFeatures();
			ibs->computeTopoFeatures();
		}
	}

	if (ibsSet.size() == 1)
	{
		pfh = ibsSet[0]->pfh;
		dirHist = ibsSet[0]->dirHist;
		distHist = ibsSet[0]->distHist;
		bettiNumbers = ibsSet[0]->bettiNumbers;
	}
	else
	{
		// compute combined features
		pfh = combinedPFH(ibsSet, reverseNormal);
		dirHist = combinedDirHist(ibsSet, reverseNormal);
		distHist = combinedDistHist(ibsSet);	
		bettiNumbers = combinedBettiNumber(ibsSet);
	}
}

QVector<double> IBS::combinedPFH(QVector<IBS*> ibsSet, QVector<bool> reverseNormal)
{
	// sampling weight
	double total = 0;
	for (auto ibs : ibsSet)
	{
		total += ibs->totalWeight;
	}

	// sample points
	for (auto ibs : ibsSet)
	{
		int num = ibsSet.size() * 1000 * ibs->totalWeight / total;
		ibs->sampling( num );
	}

	// compute distribution for each point, make sure normal points to the right direction
	QVector< QVector<double> > hist;
	for (int ibsIdx=0; ibsIdx < ibsSet.size(); ibsIdx++)
	{
		for (int i=0; i<ibsSet[ibsIdx]->samples.size(); i++)
		{
			hist << computePfhForSample(ibsIdx, i,  ibsSet, reverseNormal);
		}
	}

	// compute PFH	
	QVector<double> mean(hist[0].size(), 0);
	for(int i=0; i<mean.size(); i++)
	{
		for (auto h:hist)
		{
			mean[i] += h[i];
		}

		mean[i] /= hist.size();
	}

	QVector<double> diviation(hist[0].size(), 0);
	for(int i=0; i<diviation.size(); i++)
	{
		for (auto h:hist)
		{
			diviation[i] += pow(h[i] - mean[i], 2);
		}

		diviation[i] /= hist.size();
		diviation[i] = sqrt(diviation[i]);
	}

	QVector<double> combinedPFH;
	combinedPFH << mean << diviation;	

	return combinedPFH;
}

QVector<double> IBS::computePfhForSample(int ibsIdx, int sIdx, QVector<IBS*> ibsSet, QVector<bool> reverseNormal)
{
	int bNum = 125;
	QVector<double> samplePFH(bNum, 0);

	Vec3d n1 = ibsSet[ibsIdx]->samples[sIdx].n;
	if (reverseNormal[ibsIdx])
	{
		n1 = -n1;
	}

	int total = 0;
	for (int k=0; k<ibsSet.size(); k++)
	{
		QVector<SamplePoint> ibsSamples = ibsSet[k]->samples;
		total += ibsSamples.size();
		for (int i=0; i<ibsSamples.size(); i++)
		{
			if (i != sIdx || k != ibsIdx)
			{
				Vec3d n2 =  ibsSamples[i].n;
				if (reverseNormal[k])
				{
					n2 = -n2;
				}
				Vec3d p1p2 = ibsSamples[i].pos - ibsSet[ibsIdx]->samples[sIdx].pos;

				Vec3d v = p1p2.cross(n1).normalized();
				Vec3d w = n1.cross(v).normalized();
				Vec3d projected_n2 = Vec3d(w.dot(n2), n1.dot(n2), 0);

				double phi = acos(n1.dot(p1p2) / p1p2.norm());
				double alpha = acos(n2.dot(v));

				Vec3d local_e2(0,1,0);
				double theta = acos(local_e2.dot(projected_n2)/ projected_n2.norm());
				double cross = local_e2[0]*projected_n2[1] - local_e2[1]*projected_n2[0];
				if (cross < 0)
				{
					theta = 2*PI - theta;
				}

				double bWidth1 = PI / 5.0;
				double bWidth2 = 2 * bWidth1; 

				int phiIdx = phi / bWidth1;
				int alphaIdx = alpha / bWidth1;
				int thetaIdx = theta / bWidth2;

				int bIdx = alphaIdx * 25 + thetaIdx * 5 + phiIdx;

				samplePFH[bIdx]++;
			}
		}

	}
	
	for (auto& h:samplePFH)
	{
		h /= (total- 1);
	}

	return samplePFH;

}

QVector<double> IBS::combinedDirHist(QVector<IBS*> ibsSet, QVector<bool> reverseNormal)
{
	int n = ibsSet[0]->dirHist.size();
	QVector<double> combinedDirHist(n, 0);

	for (int i=0; i<ibsSet.size(); i++)
	{
		IBS* ibs = ibsSet[i];
		for (int j=0; j<n; j++)
		{
			if (reverseNormal[i])
			{
				combinedDirHist[j] += ibs->dirHist[n-1-j];
			}
			else
			{
				combinedDirHist[j] += ibs->dirHist[j];
			}
		}
	}

	for (int i = 0; i < combinedDirHist.size(); i++)
	{
		combinedDirHist[i] /= ibsSet.size();
	}

	return combinedDirHist;
}

QVector<double> IBS::combinedDistHist(QVector<IBS*> ibsSet)
{
	QVector<double> combinedDistHist(ibsSet[0]->dirHist.size(), 0);

	for (auto ibs : ibsSet)
	{
		for (int i=0; i<combinedDistHist.size(); i++)
		{
			combinedDistHist[i] += ibs->dirHist[i];
		}
	}

	for (int i = 0; i < combinedDistHist.size(); i++)
	{
		combinedDistHist[i] /= ibsSet.size();
	}

	return combinedDistHist;
}

QVector<int> IBS::combinedBettiNumber(QVector<IBS*> ibsSet)
{
	QVector<int> combinedBettiNumber(ibsSet[0]->bettiNumbers.size(), 0);
	return combinedBettiNumber;
}
