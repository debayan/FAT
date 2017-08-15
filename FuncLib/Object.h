#pragma once

#include "RenderObjectExt.h"
#include "SurfaceMeshModel.h"
#include "Sampler.h"
#include <QFile>

class SymGroup;
using namespace SurfaceMesh;

enum OBJECT_DRAW_MODE{
	DRAW_MESH, DRAW_SAMPLE, DRAW_BBOX_OBJ, DRAW_BBOX_COMP, 
	DRAW_NONE, DRAW_REGION_MESH, DRAW_REGION_SAMPLE, 
	DRAW_CLUSTER_MESH, DRAW_CLUSTER_SAMPLE, 
	DRAW_WIRE_ORIG, DRAW_WIRE_SUBD,
	DRAW_MESH_COLOR
};

typedef struct ObjGeoFeature
{
	double partSize;	
	QVector<double> heightDistribution;
	QVector<double> distDistribution;
	Eigen::Vector3d LPS;

	void save(QString filename)
	{
		QFile file(filename + ".gf");
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream out(&file);

		out << "partSize " << partSize << endl;

//		out << "heightDistribution";
//		for (auto hd:heightDistribution)
//			out << " " << hd;
//		out << endl;

		out << "distDistribution";
		for (auto dd:distDistribution)
			out << " " << dd;
		out << endl;

		out << "LPS";
		for (int i = 0; i < 3; i++)
				out  << " " << LPS[i];	

		file.close();
	}

	void load(QString filename)
	{
		QFile file( filename + ".gf");
		if (file.open(QIODevice::ReadOnly | QIODevice::Text)) 
		{
			QTextStream in(&file);	
			while (!in.atEnd()) 
			{
				QString line = in.readLine();
				QStringList list = line.split(" ");
				if(list.isEmpty()) continue;

				if (list[0] == "partSize")
				{
					partSize = list[1].toDouble();
				}
//				else if (list[0] == "heightDistribution")
//				{
//					heightDistribution.clear();
//					for (int i=1; i<list.size(); i++)
//					{
//						heightDistribution << list[i].toDouble();
//					}
//				}
				else if (list[0] == "distDistribution")
				{
					distDistribution.clear();
					for (int i=1; i<list.size(); i++)
					{
						distDistribution << list[i].toDouble();
					}
				}
				else if (list[0] == "LPS")
				{
					for (int i=0; i<3; i++)
					{
						LPS[i] = list[i+1].toDouble();
					}
				}
			}

			file.close();
		}
		
	}

}ObjGeoFeature;

struct Sphere  
{  
	Eigen::Vector3d center;  
    float radious;  
}; 

class Scene;
class FuncRegion;
class Interaction;

class Object
{
public:
    Object();
	Object(Scene *s);
	~Object();

public:
	void load( QString filename );
	void setMesh(SurfaceMeshModel * m);
	void combineMeshes(QVector<SurfaceMeshModel *> meshes);
	void combineObjects(QVector<Object*> objects);
	void generateNewMesh(); // generate a new mesh by combining all the clustered original interacting objects
	void combineSamples();  // collect the samples directly from the interacting objects

	void computeArea();
	void sampling(int num);
	void sampleMoreTri(QVector<int> tIdx);
	
	void draw(OBJECT_DRAW_MODE mode);
	void drawMesh(QColor c);
	void drawIR(FuncRegion* regions, QColor c);
	void drawIR(int idx, QColor color);
	void saveIR(FuncRegion* regions);
	void saveIR(int idx);
	void setCentral(bool s);
	void reverseCentralState();
		
	void prepareRegionDraw();
	void computeHeightRange();

	//////////////////////////////////////////////////////////////////////////
	// geometry features
	void computeGeoFeature(int barNum);
	void computeBS();

	void loadSymGroups();
	void identifyPatterEachRegionIn();
	int symLevel(int pIdx1, int pIdx2); // detect the symmetry of two pattern: 0: no symmetry, 1: outer symmetry, 2: inner symmetry, 3: pattern
private:
	void clearAllComponents();

	void subdivide();

	void completeRegionDeterminingScore();
	void normalizeScore();
	void propogateScore();
	void findDominateRegionForSample();
	void computeFaceScore();

	void drawBaseMesh();
	void drawMeshRegion();
	void drawMeshCluster();
	void drawMeshBBox(QColor c);
	void drawCompBBox(QColor c);
	void drawSamples(QColor c);
	void drawSampleRegion();
	void drawSampleCluster();
	void drawSegmentation();
	void drawPatterns();

	void computeLPS();
	void computeHD(int barnum);
	void computeDD(int barnum);
	
	// only used to part-level
	void determineLabelColor();
	void loadSymGroupsfromfile(QString file);

public:
	Scene * scene;

	QVector<int> origIdx;			// the consisting objects in the original scene (each object might be combined by multiple input objects)
	SurfaceMeshModel * origMesh;
	SurfaceMeshModel * subdividedMesh;  
	QVector<SurfaceMeshModel *>  components;

	QVector<int> subFaceToOrig;

	double surfaceArea;
	Eigen::AlignedBox3d bbox;
	QVector<double> heightRange;

	double sampleDensity;
	QVector<SamplePoint> samples;

	//////////////////////////////////////////////////////////////////////////
	// following properties only valid for central objects
	bool isCentral; 
	QVector<FuncRegion*> regions;	// the determining region on the object surface for each IBS
	QVector< int > dominateRegion;	// the dominating region idx for each sample: for illustration
	Surface_mesh::Face_property<QVector<Scalar>> faceDeterminingScore; // assign each face the weights for determining the IBS

	//////////////////////////////////////////////////////////////////////////
	// following properties only valid  for interacting objects
	bool isInteracting;
	int clusterIdx;     // cluster objects corresponding to the same type of interaction

	//////////////////////////////////////////////////////////////////////////
	// geometry features
	ObjGeoFeature geoFeature;
	Sphere ball;

	QColor color;
	/////////////////////////////
	QVector<Eigen::Vector3d> colors;
	QVector<Eigen::Vector3d> pattern_colors;

	//////////////////////////////////////////////////////////////////////////
	// symmetry information
	QVector<SymGroup*> symGroups; // store the symmetry formed by those repeated patterns
	QVector< QVector<int> > repeatedPatterns; // store the face idx for each pattern
	QVector<int> patternIdx;  // store the pattern idx for each triangle, -1 if not belong to any group (another way is to use face property)
	QVector< QVector<int> > patternToSym;
};