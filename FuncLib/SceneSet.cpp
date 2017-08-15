#include "SceneSet.h"
#include "UtilityGlobal.h"
#include "DistMeasure.h"
#include <QTime>
#include "writeOBJ.h"

SceneSet::SceneSet()
{

}

SceneSet::~SceneSet()
{
	clearAllScenes();
}

void SceneSet::addScenes( QStringList filenames )
{
	QVector<Scene*> newScenes(filenames.size(), NULL);
	QVector<QString> names = filenames.toVector();

#pragma omp parallel for
	for (int i = 0; i < names.size(); ++i) {
		Scene *s = new Scene();
		s->load(names[i]);
		newScenes[i] = s;
	}

	scenes << newScenes;
}

void SceneSet::addScene( Scene* s )
{
	scenes.push_back(s);
}

bool SceneSet::deleteScene( int i )
{
	if (i<0 || i>=scenes.size())
	{
		return false;
	}

	delete scenes[i];
	scenes.remove(i);

	return true;
}

void SceneSet::clearAllScenes()
{
	for (auto s : scenes) {
		if (s) delete s;
	}
	scenes.clear();
}

QVector<int> SceneSet::getScenesWithNoCentralObject()
{
	QVector<int> idx;

	for (int i=0; i<scenes.size(); i++)
	{
		if (scenes[i] && !scenes[i]->hasCentralObj())
		{
			idx.push_back(i);
		}
	}

	return idx;
}


void SceneSet::outputCentralObject(int idx)
{
	QString filename = scenes[idx]->filename + "_centric.obj";
	scenes[idx]->generateCentralObject();
	scenes[idx]->centralObj->origMesh->update_vertex_normals();
	std::string tmp = filename.toLocal8Bit().constData();
	writeOBJ::wirte(scenes[idx]->centralObj->origMesh, filename);
}

void SceneSet::outputCentralObjects()
{
	for (int i = 0; i < scenes.size(); ++i)
	{
		outputCentralObject(i);
	}
}

inline void showTime(int time)
{
	QString note = "Functional descriptor computation done! -- ";
	if (time < 1000)
	{
		note += QString::number(time) + " ms.";
	}
	else if (time/1000 < 60)
	{
		note += QString::number(time/1000) + " s.";
	}
	else
	{
		note += QString::number(time/60000) + " m " + QString::number((time%60000)/1000) + " s.";
	}
	debugBox(note);
}

void SceneSet::buildInteractionHierarchy()
{
	QTime t;
	t.start();

#pragma omp parallel for
	for (int i=0; i<scenes.size(); i++)
	{		
		scenes[i]->buildInteractionHierarchy();
	}

	showTime(t.elapsed());	
}

void SceneSet::loadInteractionHierarchy()
{
	QTime t;
	t.start();

#pragma omp parallel for
	for (int i=0; i<scenes.size(); i++)
	{		
		scenes[i]->loadInteractionHierarchy();	
	}

	showTime(t.elapsed());
}

void SceneSet::computeFeatures(QString filename, DistParameter distPara, HierParameter hierPara)
{
	QVector<QString> modelNames = getModelNames(filename);

	QTime t;
	t.start();

#pragma omp parallel for
	for (int i=0; i<modelNames.size(); i++)
	{		
		QString resultFolder = QFileInfo(modelNames[i]).dir().absolutePath() + "/../" + distPara.toString() + hierPara.toString();
		QString iconFile = resultFolder + '/' + QFileInfo(modelNames[i]).baseName();
		QFile file(iconFile+".icon");
		
		Scene * s = new Scene();
		s->load(modelNames[i]);
		s->distPara = distPara;
		s->hierPara = hierPara;
		s->buildInteractionHierarchy();	

		delete s;
	}

	showTime(t.elapsed());
}

void SceneSet::computeIBSH()
{
	for (int i=0; i<scenes.size(); i++)
	{		
		scenes[i]->computeIBSH();	
	}
}

void SceneSet::computeIBSH(QString filename, bool useTopo)
{
	QVector<QString> modelNames = getModelNames(filename);

#pragma omp parallel for
	for (int i=0; i<modelNames.size(); i++)
	{		
		QString name = QFileInfo(modelNames[i]).dir().absolutePath() + "/" + QFileInfo(modelNames[i]).baseName() + ".ibsh";
		QFile file(name);
		if (!file.exists())
		{
			Scene * s = new Scene();
			s->load(modelNames[i]);
			s->distPara.useTopo = useTopo;
			s->computeIBSH();	
			delete s;			
		}
	}
}

QVector<QString> SceneSet::getModelNames(QString filename)
{
	QFileInfo fileInfo(filename);
	QString	dirName = fileInfo.dir().absolutePath();

	// testing model names
	QVector<QString> modelNames;
	QFile file(filename);
	if (file.open(QIODevice::ReadOnly | QIODevice::Text)) 
	{
		QTextStream in(&file);	
		while (!in.atEnd())
		{
			QString line = in.readLine();
			line = line.trimmed();
			if (!line.isEmpty())
			{
				QString newModelName = dirName + '/' + line;
				if (QFile::exists(newModelName + ".txt"))
				{
					newModelName += ".txt";
				}
				else if (QFile::exists(newModelName + ".obj"))
				{
					newModelName += ".obj";
				}
				else
				{
					debugBox("Model doesn't exist: " + line);
					return QVector<QString>();
				}

				modelNames << newModelName;
			}
		}
	}

	return modelNames;
}

void SceneSet::outputInterSetFeature()
{
	QTime t;
	t.start();

#pragma omp parallel for
	for (int i=0; i<scenes.size(); i++)
	{		
		scenes[i]->outputInterSetFeature();	
	}

	showTime(t.elapsed());
}

void SceneSet::computeGeoFeature(QString filename)
{
	QVector<QString> modelNames = getModelNames(filename);

#pragma omp parallel for
	for (int i=0; i<modelNames.size(); i++)
	{
		Scene * s = new Scene();
		s->load(modelNames[i]);
		s->computeGeometryFeature(64);		
	}
}
