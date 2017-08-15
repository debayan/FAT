#include "InterSet.h"
#include <QFile>

InterSet::InterSet()
{
}

InterSet::~InterSet()
{
}

void InterSet::load(QString filename)
{
	QFile file( filename + ".iset");
	if (file.open(QIODevice::ReadOnly | QIODevice::Text)) 
	{
		int interactionNum = 0;
		QTextStream in(&file);	
		while (!in.atEnd()) {
			QString line = in.readLine();
			QStringList list = line.split(" ");
			if(list.isEmpty()) continue;

			if (list[0] == "interactionCount")
			{
				interactionNum = list[1].toInt();
			}
			else if (list[0] == "interaction")
			{
				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "origObjIdx");
				QVector<int> origObjIdx;
				for (int i=1; i<list.size(); i++)
				{
					origObjIdx << list[i].toInt();
				}

				Interaction* inter = new Interaction(NULL, origObjIdx);
				interactions << inter;

				// load the features
				inter->ibs = new IBS(NULL);
				inter->region = new FuncRegion(NULL);
				inter->region->ibs = inter->ibs;

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "ibs_pfh");
				for (int i=1; i<list.size(); i++)
				{
					inter->ibs->pfh << list[i].toDouble();
				}
				assert(inter->ibs->pfh.size() == 250);

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "ibs_dirHist");
				for (int i=1; i<list.size(); i++)
				{
					inter->ibs->dirHist << list[i].toDouble();
				}
				assert(inter->ibs->dirHist.size() == 10);

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "ibs_distHist");
				for (int i=1; i<list.size(); i++)
				{
					inter->ibs->distHist << list[i].toDouble();
				}
				assert(inter->ibs->distHist.size() == 10);

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "region_pfh");
				for (int i=1; i<list.size(); i++)
				{
					inter->region->pfh << list[i].toDouble();
				}
				assert(inter->region->pfh.size() == 250);

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "region_dirHist");
				for (int i=1; i<list.size(); i++)
				{
					inter->region->dirHist << list[i].toDouble();
				}
				assert(inter->region->dirHist.size() == 10);

				line = in.readLine();
				list = line.split(" ");
				assert(list[0] == "region_heightHist");
				for (int i=1; i<list.size(); i++)
				{
					inter->region->heightHist << list[i].toDouble();
				}
				assert(inter->region->heightHist.size() == 10);
			}
		}

		assert(interactionNum == interactions.size());
		file.close();
	}
}

void InterSet::save(QString filename)
{
	QFile file(filename + ".iset");

	if(file.open(QIODevice::WriteOnly | QIODevice::Text))
	{
		QTextStream out(&file);
		out << "interactionCount " << interactions.size() << endl;
		for (int k=0; k<interactions.size(); k++)
		{
			out << "interaction " << k << endl;

			out << "origObjIdx";
			for (int i=0; i<interactions[k]->obj->origIdx.size(); i++)
			{
				out << " " << interactions[k]->obj->origIdx[i];
			}
			out << endl;

			out << "ibs_pfh";
			for (int i=0; i<interactions[k]->ibs->pfh.size(); i++)
			{
				out << " " << interactions[k]->ibs->pfh[i];
			}
			out << endl;

			out << "ibs_dirHist";
			for (int i=0; i<interactions[k]->ibs->dirHist.size(); i++)
			{
				out << " " << interactions[k]->ibs->dirHist[i];
			}
			out << endl;

			out << "ibs_distHist";
			for (int i=0; i<interactions[k]->ibs->distHist.size(); i++)
			{
				out << " " << interactions[k]->ibs->distHist[i];
			}
			out << endl;

			out << "region_pfh";
			for (int i=0; i<interactions[k]->region->pfh.size(); i++)
			{
				out << " " << interactions[k]->region->pfh[i];
			}
			out << endl;

			out << "region_dirHist";
			for (int i=0; i<interactions[k]->region->dirHist.size(); i++)
			{
				out << " " << interactions[k]->region->dirHist[i];
			}
			out << endl;

			out << "region_heightHist";
			for (int i=0; i<interactions[k]->region->heightHist.size(); i++)
			{
				out << " " << interactions[k]->region->heightHist[i];
			}
			out << endl;
		}

		file.close();
	}
}