#include "SymGroup.h"

SymGroup::SymGroup(Object* obj)
{
	this->obj = obj;

	rotRadius = 0;;
}

SymGroup::SymGroup(Object* obj, SYM_TYPE type)
{
	this->obj = obj;
	this->type = type;

	rotRadius = 0;;
}

void SymGroup::LoadFromFile(QString file)
{
	QFile File(file);
	QVector<int> total_f;
	if(File.open(QFile::ReadOnly | QFile::Text))
	{
		QTextStream in(&File);
		QString line = in.readLine().simplified();
		int count = 0;
		while(!in.atEnd())
		{
			QStringList elements = line.split(" ");
			if(elements[0] == "end")
				break;
			if(elements[0] == "g")
			{
				count = 0;
				int typei = elements[1].toInt();
				switch (typei)
				{
				case 0:
					{
						type = GRID;
						gridT1[0] = elements[2].toDouble();
						gridT1[1] = elements[3].toDouble();
						gridT1[2] = elements[4].toDouble();
						gridT2[0] = elements[5].toDouble();
						gridT2[1] = elements[6].toDouble();
						gridT2[2] = elements[7].toDouble();
						break;
					}
				case 1:
					{
						type = TRANS;
						trans[0] = elements[2].toDouble();
						trans[1] = elements[3].toDouble();
						trans[2] = elements[4].toDouble();
						break;
					}
				case 2:
					{
						type = ROT;
						rotCenter[0] = elements[2].toDouble();
						rotCenter[1] = elements[3].toDouble();
						rotCenter[2] = elements[4].toDouble();
						rotRadius = elements[5].toDouble();
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
						pIdx1D.push_back(i);
					obj->repeatedPatterns.resize(elements[2].toInt());
				}
				else
				{
					pIdx2D.resize(elements[1].toInt());
					int pre_count = 0;
					for(int i =0; i < elements[1].toInt(); i++)
					{
						for(int j = 0; j < elements[i+2].toInt(); j++)
							pIdx2D[i].push_back(j+pre_count);
						pre_count += pIdx2D[i].size();
					}
					obj->repeatedPatterns.resize(pre_count);
				}
			}
			else
			{
				for(int i = 0; i < elements.size(); i++)
				{
					int idx = elements[i].toInt();
					if(qFind(total_f.begin(),total_f.end(),idx)==total_f.end()||total_f.size()==0)
					{
						total_f.push_back(idx);
						obj->repeatedPatterns[count].push_back(idx);
						obj->patternIdx[idx] = count;
					}
					else
					{
						//const QString w = "Repeated faces!" + QString::number(idx);
						//QMessageBox::warning(NULL,"Warning",w);
					}				
				}
				count++;
			}
			line = in.readLine().simplified();
		}
	}
	//QMessageBox::warning(NULL,"Finished","Successed!");
}

void SymGroup::addOffsetToPatternIdx(int offset)
{
	for (int i=0; i<pIdx1D.size(); i++)
	{
		pIdx1D[i] += offset;
	}

	for (int i=0; i<pIdx2D.size(); i++)
	{
		for (int j=0; j<pIdx2D[i].size(); j++)
		{
			pIdx2D[i][j] += offset;
		}
	}
}

int SymGroup::symLevel(int p1, int p2)
{
	int level = 2;

	if (type == GRID)
	{
		int s1 = -1;
		int s2 = -1;

		for (int i=0; i<pIdx2D.size(); i++)
		{
			for (int j=0; j<pIdx2D[i].size(); j++)
			{
				if (s1!=-1 && s2!=-1)
				{
					break;
				}

				if (pIdx2D[i][j] == p1)
				{
					s1 = i;
				}
				if (pIdx2D[i][j] == p2)
				{
					s2 = i;
				}
			}

			if (s1!=-1 && s2!=-1)
			{
				break;
			}
		}

		if (s1==-1 || s2==-1)
		{
			level = 0;
		}
		else if (s1!=s2)
		{
			level = 1;
		}
	}

	return level;
}
