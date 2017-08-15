#pragma once
#include "tree.h"
#include "qfile.h"
#include "Interaction.h"
#include "SceneHierarchy.h"

struct Visualization
{
	static void visualizeTree(tree<Interaction*> *hierarchy, QString path)
	{
		QFile file(path + ".dot");
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream out(&file);

		out<<"digraph G{"<<endl;

		tree<Interaction*>::sibling_iterator children;
		tree<Interaction*>::iterator iterator;
		iterator = hierarchy->begin();
		while(iterator!= hierarchy->end())
		{
			QString f = "", fl = "";
			for each(int index in (*iterator)->obj->origIdx)
			{
				f += "n" + QString::number(index);
				//fl += QString::number(index)+",";
			}
			//fl.chop(1);

			if ((*iterator)->obj->origIdx.size() > 1)
			{
				//fl += " (" + QString::number((*iterator)->mergeDist) + ")";
				fl = QString::number((*iterator)->mergeDist);
			}
			else
			{
				fl = QString::number((*iterator)->obj->origIdx[0]);
			}

			out << f << "[label=\"" << fl << "\"];";
			if(tree<Interaction*>::number_of_children(iterator)==0)
			{
				out << f << "[style=\"filled\",color=\"skyblue\"];" << endl;
			}

			children = hierarchy->begin(iterator);
			while(children != hierarchy->end(iterator))
			{
				QString  c = "", cl = "";
				for each(int index in (*children)->obj->origIdx)
				{
					c += "n" + QString::number(index);
					//cl += QString::number(index)+",";
				}
	
				//cl.chop(1);

				if ((*children)->obj->origIdx.size() > 1)
				{
					//cl += " (" + QString::number((*children)->mergeDist) + ")";
					cl = QString::number((*children)->mergeDist);
				}
				else
				{
					cl = QString::number((*children)->obj->origIdx[0]);
				}

				out << c << "[label=\"" << cl << "\"];";
				out << f << "->" << c << "[dir=back];" << endl;
				if(tree<Interaction*>::number_of_children(children)==0)
				{
					out << c << "[style=\"filled\",color=\"skyblue\"];" << endl;
				}
				++children;
			}
			++iterator;
		}

		out<<"}"<<endl;

		file.close();
		
		// convert to image
		QString commandq = QDir::currentPath() + "/../UtilityLib/graphviz/bin/dot.exe -Tpng " + path + ".dot -o " + path + ".png";
		QByteArray ba = commandq.toLocal8Bit();
		char* command = ba.data();
		system(command);
	}

	static void visualizaSubtree(tree<Interaction*> *hierarchy, QVector<tree<Interaction*>::iterator> matchedNode, QVector<double> pairSimilarity, QString path)
	{
		QFile file(path + ".dot");
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream out(&file);

		out<<"digraph G{"<<endl;

		tree<Interaction*>::sibling_iterator children;
		tree<Interaction*>::iterator iterator;
		iterator = hierarchy->begin();
		while(iterator!= hierarchy->end())
		{
			QString f = "", fl = "";
			for each(int index in (*iterator)->obj->origIdx)
			{
				f += "n" + QString::number(index);
			}
			int matchIdx1 = -1;
			for (int i=0; i<matchedNode.size(); i++)
			{
				if ( (*iterator) == (*matchedNode[i]) )
				{
					matchIdx1 = i;
					break;
				}
			}
			if ( matchIdx1 != -1 )
			{
				fl = QString::number(matchIdx1);// + ": " + QString::number(pairSimilarity[matchIdx1]);
			}


			children = hierarchy->begin(iterator);
			while(children != hierarchy->end(iterator))
			{
				QString c = "", cl = "";
				for each(int index in (*children)->obj->origIdx)
				{
					c += "n" + QString::number(index);
				}						

				int matchIdx2 = -1;
				for (int i=0; i<matchedNode.size(); i++)
				{
					if ( (*children) == (*matchedNode[i]) )
					{
						matchIdx2 = i;
						break;
					}
				}
				if ( matchIdx2 != -1 )
				{
					cl = QString::number(matchIdx2);//+ ": " + QString::number(pairSimilarity[matchIdx2]);;
				}
	
				out << f << "[label=\"" << fl << "\"];";
				out << c << "[label=\"" << cl << "\"];";
				out << f << "->" << c << "[dir=back];" << endl;
				if (matchIdx2 != -1 )
				{
					out << c << "[style=\"filled\",color=\"darkorange\"];" << endl;
				}
				++children;
			}

			if( matchIdx1 != -1 )
			{
				out << f << "[style=\"filled\",color=\"darkorange\"];" << endl;
			}

			++iterator;
		}

		out<<"}"<<endl;

		file.close();

		// convert to image
		QString commandq = QDir::currentPath() + "/../UtilityLib/graphviz/bin/dot.exe -Tpng " + path + ".dot -o " + path + ".png";
		QByteArray ba = commandq.toLocal8Bit();
		char* command = ba.data();
		system(command);
	}
	
	static void visualizeTree(tree<ObjCluster> *hierarchy, QString path)
	{
		QFile file(path + ".dot");
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream out(&file);

		out<<"digraph G{"<<endl;

		tree<ObjCluster>::sibling_iterator children;
		tree<ObjCluster>::iterator iterator;
		iterator = hierarchy->begin();
		while(iterator!= hierarchy->end())
		{
			QString f = "", fl = "";
			for each(int index in (*iterator).origObjIdx)
			{
				f += "n" + QString::number(index);
				fl += QString::number(index)+",";
			}
			fl.chop(1);

			out << f << "[label=\"" << fl << "\"];";
			if(tree<ObjCluster>::number_of_children(iterator)==0)
			{
				out << f << "[style=\"filled\",color=\"skyblue\"];" << endl;
			}

			children = hierarchy->begin(iterator);
			while(children != hierarchy->end(iterator))
			{
				QString  c = "", cl = "";
				for each(int index in (*children).origObjIdx)
				{
					c += "n" + QString::number(index);
					cl += QString::number(index)+",";
				}
				cl.chop(1);

				out << c << "[label=\"" << cl << "\"];";
				out << f << "->" << c << "[dir=back];" << endl;
				if(tree<ObjCluster>::number_of_children(children)==0)
				{
					out << c << "[style=\"filled\",color=\"skyblue\"];" << endl;
				}
				++children;
			}
			++iterator;
		}

		out<<"}"<<endl;

		file.close();

		// convert to image
		QString commandq = QDir::currentPath() + "/../UtilityLib/graphviz/bin/dot.exe -Tpng " + path + ".dot -o " + path + ".png";
		QByteArray ba = commandq.toLocal8Bit();
		char* command = ba.data();
		system(command);
	}
};

