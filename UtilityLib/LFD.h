#pragma once

#include <QtCore>
#include <Eigen\Dense>
#include "UtilityGlobal.h"

class LFD
{
public:
	static int calc(const QString filename, const QVector<QString> ObjList, const QVector<QString> NameList, Eigen::MatrixXd &DistM)
	{
		QString filepath = QDir::currentPath() + "/../UtilityLib";

		QString ListPath = filepath + "/LFD/List.txt";
		QFile file(ListPath);
		file.open(QIODevice::ReadWrite | QIODevice::Text);
		QTextStream input(&file);

		// create the cache folder if not exists
		QString cacheFolder = filepath + "/LFD/Cache/";
		QDir dir(cacheFolder);
		if (!dir.exists())
		{
			dir.mkpath(cacheFolder);
		}

		input << "L " << filename << endl;
		input << "N " << ObjList.size() << endl;
		input << "D " << cacheFolder << endl;
		for (int i = 0; i < ObjList.size(); i++)
		{
			input << ObjList[i] << endl;
			input << NameList[i] << endl;
		}

		file.close();

		QString commandq = filepath + "/LFD";
		LPCWSTR command_l = (const wchar_t*) commandq.utf16();
		wchar_t buf[1000];
		GetCurrentDirectory(1000,buf);
		SetCurrentDirectory(command_l);
		system("LFD.exe");
		SetCurrentDirectory(buf);

		///////////////////////////////////////////////////////////////////
		DistM = Eigen::MatrixXd::Zero(ObjList.size(),ObjList.size());
		QFile fileDist(filename+"_dist_mat.txt");
		fileDist.open(QIODevice::ReadWrite | QIODevice::Text);
		QTextStream output(&fileDist);

		QString line = output.readLine();
		line = output.readLine();
		line = output.readLine();
		int count = 0;
		int linewidth_pre = 0;
		int linewidth = ObjList.size() - 1;
		int linehight = 0;
		double maxDist = 0; 
		while(!line.isNull())
		{
			double dist = line.toDouble();
			if(dist > maxDist)
				maxDist = dist;
			int indexX,indexY;
			if(count < linewidth)
			{
				indexX = linehight;
				indexY = count - linewidth_pre + linehight + 1;
			}
			else
			{
				int tmp = linewidth - linewidth_pre;
				linewidth_pre = linewidth;
				linewidth += tmp - 1;
				linehight++;
				indexX = linehight;
				indexY = count - linewidth_pre + linehight + 1;
			}

			DistM(indexX,indexY) = dist;
			DistM(indexY,indexX) = dist;
			count ++;
			line = output.readLine();
		}
		DistM = DistM/maxDist;
		matrixToFile(DistM, filename+".csv");

		return 1;
	}

};