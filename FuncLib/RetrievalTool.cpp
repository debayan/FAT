#include "RetrievalTool.h"
#include "Scene.h"
RetrievalTool::RetrievalTool(QString listName, int datatype)
{
	// listName is scene list (*.sl file)

	// initialize 
	currType = ICON;
	combWeight = 0.5;
	partCombWeight = 0.5;
	combUpdate = false;
	partCombUpdate = false;
	ibshDepth = 1;
	distance.resize(8);
	prShape.resize(8);
	prClass.resize(8);
	prall.resize(8);
	
	// get the filepath
	QFileInfo fileInfo(listName);
	QString baseName = fileInfo.completeBaseName();
	dirName = fileInfo.dir().absolutePath();
	filePath = dirName + "/" + baseName;

	// load all the model names and the corresponding ICON files
	QFile file( listName );
	if (file.open(QIODevice::ReadWrite | QIODevice::Text)) 
	{
		QTextStream in(&file);
		while (!in.atEnd())
		{
			QString line = in.readLine();
			if(line=="create")
			{
				QVector<QString> modelPaths;
				if(datatype==0)
					findFiles(modelPaths,dirName,"txt");
				else
					findFiles(modelPaths,dirName,"off");
				file.resize(0);
				for each(QString s in modelPaths)
					in << s << endl;
				in.seek(0);
				line = in.readLine();
			}

			line = line.trimmed();
			if (line.isEmpty())
			{
				continue;
			}

			QString filename =  dirName + '/' + line;
			scenePaths.push_back(filename);	
			
			QFileInfo fileInfo(filename);
			QString baseName = fileInfo.completeBaseName();
			QString modelDirName = fileInfo.dir().absolutePath();
			
			QPixmap tmpimg;
			if(datatype==0)
				tmpimg = QPixmap(filename + ".png");
			else
				tmpimg = QPixmap(filename + "_pose.jpg");

			sceneImages.push_back(tmpimg);
			sceneNames.push_back(baseName);
			sceneLabels.push_back( QFileInfo(modelDirName).baseName() );
		}

		// count labelNums
		currIdx = 0;
		for each (QString s in sceneLabels)
		{
			bool flag = true;
			for each (QPair<QString,int> *p in labelNums) 
			{
				if (p->first == s)
				{
					p->second ++;
					flag = false;
					break;
				}
			}
			if(flag)
			{
				QPair<QString,int> *p = new QPair<QString,int>;
				p->first = s;
				p->second = 1;
				labelNums.push_back(p);
			}			
		}
	}
}

RetrievalTool::~RetrievalTool(void)
{
	for (int i=0; i<icons.size(); i++)
	{
		if (icons[i])
		{
			if (icons[i]->scene)
			{
				delete icons[i]->scene;
			}
			delete icons[i];
		}
	}
}

void RetrievalTool::findFiles(QVector<QString> &modelPaths, const QString &path, const QString &fileext)  
{
	QDir dir(path);  
	QStringList folders = dir.entryList(QDir::Dirs, QDir::Name);  
	for(int k=0; k<folders.size(); k++)  
	{  
		if( "." == folders[k] ||  ".." == folders[k] )  
			continue;   

		QDir subDir(path + "/" + folders[k]);
		QStringList files = subDir.entryList(QDir::Files);
		for(int i=0; i<files.size(); i++)  
		{  
			QFileInfo fileInfo(files[i]);
			QString ext = fileInfo.suffix();
			if( ext == fileext)
			{
				QString filename = fileInfo.baseName();
				modelPaths.push_back(folders[k] + "/" + filename);
			}
		}  
	}  
}

void RetrievalTool::computeDistance(FEATURE_TYPE type, QString paraStr, double weight)
{
	currType = type;

	if (!paraStr.isEmpty())
	{
		QStringList paraList = paraStr.split("(");
		assert(paraList[1].left(paraList[1].size()-1) == "ICON" || paraList[1].left(paraList[1].size()-1) == "ISET");
		paraStr = "(" + paraList[2];
		if (paraList.size() > 3)
		{
			paraStr += "(" + paraList[3];
		}
	}

	switch (type)
	{
	case ICON:
		{
			paraString = paraStr;			
			computeICONDistance();
			break;
		}
	case LFD:
		computeLFDDistance();
		break;
	case IBSH:
		computeIBSHDistance();
		break;
	case LFDICON:
		{
			if (paraString.isEmpty())
				return;
			computeLFDICONDistance(weight);
			break;
		}
	case POSE:
		computePOSEDistance();
		break;
	case ISET:
		{
			paraString = paraStr;
			computeISETDistance();
			break;
		}
	case GEO:
		computeGEODistance();
		break;
	case GEOICON:
		{
			if (paraString.isEmpty())
				return;
			computeGEOICONDistance(weight);
			break;
		}
	default:
		break;
	}
}

inline bool readFromCSVfile(QString name, Eigen::MatrixXd &matrix)
{
	std::ifstream file(name.toStdString().c_str());
	if(!file)
		return false;
	std::string value,in,line;
	QVector<double> data;
	int rows = 0;
	while(file.good())
	{
		while(std::getline(file,line))
		{
			std::istringstream stream(line);
			while(std::getline(stream,value,','))
			{
				QString qstr = QString::fromStdString(value);
				data.push_back(qstr.toDouble());
			}
			rows++;
		}
	}
	int columns = data.size()/rows;
	matrix = Eigen::MatrixXd::Zero(rows,columns);
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
		{
			matrix(i,j) = data[i*columns + j];
		}
		file.close();
		return true;
}

void RetrievalTool::computePOSEDistance()
{
	int idx = typeIdx(POSE);

	if (distance[idx].rows() > 0)
	{	
		return;
	}

	QString path = filePath + "_POSE_distance.csv";
	QFileInfo file(path);
	if(file.exists())
	{
		readFromCSVfile(path,distance[idx]);
		return;
	}

	QStringList message;
	message << "Pose file missing!";
	int pose_idx = 0;
	poses.resize(scenePaths.size());
	for(auto s: scenePaths)
	{
		QString name = s + ".pose";
		QFile file(name);
		int point_idx = 0;
		double maxX = -999999;
		double minX = 999999;
		//if (file.open(QIODevice::ReadWrite | QIODevice::Text))
		if (file.open(QIODevice::ReadOnly | QIODevice::Text))
		{
			QTextStream in(&file);
			poses[pose_idx] = Eigen::Matrix3Xd::Zero(3,27);
			while (!in.atEnd())
			{
				double x,y,z;
//				char v;
//				in >> v >> x >> y >> z;
				in >> x >> y >> z;
//				if(v!='v')
//					break;
				poses[pose_idx].col(point_idx) = Eigen::Vector3d(x,y,z);
				if(x>maxX)
					maxX = x;
				if(x<minX)
					minX = x;
				point_idx++;
				in.readLine();
			}
		}
		else
		{
			message << name;
			//break;
			continue;
		}
		poses[pose_idx] = poses[pose_idx]/(maxX-minX);
		pose_idx++;
	}

	if (message.size() > 1)
	{
		debugBoxList(message);
		return;
	}

	DistMeasure dist;
	distance[idx] = Eigen::MatrixXd::Zero(scenePaths.size(), scenePaths.size());
	double maxdist = 0;
	for (int i=0; i<scenePaths.size(); i++)
	{
		for (int j=i+1; j<scenePaths.size(); j++)
		{
			distance[idx](i,j) = dist.betweenPoses(poses[i],poses[j]);
			distance[idx](j,i) = distance[idx](i,j);
			if(distance[idx](i,j) > maxdist)
				maxdist = distance[idx](i,j);
		}
	}
	distance[idx] = distance[idx]/maxdist;
	matrixToFile(distance[idx], path);
}

void RetrievalTool::computeICONDistance()
{
	int idx = typeIdx(ICON);

	if (distance[idx].rows() > 0)
	{
		return;
	}

	QString path = filePath + "_ICON_" + paraString + "_distance.csv";
	QFileInfo file(path);
	if(file.exists())
	{
		readFromCSVfile(path, distance[idx]);
		return;
	}


	QVector<InterHierarchy*> tmpH(scenePaths.size(), NULL);
#pragma omp parallel for
	for (int i = 0; i < scenePaths.size(); ++i)
	{
		QString name = dirName + "/" + "(ICON)" + paraString + "/" + sceneNames[i] + ".icon";
		computeFeature(ICON, name);

		tmpH[i] = new InterHierarchy(scenePaths[i], paraString);
		tmpH[i]->load();
	}
	icons << tmpH;

	distance[idx] = Eigen::MatrixXd::Zero(icons.size(), icons.size());
	double maxdist = 0;
	for (int i=0; i<icons.size(); i++)
	{
		for (int j=i+1; j<icons.size(); j++)
		{
			if (icons[i] && icons[j])
			{
				DistMeasure dist(icons[i]->scene->distPara);
				CommonSubtree subtree = dist.betweenInterHierarchies(icons[i], icons[j], false, false);
				distance[idx](i,j) = subtree.distance();
				distance[idx](j,i) = distance[idx](i,j);
				if(distance[idx](i,j) > maxdist)
					maxdist = distance[idx](i,j);
			}
		}
	}
	distance[idx] = distance[idx]/maxdist;
	matrixToFile(distance[idx], path);
}

void RetrievalTool::computeLFDDistance()
{
	int idx = typeIdx(LFD);
	if (distance[idx].rows() > 0)
	{
		return;
	}

	QString path = filePath + "_LFD_distance";
	QFileInfo file(path+".csv");
	if(file.exists())
	{
		readFromCSVfile(path+".csv", distance[idx]);
		return;
	}

	QVector<QString> ObjList, NameList;
	for(auto s: scenePaths)
	{
		QString name = s + "_centric.obj";
		ObjList.push_back(name);

		computeFeature(LFD, name);
	}

	for (int i = 0; i < sceneLabels.size(); i++)
	{
		NameList.push_back(sceneLabels[i] + "_" + sceneNames[i]);
	}

	LFD::calc(path, ObjList, NameList, distance[idx]);
}

void RetrievalTool::computeLFDICONDistance(double weight)
{
	int idx = typeIdx(LFDICON);
	if (distance[idx].rows() > 0 && combWeight==weight)
	{
		combUpdate = false;
		return;
	}

	int idxl = typeIdx(LFD);
	computeLFDDistance();
	if (distance[idxl].rows() == 0)
	{
		return;
	}	

	int idxi = typeIdx(ICON);
	computeICONDistance();
	if (distance[idxi].rows() == 0)
	{
		debugBox("Please select the parameters for ICON first!");
		return;
	}	

	combWeight = weight;
	combUpdate = true;
	distance[idx] = (1.0-combWeight)*distance[idxl] + combWeight*distance[idxi];
}

void RetrievalTool::computeIBSHDistance()
{
	int idx = typeIdx(IBSH);

	if (distance[idx].rows() > 0)
	{
		return;
	}

	QString path = filePath + "_IBSH" + QString::number(ibshDepth) + "_distance.csv";
	QFileInfo file(path);
	if(file.exists())
	{
		readFromCSVfile(path,distance[idx]);
		return;
	}

	if (ibshs.isEmpty())
	{
		for (int i=0; i<scenePaths.size(); i++)
		{
			QString featureFile = scenePaths[i] + ".ibsh";
			computeFeature(IBSH, featureFile);

			ObjLevelFeature f;
			f.load(scenePaths[i]);

			if (!f.featureIBSs.isEmpty())
			{
				ibshs << f;
			}
		}
	}

	if (ibshs.size() != scenePaths.size())
	{
		debugBox("IBSH file missing!");
		return;
	}

	DistMeasure dist;
	distance[idx] = Eigen::MatrixXd::Zero(ibshs.size(), ibshs.size());
	double maxdist = 0;
	double mindist = 1e6;
	for (int i=0; i<ibshs.size(); i++)
	{
		for (int j=i+1; j<ibshs.size(); j++)
		{
			double d = dist.betweenObjLevelFeature(ibshs[i], ibshs[j], ibshDepth);
			distance[idx](i,j) = d;
			distance[idx](j,i) = d;
			if(d > maxdist)
				maxdist = d;
			if(d < mindist)
				mindist = d;
		}
	}

	distance[idx] = distance[idx] - Eigen::MatrixXd::Ones(ibshs.size(), ibshs.size()) * distance[idx].minCoeff();
	for (int i=0; i<ibshs.size(); i++)
	{
		distance[idx](i,i) = 0;
	}
	matrixToFile(distance[idx], path);
}

void RetrievalTool::computeISETDistance()
{
	int idx = typeIdx(ISET);

	if (distance[idx].rows() > 0)
	{
		return;
	}

	QString path = filePath + "_ISET_distance.csv";
	QFileInfo file(path);
	if(file.exists())
	{
		readFromCSVfile(path, distance[idx]);
		return;
	}

	sets.clear();
	for (int i=0; i<scenePaths.size(); i++)
	{
		InterSet set;
		QString name = dirName + "/" + "(ISET)" + paraString + '/' + sceneNames[i];
		computeFeature(ISET, name + ".iset");

		// load iset
		set.load(name);
		sets << set;
	}

	if (sets.size() != scenePaths.size())
	{
		debugBox("ISETs file missing!");
		return;
	}
	
	DistParameter distPara;
	QStringList list = paraString.split("(");
	distPara.readFromString(list[1].left( list[1].size()-1 ));

	DistMeasure dist(distPara);
	distance[idx] = Eigen::MatrixXd::Zero(sets.size(), sets.size());
	double maxdist = 0;
	for (int i=0; i<sets.size(); i++)
	{
		for (int j=i+1; j<sets.size(); j++)
		{
			double d = dist.betweenInterSets(sets[i], sets[j]);
			distance[idx](i,j) = d;
			distance[idx](j,i) = d;
			if(d > maxdist)
				maxdist = d;			
		}
	}
	distance[idx] = distance[idx]/maxdist;
	matrixToFile(distance[idx], path);
}

inline bool sorting(const QPair<int,double> p1, const QPair<int,double> p2)
{
	if(p1.second<p2.second)
		return true;
	return false;
}

QVector<int> RetrievalTool::returnRetrievalResult(int num)
{
	if (scenePaths.size() < num)
	{
		num = scenePaths.size();
	}

	int idx = typeIdx(currType);

	if (distance[idx].rows() == 0)
	{
		return QVector<int>();
	}

	QVector<int> result;
	QVector<QPair<int,double>> tmp_result;
	tmp_result.resize(scenePaths.size());
	result.resize(num);
	Eigen::VectorXd SimilarityLine;
	SimilarityLine = distance[idx].col(currIdx);

	for(int i = 0; i < scenePaths.size(); i++)
		tmp_result[i] = QPair<int,double>(i,SimilarityLine[i]);
	qSort(tmp_result.begin(),tmp_result.end(),sorting);
	for(int i = 0; i < num; i++)
		result[i] = tmp_result[i].first;
	return result;
}

QStringList RetrievalTool::returnRetrievalResultName( int num)
{
	if (scenePaths.size() < num)
	{
		num = scenePaths.size();
	}

	int idx = typeIdx(currType);

	QStringList filenames;

	QVector<QPair<int,double>> tmp_result;
	tmp_result.resize(scenePaths.size());
	Eigen::VectorXd SimilarityLine;
	SimilarityLine = distance[idx].col(currIdx);

	for(int i = 0; i < scenePaths.size(); i++)
		tmp_result[i] = QPair<int,double>(i, SimilarityLine[i]);

	qSort(tmp_result.begin(), tmp_result.end(), sorting);

	for(int i = 0; i < num; i++)
	{
		if ( QFile(scenePaths[tmp_result[i].first]+".obj").exists() )
		{
			filenames.append(scenePaths[tmp_result[i].first]+".obj");
		}
		else if ( QFile(scenePaths[tmp_result[i].first]+".txt").exists() )
		{
			filenames.append(scenePaths[tmp_result[i].first]+".txt");
		}
	
	}

	return filenames;
}

void RetrievalTool::analyzePR(FEATURE_TYPE type)
{
	currType = type;

	int idx = typeIdx(type);
	if (prShape[idx].rows() > 0 && (type!=LFDICON && type!=GEOICON || type==LFDICON && !combUpdate || type==GEOICON && !partCombUpdate))
	{
		return;
	}

	if (distance[idx].rows() == 0)
	{
		return;
	}

	QVector<QVector<QPair<int,double>>> tmp_result;
	tmp_result.resize(scenePaths.size());

	QVector<QVector<double>> R;
	QVector<QVector<double>> P;
	R.resize(scenePaths.size());
	P.resize(scenePaths.size());

	for(int i = 0; i < scenePaths.size(); i++)
	{
		tmp_result[i].resize(scenePaths.size());
		for(int j = 0; j < scenePaths.size(); j++)
		{
			tmp_result[i][j] = QPair<int,double>(j,distance[idx](i,j));
		}
		qSort(tmp_result[i].begin(),tmp_result[i].end(),sorting);
	}

	for(int i = 0; i < scenePaths.size(); i++)
	{
		QString gt = sceneLabels[tmp_result[i][0].first];
		int totalNum;
		for each (QPair<QString,int>* p in labelNums)
			if(p->first==gt)
			{
				totalNum = p->second - 1;
				break;
			}
			for(int j = 1; j < scenePaths.size(); j++)
			{
				if(sceneLabels[tmp_result[i][j].first] == gt)
				{
					if( j==1 )
					{
						P[i].push_back(1.0f);
						R[i].push_back(1.0f/totalNum);
					}
					else
					{
						P[i].push_back((1.0f+P[i][j-2]*(j-1))/j);
						R[i].push_back((R[i][j-2]*totalNum+1.0f)/totalNum);
					}
				}
				else
				{
					if( j==1 )
					{
						P[i].push_back(0.0f);
						R[i].push_back(0.0f);
					}
					else
					{
						P[i].push_back(P[i][j-2]*(j-1)/j);
						R[i].push_back(R[i][j-2]);
					}
				}
			}
			P[i].insert(P[i].begin(),1.0f);
			R[i].insert(R[i].begin(),0.0f);
	}

	QString prName = filePath + "_";;
	switch (type)
	{
	case ICON:
		prName += "ICON_" + paraString;
		break;
	case LFD:
		prName += "LFD";
		break;
	case IBSH:
		prName += "IBSH" + QString::number(ibshDepth);
		break;
	case LFDICON:
		prName += "LFD+ICON_" + QString::number(combWeight) + "_" + paraString;
		break;
	case POSE:
		prName += "POSE";
		break;
	case ISET:
		prName += "SET";
		break;
	case GEO:
		prName += "GEO";
		break;
	case GEOICON:
		prName += "GEO+ICON_" + QString::number(partCombWeight) + "_" + paraString;
		break;
	default:
		break;
	}
	prall[idx] = normalizePR(P,R,idx,prName);
	//	matrixToFile(prall[idx], prName + "_allPR.csv");

	prShape[idx] = Eigen::MatrixXd::Zero(2*scenePaths.size(), scenePaths.size()-1);
	for(int i = 0; i < scenePaths.size() - 1; i++)
	{
		for(int j = 0; j < scenePaths.size(); j++)
		{
			prShape[idx](2*j,i) = 1;
			if(i > R[j].size() - 1)
				continue;
			prShape[idx](2*j,i) = R[j][i];
			prShape[idx](2*j + 1,i) = P[j][i];

			if (prShape[idx](2*j,i) > 1)
			{
				prShape[idx](2*j,i) = 1;
			}
			if (prShape[idx](2*j + 1,i) > 1)
			{
				prShape[idx](2*j + 1,i) = 1;
			}
		}
	}

	prClass[idx] = Eigen::MatrixXd::Zero(2*labelNums.size(),scenePaths.size()-1);
	for(int i = 0; i < scenePaths.size() - 1; i++)
	{
		QVector<double> tmpP,tmpR;
		tmpP.resize(labelNums.size());
		tmpR.resize(labelNums.size());

		QVector<int> maxc;
		maxc.resize(labelNums.size());
		for(int j = 0; j < scenePaths.size(); j++)
		{
			if(i > R[j].size() - 1)
				continue;
			int count = 0;
			for each (QPair<QString,int>* p in labelNums)
			{
				if(p->first==sceneLabels[tmp_result[j][0].first])
					break;
				count++;
			}
			tmpP[count] += R[j][i];
			tmpR[count] += P[j][i];
			maxc[count] = R[j].size();
		}
		for(int j = 0; j < labelNums.size(); j++)
		{
			if(i > maxc[j] - 1)
			{
				prClass[idx](2*j,i) = 1;
				continue;
			}
			prClass[idx](2*j,i) = tmpP[j]/labelNums[j]->second;
			prClass[idx](2*j + 1,i) = tmpR[j]/labelNums[j]->second;
			if(prClass[idx](2*j,i) > 1 + 1e-4)
				break;
		}
	}	

	matrixToFile(prClass[idx], prName + "_PR_category.csv");
	matrixToFile(prall[idx], prName + "_PR_dataset.csv");
}

void addNonReapeatedelement(QVector<double> &v, double ele)
{
	if(ele>1)
		return;
	if(v.size()==0)
	{
		v.push_back(ele);
		return;
	}
	if(qFind(v.begin(),v.end(),ele)==v.end())
		v.push_back(ele);
}

Eigen::MatrixXd RetrievalTool::normalizePR(QVector<QVector<double>> &P,QVector<QVector<double>> &R,int idx,QString fn)
{
	QVector<QVector<double>> normalizedP;
	QVector<QVector<double>> normalizedR;

	normalizedR.resize(scenePaths.size());
	normalizedP.resize(scenePaths.size());
	QVector<double> totalR;
	for (int i = 0; i < scenePaths.size(); i++)
	{
		QVector<int> count_P;
		for (int j = 0; j < R[i].size(); j++)
		{
			if(normalizedR[i].size()==0)
			{
				normalizedR[i].push_back(R[i][j]);
				normalizedP[i].push_back(P[i][j]);
				count_P.push_back(1);
				addNonReapeatedelement(totalR,R[i][j]);
				continue;
			}
			if(qFind(normalizedR[i].begin(),normalizedR[i].end(),R[i][j])!=normalizedR[i].end())
			{
				int idx = qFind(normalizedR[i].begin(),normalizedR[i].end(),R[i][j]) - normalizedR[i].begin();
				normalizedP[i][idx] += P[i][j];
//				if(P[i][j] > normalizedP[i][idx])
//					normalizedP[i][idx] = P[i][j];
				count_P[idx] ++;
			}
			else
			{
				normalizedR[i].push_back(R[i][j]);
				normalizedP[i].push_back(P[i][j]);
				count_P.push_back(1);
				addNonReapeatedelement(totalR,R[i][j]);
			}
		}
		for (int j = 0; j < normalizedR[i].size(); j++)
		{
			normalizedP[i][j] = normalizedP[i][j]/count_P[j];
		}
	}
	//////////////////////////////////////////////////////////////////////////
	qSort(totalR.begin(),totalR.end());
	QVector<double> totalP;
	QVector<int> totalP_num;
	totalP.resize(totalR.size());
	totalP_num.resize(totalR.size());
	QVector<QVector<double>> interP;
	interP.resize(scenePaths.size());

	for (int i = 0; i < scenePaths.size(); i++)
	{
		for (int j = 0; j < totalR.size(); j++)
		{
			for (int k = 0; k < normalizedR[i].size()-1; k++)
			{
				if(totalR[j] >= normalizedR[i][k] && totalR[j] < normalizedR[i][k+1])
				{
					totalP[j] += normalizedP[i][k] + 
						(normalizedP[i][k+1]-normalizedP[i][k])*(totalR[j]-normalizedR[i][k])/(normalizedR[i][k+1]-normalizedR[i][k]);
					break;
				}
			}
		}
	}

	for (int i = 0; i < totalP.size(); i++)
		totalP[i]/=scenePaths.size();

	int max_num = 0;
	for (int i = 0; i < normalizedR.size(); i++)
		if(normalizedR[i].size()>max_num)
			max_num = normalizedR[i].size();
	Eigen::MatrixXd all_M = Eigen::MatrixXd::Ones(2*scenePaths.size(),max_num)*-1;
	for (int i = 0; i < normalizedR.size(); i++)
	{
		for (int j = 0; j < normalizedR[i].size(); j++)
		{
			all_M(2*i,j) = normalizedP[i][j];
			all_M(2*i+1,j) = normalizedR[i][j];
		}
	}
	Eigen::MatrixXd class_M = Eigen::MatrixXd::Ones(2*scenePaths.size(),max_num)*-1;
	for (int i = 0; i < normalizedR.size(); i++)
	{

	}
	matrixToFile(all_M, fn + "_PR_shape.csv");
	//////////////////////////////////////////////////////////////////////////
	Eigen::MatrixXd avg_pr = Eigen::MatrixXd::Zero(2, totalR.size());
	for (int i = 0; i < totalP.size(); i++)
	{
		avg_pr(1,i) = totalP[i];
		avg_pr(0,i) = totalR[i];
	}
	P.clear();
	R.clear();
	P = normalizedP;
	R = normalizedR;
	return avg_pr;
}

Eigen::MatrixXd RetrievalTool::getPRShape()
{
	int idx = typeIdx(currType);

	Eigen::MatrixXd pr = Eigen::MatrixXd::Zero(2, scenePaths.size()-1);
	
	for (int j=0; j<scenePaths.size()-1; j++)
	{
		pr(0,j) = prShape[idx](2*currIdx, j);
		pr(1,j) = prShape[idx](2*currIdx+1, j);
	}

	return pr;
}

Eigen::MatrixXd RetrievalTool::getPRall()
{
	int idx = typeIdx(currType);
	return prall[idx];
}

Eigen::MatrixXd RetrievalTool::getPRClass()
{
	int idx = typeIdx(currType);

	int classIdx = 0;
	for each (QPair<QString,int>* p in labelNums)
	{
		if( p->first == sceneLabels[currIdx] )
		{
			break;
		}
		classIdx++;
	}

	Eigen::MatrixXd pr = Eigen::MatrixXd::Zero(2*(1+labelNums[classIdx]->second), scenePaths.size()-1);

	for (int j=0; j<scenePaths.size()-1; j++)
	{
		pr(0,j) = prClass[idx](2*classIdx, j);
		pr(1,j) = prClass[idx](2*classIdx+1, j);
	}

	int k = 1;
	for (int i=0; i<scenePaths.size(); i++)
	{
		if (sceneLabels[i] == sceneLabels[currIdx])
		{
			for (int j=0; j<scenePaths.size()-1; j++)
			{
				pr(2*k,j) = prShape[idx](2*i, j);
				pr(2*k+1,j) = prShape[idx](2*i+1, j);
			}

			k++;
		}
	}

	return pr;
}

bool RetrievalTool::prReady()
{
	int idx = typeIdx(currType);

	if (prShape[idx].rows() > 0 )
	{
		return true;
	}
	else
	{
		return false;
	}
}

int RetrievalTool::typeIdx( FEATURE_TYPE type )
{
	int idx;
	switch (type)
	{
	case ICON:
		idx = 0;
		break;
	case LFD:
		idx = 1;
		break;
	case IBSH:
		idx = 2;
		break;
	case LFDICON:
		idx = 3;
		break;
	case POSE:
		idx = 4;
		break;
	case ISET:
		idx = 5;
		break;
	case GEO:
		idx = 6;
		break;
	case GEOICON:
		idx = 7;
		break;
	default:
		break;
	}

	return idx;
}

void RetrievalTool::computeGEOICONDistance(double weight)
{
	int idx = typeIdx(GEOICON);
	if (distance[idx].rows() > 0 && partCombWeight==weight)
	{
		partCombUpdate = false;
		return;
	}

	int idxg = typeIdx(GEO);
	computeGEODistance();
	if (distance[idxg].rows() == 0)
	{
		return;
	}	

	int idxi = typeIdx(ICON);
	computeICONDistance();
	if (distance[idxi].rows() == 0)
	{
		debugBox("Please select the parameters for ICON first!");
		return;
	}	

	partCombWeight = weight;
	partCombUpdate = true;
	distance[idx] = (1.0-partCombWeight)* distance[idxg] + partCombWeight*distance[idxi];
}

void RetrievalTool::computeGEODistance()
{
	int idx = typeIdx(GEO);
	if (distance[idx].rows() > 0)
	{
		return;
	}

	QString path = filePath + "_GEO_distance.csv";
	QFileInfo file(path);
	if(file.exists())
	{
		readFromCSVfile(path, distance[idx]);
	}
	else
	{
		QVector<ObjGeoFeature> geoFeatures;
		for (int i=0; i<scenePaths.size(); i++)
		{
			QString featureFile = scenePaths[i] + ".gf";		
			computeFeature(GEO, featureFile);

			ObjGeoFeature f;
			f.load(scenePaths[i]);
			geoFeatures << f;
		}

		DistMeasure dist;
		distance[idx] = Eigen::MatrixXd::Zero(geoFeatures.size(), geoFeatures.size());
		double maxdist = 0;
		for (int i=0; i<geoFeatures.size(); i++)
		{
			for (int j=i+1; j<geoFeatures.size(); j++)
			{
				double d = dist.betweenGeoFeatures(geoFeatures[i], geoFeatures[j]);
				distance[idx](i,j) = d;
				distance[idx](j,i) = d;
				if(d > maxdist)
					maxdist = d;
			}
		}

		distance[idx] /= maxdist;
//		distance[idx] = dist.betweenGeoFeatures(geoFeatures);

		matrixToFile(distance[idx], path);
	}
}

void RetrievalTool::computeFeature(FEATURE_TYPE type, QString featureFile)
{
	QFile file(featureFile);
	if (!QFile::exists(featureFile))
	{
		QFileInfo info(file);
		QString sceneName, baseFileName;

		if (type == LFD)								// for LFD
		{
			QStringList centralBaeeFile;
			QStringList temp = info.baseName().split("_");
			for (int i = 0; i < temp.size() - 1; ++i)
			{
				centralBaeeFile << temp[i];
			}
			baseFileName = centralBaeeFile.join("_");
		}
		else
		{
			baseFileName = info.baseName();
		}

		int idx = -1;
		for (int i = 0; i < scenePaths.size(); ++i)
		{
			QFileInfo tempInfo(scenePaths[i]);
			if (baseFileName == tempInfo.baseName())
			{
				idx = i;
				break;
			}
		}
		assert(idx != -1);
		sceneName = scenePaths[idx] + ".txt";

		Scene *s = new Scene();
		s->load(sceneName);

		switch (type)
		{
		case ICON:
			{
				QStringList list = paraString.split("(");
				s->distPara.readFromString(list[1].left( list[1].size()-1 ));
				if (list.size() > 2)
				{
					s->hierPara.readFromString(list[2].left( list[2].size()-1 ));
				}
				else
				{
					s->hierPara.readFromString(QString());
				}
				s->buildInteractionHierarchy();
				break;
			}
		case LFD:
			{
				s->outputCentricObject();
				break;
			}
		case IBSH:
			{
				s->visualizeHierarchy = true;
				s->computeIBSH();
				break;
			}
		case LFDICON:
			{
				break;
			}
		case POSE:
			{
				// need to load .pose file manually
				break;
			}
		case ISET:
			{
				QStringList list = paraString.split("(");
				s->distPara.readFromString(list[1].left( list[1].size()-1 ));
				s->outputInterSetFeature();
				break;
			}
		case GEO:
			{
				s->computeGeometryFeature(64);
				break;
			}
		case GEOICON:
			{
				break;
			}
		default:
			break;
		}

		delete s;
	}
}

void RetrievalTool::computeICONFeatures(DistParameter distPara, HierParameter hierPara)
{
	paraString = distPara.toString() + hierPara.toString();
#pragma omp parallel for
	for (int i = 0; i < scenePaths.size(); ++i)
	{
		QString name = dirName + "/" + "(ICON)" + paraString + "/" + sceneNames[i] + ".icon";
		computeFeature(ICON, name);
	}
}


void RetrievalTool::computeISETFeatures(DistParameter distPara)
{
	paraString = distPara.toString();
#pragma omp parallel for
	for (int i=0; i<scenePaths.size(); i++)
	{
		QString name = dirName + "/" + "(ISET)" + paraString + '/' + sceneNames[i] + ".iset";
		computeFeature(ISET, name);
	}
}

void RetrievalTool::setCurrentRetrievalMode(RETRIEVAL_MODE mode)
{
	currentRetrievalMode = mode;
}