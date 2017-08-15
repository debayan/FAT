
#include "DistMeasure.h"
#include "Scene.h"
#include "histogramDistance.h"
#include "Visualization.h"
#include "UtilityGlobal.h"

double DistMeasure::betweenPoses(Eigen::Matrix3Xd X,Eigen::Matrix3Xd Y)
{
	Eigen::Vector3d src_mean, dst_mean;
	for(unsigned int i=0; i<3; ++i) {
		src_mean(i) = (X.row(i).array().transpose().array()).sum();
		dst_mean(i) = (Y.row(i).array().transpose().array()).sum();
	}
	X.colwise() -= src_mean;
    Y.colwise() -= dst_mean;
	Eigen::Matrix3d sigma = X * Y.transpose();
	Eigen::JacobiSVD<Eigen::Matrix3d,Eigen::FullPivHouseholderQRPreconditioner> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Affine3d transformation;
	if(svd.matrixU().determinant()*svd.matrixV().determinant()<0.0) {
		Eigen::Vector3d S = Eigen::Vector3d::Ones(); S(2) = -1.0;
		transformation.linear().noalias() = svd.matrixV()*S.asDiagonal()*svd.matrixU().transpose();
	} else {
		transformation.linear().noalias() = svd.matrixV()*svd.matrixU().transpose();
	}
	transformation.translation().noalias() = dst_mean - transformation.linear()*src_mean;
	X = transformation.linear()*X;
	X.colwise() += src_mean + transformation.translation();
	Y.colwise() += dst_mean;

	Eigen::Matrix3Xd meanDist = X - Y;
	double dist = 0;
	for(int i = 0; i < meanDist.cols(); i++)
	{
		dist += meanDist(0,i)*meanDist(0,i) + meanDist(1,i)*meanDist(1,i) + meanDist(2,i)*meanDist(2,i);
	}
	return dist;
}

double DistMeasure::betweenInteractions( Interaction* iter1, Interaction* iter2, bool insideScene)
{
	double ibs_geo = computeIbsGeomDistance(iter1->ibs, iter2->ibs);
	int ibs_top = computeIbsTopoDistance(iter1->ibs, iter2->ibs);
	double region = computeRegionDistance(iter1->region, iter2->region);

	double dist;

	if (distPara.useTopo && ibs_top == 1)
	{
		ibs_geo = 1;
	}
		
	dist = distPara.weight[0] * ibs_geo + distPara.weight[1] * region;

	if (insideScene && distPara.symRatio != 1) 
	{
		int level = symLevel(iter1, iter2);

		dist /= pow(distPara.symRatio, level);
	}

	return dist;
}

int DistMeasure::symLevel(Interaction* iter1, Interaction* iter2) // for now, assume each interaction only have one interacting object: only leaf node in the interHierarchy
{
	int level = 0;

	if (iter1->region->object->symGroups.isEmpty() || iter2->region->object->symGroups.isEmpty())
	{
		return level;
	}

	if (iter1->region->object == iter2->region->object)
	{
		int pIdx1 = iter1->region->patternIdx;
		int pIdx2 = iter2->region->patternIdx;

		level = iter1->region->object->symLevel(pIdx1, pIdx2);
	}

	return level;
}

double DistMeasure::computeIbsGeomDistance( IBS* ibs1, IBS* ibs2 )
{
	double d_pfh = histogramDistance(ibs1->pfh, ibs2->pfh, distPara.distType);
	double d_dir = histogramDistance(ibs1->dirHist, ibs2->dirHist, distPara.distType);
	double d_dist = histogramDistance(ibs1->distHist, ibs2->distHist, distPara.distType);

	double d_geo = distPara.ibsWeight[0] * d_pfh + distPara.ibsWeight[1] * d_dir + distPara.ibsWeight[2] * d_dist; 	

	return d_geo;
}

int DistMeasure::computeIbsTopoDistance( IBS* ibs1, IBS* ibs2 )
{
	int d_topo = 0;

	for (int k=1; k<3; k++)
	{
		if (ibs1->bettiNumbers[k] != ibs2->bettiNumbers[k])
		{
			d_topo = 1;
			break;
		}
	}

	return d_topo;
}

double DistMeasure::computeRegionDistance( FuncRegion* r1, FuncRegion* r2 )
{
	double d_pfh = histogramDistance(r1->pfh, r2->pfh, distPara.distType);
	double d_dir = histogramDistance(r1->dirHist, r2->dirHist, distPara.distType);
	double d_height = histogramDistance(r1->heightHist, r2->heightHist, distPara.distType);

	double d_geo = distPara.regionWeight[0] * d_pfh + distPara.regionWeight[1] * d_dir + distPara.regionWeight[2] * d_height;

	return d_geo;
}

CommonSubtree DistMeasure::betweenInterHierarchies( InterHierarchy* hier1, InterHierarchy* hier2, bool visualize, bool updateMatching)
{
	int idx1 = -1;
	int idx2 = -1;
	double dist = 1;
	CommonSubtree bestMatch(hier1->interTree[0], hier2->interTree[0]);
	for (int i=0; i<hier1->interTree.size(); i++)
	{
		for (int j=0; j<hier2->interTree.size(); j++)
		{
			CommonSubtree subtree(hier1->interTree[i], hier2->interTree[j]);
			subtree.findMaxMatch();
			double d = subtree.distance();

			if (d < dist)
			{
				dist = d;
				idx1 = i;
				idx2 = j;
				bestMatch = subtree;
			}
		}
	}


	if (visualize)
	{
		QString name = hier1->scene->filename + "_" + QString::number(idx1) + "-" + hier2->scene->name + "_" + QString::number(idx2);
		if (!QFile(name + ".png").exists() || updateMatching)
		{
			QVector<tree<Interaction*>::iterator> matchedNode;
			for (int i=0; i<bestMatch.mapping.size(); i++)
			{
				matchedNode << bestMatch.mapping[i].first;
			}
			Visualization::visualizaSubtree(hier1->interTree[idx1], matchedNode, bestMatch.pairSimilarity, name);
		}		

		name = hier2->scene->filename + "_" + QString::number(idx2) + "-" + hier1->scene->name + "_" + QString::number(idx1);
		if (!QFile(name + ".png").exists() || updateMatching)
		{
			QVector<tree<Interaction*>::iterator> matchedNode;
			for (int i=0; i<bestMatch.mapping.size(); i++)
			{
				matchedNode << bestMatch.mapping[i].second;
			}			
			Visualization::visualizaSubtree(hier2->interTree[idx2], matchedNode, bestMatch.pairSimilarity, name);
		}
	}

	return bestMatch;
}

double DistMeasure::betweenObjLevelFeature(ObjLevelFeature f1, ObjLevelFeature f2, int level)
{
	if (f1.featureIBSs.size() < level)
	{
		level = f1.featureIBSs.size();
	}
	if (f2.featureIBSs.size() < level)
	{
		level = f2.featureIBSs.size();
	}

	double s = 0;
	for (int i=0; i<level; i++)
	{
		double s12 = 0;
		for (int j=0; j<f1.featureIBSs[i].size(); j++)
		{
			for (int k=0; k<f2.featureIBSs[i].size(); k++)
			{
				double d_pfh = histogramDistance(f1.featureIBSs[i][j].pfh, f2.featureIBSs[i][k].pfh, L1);
				double d_dir = histogramDistance(f1.featureIBSs[i][j].dirHist, f2.featureIBSs[i][k].dirHist, L1);
				double d_dist = histogramDistance(f1.featureIBSs[i][j].distHist, f2.featureIBSs[i][k].distHist, L1);

				s12 += 0.1 * (1-d_pfh) + 0.5 * (1-d_dir) + 0.4 * (1-d_dist);
			}
		}

		double s11 = 0;
		for (int j=0; j<f1.featureIBSs[i].size(); j++)
		{
			for (int k=0; k<f1.featureIBSs[i].size(); k++)
			{
				double d_pfh = histogramDistance(f1.featureIBSs[i][j].pfh, f1.featureIBSs[i][k].pfh, L1);
				double d_dir = histogramDistance(f1.featureIBSs[i][j].dirHist, f1.featureIBSs[i][k].dirHist, L1);
				double d_dist = histogramDistance(f1.featureIBSs[i][j].distHist, f1.featureIBSs[i][k].distHist, L1);

				s11 += 0.1 * (1-d_pfh) + 0.5 * (1-d_dir) + 0.4 * (1-d_dist);
			}
		}

		double s22 = 0;
		for (int j=0; j<f2.featureIBSs[i].size(); j++)
		{
			for (int k=0; k<f2.featureIBSs[i].size(); k++)
			{
				double d_pfh = histogramDistance(f2.featureIBSs[i][j].pfh, f2.featureIBSs[i][k].pfh, L1);
				double d_dir = histogramDistance(f2.featureIBSs[i][j].dirHist, f2.featureIBSs[i][k].dirHist, L1);
				double d_dist = histogramDistance(f2.featureIBSs[i][j].distHist, f2.featureIBSs[i][k].distHist, L1);

				s22 += 0.1 * (1-d_pfh) + 0.5 * (1-d_dir) + 0.4 * (1-d_dist);
			}
		}

		double s_i = 0;
		if (s11>=s22 && s11!=0)
		{
			s_i = s12 / s11;
		}
		else
		{
			s_i = s12 / s22;
		}

		s += pow(0.5, i) * s_i;

		if (s>1)
		{
			qDebug() << s12 << " " << s11 << " " << s22 << endl;
		}
	}

	double d = 1 - s;

	return d;
}

double DistMeasure::betweenInterSets(InterSet set1, InterSet set2)
{
	double s12 = 0;
	for (int j=0; j<set1.interactions.size(); j++)
	{
		for (int k=0; k<set2.interactions.size(); k++)
		{
			double d = betweenInteractions(set1.interactions[j], set2.interactions[k]);
			s12 += 1 - d;
		}
	}

	double s11 = 0;
	for (int j=0; j<set1.interactions.size(); j++)
	{
		for (int k=0; k<set1.interactions.size(); k++)
		{
			double d = betweenInteractions(set1.interactions[j], set1.interactions[k]);

			s11 += 1 - d;
		}
	}

	double s22 = 0;
	for (int j=0; j<set2.interactions.size(); j++)
	{
		for (int k=0; k<set2.interactions.size(); k++)
		{
			double d = betweenInteractions(set2.interactions[j], set2.interactions[k]);

			s22 += 1 - d;
		}
	}

	double s = 0;
	if (s11>=s22 && s11!=0)
	{
		s = s12 / s11;
	}
	else
	{
		s = s12 / s22;
	}


	double d =  1 - s;

	return d;
}

double DistMeasure::betweenGeoFeatures(ObjGeoFeature f1, ObjGeoFeature f2)
{
	double d1 = abs(f1.partSize - f2. partSize);
	//double d2 = histogramDistance(f1.heightDistribution, f2.heightDistribution, L1);
	double d3 = histogramDistance(f1.distDistribution, f2.distDistribution, L1);
	double d4 = (f1.LPS - f2.LPS).norm();

	//double d = (d1 + d2 + d3 + d4) / 4;
	double d = (d1 + d3 + d4) / 3;

	return d;
}
Eigen::MatrixXd DistMeasure::betweenGeoFeatures(QVector<ObjGeoFeature> features)
{
	Eigen::MatrixXd dist1 = Eigen::MatrixXd::Zero(features.size(), features.size());
	Eigen::MatrixXd dist2 = Eigen::MatrixXd::Zero(features.size(), features.size());
	Eigen::MatrixXd dist3 = Eigen::MatrixXd::Zero(features.size(), features.size());
	Eigen::MatrixXd dist = Eigen::MatrixXd::Zero(features.size(), features.size());
	double maxdist1 = 0,maxdist2 = 0,maxdist3 = 0;
	for (int i=0; i<features.size(); i++)
	{
		for (int j=i+1; j<features.size(); j++)
		{
			double d1 = abs(features[i].partSize - features[j]. partSize);
			double d2 = histogramDistance(features[i].distDistribution, features[j].distDistribution, L1);
			double d3 = (features[i].LPS - features[j].LPS).norm();

			dist1(i,j) = d1;
			dist2(i,j) = d2;
			dist3(i,j) = d3;			
			if(d1 > maxdist1)
				maxdist1 = d1;
			if(d2 > maxdist2)
				maxdist2 = d2;
			if(d3 > maxdist3)
				maxdist3 = d3;
		}
	}
	maxdist1 = maxdist1/2;
	maxdist2 = maxdist2/2;
	maxdist3 = maxdist3/2;
	double maxdist = 0;
	for (int i=0; i<features.size(); i++)
	{
		for (int j=i+1; j<features.size(); j++)
		{
			double d = (std::exp(-dist1(i,j)*dist1(i,j)/maxdist1) + std::exp(-dist2(i,j)*dist2(i,j)/maxdist2) + std::exp(-dist3(i,j)*dist3(i,j)/maxdist3))/3.0f;
			dist(i,j) = d;
			dist(j,i) = d;
			if(dist(i,j)>maxdist)
				maxdist = dist(i,j);
		}
	}
	dist = dist/maxdist;
	dist = Eigen::MatrixXd::Ones(features.size(), features.size()) - dist;
	for (int i=0; i<features.size(); i++)
		dist(i,i) = 0;
	return dist;
}

double DistMeasure::betweenSymGroups(SymGroup sg1, SymGroup sg2)
{
	double d = 0;

	// to do

	return d;
}


