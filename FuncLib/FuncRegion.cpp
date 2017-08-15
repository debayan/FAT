#include "FuncRegion.h"
#include "Scene.h"

#define PI 3.1415926

FuncRegion::FuncRegion( Object* obj )
{
	object = obj;
	ibs = NULL;
	patternIdx = -1;
}

void FuncRegion::mapFromIbs(IBS* ibs)
{
	this->ibs = ibs;

	ScalarFaceProperty fweight = ibs->mesh->get_face_property<Scalar>("f:weight");
	ScalarFaceProperty fdist = ibs->mesh->get_face_property<Scalar>("f:dist");
	Vector3FaceProperty fnormal = ibs->mesh->get_face_property<Vector3>("f:normal");	// FNORMAL

	determiningScore = QVector<double>(object->samples.size(), 0);
	reverseNormal = QVector<bool>(object->samples.size(), false);
	QVector<bool> isContacting = QVector<bool>(object->samples.size(), false);
	double distThreshold = object->scene->bbox.diagonal().norm() / 100;

	for (int i=0; i<ibs->samplePairs.size(); i++)
	{
		QPair<int,int> pair = ibs->samplePairs[i];

		int sIdx = pair.first;
		if (ibs->pointToCentralObject)
		{
			sIdx = pair.second;
		}

		determiningScore[sIdx] += fweight[Face(i)];

		// check whether this sample is contacting with the interacting object
		if (fdist[Face(i)] < distThreshold)
		{
			isContacting[sIdx] = true;
		}

		// check whether the normal of this sample is pointing towards the interacting object
		bool sameDir = fnormal[Face(i)].dot(object->samples[sIdx].n) > 0;		
		if (ibs->pointToCentralObject)
		{
			reverseNormal[sIdx] =  sameDir ;		
		}	
		else
		{
			reverseNormal[sIdx] = !sameDir;
		}
	}

	sampleIdxs = QVector<int>();
	maxScore = 0;
	for (int i=0; i<determiningScore.size(); i++)
	{
		if (determiningScore[i] > 0)
		{
			//if (isContacting[i])
			//{
			//	determiningScore[i] *= 2;
			//}
			sampleIdxs.push_back(i);
		}

		if (determiningScore[i] > maxScore)
		{
			maxScore = determiningScore[i];
		}
	}
}

void FuncRegion::computeFeature()
{
	computePFH();
	computeDirHist();
	computeHeightHist();
}

void FuncRegion::computePFH()
{
	QVector< QVector<double> > hist;
	double scoreSum = 0;
	for (int j=0; j<sampleIdxs.size(); j++)
	{
		hist << computePfhForSample(j);
		scoreSum += determiningScore[sampleIdxs[j]];
	}	

	QVector<double> mean(hist[0].size(), 0);
	for(int k=0; k<mean.size(); k++)
	{
		for (int j=0; j<hist.size(); j++)
		{
			mean[k] += hist[j][k] * determiningScore[sampleIdxs[j]];
		}

		mean[k] /= scoreSum;
		//mean[k] /= hist.size();
	}

	QVector<double> deviation(hist[0].size(), 0);
	for(int k=0; k<deviation.size(); k++)
	{
		for (int j=0; j<hist.size(); j++)
		{
			deviation[k] += pow(hist[j][k] - mean[k], 2) * determiningScore[sampleIdxs[j]];
		}

		deviation[k] /= scoreSum;
		//diviation[k] /= hist.size();
		deviation[k] = sqrt(deviation[k]);
	}

	pfh.clear();
	pfh << mean << deviation;		
}

QVector<double> FuncRegion::computePfhForSample( int sIdx )
{
	int bNum = 125;
	QVector<double> samplePFH(bNum, 0);
	Vec3d n1 = object->samples[sampleIdxs[sIdx]].n;
	if (reverseNormal[sampleIdxs[sIdx]])
	{
		n1 = -n1;
	}

	double scoreSum = 0;
	for (int j=0; j<sampleIdxs.size(); j++)
	{
		if (j != sIdx)
		{
			int s = sampleIdxs[j];

			Vec3d n2 =  object->samples[s].n;
			if (reverseNormal[s])
			{
				n2 = -n2;
			}
			Vec3d p1p2 = object->samples[s].pos - object->samples[sampleIdxs[sIdx]].pos;
			Vec3d v = p1p2.cross(n1).normalized();
			Vec3d w = n1.cross(v).normalized();
			Vec3d projected_n2 = Vec3d(w.dot(n2), n1.dot(n2), 0);

			double phi = acos(n1.dot(p1p2) / p1p2.norm());
			double dot = n2.dot(v);
			if (dot < -1)
			{
				dot = -1;
			}
			else if (dot > 1)
			{
				dot = 1;
			}
			double alpha = acos(dot);

			Vec3d local_e2(0,1,0);
			double theta = acos(local_e2.dot(projected_n2)/ projected_n2.norm());
			if (projected_n2.norm() < 1.0e-10)
			{
				theta = 0;
			}

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

			if (phiIdx == 5)
			{
				phiIdx = 4;
			}

			if (alphaIdx == 5)
			{
				alphaIdx = 4;
			}

			if (thetaIdx == 5)
			{
				thetaIdx = 4;
			}

			int bIdx = alphaIdx * 25 + thetaIdx * 5 + phiIdx;

			samplePFH[bIdx] += determiningScore[s];
			scoreSum += determiningScore[s];
		}
	}

	for (auto& h:samplePFH)
	{
		h /= scoreSum;
	}

	return samplePFH;
}

void FuncRegion::computeDirHist()
{
	int bNum = 10;
	QVector<double> h(bNum, 0);
	double scoreSum = 0;
	for (int j=0; j<sampleIdxs.size(); j++)
	{
		int s = sampleIdxs[j];
		Vec3d n = object->samples[s].n;
		if (reverseNormal[s])
		{
			n = -n;
		}

		double angle = acos(object->scene->upright.dot(n));

		int bIdx = angle / PI * bNum;
		bIdx = (bIdx > bNum-1)? (bNum - 1):bIdx;

		h[bIdx] += determiningScore[s];
		scoreSum += determiningScore[s];
	}

	for (auto& v:h)
	{
		v /= scoreSum;
	}	

	dirHist = h;
}

void FuncRegion::computeHeightHist()
{
	int bNum = 10;
	QVector<double> h(bNum, 0);
	double scoreSum = 0;
	double range = object->heightRange[1] - object->heightRange[0];		// heightRange[0]: min
	for (int j=0; j<sampleIdxs.size(); j++)
	{
		int s = sampleIdxs[j];
		Vec3d p = object->samples[s].pos;
		double d = p.dot(object->scene->upright);

		int bIdx =  (d - object->heightRange[0]) / range * bNum;
		bIdx = (bIdx > bNum-1)? (bNum-1):bIdx;

		h[bIdx] += determiningScore[s];
		scoreSum += determiningScore[s];
	}

	for (auto& v:h)
	{
		v /= scoreSum;
	}	

	heightHist = h;
}

void FuncRegion::completeDetermingScore()
{
	if (determiningScore.size() < object->samples.size())
	{
		for (int i=determiningScore.size(); i<object->samples.size(); i++)
		{
			determiningScore.push_back(0);	// just for convenience, it should be the weighted score of the new sample's neighbor
			reverseNormal.push_back(false); // should refine accordingly
		}
	}
}