#pragma once

#include "emd.h"
#include <QVector>


//enum DIST_TYPE {L1, EMD};

float dist(feature_t *F1, feature_t *F2) 
{ 
	return abs(*F1-*F2); 
}

signature_t convertSignal(QVector<double> h)
{
	QVector<feature_t> f;
	QVector<float> w;

	for (int i=0; i<h.size(); i++)
	{
		if (h[i] > 0)
		{
			f << i;
			w << h[i];
		}
	}

	feature_t* feature= new feature_t[f.size()];
	float* weight = new float[f.size()];
	for (int i=0; i<f.size(); i++)
	{
		feature[i] = f[i];
		weight[i] = w[i];
	}

	signature_t s = {f.size(), feature, weight};

	return s;
}

double computeL1Distance(QVector<double> h1, QVector<double> h2)
{
	assert(h1.size() == h2.size());

	double d = 0;
	for (int i=0; i<h1.size(); i++)
	{
		d += abs(h1[i] - h2[i]);
	}
	d /= (2*sqrt(2)); // normalize the distance into [0,1]; for normalized vectors, the maximal L1 distance is 2*sqrt(2)

	return d;
}

double computeEMD(QVector<double> h1, QVector<double> h2)
{
	signature_t s1 = convertSignal(h1);
	signature_t s2 = convertSignal(h2);

	double d = emd(&s1, &s2, dist, 0, 0);
	d /= (h1.size()-1); // normalize the distance into [0,1]; for normalized vectors, the maximal EMD distance is the maximal distance between bins

	delete s1.Features;
	delete s1.Weights;
	delete s2.Features;
	delete s2.Weights;
		
	return d;
}

double histogramDistance(QVector<double> h1, QVector<double> h2, DIST_TYPE type)
{
	double d = -1;
	switch (type)
	{
	case L1:
		d = computeL1Distance(h1,h2);
		break;
	case EMD:
		d = computeEMD(h1, h2);
		break;
	default:
		break;
	}

	return d;
}