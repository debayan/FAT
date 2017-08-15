#include "HierarchicalCluster.h"
#include <QString>
#include <QDebug>
#include "UtilityGlobal.h"

HierarchicalCluster::HierarchicalCluster()
{
}

HierarchicalCluster::~HierarchicalCluster()
{
}

inline QString convertMatrixToString(Eigen::MatrixXd M)
{
	QString string;
	string += "[";
	for (int i=0; i<M.rows(); i++)
	{
		string += "[";
		for (int j=0; j<M.cols(); j++)
		{
			if ( j < M.cols() - 1 )
			{
				string += QString::number( M(i, j) ) + ",";
			}
			else
			{
				string += QString::number( M(i, M.cols()-1) ) + "]";
			}
		}		

		if ( i < M.rows() - 1 )
		{
			string += ",";
		}
	}
	string += "]";

	return string;
}

alglib::ahcreport* HierarchicalCluster::clustering( Eigen::MatrixXd distance )
{
	alglib::clusterizerstate s;
	alglib::ahcreport* rep = new alglib::ahcreport();
	clusterizercreate(s);

	QString string = convertMatrixToString(distance);
	alglib::real_2d_array d(string.toStdString().c_str());
	clusterizersetdistances(s, d, true);
	clusterizerrunahc(s, *rep);

	//QStringList ahcInfor;
	//ahcInfor << QString::fromStdString(rep->p.tostring());
	//ahcInfor << QString::fromStdString(rep->pz.tostring());
	//ahcInfor << QString::fromStdString(rep->pm.tostring());
	//ahcInfor << QString::fromStdString(rep->mergedist.tostring(5));
	//debugBoxList(ahcInfor);

	return rep;
}
