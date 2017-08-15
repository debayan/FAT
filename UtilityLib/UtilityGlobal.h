#pragma once

#include <QtXml/QDomDocument>

#include <QSharedPointer>
#include <QString>
#include <QVector>
#include <QQueue>
#include <QMap>
#include <QSet>
#include <QFile>

#include "SurfaceMeshModel.h"
using namespace SurfaceMesh;

namespace SurfaceMesh{
typedef Eigen::Vector2d Vector2;
typedef Eigen::Vector4d Vector4;
}

typedef QSharedPointer<SurfaceMeshModel> MeshPtr;

QString qStr(const Vector2 &v, char sep = ' ');
QString qStr(const Vector3 &v, char sep = ' ');
QString qStr(const Vector4 &v, char sep = ' ');

Vector3 toVector3(QString string);

typedef QVector< QVector<QString> > StrArray2D;

typedef QMap< QString, QVariant > PropertyMap;

#ifdef Q_OS_WIN
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif
QString getcwd();

// Simplify runtime debugging
#include <QMessageBox>
template<typename DataType>
static inline void debugBox( DataType message ){
	QMessageBox msgBox;
	msgBox.setTextFormat(Qt::RichText);
	msgBox.setText( QString("%1").arg(message) );
	msgBox.setStandardButtons(QMessageBox::Ok);
	msgBox.exec();
}
static inline void debugBoxList( QStringList messages ){
	debugBox( messages.join("<br>") );
}
template<typename Container>
static inline void debugBoxVec( Container data ){
	QStringList l;
	for(auto d : data) l << QString("%1").arg(d);
	debugBoxList(l);
}
template<typename Container2D>
static inline void debugBoxVec2( Container2D data, int limit = -1 ){
	QStringList l;
	for(auto row : data){
		QStringList line;
		for(auto d : row) line << QString("%1").arg( d );
		l << QString("%1").arg( line.join(", ") );
		if(limit > 0 && l.size() == limit) break;
	}
	if(limit > 0 && data.size() - l.size() > 0) l << QString("... (%1) more").arg(data.size() - l.size());
	debugBoxList(l);
}

static inline void matrixToFile(const Eigen::MatrixXd & M, QString filename)
{
	QFile file( filename );
	if(!file.open(QFile::WriteOnly | QFile::Text)) return;
	QTextStream out(&file);
	for(unsigned int i = 0; i < M.rows(); i++)
	{
		QStringList row;
		for(unsigned int j = 0; j < M.cols(); j++)
			row << QString::number(M(i,j));
		out << (row.join(",") + "\n");
	}
}

#include "BasicTable.h"
//#include "Colormap.h"