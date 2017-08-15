#pragma once

#include <QWidget>
#include "FuncPlugin.h"
#include <qwt_plot.h>
#include <qwt_scale_draw.h>
#include <QListWidget>

namespace Ui {
class IconWidget;
}

class IconWidget : public QWidget
{
    Q_OBJECT

public:
    explicit IconWidget(FuncPlugin *f, QWidget *parent = 0);
    ~IconWidget();

public:
	void updateObjectList();
	void updateIbsList();
	void updateSceneListTable(int width = -1);
	
	void updateHierarchyList();
	void updateIbsListCheckState(QVector<bool> sel);

	void setSelectedScene(int idx);
	void deleteSelectedScene();
	void updateSelectedScene();
	
	void buildICON();
	bool visualizeHierarchy();

	QVector<bool> getObjectCheckState();
	QVector<bool> getIbsCheckState();

	void addInteractionFeatureWidget();
	void drawHistogram(QwtPlot* plot, QVector<double> distribution);

	int getSeletedInteractionIdx();
	int getSelectedHierarchyIdx();

	void resetIbsFeatures();
	void setInteractionFeatureYScale(QVector<double> yMax);
	void setClusterFeatureYScale(QVector<double> yMax);

	Vec3d getUprightDirection();
	void setUprightDirection(Vec3d upright);
	void setObjectDrawMode(OBJECT_DRAW_MODE m);

	void setCurrTabToDescriptor();
	void assignParaToScene(Scene* scene);
	void getParaFromScene(Scene* scene);

	void assignParameters(DistParameter &distPara, HierParameter &hierPara);

	void setCurrentTable(int idx);
	void enableComputeFeatures(int type);		// type = 0 for ICON, = 1 for ISET

public slots:
	void checkObjectDrawMode();
	void checkIbsDrawMode();
	void checkIbsSampleShow();
	void checkIbsWeightShow();
	void checkInteractionShow();
	void checkInteractionIBSShow();
	void checkInteractionIRShow();
	void checkInteractionObjShow();
	void checkIBSWeightAndSamplesDrawMode();

	void updateInteractionFeatures();

	void normalizeWeight();
	void normalizeIbsWeightBasedOnPfh();
	void normalizeIbsWeightBasedOnDir();
	void normalizeRegionWeightBasedOnPfh(); 
	void normalizeRegionWeightBasedOnDir(); 

	void updateInteractionColor(QTreeWidgetItem * item);
	void changeBaseColor();

signals:
	void objectDrawModeChanged(OBJECT_DRAW_MODE m);	
	void ibsDrawModeChanged(IBS_DRAW_MODE m);
	void ibsSampleShow(bool show);
	void ibsWeightShow(bool show);
	void interactionShow(bool show);
	void interactionIBSShow(bool show);
	void interactionIRShow(bool show);
	void interactionObjShow(bool show);
	void uprightChanged(Vec3d up);

private:
    Ui::IconWidget *ui;
	FuncPlugin * func;

	QVector<QwtPlot*>  interactionFeaturePlot;  // 5 histograms: pfh, dir, dist, region:pfh, region, dir
	QVector<QwtPlot*>  clusterFeaturePlot;

public:
	bool isDeleting;
};

class MyScaleDraw: public QwtScaleDraw
{
public:
	MyScaleDraw(Scene *scene):QwtScaleDraw(){s = scene; prepareMapping();};

protected:
	virtual QwtText label(double value) const
	{
		int idx = value;

		QString str;
		if (abs(value-idx) < 1e-15 && idx>=0 && idx<objIdx.size())
		{
			str = QString::number(objIdx[idx]);
		}

		return str;
	}

private:
	void prepareMapping() 
	{
		int n = s->ahcreport->npoints;

		objIdx = QVector<int>(n, -1);
		for (int i=0; i<n; i++)
		{
			objIdx[s->ahcreport->p[i]] = s->interObjIdx[i];
		}
	};

public:
	Scene* s;
	QVector<int> objIdx;
};
