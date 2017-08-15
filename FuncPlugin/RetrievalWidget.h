#pragma once

#include <QWidget>
#include "FuncPlugin.h"
#include <qwt_plot.h>
#include <qwt_scale_draw.h>

namespace Ui {
class RetrievalWidget;
}

class RetrievalWidget : public QWidget
{
    Q_OBJECT

public:
    explicit RetrievalWidget(FuncPlugin * f, QWidget *parent = 0);
    ~RetrievalWidget();

public:
	void updateSceneTable_Retrieval(int width = -1);
	void updateSceneTable_Result(QVector<int> RRR);

	void addPRWidget();
	void setCurrTabToScenes();

	void updateFeatureType(FEATURE_TYPE type);
	double getCombWeight();

	void enableICON(bool mode = true);
	void enableISET(bool mode = true);
	void enableIBSH(bool mode = true);
	void enablePOSE(bool mode = true);
	void enableLFD(bool mode = true);
	void enableGEO(bool mode = true);
	void enableLFDICON(bool mode = true);
	void enableGEOICON(bool mode = true);
	void enablePR();

	void enableOurDataRetriComp(bool);				// true is enabled, false is disabled
	void enableShape2PoseDataRetriComp(bool);
	void enablePartDataRetriComp(bool);

	void enablePoseRetriMode(bool);
	void enablePoseRetriModeVisiable(bool);

	void updatePoseRetrievalMode(int mode);			// 0 normal, 1 shape2pose

	int getDepth();
	int getReturnRetrievalNum();
	double getPartCombWeight();

public slots:
	void updatePRWidget();

private slots:
	void updateFeatureType();
	void updateSelectedScene_Retrieval();

private:
    Ui::RetrievalWidget *ui;
	FuncPlugin * func;

public:
	QwtPlot* prPlot;
};