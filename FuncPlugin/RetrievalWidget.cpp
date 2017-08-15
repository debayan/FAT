#include "RetrievalWidget.h"
#include "ui_RetrievalWidget.h"
#include <qwt_plot.h>
#include <qwt_plot_histogram.h>
#include <QColorDialog>
#include <qwt_plot_curve.h>

RetrievalWidget::RetrievalWidget(FuncPlugin * f, QWidget *parent) : QWidget(parent), ui(new Ui::RetrievalWidget)
{
    ui->setupUi(this);
	func = f;

	addPRWidget();

	/* retrieval */
	func->connect(ui->loadSceneListPushButton, SIGNAL(clicked()), SLOT(loadSceneListForRetrieval()));
	func->connect(ui->loadPartDataPushButton, SIGNAL(clicked()), SLOT(loadPartDataForRetrieval()));
	func->connect(ui->loadPoseDataPushButton, SIGNAL(clicked()), SLOT(loadShape2PoseDataForRetrieval()));
	this->connect(ui->tableWidget, SIGNAL(itemSelectionChanged()), SLOT(updateSelectedScene_Retrieval()));

	// different feature type
	func->connect(ui->DoRetrieval_ICON, SIGNAL(clicked()), SLOT(doRetrieval_icon()));
	func->connect(ui->DoRetrieval_LFD, SIGNAL(clicked()), SLOT(doRetrieval_lfd()));
	func->connect(ui->DoRetrieval_LFDICON, SIGNAL(clicked()), SLOT(doRetrieval_lfd_icon()));
	func->connect(ui->DoRetrieval_POSE, SIGNAL(clicked()), SLOT(doRetrieval_pose()));
	func->connect(ui->DoRetrieval_IBSH, SIGNAL(clicked()), SLOT(doRetrieval_ibsh()));
	func->connect(ui->DoRetrieval_ISET, SIGNAL(clicked()), SLOT(doRetrieval_iset()));
	func->connect(ui->DoRetrieval_GEO, SIGNAL(clicked()), SLOT(doRetrieval_geo()));
	func->connect(ui->DoRetrieval_GEOICON, SIGNAL(clicked()), SLOT(doRetrieval_geo_icon()));

	// for shape2pose retrieval mode
	connect(ui->modelSnapRadioButton, &QRadioButton::clicked, func, &FuncPlugin::loadShape2PoseDataForICON);
	connect(ui->poseSnapRadioButton, &QRadioButton::clicked, func, &FuncPlugin::loadShape2PoseDataForPOSE);

	ui->modelSnapRadioButton->setVisible(false);
	ui->poseSnapRadioButton->setVisible(false);

	// for pr-curve illustration
	this->connect(ui->iconRadioButton, SIGNAL(clicked()), SLOT(updateFeatureType()));
	this->connect(ui->isetRadioButton, SIGNAL(clicked()), SLOT(updateFeatureType()));
	this->connect(ui->lfdRadioButton, SIGNAL(clicked()), SLOT(updateFeatureType()));
	this->connect(ui->lfd_iconRadioButton, SIGNAL(clicked()), SLOT(updateFeatureType()));
	this->connect(ui->ibshRadioButton, SIGNAL(clicked()), SLOT(updateFeatureType()));
	this->connect(ui->poseRadioButton, SIGNAL(clicked()), SLOT(updateFeatureType()));
	this->connect(ui->geoRadioButton, SIGNAL(clicked()), SLOT(updateFeatureType()));
	this->connect(ui->geo_iconRadioButton, SIGNAL(clicked()), SLOT(updateFeatureType()));
	this->connect(ui->prShapeRadioButton, SIGNAL(clicked()), SLOT(updatePRWidget()));
	this->connect(ui->prClassRadioButton, SIGNAL(clicked()), SLOT(updatePRWidget()));

	/* results */
	func->connect(ui->loadRetrievalPushButton, SIGNAL(clicked()), SLOT(loadRetrievalResults()));
}

RetrievalWidget::~RetrievalWidget()
{
    delete ui;

	if (prPlot)
	{
		delete prPlot;
	}
}

void RetrievalWidget::updateSceneTable_Result(QVector<int> retrievalResults)
{
	int sceneNum = retrievalResults.size();
	int width = ui->retrievalWidget->width() / 2 - 1;

	ui->retrievalWidget->setColumnCount(2);
	ui->retrievalWidget->setRowCount(sceneNum);
	ui->retrievalWidget->setColumnWidth(0, width);
	ui->retrievalWidget->setColumnWidth(1, width);
	ui->retrievalWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);

	for (int i=0; i<sceneNum; i++)
	{
		QString label = "Scene #" + QString::number(retrievalResults[i]) + ":\n\n" + func->retrievalTool->sceneNames[retrievalResults[i]];
		QTableWidgetItem* itemText = new QTableWidgetItem(label);
		itemText->setTextAlignment(Qt::AlignCenter);
		ui->retrievalWidget->setItem(i, 0, itemText);

		QTableWidgetItem* itemImage = new QTableWidgetItem();
		itemImage->setData(Qt::DecorationRole, func->retrievalTool->sceneImages[retrievalResults[i]].scaled(width, width, Qt::KeepAspectRatio, Qt::SmoothTransformation));
		ui->retrievalWidget->setItem(i, 1, itemImage);

		ui->retrievalWidget->setRowHeight(i, width);
	}
}

void RetrievalWidget::updateSceneTable_Retrieval(int width)
{
	if (!func->retrievalTool)
	{
		return;
	}

	int sceneNum = func->retrievalTool->sceneNames.size();	
	func->retrievalTool->currIdx = 0;
	if (width == -1)
	{
		width = ui->tableWidget->width() / 2 - 10;
	}

	ui->tableWidget->setColumnCount(2);
	ui->tableWidget->setRowCount(sceneNum);
	ui->tableWidget->setColumnWidth(0, width);
	ui->tableWidget->setColumnWidth(1, width);
	ui->tableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);

	for (int i=0; i<sceneNum; i++)
	{
		QString label = "Scene #" + QString::number(i) + ":\n\n" + func->retrievalTool->sceneNames[i];
		QTableWidgetItem* itemText = new QTableWidgetItem(label);
		itemText->setTextAlignment(Qt::AlignCenter);
		ui->tableWidget->setItem(i, 0, itemText);

		QTableWidgetItem* itemImage = new QTableWidgetItem();
		itemImage->setData(Qt::DecorationRole, func->retrievalTool->sceneImages[i].scaled(width, width, Qt::KeepAspectRatio, Qt::SmoothTransformation));
		ui->tableWidget->setItem(i, 1, itemImage);

		ui->tableWidget->setRowHeight(i, width);
	}

	ui->tableWidget->setCurrentCell(0, 0);
}

void RetrievalWidget::updateSelectedScene_Retrieval()
{
	int idx = ui->tableWidget->currentRow();

	func->setCurrScene_Ratrieval(idx);
	func->updateRetrieval();

	updatePRWidget();
}

void RetrievalWidget::addPRWidget()
{
	prPlot = new QwtPlot(ui->retrievalTab);
	prPlot->setCanvasBackground( Qt::white );
	prPlot->enableAxis(QwtPlot::xBottom, false);
	prPlot->enableAxis(QwtPlot::yLeft, false);
	prPlot->setFixedHeight(200);
	ui->prLayout->addWidget(prPlot);	
}

void RetrievalWidget::updatePRWidget()
{
	if (!prPlot)
	{
		return;
	}

	prPlot->detachItems();
	prPlot->enableAxis(QwtPlot::xBottom, false);
	prPlot->enableAxis(QwtPlot::yLeft, false);

	if ( !func->retrievalTool || !func->retrievalTool->prReady())
	{
		return;
	}

	Eigen::MatrixXd pr;
	if (ui->prShapeRadioButton->isChecked())
	{
		pr = func->retrievalTool->getPRShape();
	}
	else
	{
		pr = func->retrievalTool->getPRClass();
	}
	int curvNum = pr.rows()/2;
	for (int i=curvNum-1; i>=0; i--)
	{
		QVector<QPointF> samples;	

		for (int j=0; j<pr.cols(); j++)
		{
			double x = pr(2*i, j);	
			double y = pr(2*i+1, j);
			
			samples << QPointF(x,y);	
		}

		QwtPlotCurve * curve = new QwtPlotCurve();
		curve->setSamples(samples);

		if (i==0 && ui->prClassRadioButton->isChecked())
		{
			curve->setPen( QPen(Qt::red, 2) );
		}
		else
		{
			curve->setPen( QPen(Qt::blue, 2) );
		}

		curve->setStyle(QwtPlotCurve::Lines);
		curve->attach(prPlot);
	}

	prPlot->setAxisScale(QwtPlot::xBottom, 0, 1);
	prPlot->setAxisScale(QwtPlot::yLeft, 0, 1); 
	prPlot->enableAxis(QwtPlot::xBottom, true);
	prPlot->enableAxis(QwtPlot::yLeft, true);
	prPlot->replot();
}

void RetrievalWidget::setCurrTabToScenes()
{
	ui->tabWidget->setCurrentIndex(0);
}

void RetrievalWidget::updateFeatureType( FEATURE_TYPE type )
{
	switch (type)
	{
	case ICON:
		ui->iconRadioButton->setChecked(true);
		break;
	case LFDICON:
		ui->lfd_iconRadioButton->setChecked(true);
		break;
	case LFD:
		ui->lfdRadioButton->setChecked(true);
		break;
	case IBSH:
		ui->ibshRadioButton->setChecked(true);
		break;
	case POSE:
		ui->poseRadioButton->setChecked(true);
		break;
	case ISET:
		ui->isetRadioButton->setChecked(true);
		break;
	case GEO:
		ui->geoRadioButton->setChecked(true);
		break;
	case GEOICON:
		ui->geo_iconRadioButton->setChecked(true);
		break;
	default:
		break;
	}
}

void RetrievalWidget::updateFeatureType()
{
	FEATURE_TYPE type;
	if (ui->iconRadioButton->isChecked())
	{
		type = ICON;
	}
	else if (ui->lfdRadioButton->isChecked())
	{
		type = LFD;
	}
	else if (ui->ibshRadioButton->isChecked())
	{
		type = IBSH;
	}
	else if (ui->lfd_iconRadioButton->isChecked())
	{
		type = LFDICON;
	}
	else if (ui->poseRadioButton->isChecked())
	{
		type = POSE;
	}
	else if (ui->isetRadioButton->isChecked())
	{
		type = ISET;
	}
	else if (ui->geoRadioButton->isChecked())
	{
		type = GEO;
	}
	else if (ui->geo_iconRadioButton->isChecked())
	{
		type = GEOICON;
	}

	if (func->retrievalTool && func->retrievalTool->currType != type)
	{
		func->retrievalTool->currType = type;
		updatePRWidget();
	}
}

double RetrievalWidget::getCombWeight()
{
	return ui->combWeightSpinBox->value();
}

void RetrievalWidget::enableICON(bool mode)
{
	if (mode)
	{
		ui->DoRetrieval_ICON->setEnabled(true);
		ui->iconRadioButton->setEnabled(true);
	}
	else
	{
		ui->DoRetrieval_ICON->setEnabled(false);
		ui->iconRadioButton->setEnabled(false);
	}
}

void RetrievalWidget::enableISET(bool mode)
{
	if (mode)
	{
		ui->DoRetrieval_ISET->setEnabled(true);
		ui->isetRadioButton->setEnabled(true);
	}
	else
	{
		ui->DoRetrieval_ISET->setEnabled(false);
		ui->isetRadioButton->setEnabled(false);
	}
}

void RetrievalWidget::enableIBSH(bool mode)
{
	if (mode)
	{
		ui->DoRetrieval_IBSH->setEnabled(true);
		ui->ibshRadioButton->setEnabled(true);
		ui->depthSpinBox->setEnabled(true);
	}
	else
	{
		ui->DoRetrieval_IBSH->setEnabled(false);
		ui->ibshRadioButton->setEnabled(false);
		ui->depthSpinBox->setEnabled(false);
	}
}

void RetrievalWidget::enablePOSE(bool mode)
{
	if (mode)
	{
		ui->DoRetrieval_POSE->setEnabled(true);
		ui->poseRadioButton->setEnabled(true);
	}
	else
	{
		ui->DoRetrieval_POSE->setEnabled(false);
		ui->poseRadioButton->setEnabled(false);
	}
}

void RetrievalWidget::enableLFD(bool mode)
{
	if (mode)
	{
		ui->DoRetrieval_LFD->setEnabled(true);
		ui->lfdRadioButton->setEnabled(true);
	}
	else
	{
		ui->DoRetrieval_LFD->setEnabled(false);
		ui->lfdRadioButton->setEnabled(false);
	}
}

void RetrievalWidget::enableGEO(bool mode)
{
	if (mode)
	{
		ui->DoRetrieval_GEO->setEnabled(true);
		ui->geoRadioButton->setEnabled(true);
	}
	else
	{
		ui->DoRetrieval_GEO->setEnabled(false);
		ui->geoRadioButton->setEnabled(false);
	}
}

void RetrievalWidget::enableLFDICON(bool mode)
{
	if (mode)
	{
		ui->DoRetrieval_LFDICON->setEnabled(true);
		ui->lfd_iconRadioButton->setEnabled(true);
		ui->combWeightSpinBox->setEnabled(true);
	}
	else
	{
		ui->DoRetrieval_LFDICON->setEnabled(false);
		ui->lfd_iconRadioButton->setEnabled(false);
		ui->combWeightSpinBox->setEnabled(false);
	}
}

void RetrievalWidget::enableGEOICON(bool mode)
{
	if (mode)
	{
		ui->DoRetrieval_GEOICON->setEnabled(true);
		ui->geo_iconRadioButton->setEnabled(true);
		ui->partCombWeightSpinBox->setEnabled(true);
	}
	else
	{
		ui->DoRetrieval_GEOICON->setEnabled(false);
		ui->geo_iconRadioButton->setEnabled(false);
		ui->partCombWeightSpinBox->setEnabled(false);
	}
}

void RetrievalWidget::enableOurDataRetriComp(bool mode)
{
	enableICON(mode);
	enableISET(mode);
	enableIBSH(mode);
	enableLFD(mode);
	enableLFDICON(mode);
	enableGEOICON(false);
	enablePR();
}

void RetrievalWidget::enableShape2PoseDataRetriComp(bool mode)
{
	enableICON(mode);
	enablePOSE(mode);
	enableGEOICON(false);
	enablePoseRetriMode(mode);
}

void RetrievalWidget::enablePartDataRetriComp(bool mode)
{
	enableICON(true);
	enableGEO(mode);
	//enableGEOICON(mode);
}

void RetrievalWidget::enablePoseRetriMode(bool mode)
{
	ui->modelSnapRadioButton->setEnabled(mode);
	ui->poseSnapRadioButton->setEnabled(mode);
	ui->modelSnapRadioButton->setChecked(mode);
}

void RetrievalWidget::enablePoseRetriModeVisiable(bool mode)
{
	ui->modelSnapRadioButton->setVisible(mode);
	ui->poseSnapRadioButton->setVisible(mode);
}

void RetrievalWidget::enablePR()
{
	ui->prShapeRadioButton->setEnabled(true);
	ui->prClassRadioButton->setEnabled(true);
}

void RetrievalWidget::updatePoseRetrievalMode(int mode)
{
	switch (mode)
	{
	case 0:		// normal
		{
			ui->modelSnapRadioButton->setChecked(true);
			break;
		}
	case 1:		// shape2pose
		{
			ui->poseSnapRadioButton->setChecked(true);
			break;
		}
	default:
		break;
	}
}

int RetrievalWidget::getDepth()
{
	return ui->depthSpinBox->value();
}

int RetrievalWidget::getReturnRetrievalNum()
{
	return ui->returnRetrievalSpinBox->value();
}

double RetrievalWidget::getPartCombWeight()
{
	return ui->partCombWeightSpinBox->value();
}