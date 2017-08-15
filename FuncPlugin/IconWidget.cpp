#include "IconWidget.h"
#include "ui_IconWidget.h"
#include <qwt_plot.h>
#include <qwt_plot_histogram.h>
#include <QColorDialog>
#include <qwt_plot_curve.h>

IconWidget::IconWidget(FuncPlugin * f, QWidget *parent) :
	QWidget(parent), ui(new Ui::IconWidget)
{
    ui->setupUi(this);
	this->func = f;
	isDeleting = false;
	
	addInteractionFeatureWidget();

	/* scene panel */
	// scene list
	func->connect(ui->loadScenes, SIGNAL(clicked()), SLOT(loadScenes()));
	func->connect(ui->addScenes, SIGNAL(clicked()), SLOT(addScenes()));
	func->connect(ui->updateSnapshot, SIGNAL(clicked()), SLOT(updateSnapshot()));
	this->connect(ui->deleteScene, &QPushButton::clicked, this, &IconWidget::deleteSelectedScene);		// problem???
	this->connect(ui->sceneTableWidget, &QTableWidget::itemSelectionChanged, this, &IconWidget::updateSelectedScene);
	// construction
	this->connect(ui->buildIconButton, &QPushButton::clicked, this, &IconWidget::buildICON);		// problem???

	/* settings panel */
	func->connect(ui->objectListWidget, SIGNAL(itemChanged(QListWidgetItem *)), SLOT(updateSelectedObjects()));
	// upright direction
	func->connect(ui->uprightCurrPushButton, SIGNAL(clicked()), SLOT(setUprightDirection()));
	func->connect(ui->uprightAllPushButton, SIGNAL(clicked()), SLOT(setUprightDirectionForAll()));
	// parameters
	this->connect(ui->ibsWeightSpinBox, SIGNAL(valueChanged(double)), SLOT(normalizeWeight()));
	this->connect(ui->ibsPfhSpinBox, SIGNAL(valueChanged(double)), SLOT(normalizeIbsWeightBasedOnPfh()));
	this->connect(ui->ibsDirSpinBox, SIGNAL(valueChanged(double)), SLOT(normalizeIbsWeightBasedOnDir()));
	this->connect(ui->regionPfhSpinBox, SIGNAL(valueChanged(double)), SLOT(normalizeRegionWeightBasedOnPfh()));
	this->connect(ui->regionDirSpinBox, SIGNAL(valueChanged(double)), SLOT(normalizeRegionWeightBasedOnDir()));
	// compute feature for ICON and ISET
	func->connect(ui->computeFeatureButton, SIGNAL(clicked()), SLOT(computeFeatures()));

	/* render panel */	
	// object mode
	this->connect(ui->meshRadioButton, SIGNAL(clicked()), SLOT(checkObjectDrawMode()));
	this->connect(ui->colorRadioButton, SIGNAL(clicked()), SLOT(checkObjectDrawMode()));
	this->connect(ui->sampleRadioButton, SIGNAL(clicked()), SLOT(checkObjectDrawMode()));
	this->connect(ui->noneRadioButton, SIGNAL(clicked()), SLOT(checkObjectDrawMode()));
	this->connect(ui->wireOrigRadioButton, SIGNAL(clicked()), SLOT(checkObjectDrawMode()));
	this->connect(ui->boxObjRadioButton, SIGNAL(clicked()), SLOT(checkObjectDrawMode()));
	func->connect(this, SIGNAL(objectDrawModeChanged(OBJECT_DRAW_MODE)), SLOT(setObjectDrawMode(OBJECT_DRAW_MODE)));
	// candidate hierarchies
	this->connect(ui->hierarchyListWidget, &QListWidget::itemSelectionChanged, this, &IconWidget::updateIbsList);
	// interactions
	func->connect(ui->ibsTreeWidget, SIGNAL(itemChanged(QTreeWidgetItem*, int)), SLOT(updateSelectedIBS()));
	this->connect(ui->ibsTreeWidget, SIGNAL(itemSelectionChanged()), SLOT(updateInteractionFeatures()));
	this->connect(ui->ibsTreeWidget, SIGNAL(itemDoubleClicked(QTreeWidgetItem*, int)), SLOT(updateInteractionColor(QTreeWidgetItem *)));
	// modes
	this->connect(ui->partialIbsRadioButton, SIGNAL(clicked()), SLOT(checkIbsDrawMode()));

	func->connect(this, SIGNAL(ibsDrawModeChanged(IBS_DRAW_MODE)), SLOT(setIbsDrawMode(IBS_DRAW_MODE)));
	this->connect(ui->ibsSampleCheckBox, SIGNAL(clicked()), SLOT(checkIbsSampleShow()));
	func->connect(this, SIGNAL(ibsSampleShow(bool)), SLOT(setIbsSampleShow(bool)));
	this->connect(ui->ibsWeightCheckBox, SIGNAL(clicked()), SLOT(checkIbsWeightShow()));
	func->connect(this,SIGNAL(ibsWeightShow(bool)), SLOT(setIbsWeightShow(bool)));
	this->connect(ui->drawInteractionRadioBtn, SIGNAL(clicked()), SLOT(checkInteractionShow()));
	func->connect(this, SIGNAL(interactionShow(bool)), SLOT(setInteractionShow(bool)));
	this->connect(ui->interIBSCheckBox, SIGNAL(clicked()), SLOT(checkInteractionIBSShow()));
	func->connect(this, SIGNAL(interactionIBSShow(bool)), SLOT(setInteractionIBSShow(bool)));
	this->connect(ui->interIRCheckBox, SIGNAL(clicked()), SLOT(checkInteractionIRShow()));
	func->connect(this, SIGNAL(interactionIRShow(bool)), SLOT(setInteractionIRShow(bool)));
	this->connect(ui->interObjCheckBox, SIGNAL(clicked()), SLOT(checkInteractionObjShow()));
	func->connect(this, SIGNAL(interactionObjShow(bool)), SLOT(setInteractionObjShow(bool)));
	this->connect(ui->baseColorPushButton, SIGNAL(clicked()), SLOT(changeBaseColor()));
	func->connect(ui->outputIbsButton, SIGNAL(clicked()), SLOT(outputIBS()));
	func->connect(ui->outputIRButton, SIGNAL(clicked()), SLOT(outputIR()));
}

IconWidget::~IconWidget()
{
    delete ui;

	for (auto plot:interactionFeaturePlot)
	{
		if (plot)
		{
			delete plot;
		}
	}

	for (auto plot:clusterFeaturePlot)
	{
		if (plot)
		{
			delete plot;
		}
	}
}

void IconWidget::updateHierarchyList()
{
	ui->hierarchyListWidget->clear();

	if (func->currScene) {
		for (int i = 0; i < func->currScene->getCandiHierNum(); ++i) {
			QListWidgetItem* item = new QListWidgetItem("Hierarchy #" + QString::number(i));
			item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
			ui->hierarchyListWidget->addItem(item);
		}
		ui->hierarchyListWidget->setCurrentRow(0);
	}
}

void IconWidget::updateObjectList()
{
	ui->objectListWidget->clear();

	if (func->currScene)
	{
		for (int i=0; i<func->currScene->objects.size(); i++)
		{
			QListWidgetItem* item = new QListWidgetItem("Object #" + QString::number(i));
			item->setFlags(item->flags() | Qt::ItemIsUserCheckable);	// set checkable flag

			if (func->currScene->objects[i]->isCentral)
			{
				item->setCheckState(Qt::Checked); 
			}
			else
			{
				item->setCheckState(Qt::Unchecked); 
			}		

			ui->objectListWidget->addItem(item);
		}	
	}

	updateHierarchyList();
	updateIbsList();
}

void IconWidget::updateIbsList()
{
	ui->ibsTreeWidget->clear();

	if (func->currScene)
	{
		int hierarchyIdx = ui->hierarchyListWidget->currentRow();
		QVector<int> interactionsIdx = func->currScene->getInteractions(hierarchyIdx);
		QVector<int> interHierParentIdx = func->currScene->getInterHierParentNodeIdx(hierarchyIdx);

		func->currScene->interactionsIdx = func->currScene->getInteractionsPreOrder(hierarchyIdx);
		func->currScene->intersParentIdx = func->currScene->getInterHierParentNodeIdxPreOrder(hierarchyIdx);
		//func->currScene->intersParentIdx = interHierParentIdx;

		func->currScene->lastSelectedInters = QVector<bool>(interactionsIdx.size(), false);
		assert(interactionsIdx.size() == interHierParentIdx.size());

		if (interactionsIdx.size() > 0)
		{
			QVector<QTreeWidgetItem*> treeNodeList;
			for (int i = 0; i < interactionsIdx.size(); ++i)
			{
				int interIdx = interactionsIdx[i];
				QVector<int> interObjIndex = func->currScene->allInteractions[interIdx]->obj->origIdx;
			
				// push back tree item
				QString itemLabel;
				for (int j = 0; j < interObjIndex.size(); ++j) {
					if (j == interObjIndex.size() - 1)
						itemLabel = itemLabel + QString::number(interObjIndex[j]);
					else
						itemLabel = itemLabel + QString::number(interObjIndex[j]) + ",";
				}

				QTreeWidgetItem *item = new QTreeWidgetItem();
				item->setText(0, itemLabel);
				item->setCheckState(0, Qt::Unchecked);
				item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
				treeNodeList.push_back(item);
			}

			int rootNodeIdx = -1;
			int currParentIdx = -1;
			int nLevel = 0;
			for (int i = 0; i< interHierParentIdx.size(); ++i)
			{
				int parentIdx = interHierParentIdx[i];
				if (parentIdx != currParentIdx)
				{
					++nLevel;
					currParentIdx = parentIdx;
				}
				//else
				{
					QString itemLabel = treeNodeList[i]->text(0);
					itemLabel.push_front("Level #" + QString::number(nLevel) + ": ");
					treeNodeList[i]->setText(0, itemLabel);
				}
				
				if (parentIdx == -1)	// not root node
				{
					rootNodeIdx = i;
				}
				else
				{
					int childIdx = i;
					treeNodeList[parentIdx]->addChild(treeNodeList[childIdx]);
				}
			}

			ui->ibsTreeWidget->addTopLevelItem(treeNodeList[rootNodeIdx]);

			// expand IBS tree widget
			QTreeWidgetItemIterator it(ui->ibsTreeWidget);
			while (*it) {
				ui->ibsTreeWidget->expandItem(*it);
				++it;
			}
		}
	}
	
	resetIbsFeatures();
}

void IconWidget::updateIbsListCheckState(QVector<bool> sel)
{
	QTreeWidgetItemIterator it(ui->ibsTreeWidget);
	int idx = 0;
	isDeleting = true;
	while ((*it)) {
		if (sel[idx])
			(*it)->setCheckState(0, Qt::Checked);
		else
			(*it)->setCheckState(0, Qt::Unchecked);
		++idx;
		++it;
	}
	isDeleting = false;
}

void IconWidget::updateSceneListTable(int width /* = -1 */)
{
	isDeleting = true;
	while (ui->sceneTableWidget->rowCount() > 0)
	{
		ui->sceneTableWidget->removeRow(0);
	}
	isDeleting = false;

	if (!func->sceneSet)
	{
		return;
	}

	int sceneNum = func->sceneSet->scenes.size();		 
	int currRow = func->currSceneIdx;

	if (width == -1)
	{
		width = ui->sceneTableWidget->width() / 2 - 10;
	}

	ui->sceneTableWidget->setColumnCount(2);
	ui->sceneTableWidget->setRowCount(sceneNum);
	ui->sceneTableWidget->setColumnWidth(0, width);
	ui->sceneTableWidget->setColumnWidth(1, width);
	ui->sceneTableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);

	for (int i=0; i<sceneNum; i++)
	{
		Scene * scene = func->sceneSet->scenes[i];
		QString label = "Scene #" + QString::number(i) + ":\n\n" + scene->name;
		QTableWidgetItem* itemText = new QTableWidgetItem(label);
		itemText->setTextAlignment(Qt::AlignCenter);
		ui->sceneTableWidget->setItem(i, 0, itemText);

		QTableWidgetItem* itemImage = new QTableWidgetItem();
		itemImage->setData(Qt::DecorationRole, scene->snapshot.scaled(width, width, Qt::KeepAspectRatio, Qt::SmoothTransformation));
		ui->sceneTableWidget->setItem(i, 1, itemImage);

		ui->sceneTableWidget->setRowHeight(i, width);
	}

	if (currRow != -1)
	{
		ui->sceneTableWidget->setCurrentCell(currRow, 0);
	}
}

void IconWidget::setSelectedScene(int idx)
{
	updateSceneListTable();
	if (idx != -1)
	{
		ui->sceneTableWidget->setCurrentCell(idx, 0);
	}
}

void IconWidget::deleteSelectedScene()
{
	int idx = ui->sceneTableWidget->currentRow();
	func->deleteScene(idx);
}

void IconWidget::updateSelectedScene()
{
	if (!isDeleting) 
	{
		int idx = ui->sceneTableWidget->currentRow();
		func->setCurrScene(idx);		
	}
}

void IconWidget::buildICON()
{
	if (ui->buildForSelectedRadioBtn->isChecked())
	{
		func->buildInteractionHierarchy();
	}

	if (ui->buildForAllRadioBtn->isChecked())
	{
		func->buildInteractionHierarchyForAllScenes();
	}

	if (ui->buildForSceneListRadioBtn->isChecked())
	{
		func->buildInteractionHierarchyForSceneList();
	}
}

bool IconWidget::visualizeHierarchy()
{
	return ui->visualizationCheckBox->isChecked();
}

QVector<bool> IconWidget::getObjectCheckState()
{
	QVector<bool> state;
	for (int i=0; i<ui->objectListWidget->count(); i++)	{
		state.push_back(ui->objectListWidget->item(i)->checkState());		
	}

	return state;
}

QVector<bool> IconWidget::getIbsCheckState()
{
	QVector<bool> state;
	
	QTreeWidgetItemIterator it(ui->ibsTreeWidget);
	while (*it) 
	{
		state.push_back((*it)->checkState(0));
		++it;
	}
	
	return state;
}

void IconWidget::checkObjectDrawMode()
{
	OBJECT_DRAW_MODE m;
	if (ui->meshRadioButton->isChecked())
	{
		m = DRAW_MESH;
	}
	else if (ui->sampleRadioButton->isChecked())
	{
		m = DRAW_SAMPLE;
	}
	else if (ui->boxObjRadioButton->isChecked())
	{
		m = DRAW_BBOX_OBJ;
	}
	else if (ui->noneRadioButton->isChecked())
	{
		m = DRAW_NONE;
	}
	else if (ui->wireOrigRadioButton->isChecked())
	{
		m = DRAW_WIRE_ORIG;
	}
	else if (ui->colorRadioButton->isChecked())
	{
		m = DRAW_MESH_COLOR;
	}
	emit(objectDrawModeChanged(m));
}

void IconWidget::checkIbsDrawMode()
{
	IBS_DRAW_MODE m;
	m = DRAW_PARTIAL;

	checkInteractionShow();
	emit(ibsDrawModeChanged(m));
}

void IconWidget::checkIbsSampleShow()
{
	bool show = ui->ibsSampleCheckBox->checkState();

	emit(ibsSampleShow(show));
}

void IconWidget::checkIbsWeightShow()
{
	bool show = ui->ibsWeightCheckBox->checkState();

	emit(ibsWeightShow(show));
}

void IconWidget::checkInteractionShow()
{
	bool show = ui->drawInteractionRadioBtn->isChecked();
	if (show) {
		ui->interIBSCheckBox->setEnabled(true);
		ui->interIRCheckBox->setEnabled(true);
		ui->interObjCheckBox->setEnabled(true);
		ui->interIBSCheckBox->setChecked(true);
		ui->interIRCheckBox->setChecked(true);
		ui->interObjCheckBox->setChecked(true);
		ui->baseColorPushButton->setEnabled(true);
	} else {
		ui->interIBSCheckBox->setChecked(false);
		ui->interIRCheckBox->setChecked(false);
		ui->interObjCheckBox->setChecked(false);
		ui->interIBSCheckBox->setDisabled(true);
		ui->interIRCheckBox->setDisabled(true);
		ui->interObjCheckBox->setDisabled(true);
		ui->baseColorPushButton->setEnabled(false);
	}

	emit(interactionShow(show));

	checkInteractionIBSShow();
	checkInteractionIRShow();
	checkInteractionObjShow();
	checkIBSWeightAndSamplesDrawMode();
}

void IconWidget::checkInteractionIBSShow()
{
	bool show = ui->interIBSCheckBox->checkState();
	checkIBSWeightAndSamplesDrawMode();

	emit(interactionIBSShow(show));
}

void IconWidget::checkInteractionIRShow()
{
	bool show = ui->interIRCheckBox->checkState();

	emit(interactionIRShow(show));
}

void IconWidget::checkInteractionObjShow()
{
	bool show = ui->interObjCheckBox->checkState();

	emit(interactionObjShow(show));
}

void IconWidget::checkIBSWeightAndSamplesDrawMode()
{
	bool show = ui->interIBSCheckBox->checkState();

	if (show) {
		ui->ibsWeightCheckBox->setEnabled(true);
		ui->ibsSampleCheckBox->setEnabled(true);
		ui->ibsWeightCheckBox->setChecked(true);
		ui->ibsSampleCheckBox->setChecked(false);
	} else {
		ui->ibsWeightCheckBox->setChecked(false);
		ui->ibsSampleCheckBox->setChecked(false);
		ui->ibsWeightCheckBox->setEnabled(false);
		ui->ibsSampleCheckBox->setEnabled(false);
	}
}

void IconWidget::addInteractionFeatureWidget()
{
	// create 5 new qwt plot widgets for 3 IBS features and 2 region feature, one empty
	for (int i=0; i<6; i++)
	{
		QwtPlot* plot = new QwtPlot(ui->ibsGroup);
		plot->setCanvasBackground( Qt::white );
		plot->setAxisScale( QwtPlot::yLeft, 0.0, 1.0 );
		plot->enableAxis(QwtPlot::xBottom, false);
		plot->enableAxis(QwtPlot::yLeft, false);
		plot->setFixedHeight(50);

		interactionFeaturePlot.push_back(plot);		
	}	
	interactionFeaturePlot[0]->setFooter("IBS:PFH");
	interactionFeaturePlot[1]->setFooter("IBS:Dir");
	interactionFeaturePlot[2]->setFooter("IBS:Dist");
	interactionFeaturePlot[3]->setFooter("IR:PFH");
	interactionFeaturePlot[4]->setFooter("IR:Dir");
	interactionFeaturePlot[5]->setFooter("IR:Height");

	QHBoxLayout *ibsFeatureLayout = new QHBoxLayout(ui->ibsGroup);
	ibsFeatureLayout->addWidget(interactionFeaturePlot[0]);
	ibsFeatureLayout->addWidget(interactionFeaturePlot[1]);
	ibsFeatureLayout->addWidget(interactionFeaturePlot[2]);
	
	QHBoxLayout *regionFeatureLayout = new QHBoxLayout(ui->ibsGroup);
	regionFeatureLayout->addWidget(interactionFeaturePlot[3]);
	regionFeatureLayout->addWidget(interactionFeaturePlot[4]);
	regionFeatureLayout->addWidget(interactionFeaturePlot[5]);

	ui->ibsFeatureLayout->addLayout(ibsFeatureLayout);
	ui->irFeatureLayout->addLayout(regionFeatureLayout);
}

void IconWidget::updateInteractionFeatures()
{
	// Geometrical features
	for (int i=0; i<6; i++)
	{
		QVector<double> distribution = func->getInteractionFeature(i);
		drawHistogram(interactionFeaturePlot[i], distribution);
	}

	// Topological features
	QVector<int> bettiNumbers = func->getBettiNumber();
}

void IconWidget::drawHistogram( QwtPlot* plot, QVector<double> distribution )
{
	if (distribution.isEmpty())
	{
		plot->detachItems();
		return;
	}

	QwtPlotHistogram *histogram = new QwtPlotHistogram();	
	histogram->setStyle(QwtPlotHistogram::Columns);  

	int numValues = distribution.size();
	plot->setAxisScale( QwtPlot::xBottom, 0.0, numValues+1);

	QVector<QwtIntervalSample> samples(numValues);  
	for (int i=0; i<distribution.size(); i++)
	{
		QwtInterval interval(i, i+1); 
		interval.setBorderFlags(QwtInterval::ExcludeMaximum);  
		samples[i] = QwtIntervalSample(distribution[i], interval);  
	}
	histogram->setData(new QwtIntervalSeriesData(samples));

	plot->detachItems();
	histogram->attach(plot);
	plot->replot();
}

inline QColor qtJetColorMap(double value, double min = 0.0, double max = 1.0){
	unsigned char rgb[3];
	unsigned char c1=144;
	float max4=(max-min)/4;
	value-=min;
	if(value==HUGE_VAL)
	{rgb[0]=rgb[1]=rgb[2]=255;}
	else if(value<0)
	{rgb[0]=rgb[1]=rgb[2]=0;}
	else if(value<max4)
	{rgb[0]=0;rgb[1]=0;rgb[2]=c1+(unsigned char)((255-c1)*value/max4);}
	else if(value<2*max4)
	{rgb[0]=0;rgb[1]=(unsigned char)(255*(value-max4)/max4);rgb[2]=255;}
	else if(value<3*max4)
	{rgb[0]=(unsigned char)(255*(value-2*max4)/max4);rgb[1]=255;rgb[2]=255-rgb[0];}
	else if(value<max)
	{rgb[0]=255;rgb[1]=(unsigned char)(255-255*(value-3*max4)/max4);rgb[2]=0;}
	else {rgb[0]=255;rgb[1]=rgb[2]=0;}
	return QColor(rgb[0],rgb[1],rgb[2]);
}

int IconWidget::getSeletedInteractionIdx()
{
	int idx = -1, k =0;
	QTreeWidgetItemIterator it(ui->ibsTreeWidget);

	if (ui->ibsTreeWidget->selectedItems().isEmpty())
		return -1;

	while (*it)
	{
		if (ui->ibsTreeWidget->currentItem()->text(0) == (*it)->text(0))
		{
			idx = k;
			break;
		}
		++k;
		++it;
	}
	
	if (idx == -1)
		return idx;
	else
		return func->currScene->interactionsIdx[idx];
}

int IconWidget::getSelectedHierarchyIdx()
{
	return ui->hierarchyListWidget->currentRow();
}

void IconWidget::resetIbsFeatures()
{
	// Geometrical features
	for (int i=0; i<6; i++)
	{
		interactionFeaturePlot[i]->detachItems();
		interactionFeaturePlot[i]->replot();
	}
}

void IconWidget::setInteractionFeatureYScale( QVector<double> yMax )
{
	if (yMax.size() != interactionFeaturePlot.size()-1)
	{
		return;
	}

	for (int i=0; i<interactionFeaturePlot.size()-1; i++)
	{
		interactionFeaturePlot[i]->setAxisScale( QwtPlot::yLeft, 0.0, yMax[i]);
	}
}

void IconWidget::setClusterFeatureYScale( QVector<double> yMax )
{
	if (yMax.size() != clusterFeaturePlot.size()-1)
	{
		return;
	}

	for (int i=0; i<clusterFeaturePlot.size()-1; i++)
	{
		clusterFeaturePlot[i]->setAxisScale( QwtPlot::yLeft, 0.0, yMax[i]);
	}
}

Vec3d IconWidget::getUprightDirection()
{
	Vec3d upright;

	if (ui->xRadioButton->isChecked())
	{
		upright = Vec3d(1, 0, 0);
	}
	else if (ui->yRadioButton->isChecked())
	{
		upright = Vec3d(0, 1, 0);
	}
	else if (ui->zRadioButton->isChecked())
	{
		upright = Vec3d(0, 0, 1);
	}

	return upright;
}

void IconWidget::setUprightDirection(Vec3d upright)
{
	if (upright == Vec3d(1, 0, 0))
	{
		ui->xRadioButton->setChecked(true);
	}
	else if (upright == Vec3d(0, 1, 0))
	{
		ui->yRadioButton->setChecked(true);
	}
	else if (upright == Vec3d(0, 0, 1))
	{
		ui->zRadioButton->setChecked(true);
	}
}

void IconWidget::setObjectDrawMode( OBJECT_DRAW_MODE m )
{
	if (m == DRAW_MESH)
	{
		ui->meshRadioButton->setChecked(true);
	}
	else if (m == DRAW_SAMPLE)
	{
		ui->sampleRadioButton->setChecked(true);
	}
	else if (m == DRAW_BBOX_OBJ)
	{
		ui->boxObjRadioButton->setChecked(true);
	}
	else if (m == DRAW_BBOX_COMP)
	{
		ui->meshRadioButton->setChecked(true);
	}
	else if (m == DRAW_NONE)
	{		
		ui->noneRadioButton->setChecked(true);
	}
	else if (m == DRAW_WIRE_ORIG)
	{
		ui->wireOrigRadioButton->setChecked(true);
	}
	emit(objectDrawModeChanged(m));
}

void IconWidget::setCurrTabToDescriptor()
{
	ui->tabWidget->setCurrentIndex(1);
}

void IconWidget::getParaFromScene( Scene* scene )
{
	if (scene)
	{
		DistParameter distPara = scene->distPara;

		if (distPara.distType == L1)
		{
			ui->L1RadioButton->setChecked(true);
		}
		else
		{
			ui->emdRadioButton->setChecked(true);
		}

		ui->symRatioSpinBox->setValue(distPara.symRatio);

		ui->ibsWeightSpinBox->setValue(distPara.weight[0]);
		ui->regionWeightSpinBox->setValue(distPara.weight[1]);

		ui->ibsPfhSpinBox->setValue(distPara.ibsWeight[0]);
		ui->ibsDirSpinBox->setValue(distPara.ibsWeight[1]);
		ui->ibsDistSpinBox->setValue(distPara.ibsWeight[2]);

		ui->regionPfhSpinBox->setValue(distPara.regionWeight[0]);
		ui->regionDirSpinBox->setValue(distPara.regionWeight[1]);
		ui->regionHeightSpinBox->setValue(distPara.regionWeight[2]);

		// parameters related to interaction hierarchy construction
		HierParameter hierPara = scene->hierPara;

		ui->multiTreeRadioButton->setChecked(hierPara.useMultipleTree);		
		ui->singleTreeRadioButton->setChecked(!hierPara.useMultipleTree);
		
		ui->similarSpinBox->setValue(hierPara.similarThreshold);
		ui->lowSpinBox->setValue(hierPara.lowerBound);
		ui->highSpinBox->setValue(hierPara.higherBound);
	}
}

void IconWidget::assignParaToScene( Scene* scene )
{
	if (scene)
	{
		scene->clusterType = AHC;

		// parameters related to distance between interactions
		if (ui->L1RadioButton->isChecked())
		{
			scene->distPara.distType = L1;
		}
		else
		{
			scene->distPara.distType = EMD;
		}

		//scene->distPara.useTopo = ui->topoCheckBox->isChecked();	
		scene->distPara.useTopo = false;
		scene->distPara.symRatio = ui->symRatioSpinBox->value();

		scene->distPara.weight[0] = ui->ibsWeightSpinBox->value();
		scene->distPara.weight[1] = ui->regionWeightSpinBox->value();

		scene->distPara.ibsWeight[0] = ui->ibsPfhSpinBox->value();
		scene->distPara.ibsWeight[1] = ui->ibsDirSpinBox->value();
		scene->distPara.ibsWeight[2] = ui->ibsDistSpinBox->value();

		scene->distPara.regionWeight[0] = ui->regionPfhSpinBox->value();
		scene->distPara.regionWeight[1] = ui->regionDirSpinBox->value();
		scene->distPara.regionWeight[2] = ui->regionHeightSpinBox->value();

		// parameters related to interaction hierarchy construction
		scene->hierPara.useMultipleTree = ui->multiTreeRadioButton->isChecked();
		scene->hierPara.topDown = true;
		scene->hierPara.similarThreshold = ui->similarSpinBox->value();
		scene->hierPara.lowerBound = ui->lowSpinBox->value();
		scene->hierPara.higherBound = ui->highSpinBox->value();
	}
}

void IconWidget::normalizeWeight()
{
	ui->regionWeightSpinBox->setValue( 1.0 - ui->ibsWeightSpinBox->value());
}

void IconWidget::normalizeIbsWeightBasedOnPfh()
{
	double pfh = ui->ibsPfhSpinBox->value();
	double dir = ui->ibsDirSpinBox->value();
	if (pfh+dir >= 1)
	{
		ui->ibsDistSpinBox->setValue(0);
		ui->ibsDirSpinBox->setValue(1.0 - pfh);
	}
	else
	{
		ui->ibsDistSpinBox->setValue(1.0 - pfh - dir);
	}
}

void IconWidget::normalizeIbsWeightBasedOnDir()
{
	double pfh = ui->ibsPfhSpinBox->value();
	double dir = ui->ibsDirSpinBox->value();
	if (pfh+dir >= 1)
	{
		ui->ibsDistSpinBox->setValue(0);
		ui->ibsPfhSpinBox->setValue(1.0 - dir);
	}
	else
	{
		ui->ibsDistSpinBox->setValue(1.0 - pfh - dir);
	}
}

void IconWidget::normalizeRegionWeightBasedOnPfh()
{
	double pfh = ui->regionPfhSpinBox->value();
	double dir = ui->regionDirSpinBox->value();
	if (pfh+dir >= 1)
	{
		ui->regionHeightSpinBox->setValue(0);
		ui->regionDirSpinBox->setValue(1.0 - pfh);
	}
	else
	{
		ui->regionHeightSpinBox->setValue(1.0 - pfh - dir);
	}
}

void IconWidget::normalizeRegionWeightBasedOnDir()
{
	double pfh = ui->regionPfhSpinBox->value();
	double dir = ui->regionDirSpinBox->value();
	if (pfh+dir >= 1)
	{
		ui->regionHeightSpinBox->setValue(0);
		ui->regionPfhSpinBox->setValue(1.0 - dir);
	}
	else
	{
		ui->regionHeightSpinBox->setValue(1.0 - pfh - dir);
	}
}

void IconWidget::assignParameters(DistParameter &distPara, HierParameter &hierPara)
{
	// parameters related to distance between interactions
	if (ui->L1RadioButton->isChecked())
	{
		distPara.distType = L1;
	}
	else
	{
		distPara.distType = EMD;
	}

	distPara.useTopo = false;
	distPara.symRatio = ui->symRatioSpinBox->value();

	distPara.weight[0] = ui->ibsWeightSpinBox->value();
	if (distPara.weight[0] < 10e-10)
		distPara.weight[0] = 0;
	distPara.weight[1] = 1 - distPara.weight[0];

	distPara.ibsWeight[0] = ui->ibsPfhSpinBox->value();
	distPara.ibsWeight[1] = ui->ibsDirSpinBox->value();
	distPara.ibsWeight[2] = ui->ibsDistSpinBox->value();

	distPara.regionWeight[0] = ui->regionPfhSpinBox->value();
	distPara.regionWeight[1] = ui->regionDirSpinBox->value();
	distPara.regionWeight[2] = ui->regionHeightSpinBox->value();

	// parameters related to interaction hierarchy construction
	hierPara.methodtype = true;
	hierPara.useMultipleTree = ui->multiTreeRadioButton->isChecked();
	hierPara.topDown = true;
	hierPara.similarThreshold = ui->similarSpinBox->value();
	hierPara.lowerBound = ui->lowSpinBox->value();
	hierPara.higherBound = ui->highSpinBox->value();
}

void IconWidget::updateInteractionColor(QTreeWidgetItem * item)
{
	if (func->currScene)
	{
		int idx = -1;
		QTreeWidgetItemIterator it(ui->ibsTreeWidget);
		while (*it) 
		{
			++idx;
			if ((*it)->text(0) == item->text(0))
				break;
			++it;			
		}
	
		QColor c = func->currScene->colorMap[idx];
		QColorDialog dialog(func->mainWindow());
		dialog.setCurrentColor(c);

		if (dialog.exec())
		{
			c = dialog.selectedColor();
			func->setInteractionColor(idx, c);
		}
	}
}

void IconWidget::changeBaseColor()
{
	if (func->currScene)
	{
		QColor c = func->currScene->baseColor;
		QColorDialog dialog(func->mainWindow());
		dialog.setCurrentColor(c);
		if (dialog.exec())
		{
			c = dialog.selectedColor();
			func->currScene->baseColor = c;
			func->drawArea()->updateGL();
		}
	}
}


void IconWidget::setCurrentTable(int idx)
{
	ui->tabWidget->setCurrentIndex(idx);
}

void IconWidget::enableComputeFeatures(int type)
{
	switch (type)
	{
	case 0:
		{
			ui->compFeatureForICON->setCheckable(true);
			ui->compFeatureForICON->setChecked(true);
			ui->compFeatureForISET->setCheckable(false);
			break;
		}
	case 1:
		{
			ui->compFeatureForISET->setCheckable(true);
			ui->compFeatureForISET->setChecked(true);
			ui->compFeatureForICON->setCheckable(false);
			break;
		}
	default:
		{
			ui->compFeatureForICON->setCheckable(false);
			ui->compFeatureForISET->setCheckable(false);
			ui->compFeatureForICON->setChecked(false);
			ui->compFeatureForISET->setChecked(false);
			break;
		}
	}

	if (ui->compFeatureForICON->isChecked() || ui->compFeatureForISET->isChecked())
	{
		ui->computeFeatureButton->setEnabled(true);
	}
	else
	{
		ui->computeFeatureButton->setEnabled(false);
	}
}