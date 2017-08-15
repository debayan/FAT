#include "FuncPlugin.h"
#include "IconWidget.h"
#include "RetrievalWidget.h"
#include <qwt_plot_renderer.h>
#include <QFileDialog>

FuncPlugin::FuncPlugin() : 
	iconWidget(NULL), 
	retrievalWidget(NULL),
	iconDockWidget(NULL), 
	retrievalDockWidget(NULL)
{
	sceneSet = NULL;
	currScene = NULL;
	currSceneIdx = -1;
	currentFeatureType = -1;

	para.objMode = DRAW_MESH;
	para.ibsMode = DRAW_ENTIRE;
	para.drawIbsSample = false;
	para.drawIbsWeight = true;
	para.drawClusterIbs = true;
	para.drawClusterIbsWeight = true;
	para.drawClusterIbsSample = false;
	para.drawInteraction = false;
	para.drawInteractionIBS = true;
	para.drawInteractionIR = true;
	para.drawInteractionObj = true;

	retrievalTool = NULL;
}

FuncPlugin::~FuncPlugin()
{
	if (sceneSet) 
	{
		delete sceneSet;
	}
	if (retrievalTool) 
	{
		delete retrievalTool;
	}
}

void FuncPlugin::create()
{
	if (this->iconWidget || this->retrievalWidget)
		return;

	iconDockWidget = new ModePluginDockWidget("ICON", mainWindow());
	iconWidget = new IconWidget(this);
	iconDockWidget->setWidget(iconWidget);
	iconDockWidget->setFixedWidth(400);
	mainWindow()->addDockWidget(Qt::LeftDockWidgetArea, iconDockWidget);

	retrievalDockWidget = new ModePluginDockWidget("Retrieval", mainWindow());
	retrievalWidget = new RetrievalWidget(this);
	retrievalDockWidget->setWidget(retrievalWidget);
	retrievalDockWidget->setFixedWidth(400);
	mainWindow()->addDockWidget(Qt::RightDockWidgetArea, retrievalDockWidget);

	drawArea()->setFPSIsDisplayed(false);
}

void FuncPlugin::decorate()
{
	if (currScene)
	{
		currScene->draw(para);
	}
}

bool FuncPlugin::keyPressEvent(QKeyEvent *keyEvent)
{
	if (sceneSet) {
		int sceneCount = sceneSet->scenes.size();

		if (keyEvent->key() == Qt::Key_PageUp) {
			currSceneIdx = ((currSceneIdx - 1) + sceneCount) % sceneCount;
			setCurrScene(currSceneIdx);

			updateSelectedIBS();

			iconWidget->setSelectedScene(currSceneIdx);
			return true;
		}

		if (keyEvent->key() == Qt::Key_PageDown) {
			currSceneIdx = (currSceneIdx + 1) % sceneCount;
			setCurrScene(currSceneIdx);

			updateSelectedIBS();
			
			iconWidget->setSelectedScene(currSceneIdx);
			return true;
		}
	}
	return false;
}

void FuncPlugin::loadCurrentCameraSetting()
{
	if (currScene)
	{
		QDomDocument document;
		QFile f(currScene->filename + ".xml");
		if (f.open(QIODevice::ReadOnly))
		{
			document.setContent(&f);
			f.close();
		}
		// Parse the DOM tree
		QDomElement main = document.documentElement();
		drawArea()->camera()->initFromDOMElement(main);
		drawArea()->camera()->initFromDOMElement(main); // bug; use to set the right camera position
	}
}

void FuncPlugin::saveCurrentCameraSetting()
{
	if (currScene)
	{
		QDomDocument document("myCamera");
		document.appendChild( drawArea()->camera()->domElement("Camera", document) );
		QFile f(currScene->filename + ".xml");
		if (f.open(QIODevice::WriteOnly))
		{
			QTextStream out(&f);
			document.save(out, 2);
		}
	}
}

void FuncPlugin::resetCamera()
{
	if (currScene)
	{
		if (QFile(currScene->filename + ".xml").exists())
		{
			loadCurrentCameraSetting();
		}
		else
		{
			Eigen::AlignedBox3d bbox = currScene->bbox;

			drawArea()->camera()->setSceneRadius(bbox.sizes().norm() * 2);
			if (currScene->upright == Vec3d(0,1,0))
			{
				drawArea()->camera()->setUpVector(qglviewer::Vec(0,1,0));
				drawArea()->camera()->setPosition(qglviewer::Vec(2,1.5,2));
			}
			else
			{
				drawArea()->camera()->setUpVector(qglviewer::Vec(0,0,1));
				drawArea()->camera()->setPosition(qglviewer::Vec(2,-2,1.5));
			}
			drawArea()->camera()->lookAt(qglviewer::Vec());
			drawArea()->camera()->showEntireScene();
			double S = bbox.sizes().norm() * 0.2;
			bbox.extend(bbox.min() + (Vector3(-1,-1,-1) * S));
			bbox.extend(bbox.max() + (Vector3(1,1,1) * S));
			drawArea()->camera()->setRevolveAroundPoint(qglviewer::Vec(bbox.center()));
			drawArea()->camera()->setSceneCenter(qglviewer::Vec(bbox.center()));
			drawArea()->camera()->fitBoundingBox(qglviewer::Vec(bbox.min().data()),qglviewer::Vec(bbox.max().data()));

			saveCurrentCameraSetting();
		}

		drawArea()->updateGL();
	}
}

void FuncPlugin::updateSelectedObjects()
{
	if (!currScene)
	{
		return;
	}

	QVector<bool> selected = iconWidget->getObjectCheckState();

	if (selected.size() == currScene->objects.size())
	{
		currScene->setObjCentralState(selected);
		iconWidget->updateHierarchyList();
		iconWidget->updateIbsList();

		drawArea()->updateGL();
	}
}

void FuncPlugin::updateSelectedIBS()
{
	if (!currScene || iconWidget->isDeleting)
	{
		return;
	}
	
	QVector<bool> selected = iconWidget->getIbsCheckState();

	if (selected != currScene->lastSelectedInters)
	{
		// for hierarchy
		if (selected.size() == currScene->interactionsIdx.size())
		{
			currScene->setInteractionHierSelectState(selected);
			iconWidget->updateIbsListCheckState(currScene->isSelectedInterGroup);

			drawArea()->updateGL();
		}

		if (selected.size() == currScene->interactionsIdx.size())
		{
			currScene->setInteractionSelectState(currScene->isSelectedInterGroup);
			drawArea()->updateGL();
		}
	}

	// default
	drawArea()->updateGL();
}

void FuncPlugin::updateSelectedHierarchy()
{
	if (!currScene)
	{
		return;
	}
}

bool FuncPlugin::mousePressEvent( QMouseEvent*  event)
{
	if(event->modifiers() & Qt::SHIFT)
	{
		auto pos = event->pos();
		qglviewer::Vec orig, dir;
		drawArea()->camera()->convertClickToLine( pos, orig, dir );

		// Do intersection
		int seletedID = currScene->getSelectedObjectID(Eigen::Vector3d(orig), Eigen::Vector3d(dir));
		if (seletedID != -1)
		{			
			currScene->reverseObjCentralState(seletedID);
			iconWidget->updateObjectList();
		}

		drawArea()->updateGL();
	}

	return false;
}

void FuncPlugin::setObjectDrawMode( OBJECT_DRAW_MODE m )
{
	para.objMode = m;

	if (currScene)
	{
		drawArea()->updateGL();
	}
}

void FuncPlugin::setIbsDrawMode( IBS_DRAW_MODE m )
{
	para.ibsMode = m;

	if (currScene)
	{		
		drawArea()->updateGL();
	}
}

void FuncPlugin::setIbsSampleShow( bool show )
{
	para.drawIbsSample = show;

	if (currScene)
	{
		drawArea()->updateGL();
	}
}

void FuncPlugin::setIbsWeightShow( bool show )
{	
	para.drawIbsWeight = show;

	if (currScene)
	{
		drawArea()->updateGL();
	}
}

void FuncPlugin::setInteractionShow(bool show)
{
	para.drawInteraction = show;

	if (currScene)
	{
		drawArea()->updateGL();
	}
}

void FuncPlugin::setInteractionIBSShow(bool show)
{
	para.drawInteractionIBS = show;

	if (currScene)
	{
		drawArea()->updateGL();
	}
}

void FuncPlugin::setInteractionIRShow(bool show)
{
	para.drawInteractionIR = show;

	if (currScene)
	{
		drawArea()->updateGL();
	}
}

void FuncPlugin::setInteractionObjShow(bool show)
{
	para.drawInteractionObj = show;

	if (currScene)
	{
		drawArea()->updateGL();
	}
}

QVector<double> FuncPlugin::getInteractionFeature( int idx )
{
	QVector<double> distribution;

	int seletedIdx = iconWidget->getSeletedInteractionIdx();

	if (seletedIdx >= 0 && seletedIdx < currScene->allInteractions.size())
	{
		distribution = currScene->allInteractions[seletedIdx]->getFeature(idx);
	}
	
	return distribution;
}

QVector<int> FuncPlugin::getBettiNumber()
{
	QVector<int> bettiNumbers;

	if (currScene)
	{
		int selectedIdx = iconWidget->getSeletedInteractionIdx();
		if (selectedIdx >= 0 && selectedIdx < currScene->allInteractions.size())
		{
			IBS* ibs = currScene->allInteractions[selectedIdx]->ibs;

			if (ibs)
			{
				bettiNumbers = ibs->bettiNumbers;
			}
		}

	}
	
	return bettiNumbers;
}

void FuncPlugin::outputIBS()
{
	if (currScene)
	{
		QString ibsfile = currScene->filename + "_ibs" + ".off";
		int selectedIdx = iconWidget->getSeletedInteractionIdx();

		if (selectedIdx < 0 || selectedIdx > currScene->interIbsSet.size())
		{
			return;
		}

		if (para.drawInteraction) 
		{
			currScene->allInteractions[selectedIdx]->ibs->mesh->write(ibsfile.toStdString());
		}
		else 
		{
			QVector<int> origObj = currScene->allInteractions[selectedIdx]->obj->origIdx;
			QVector<int> selectedIBS(currScene->interObjIdx.size(), false);
			for (int i = 0; i < currScene->interObjIdx.size(); ++i)
			{
				if (origObj.contains(currScene->interObjIdx[i]))
					selectedIBS[i] = true;
			}

			for (int i = 0; i < currScene->interObjIdx.size(); ++i)
			{
				QString ibsfile = currScene->filename + "_ibs_partial_" + QString::number(i) + ".off";

				if (selectedIBS[i])
				{
					for (auto ibsIdx : currScene->interIbsSet[i])
					{
						if (ibsIdx < currScene->ibsSetScene.size())
							currScene->ibsSetScene[ibsIdx]->mesh->write(ibsfile.toStdString());
					}
				}
			}
		}		
	}
}

void FuncPlugin::outputIR()
{
	if (currScene)
	{
		if (para.drawInteraction)
		{
			int seletedIdx = iconWidget->getSeletedInteractionIdx();
			if (seletedIdx < 0 || seletedIdx > currScene->interIbsSet.size())
			{
				return;
			}

			FuncRegion *region = currScene->allInteractions[seletedIdx]->region;
			currScene->centralObj->saveIR(region);
		}
	}
}

void FuncPlugin::setUprightDirection()
{
	if (currScene)
	{
		Vec3d up = iconWidget->getUprightDirection();
		currScene->upright = up;
		currScene->saveCentralIdx();
	}
}

void FuncPlugin::setUprightDirectionForAll()
{
	if (sceneSet)
	{
		Vec3d up = iconWidget->getUprightDirection();
		for (int i=0; i<sceneSet->scenes.size(); i++)
		{
			if (sceneSet->scenes[i])
			{
				sceneSet->scenes[i]->upright = up;
				sceneSet->scenes[i]->saveCentralIdx();
			}

		}
	}
}

void FuncPlugin::buildInteractionHierarchy()
{
	if (currScene)
	{
		showMsgOnStatusBar("Computing ICON for current scene...");

		iconWidget->assignParaToScene(currScene);
		currScene->visualizeHierarchy = iconWidget->visualizeHierarchy();

		currScene->buildInteractionHierarchy();

		iconWidget->updateHierarchyList();
		iconWidget->updateIbsList();
		iconWidget->setInteractionFeatureYScale(currScene->yMaxHistInteraction);

		drawArea()->updateGL();
	}
}

void FuncPlugin::loadInteractionHierarchy()
{
	if (currScene)
	{
		iconWidget->assignParaToScene(currScene);	
		currScene->loadInteractionHierarchy();

		drawArea()->updateGL();
	}
}

void FuncPlugin::loadScenes()
{
 	QFileDialog dialog(mainWindow());
 	dialog.setDirectory(settings()->getString("lastUsedDirectory"));
 	dialog.setFileMode(QFileDialog::ExistingFiles);
 	dialog.setNameFilter(tr("Text Files (*.txt);; Obj Files (*.obj)"));	

	if (dialog.exec()) {
		if (!sceneSet)
			sceneSet = new SceneSet();
		else
			sceneSet->clearAllScenes();

		QStringList filenames = dialog.selectedFiles();
		sceneSet->addScenes(filenames);

		prepareSnapshot();

		setCurrScene(0);
		iconWidget->updateSceneListTable();

		QFileInfo info(filenames.back());
		QString lastUsedDir = info.absolutePath();
		settings()->set("lastUsedDirectory", lastUsedDir);

		updateCurrentSetting();
		drawArea()->updateGL();			// after loading SceneSet, re-draw the whole scenes
	}
}

void FuncPlugin::addScenes()
{
	QFileDialog dialog(mainWindow());
	dialog.setDirectory(settings()->getString("lastUsedDirectory"));
	dialog.setFileMode(QFileDialog::ExistingFiles);
	dialog.setNameFilter(tr("Text Files (*.txt);; Obj Files (*.obj)"));	

	if (dialog.exec())
	{
		QStringList filenames = dialog.selectedFiles();

		if (!sceneSet)
		{
			sceneSet = new SceneSet();
		}

		sceneSet->addScenes(filenames);

		prepareSnapshot();
		
		setCurrScene( sceneSet->scenes.size()-1 );
		
		iconWidget->updateSceneListTable();

		// update lastUsedDirectory
		QFileInfo info(filenames.back());
		QString lastUsedDir = info.absolutePath();
		settings()->set("lastUsedDirectory", lastUsedDir);

		drawArea()->updateGL();
	}
}

void FuncPlugin::setCurrScene( int idx )
{
	currSceneIdx = idx;
	if (idx == -1)
	{
		currScene = NULL;
		if (sceneSet)
		{
			delete sceneSet;
			sceneSet = NULL;
		}
	}
	else
	{
		currScene = sceneSet->scenes[idx];
		currScene->resetInteractionSelectState();
	}
	
	updateCurrentSetting();
	drawArea()->updateGL();
}

void FuncPlugin::setCurrScene_Ratrieval( int idx )
{
	retrievalTool->currIdx = idx;
}

void FuncPlugin::deleteScene( int idx )
{
	if (sceneSet->deleteScene(idx))
	{
		if (idx == sceneSet->scenes.size())
		{
			idx = sceneSet->scenes.size() - 1;
		}
	
		setCurrScene(idx);
		iconWidget->setSelectedScene(idx);

		drawArea()->updateGL();
	}
}

void FuncPlugin::updateCurrentSetting()
{
	if (currScene)
	{
		// the upright direction is initialized when the object is created
		iconWidget->setUprightDirection(currScene->upright);
		iconWidget->getParaFromScene(currScene);
	}
	
	iconWidget->updateObjectList();
	resetCamera();
}

void FuncPlugin::buildInteractionHierarchyForAllScenes()
{
	if (sceneSet)
	{
		showMsgOnStatusBar("Computing ICON for scene set...");

		QVector<int> idx = sceneSet->getScenesWithNoCentralObject();		// scenes w/o central objects

		if (!idx.isEmpty())
		{
			QString sceneNames;
			for (auto i:idx)
			{
				sceneNames += " #" + QString::number(i);
			}

			QMessageBox::warning(mainWindow(), "Warning", "Please select the functional object for scene:" + sceneNames) ;
		}
		else
		{
			for (auto s:sceneSet->scenes)
			{
				iconWidget->assignParaToScene(s);
				s->visualizeHierarchy = iconWidget->visualizeHierarchy();
			}
			sceneSet->buildInteractionHierarchy();

			iconWidget->updateHierarchyList();
			iconWidget->updateIbsList();
			iconWidget->setInteractionFeatureYScale(currScene->yMaxHistInteraction);

			drawArea()->updateGL();
		}
	}
}

inline QImage cropImage(QImage original, bool rectangle)
{
	int left = 0;
	bool found = false;
	for (int i=0; i<original.width(); i++)
	{
		for (int j=0; j<original.height(); j++)
		{
			QColor c = original.pixel(i,j);
			if (c.red() != 255 || c.green() != 255 || c.blue() != 255)
			{
				left = i-1;
				found = true;
				break;
			}
		}

		if (found)
		{
			break;
		}
	}
	
	int right = original.width();
	found = false;
	for (int i=original.width()-1; i>=0; i--)
	{
		for (int j=0; j<original.height(); j++)
		{
			QColor c = original.pixel(i,j);
			if (c.red() != 255 || c.green() != 255 || c.blue() != 255)
			{
				right = i+1;
				found = true;
				break;
			}
		}

		if (found)
		{
			break;
		}
	}

	int up = original.height();
	found = false;
	for (int j=0; j<original.height(); j++)
	{
		for (int i=0; i<original.width(); i++)
		{
			QColor c = original.pixel(i,j);
			if (c.red() != 255 || c.green() != 255 || c.blue() != 255)
			{
				up = j-1;
				found = true;
				break;
			}
		}

		if (found)
		{
			break;
		}
	}
	
	int down = original.height();
	found = false;
	for (int j=original.height()-1; j>=0; j--)
	{
		for (int i=0; i<original.width(); i++)
		{
			QColor c = original.pixel(i,j);
			if (c.red() != 255 || c.green() != 255 || c.blue() != 255)
			{
				down = j+1;
				found = true;
				break;
			}
		}

		if (found)
		{
			break;
		}
	}
	int width = right-left+1;
	int height = down-up+1;

	if (rectangle)
	{
		if (width > height)
		{
			int d = (width-height) / 2;
			up -= d;
			height = width;
		}
		else if (height > width)
		{
			int d = (height-width) / 2;
			left -= d;
			width = height;
		}

		left = (left>0) ? left:0;
		up = (up>0) ? up:0;
		width = ((width+left)<original.width()) ? width:original.width()-left;
		height = ((height+up)<original.height()) ? height:original.height()-up;
	}

	QRect rect(left, up, width, height);
	QImage cropped = original.copy(rect);

	return cropped;
}

QImage FuncPlugin::getCurrentSnapshot(bool rectangle = true)
{
	QImage snapshot = drawArea()->grabFrameBuffer();
	QImage cropped = cropImage(snapshot, rectangle);

	return cropped;
}

void FuncPlugin::prepareSnapshot()
{
	if (sceneSet)
	{
		for (auto s:sceneSet->scenes)
		{
			if (s->snapshot.isNull())
			{
				currScene = s;

				resetCamera();
				drawArea()->updateGL();
				
				QImage image = getCurrentSnapshot();
				image.save(s->filename + ".png");
				s->loadSnapshot();
			}
		}
	}
}

void FuncPlugin::updateSnapshot()
{
	if (currScene)
	{
		QImage image = getCurrentSnapshot();
		image.save(currScene->filename + ".png");
		currScene->loadSnapshot();

		// save the viewport
		saveCurrentCameraSetting();

		iconWidget->updateSceneListTable();
	}
}

//////////////////////////////////////////////////////////////////////////
// for retrieval
void FuncPlugin::loadSceneListForRetrieval()
{
	QString filename = QFileDialog::getOpenFileName(mainWindow(),tr("Open Scene List File"), "/home", tr("Scene List Files (*.sl)"));

	if(filename.isNull()) return;

	updateRetrieval();

	retrievalWidget->enablePoseRetriModeVisiable(false);
	retrievalWidget->enableShape2PoseDataRetriComp(false);
	retrievalWidget->enablePartDataRetriComp(false);
	retrievalWidget->enableOurDataRetriComp(true);

	if (retrievalTool)
	{
		delete retrievalTool;	
	}
	
	retrievalTool = new RetrievalTool(filename);
	retrievalTool->setCurrentRetrievalMode(NORMAL);
	setCurrScene_Ratrieval( 0 );
	retrievalWidget->updateSceneTable_Retrieval();
	retrievalWidget->updatePRWidget();			// update PR curve (re-plot)
	retrievalWidget->updatePoseRetrievalMode(-1);
}

void FuncPlugin::loadPartDataForRetrieval()
{
	QString filename = QFileDialog::getOpenFileName(mainWindow(),tr("Open Scene List File"), "/home", tr("Scene List Files (*.sl)"));

	if(filename.isNull()) return;

	retrievalWidget->enablePoseRetriModeVisiable(false);
	retrievalWidget->enableShape2PoseDataRetriComp(false);
	retrievalWidget->enableOurDataRetriComp(false);
	retrievalWidget->enablePartDataRetriComp(true);

	if (retrievalTool)
	{
		delete retrievalTool;	
	}

	retrievalTool = new RetrievalTool(filename);
	retrievalTool->setCurrentRetrievalMode(PART);
	setCurrScene_Ratrieval( 0 );
	retrievalWidget->updateSceneTable_Retrieval();
	retrievalWidget->updatePRWidget();			// update PR curve (re-plot)
	retrievalWidget->updatePoseRetrievalMode(-1);
	
}

void FuncPlugin::loadShape2PoseDataForRetrieval()
{
	QString filename = QFileDialog::getOpenFileName(mainWindow(),tr("Open Scene List File"), "/home", tr("Scene List Files (*.sl)"));

	if(filename.isNull()) return;
	
	retrievalWidget->enablePoseRetriModeVisiable(true);
	retrievalWidget->enableOurDataRetriComp(false);
	retrievalWidget->enablePartDataRetriComp(false);
	retrievalWidget->enableShape2PoseDataRetriComp(true);
	
	if (retrievalTool)
	{
		delete retrievalTool;
	}

	retrievalTool = new RetrievalTool(filename, 1);
	retrievalTool->setCurrentRetrievalMode(SHAPE2POSE);
	setCurrScene_Ratrieval( 0 );
	retrievalWidget->updateSceneTable_Retrieval();
	retrievalWidget->updatePRWidget();
	retrievalWidget->updatePoseRetrievalMode(1);
}

void FuncPlugin::loadShape2PoseDataForICON()
{
	QString filename;

	if (retrievalTool)
	{
		filename = retrievalTool->filePath + ".sl";
		delete retrievalTool;
	}
	
	retrievalTool = new RetrievalTool(filename);
	retrievalTool->setCurrentRetrievalMode(SHAPE2POSE);
	setCurrScene_Ratrieval( 0 );
	retrievalWidget->updateSceneTable_Retrieval();
	retrievalWidget->updatePRWidget();
	retrievalWidget->updatePoseRetrievalMode(0);
}

void FuncPlugin::loadShape2PoseDataForPOSE()
{
	QString filename;

	if (retrievalTool)
	{
		filename = retrievalTool->filePath + ".sl";
		delete retrievalTool;
	}

	retrievalTool = new RetrievalTool(filename, 1);
	retrievalTool->setCurrentRetrievalMode(SHAPE2POSE);
	setCurrScene_Ratrieval( 0 );
	retrievalWidget->updateSceneTable_Retrieval();
	retrievalWidget->updatePRWidget();
	retrievalWidget->updatePoseRetrievalMode(1);
}

// do retrieval - different method
void FuncPlugin::doRetrieval_pose()
{
	if (retrievalTool)
	{
		if (retrievalTool->currentRetrievalMode == SHAPE2POSE)				// shape2pose mode
		{
			QString filename = retrievalTool->filePath + ".sl";
			retrievalWidget->updatePoseRetrievalMode(1);
			loadShape2PoseDataForPOSE();
		}

		showMsgOnStatusBar("Computing POSE...");

		retrievalTool->computeDistance(POSE);
		retrievalTool->analyzePR(POSE);
		
		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(POSE);
		retrievalWidget->updatePRWidget();
	}
}

void FuncPlugin::doRetrieval_icon()
{
	if (retrievalTool)
	{
		if (retrievalTool->currentRetrievalMode == SHAPE2POSE)				// shape2pose mode
		{
			QString filename = retrievalTool->filePath + ".sl";
			retrievalWidget->updatePoseRetrievalMode(0);
			loadShape2PoseDataForICON();
		}

		QString paraString;
		QMessageBox loadMsgBox(QMessageBox::Question, tr(""), tr("Loading existing features?"), QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
		int loadRet = loadMsgBox.exec();

		if (loadRet == QMessageBox::Yes)
		{
			// if Yes, load features
			QFileDialog dialog(mainWindow());
			dialog.setFileMode(QFileDialog::Directory);
			if (!dialog.exec()) return;
			paraString = dialog.directory().dirName();

			if (!paraString.startsWith("(ICON)"))
			{
				debugBox("Please select folder whose name starts with \"(ICON)\"");
				return;
			}
		}
		else if (loadRet == QMessageBox::No)
		{
			iconWidget->setCurrentTable(1);

			QMessageBox msgBox(QMessageBox::Question, tr("Computing features"), tr("Please first set the parameters and then click \"Compute\" button on the ICON widget."), QMessageBox::Yes | QMessageBox::No);
			msgBox.setButtonText(QMessageBox::Yes, tr("OK"));
			msgBox.setButtonText(QMessageBox::No, tr("Cancel"));
			int ret = msgBox.exec();

			if (ret == QMessageBox::Yes)
			{
				// if Yes, enable computing ICON
				currentFeatureType = 0;
			}
			else
			{
				currentFeatureType = -1;
			}
			iconWidget->enableComputeFeatures(currentFeatureType);

			return;
		}
		else //if (loadRet == QMessageBox::Cancel)
		{
			return;
		}

		showMsgOnStatusBar("Computing ICON...");

		retrievalTool->computeDistance(ICON, paraString);
		retrievalTool->analyzePR(ICON);

		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(ICON);
		retrievalWidget->updatePRWidget();

		if (retrievalTool->currentRetrievalMode == NORMAL)
		{
			retrievalWidget->enableLFDICON();
		}

		if (retrievalTool->currentRetrievalMode == PART)
		{
			retrievalWidget->enableGEOICON();
		}
	}
}

void FuncPlugin::doRetrieval_lfd()
{
	if(retrievalTool)
	{
		showMsgOnStatusBar("Computing LFD...");

		retrievalTool->computeDistance(LFD);
		retrievalTool->analyzePR(LFD);

		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(LFD);
		retrievalWidget->updatePRWidget();
	}
}

void FuncPlugin::doRetrieval_lfd_icon()
{
	if(retrievalTool)
	{
		showMsgOnStatusBar("Computing LFD+ICON...");

		double weight = retrievalWidget->getCombWeight();
		retrievalTool->computeDistance(LFDICON, "", weight);
		retrievalTool->analyzePR(LFDICON);

		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(LFDICON);
		retrievalWidget->updatePRWidget();
	}
}

void FuncPlugin::doRetrieval_ibsh()
{
	if (retrievalTool)
	{
		showMsgOnStatusBar("Computing IBSH...");

		retrievalTool->ibshDepth = retrievalWidget->getDepth();
		retrievalTool->computeDistance(IBSH);
		retrievalTool->analyzePR(IBSH);

		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(IBSH);
		retrievalWidget->updatePRWidget();
	}	
}

void FuncPlugin::doRetrieval_iset()
{	
	if (retrievalTool)
	{
		QString paraString;
		if (QMessageBox::Yes == QMessageBox::question(mainWindow(), tr(""), tr("Loading existing features?"), QMessageBox::Yes | QMessageBox::No))
		{
			// if Yes, load features
			QFileDialog dialog(mainWindow());
			dialog.setFileMode(QFileDialog::Directory);
			if (!dialog.exec()) return;
			paraString = dialog.directory().dirName();
			
			if (!paraString.startsWith("(ISET)"))
			{
				debugBox("Please select folder whose name starts with \"(ISET)\"");
				return;
			}
		}
		else
		{
			iconWidget->setCurrentTable(1);

			QMessageBox msgBox(QMessageBox::Question, tr("Computing features"), tr("Please first set the parameters and then click \"Compute\" button on the ICON widget."), QMessageBox::Yes | QMessageBox::No);
			msgBox.setButtonText(QMessageBox::Yes, tr("OK"));
			msgBox.setButtonText(QMessageBox::No, tr("Cancel"));
			int ret = msgBox.exec();

			if (ret == QMessageBox::Yes)
			{
				// if Yes, enable computing ISET
				currentFeatureType = 1;
			}
			else
			{
				currentFeatureType = -1;		// nothing
			}
			iconWidget->enableComputeFeatures(currentFeatureType);

			return;
		}

		showMsgOnStatusBar("Computing ISET...");

		retrievalTool->computeDistance(ISET, paraString);
		retrievalTool->analyzePR(ISET);

		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(ISET);
		retrievalWidget->updatePRWidget();
	}
}

void FuncPlugin::loadRetrievalResults()
{
	if (sceneSet)
	{
		sceneSet->clearAllScenes();
	}
	else
	{
		sceneSet = new SceneSet();
	}

	if (retrievalTool && retrievalTool->prReady() )
	{
		int returnRetrievalNum = retrievalWidget->getReturnRetrievalNum() + 1;

		QStringList filenames = retrievalTool->returnRetrievalResultName(returnRetrievalNum);
		sceneSet->addScenes(filenames);

		setCurrScene( sceneSet->scenes.size()-1 );
		iconWidget->updateSceneListTable();				// update icon widget

		retrievalWidget->setCurrTabToScenes();
		iconWidget->setCurrentTable(0);

		updateCurrentSetting(); 
		drawArea()->updateGL();
	}
}

void FuncPlugin::updateRetrieval()
{
	if (retrievalTool)
	{
		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(retrievalTool->currType);
		retrievalWidget->updatePRWidget();
	}
}

void FuncPlugin::getParameters(DistParameter &distPara, HierParameter &hierPara)
{
	iconWidget->assignParameters(distPara, hierPara);
}

void FuncPlugin::buildInteractionHierarchyForSceneList()
{
	QString filename = QFileDialog::getOpenFileName(mainWindow(),tr("Open Scene List File"), "/home", tr("Scene List Files (*.sl)"));
	if(filename.isNull()) return;

	showMsgOnStatusBar("Computing ICON for scene list...");

	// generate parameters
	DistParameter distPara;
	HierParameter hierPara;
	getParameters(distPara, hierPara);

	QString resultFolder = QFileInfo(filename).dir().absolutePath() + "/(ICON)" + distPara.toString() + hierPara.toString();
	QDir dir(resultFolder);
	if (!dir.exists())
	{
		dir.mkpath(resultFolder);
	}

	// compute the corresponding features
	if (sceneSet)
	{
		delete sceneSet;
	}
	sceneSet = new SceneSet();
	sceneSet->computeFeatures(filename, distPara, hierPara);

	debugBox("Finish computing ICON");
}

void FuncPlugin::doRetrieval_geo()
{
	if(retrievalTool)
	{
		showMsgOnStatusBar("Computing GEO...");

		retrievalTool->computeDistance(GEO);
		retrievalTool->analyzePR(GEO);

		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(GEO);
		retrievalWidget->updatePRWidget();
	}
}

void FuncPlugin::doRetrieval_geo_icon()
{
	if(retrievalTool)
	{
		showMsgOnStatusBar("Computing GEO+ICON...");

		double weight = retrievalWidget->getPartCombWeight();
		retrievalTool->computeDistance(GEOICON, "", weight);
		retrievalTool->analyzePR(GEOICON);

		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(GEOICON);
		retrievalWidget->updatePRWidget();
	}
}

void FuncPlugin::setInteractionColor(int idx, QColor color)
{
	if (currScene && idx >= 0 && idx < currScene->allInteractions.size())
	{
		currScene->colorMap[idx] = color;
		drawArea()->updateGL();
	}
}

void FuncPlugin::showMsgOnStatusBar(QString message)
{
	mainWindow()->statusBar()->showMessage(message, 5000);		// 5000 ms = 5 s
}

void FuncPlugin::computeFeatures()
{
	DistParameter distPara;
	HierParameter hierPara;
	getParameters(distPara, hierPara);

	if (currentFeatureType == 0)
	{
		showMsgOnStatusBar("Computing ICON...");
		retrievalTool->computeICONFeatures(distPara, hierPara);

		// compute distance
		QString paraString = "(ICON)" + distPara.toString() + hierPara.toString();
		retrievalTool->computeDistance(ICON, paraString);
		retrievalTool->analyzePR(ICON);

		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(ICON);
		retrievalWidget->updatePRWidget();
		retrievalWidget->enableLFDICON();
		retrievalWidget->enableGEOICON();
	}
	else if (currentFeatureType == 1)
	{
		showMsgOnStatusBar("Computing ISET...");
		retrievalTool->computeISETFeatures(distPara);

		// compute distance
		QString paraString = "(ISET)" + distPara.toString();
		retrievalTool->computeDistance(ISET, paraString);
		retrievalTool->analyzePR(ISET);

		QVector<int> retrievalResults = retrievalTool->returnRetrievalResult(11);
		retrievalWidget->updateSceneTable_Result(retrievalResults);
		retrievalWidget->updateFeatureType(ISET);
		retrievalWidget->updatePRWidget();
	}
	else
		return;
}