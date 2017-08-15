#pragma once
#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"
#include "interfaces/ModePluginDockWidget.h"
#include "SceneSet.h"
#include "RetrievalTool.h"

class IconWidget;
class RetrievalWidget;

class FuncPlugin: public SurfaceMeshModePlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "function.plugin.starlab")
    Q_INTERFACES(ModePlugin)

    QIcon icon(){ return QIcon(":/images/function.png"); }

public:
    FuncPlugin();
	~FuncPlugin();

    // Main functionality
    void create();
    void destroy(){}
    void decorate();

    bool keyPressEvent(QKeyEvent*);

    IconWidget* iconWidget;
	RetrievalWidget* retrievalWidget;
	ModePluginDockWidget* iconDockWidget;
	ModePluginDockWidget* retrievalDockWidget;

    // Always usable
    bool isApplicable() { return true; }

private:
	void resetCamera();
	void loadCurrentCameraSetting();
	void saveCurrentCameraSetting();

	void updateCurrentSetting();

	QImage getCurrentSnapshot(bool rectangle);
	void prepareSnapshot();
	void updateParameters(Scene * scene);

	void getParameters(DistParameter &distPara, HierParameter &hierPara);

public:
	SceneSet * sceneSet;
	Scene *currScene;
	int currSceneIdx;

	DrawParameter para;

	int currentFeatureType;			// 0 = ICON, 1 = ISET

	//For retrieval
	RetrievalTool * retrievalTool;
public:
	bool mousePressEvent(QMouseEvent*);

public slots:
	void loadScenes();
	void addScenes();

	void buildInteractionHierarchy();
	void loadInteractionHierarchy();

	void updateSelectedObjects();
	void updateSelectedIBS();
	void updateSelectedHierarchy();

	void setObjectDrawMode(OBJECT_DRAW_MODE m);
	void setIbsDrawMode(IBS_DRAW_MODE m);
	void setIbsSampleShow(bool show);
	void setIbsWeightShow(bool show);
	void setUprightDirection();	
	void setUprightDirectionForAll();

	void setInteractionShow(bool show);
	void setInteractionIBSShow(bool show);
	void setInteractionIRShow(bool show);
	void setInteractionObjShow(bool show);

	void outputIBS();
	void outputIR();

	void setCurrScene(int idx);
	void updateSnapshot();
	void deleteScene(int idx);


	void setCurrScene_Ratrieval(int idx);

	//For retrieval
	void loadSceneListForRetrieval();
	void loadPartDataForRetrieval();
	void loadShape2PoseDataForRetrieval();
	void loadShape2PoseDataForICON();
	void loadShape2PoseDataForPOSE();

	void doRetrieval_icon();
	void doRetrieval_lfd();
	void doRetrieval_lfd_icon();
	void doRetrieval_pose();
	void doRetrieval_ibsh();
	void doRetrieval_iset();
	void doRetrieval_geo();
	void doRetrieval_geo_icon();

	void computeFeatures();			// for ICON and ISET

	void loadRetrievalResults();
	void buildInteractionHierarchyForAllScenes();
	void buildInteractionHierarchyForSceneList();

public:
	QVector<double> getInteractionFeature(int idx);
	QVector<int> getBettiNumber();
	void updateRetrieval();
	void setInteractionColor(int idx, QColor color);

	void showMsgOnStatusBar(QString message);
};