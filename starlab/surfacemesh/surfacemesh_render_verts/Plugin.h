#pragma once
#include "SurfaceMeshPlugins.h"
class Plugin : public SurfaceMeshRenderPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "surfacemesh_render_verts.plugin.starlab")
    Q_INTERFACES(RenderPlugin)
  
public: 
    QString name() { return "Vertices as Dots"; }
    QIcon icon(){ return QIcon(":/icons/cloud.png"); }
    Renderer* instance();
};
