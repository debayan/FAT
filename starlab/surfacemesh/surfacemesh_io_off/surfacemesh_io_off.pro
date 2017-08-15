include($$[STARLAB])
include($$[SURFACEMESH])
QT += gui widgets xml opengl
StarlabTemplate(plugin)

HEADERS += surfacemesh_io_off.h
SOURCES += surfacemesh_io_off.cpp
INCLUDEPATH += ../surfacemesh/ ../../core/starlib/interfaces/ ../../core/starlib/



