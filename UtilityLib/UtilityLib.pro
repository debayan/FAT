include($$[STARLAB])
include($$[SURFACEMESH])

QT += gui opengl xml svg

TEMPLATE = lib
CONFIG += staticlib

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

# Library name and destination
TARGET = UtilityLib
DESTDIR = $$PWD/$$CFG/lib

# qwt-6.1.0
LIBS += -L$$PWD/qwt-6.1.0/lib/$$CFG/ -lqwt
INCLUDEPATH += qwt-6.1.0/include

# rapid
LIBS += -L$$PWD/rapid-mine/lib/$$CFG/ -lRAPID
INCLUDEPATH += rapid-mine/inc

# pqp
INCLUDEPATH += pqp-2.0

# ALGLIB
LIBS += -L$$PWD/ALGLIB/lib/$$CFG/ -lALGLIB
INCLUDEPATH += ALGLIB/src


HEADERS += \
    DetectCollision.h \
    QuickMeshDraw.h \
    Sampler.h \
    primitives.h \
    OrientHelper.h \
    SegMeshLoader.h \
    UtilityGlobal.h \
	Colormap.h \
	BasicTable.h \
	histogramDistance.h \
	HierarchicalCluster.h \
	tree.h \
	hungarian/hungarian.h \
	LFD.h \
	writeOBJ.h
	
SOURCES += \
    DetectCollision.cpp \
    Sampler.cpp \	
    OrientHelper.cpp \
    SegMeshLoader.cpp \
    UtilityGlobal.cpp \
	HierarchicalCluster.cpp \
	hungarian/hungarian.c
 
FORMS += \

RESOURCES += \
