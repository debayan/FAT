include($$[STARLAB])
include($$[SURFACEMESH])
include($$[QHULL])

QT += gui opengl xml svg

TEMPLATE = lib
CONFIG += staticlib

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

# Library name and destination
TARGET = FuncLib
DESTDIR = $$PWD/$$CFG/lib

# Utility library
LIBS += -L$$PWD/../UtilityLib/$$CFG/lib -lUtilityLib
INCLUDEPATH += ../UtilityLib

# qwt-6.1.0
LIBS += -L$$PWD/../UtilityLib/qwt-6.1.0/lib/$$CFG/ -lqwt
INCLUDEPATH += ../UtilityLib/qwt-6.1.0/include

# rapid
LIBS += -L$$PWD/../UtilityLib/rapid-mine/lib/$$CFG/ -lRAPID
INCLUDEPATH += ../UtilityLib/rapid-mine/inc

# pqp
INCLUDEPATH += ../UtilityLib/pqp-2.0

# ALGLIB
LIBS += -L$$PWD/../UtilityLib/ALGLIB/lib/$$CFG/ -lALGLIB
INCLUDEPATH += ../UtilityLib/ALGLIB/src

HEADERS += \
    Object.h \
    Scene.h \
    IBS.h \
    IbsGenerator.h \
	FuncRegion.h \ 
    SceneSet.h \
	Interaction.h \
	InterHierarchy.h \
	Visualization.h \
	DistMeasure.h \
	CommonSubtree.h \
	RetrievalTool.h \
	SceneHierarchy.h \
	InterSet.h \
	SymGroup.h \

SOURCES += \
    Object.cpp \
    Scene.cpp \
    IBS.cpp \
    IbsGenerator.cpp \ 
	FuncRegion.cpp \
    SceneSet.cpp \
	Interaction.cpp \
	InterHierarchy.cpp \
	DistMeasure.cpp \
	CommonSubtree.cpp \
	RetrievalTool.cpp \
	SceneHierarchy.cpp \
	InterSet.cpp \
	SymGroup.cpp

FORMS += \

RESOURCES += \
