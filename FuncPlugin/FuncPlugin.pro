include($$[STARLAB])
include($$[SURFACEMESH])
include($$[QHULL])
StarlabTemplate(plugin)

QT += gui opengl xml svg printsupport

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

# Foldabilizer library
LIBS += -L$$PWD/../FuncLib/$$CFG/lib -lFuncLib
INCLUDEPATH += ../FuncLib

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
    FuncPlugin.h \
    IconWidget.h \
    RetrievalWidget.h

SOURCES += \
    FuncPlugin.cpp \
    IconWidget.cpp \
    RetrievalWidget.cpp
 
FORMS += \
    IconWidget.ui \
    RetrievalWidget.ui

RESOURCES += \
    FuncPlugin.qrc
 
