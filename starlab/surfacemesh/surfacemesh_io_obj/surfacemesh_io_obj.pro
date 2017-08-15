include($$[STARLAB])
include($$[SURFACEMESH])
#StarlabTemplate(plugin)
QT += widgets
INCLUDEPATH += ../surfacemesh/ ../../core/starlib/interfaces/
HEADERS += \
    plugin.h
SOURCES += \
    plugin.cpp
