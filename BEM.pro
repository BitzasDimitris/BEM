#-------------------------------------------------
#
# Project created by QtCreator 2017-07-29T19:50:57
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = BEM
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
        main.cpp \
        mainwindow.cpp \
    point.cpp \
    glview.cpp \
    problem.cpp \
    pointinput.cpp \
    node.cpp \
    lu_decomposition.cpp

HEADERS += \
        mainwindow.h \
    point.h \
    glview.h \
    problem.h \
    pointinput.h \
    node.h \
    lu_decomposition.hpp

FORMS += \
        mainwindow.ui \
    pointinput.ui

LIBS +=-LC:\Qt\Qt5.9.1\5.9.1\mingw53_32\lib\libQt5OpenGL.a -lopengl32
