TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += Ising_wl_int_gsl2d.cpp

TARGET = 2D_WL_ASI

#LIBS += -L/usr/local/lib -lgsl -lgslcblas
LIBS += -lgsl -lgslcblas
