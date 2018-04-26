#ifndef GLVIEW_H
#define GLVIEW_H
#include <QOpenGLWidget>
#include "point.h"
#include <QPainter>
#include <QFont>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QCursor>
#include <math.h>

class GLView: public QOpenGLWidget
{
    Q_OBJECT
public:
    explicit GLView(QWidget *parent = 0);
    void draw();
    void SetVecs(std::vector<Point> p,std::vector<Point> v);
    void Clear();
    void SetGridShow(bool);
    void Restore();
protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
private:
    std::vector<Point> P,PIN,VP,V;
    bool pshow=false;
    Point point2show;
    bool bpoint2show=false;
    int Width,Height;
    float left=-2,right=2,top=2,bottom=-2;
    float viewfactor=4;
    float aspect;
    bool gridshow=true;
    float xstepper=0.5,ystepper=0.5;
    int xsteps,ysteps;
    Point LastPos;
    float movespeed=0.3;

    void UpdateBoundaries();
    Point FromScreen2Space(Point p);
    Point FromSpace2Screen(Point p);
    void DrawGrid();
    void UpdateGrid();
signals:
    void PointClicked(float x, float y);
private slots:
    void SetPoints(std::vector<Point> p);
    void SetInnerPoints(std::vector<Point> pin);
};

#endif // GLVIEW_H
