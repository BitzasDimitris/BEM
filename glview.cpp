#include "glview.h"
#include <QDebug>

GLView::GLView(QWidget *parent): QOpenGLWidget(parent)
{
}


void GLView::initializeGL()
{
    glClearColor(0,0,0,0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    UpdateGrid();

}

void GLView::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -10.0);
    glRotatef(0, 1.0, 0.0, 0.0);
    glRotatef(0, 0.0, 1.0, 0.0);
    glRotatef(0, 0.0, 0.0, 1.0);
    draw();
}

void GLView::resizeGL(int width, int height)
{
    int side = qMin(width, height);
    glViewport((width - side) / 2, (height - side) / 2, side, side);
    Width=width;
    Height=height;
    aspect=((float)Height)/Width;
    right=left+viewfactor;
    bottom=left*aspect;
    top=bottom+viewfactor*aspect;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(left,right,bottom,top,1.0,15.0);
    glMatrixMode(GL_MODELVIEW);

}

void GLView::mouseReleaseEvent(QMouseEvent *event){
    Point p=Point(event->pos().x(),event->pos().y());
    if(p.distance(LastPos)<0.01){
        p=FromScreen2Space(p);
        emit PointClicked(p.x,p.y);
        float minp=-1;
        int minindexp=-1;

        for(int i=0;i<P.size();i++){
            if(minp<0||minp>p.distance(P[i])){
                minp=p.distance(P[i]);
                minindexp=i;
            }
        }
        float minpin=-1;
        int minindexpin=-1;

        for(int i=0;i<PIN.size();i++){
            if(minpin<0||minpin>p.distance(PIN[i])){
                minpin=p.distance(PIN[i]);
                minindexpin=i;
            }
        }
        float delta=viewfactor*0.02;
        if(minp<minpin&&minp>0&&minpin>0){
            if(minp<delta){
                point2show=P[minindexp];
                pshow=true;
                bpoint2show=true;
            }
            else{
                pshow=false;
            }
        }
        else if(minpin>0&&minp>0){
            if(minpin<delta){
                point2show=PIN[minindexpin];
                pshow=true;
                bpoint2show=false;
            }
            else{
                pshow=false;
            }
        }
        else if(minpin>=0){
            if(minpin<delta){
                point2show=PIN[minindexpin];
                pshow=true;
                bpoint2show=false;
            }
            else{
                pshow=false;
            }
        }
        else if(minp>=0){
            if(minp<delta){
                point2show=P[minindexp];
                pshow=true;
                bpoint2show=true;
            }
            else{
                pshow=false;
            }
        }
        update();
    }
}

void GLView::mousePressEvent(QMouseEvent *event)
{
    QPoint mp = event->pos();
    Point p=Point(mp.x(),mp.y());
    LastPos=p;

}

void GLView::mouseMoveEvent(QMouseEvent *event)
{
    QPoint mp = event->pos();
    Point p=Point(mp.x(),mp.y());
    if(event->buttons()& Qt::LeftButton){
        p=FromScreen2Space(LastPos)-FromScreen2Space(p);
        left+=p.x;
        bottom+=p.y;
    }
    right=left+viewfactor;
    top=bottom+viewfactor*aspect;
    LastPos=Point(mp.x(),mp.y());
    UpdateGrid();
    update();
}

void GLView::wheelEvent(QWheelEvent *event){
    float zoom=event->delta()/8;
    float newviewfactor=viewfactor/(1+zoom/360.0);
    Point p=FromScreen2Space(Point(event->pos().x(),event->pos().y()));
    if(zoom>0){
        left=p.x-newviewfactor/2.0;
        right=left+newviewfactor;
        bottom=p.y-newviewfactor*aspect/2.0;
        top=bottom+newviewfactor*aspect;
    }
    else{
        left=left-(newviewfactor-viewfactor)/2.0;
        right=left+viewfactor;
        bottom=bottom-(newviewfactor-viewfactor)*aspect/2.0;
        top=bottom+viewfactor*aspect;
    }
    viewfactor=newviewfactor;
    QCursor::setPos(this->mapToGlobal(QPoint(Width/2,Height/2)));
    UpdateGrid();
    update();
}

void GLView::UpdateBoundaries(){
    float minx,miny,maxx,maxy;
    if(P.size()>0){
        minx=P[0].x;
        maxx=minx;
        miny=P[0].y;
        maxy=miny;
        for(int i=0;i<P.size();i++){
            if(minx>P[i].x){
                minx=P[i].x;
            }
            if(maxx<P[i].x){
                maxx=P[i].x;
            }
            if(miny>P[i].y){
                miny=P[i].y;
            }
            if(maxy<P[i].y){
                maxy=P[i].y;
            }
        }
    }
    if(PIN.size()>0){
        if(P.size()==0){
            minx=PIN[0].x;
            maxx=minx;
            miny=PIN[0].y;
            maxy=miny;
        }
        for(int i=0;i<PIN.size();i++){
            if(minx>PIN[i].x){
                minx=PIN[i].x;
            }
            if(maxx<PIN[i].x){
                maxx=PIN[i].x;
            }
            if(miny>PIN[i].y){
                miny=PIN[i].y;
            }
            if(maxy<PIN[i].y){
                maxy=PIN[i].y;
            }
        }
    }
    if(P.size()>0||PIN.size()>0){
        left=minx-0.3;
        right=maxx+0.3;
        bottom=miny-0.3;
        top=maxy+0.3;
        if(right-left>top-bottom){
            viewfactor=right-left;
            float diff=(viewfactor-(top-bottom))/2.0;
            top+=diff;
            bottom-=diff;
        }
        else{
            viewfactor=(top-bottom)/aspect;
            float diff=(viewfactor-(right-left))/2.0;
            right+=diff;
            left-=diff;
        }

    }
    UpdateGrid();
}

void GLView::UpdateGrid(){
    xsteps=(int)((ceil(right)-floor(left))/xstepper);
    while(xsteps<=10){
        xstepper/=2.0;
        xsteps=(int)((ceil(right)-floor(left))/xstepper);
    }
    while(xsteps>=14){
        xstepper*=2.0;
        xsteps=(int)((ceil(right)-floor(left))/xstepper);
    }
    ysteps=(int)((ceil(top)-floor(bottom))/ystepper);
    while(ysteps<=10){
        ystepper/=2.0;
        ysteps=(int)((ceil(top)-floor(bottom))/ystepper);
    }
    while(ysteps>=14){
        ystepper*=2.0;
        ysteps=(int)((ceil(top)-floor(bottom))/ystepper);
    }
    if(xstepper<ystepper){
        ystepper=xstepper;
        ysteps=(int)((ceil(top)-floor(bottom))/ystepper);
    }
    else{
        xstepper=ystepper;
        xsteps=(int)((ceil(right)-floor(left))/xstepper);
    }
}

void GLView::Clear(){
    P.clear();
    PIN.clear();
    VP.clear();
    V.clear();
    pshow=false;
    update();
}

void GLView::SetPoints(std::vector<Point>p){
    P=p;
    UpdateBoundaries();
    update();
}

void GLView::SetInnerPoints(std::vector<Point>pin){
    PIN=pin;
    UpdateBoundaries();
    update();
}

void GLView::SetVecs(std::vector<Point>p, std::vector<Point> v){
    VP=p;
    V=v;
    update();
}

void GLView::SetGridShow(bool b){
    gridshow=b;
    update();
}

void GLView::Restore(){
    UpdateBoundaries();
    pshow=false;
    update();
}

void GLView::DrawGrid(){


    if(gridshow){
        for(int i=0;i<=xsteps;i++){
            float curx=floor(left)+i*xstepper;
            glBegin(GL_LINES);
                glLineWidth(1);
                glColor3f(0.25,0.25,0.25);
                glVertex2f(curx,top);
                glVertex2f(curx,bottom);
            glEnd();
        }
        for(int i=0;i<=ysteps;i++){
            float cury=floor(bottom)+i*ystepper;
            glBegin(GL_LINES);
                glLineWidth(1);
                glColor3f(0.25,0.25,0.25);
                glVertex2f(left,cury);
                glVertex2f(right,cury);
            glEnd();
        }
    }
    QPainter painter(this);
    QPen qpen=QPen(Qt::PenStyle::SolidLine);
    qpen.setColor(QColor(255,255,255));
    painter.setPen(qpen);
    painter.setFont(QFont("Arial",8,5));
    Point p;
    int diff=10;
    if(!gridshow){
        diff=0;
    }
    for(int i=0;i<=xsteps;i++){
        p=FromSpace2Screen(Point(floor(left)+i*xstepper,top));
        painter.drawText((int)(p.x+diff),(int) (p.y+20), QString::number(floor(left)+i*xstepper));
        p=FromSpace2Screen(Point(floor(left)+i*xstepper,bottom));
        painter.drawText((int)(p.x+diff),(int) (p.y-20), QString::number(floor(left)+i*xstepper));
    }
    for(int i=0;i<=ysteps;i++){
        p=FromSpace2Screen(Point(left,floor(bottom)+i*ystepper));
        painter.drawText((int)(p.x+20),(int) (p.y-diff), QString::number(floor(bottom)+i*ystepper));
        p=FromSpace2Screen(Point(right,floor(bottom)+i*ystepper));
        painter.drawText((int)(p.x-20),(int) (p.y-diff), QString::number(floor(bottom)+i*ystepper));
    }
    painter.end();
}


void GLView::draw()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(left,right,bottom,top,1.0,15.0);
    glMatrixMode(GL_MODELVIEW);

    DrawGrid();

    glPointSize(8);
    glBegin(GL_POINTS);
        glColor3f(1.0,0.0,0.0);
        for(int i=0;i<P.size();i++){
            glVertex2f(P[i].x,P[i].y);
        }
    glEnd();
    if(P.size()>0){
        glLineWidth(1);
        glBegin(GL_LINES);
            glColor3f(1.0,0.0,0.0);
            glVertex2f(P[0].x,P[0].y);
            for(int i=1;i<P.size();i++){
                glVertex2f(P[i].x,P[i].y);
                glVertex2f(P[i].x,P[i].y);
            }
            glVertex2f(P[0].x,P[0].y);
        glEnd();
    }
    glPointSize(8);
    glBegin(GL_POINTS);
        glColor3f(0.0,1.0,0.0);
        for(int i=0;i<PIN.size();i++){
            glVertex2f(PIN[i].x,PIN[i].y);
        }
    glEnd();
    if(pshow){//chekarw gt den to deixnei
        QPainter painter(this);
        QPen qpen=QPen(Qt::PenStyle::SolidLine);
        if(bpoint2show){
            qpen.setColor(QColor(255,0,0));
        }
        else{
            qpen.setColor(QColor(0,255,0));
        }
        painter.setPen(qpen);
        painter.setFont(QFont("Arial",12,5));
        Point p=FromSpace2Screen(point2show);
        painter.drawText((int)(p.x+10),(int) (p.y-6), QString::number(point2show.x)+" , "+QString::number(point2show.y));
        painter.end();
    }
}

Point GLView::FromScreen2Space(Point p){
    float x,y;
    x=(right-left)*p.x/Width+left;
    y=(top-bottom)*(1-p.y/Height)+bottom;
    return Point(x,y);
}

Point GLView::FromSpace2Screen(Point p){
    float x,y;
    x=(p.x-left)*Width/(right-left);
    y=((top-bottom)-(p.y-bottom))*Height/(top-bottom);
    return Point(x,y);
}
