#ifndef PROBLEM_H
#define PROBLEM_H
#define PI 3.14159265358979323846
#include <vector>
#include <QTableWidget>
#include <QDebug>
#include "point.h"
#include "math.h"
#include "node.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lu_decomposition.hpp"

class Problem : public QObject
{
    Q_OBJECT
public:
    explicit Problem(QWidget *parent=0);
    ~Problem();
    void LoadFromFile(QString Filename);
    void SaveToFile(QString Filename);
    void AddBoundaryTable(QTableWidget *table);
    void AddInnerTable(QTableWidget *table);
    void DeleteBoundaryPoint(int index);
    void DeleteInnerPoint(int index);
    void AddBoundaryPoint(float x,float y);
    void AddInnerPoint(float x,float y);
    void UpdateBoundaryPoint(int index, float x,float y,int Known,int type,float ub);
    void UpdateInnerPoint(int index,float x, float y);
    void Solve();
    void Clear();
private:
    void gmatr();
    void hmatr();
    void fmatr();
    float anode(float U);
    float cathode(float U);
    float anodeDF(float U);
    float cathodeDF(float U);
    bool iterationCheck();
    void abmatr();
    void solveq();
    void solveq_lu();
    bool legs();
    bool ludcmp();
    void reorder();
    void uinter();
    void deriv();
    void UpdateResults();
    float rlintc(Point p1,Point p2,Point p3);
    float slintc(Point p1,Point p2);
    float dalpha(Point p1,Point p2,Point p3);

    int N=0,InnerN=0;
    float initialPhiAnode=0.1f;
    float initialPhiCathode=-0.1f;
    float Ia=1;
    float Ic=1;
    float Phia=-0.5f;
    float Phic=0.5f;
    float Ba=0.05f;
    float Bc=0.05f;
    float Aa=0.05f;
    float Ac=0.05f;
    float C=10.0f;
    float e=0.0001;
    std::vector<int> INDEX;
    std::vector<int> TYPE;
    std::vector<node> nodes;
    std::vector<float> UB,B,UIN,F,DF,X;
    std::vector<std::vector<float>> G,H,A;
    std::vector<Point> PM,PL,PIN,UPIN;
    bool singular=false,initialized=false;
    std::vector<QTableWidget *> btables;
    QTableWidget *innertable;
signals:
    void UpdateBoundaryPoints(std::vector<Point>);
    void UpdateInnerPoints(std::vector<Point>);
    void Error(QString error);
};

#endif // PROBLEM_H
