#ifndef MAINWINDOW_H
#define MAINWINDOW_H


#include <QMainWindow>
#include <QFileDialog>
#include <QMessageBox>
#include "problem.h"
#include "glview.h"
#include "point.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

    void on_actionLoad_triggered();

    void on_actionSave_triggered();

    void on_actionClose_triggered();

    void on_boundary_point(float x,float y);

    void on_inner_point(float x,float y);

    void on_GL_clicked(float x, float y);

    void on_pushButton_2_clicked();

    void on_Context_menu(QPoint);

    void removeDataPoint();

    void on_checkBox_clicked(bool checked);

    void on_pushButton_3_clicked();

    void on_Cell_Changed(int ,int);

    void on_Problem_Error(QString error);

private:
    Ui::MainWindow *ui;
    Problem *problem;
    QAbstractItemModel *model;
    bool loaded=false;
};

#endif // MAINWINDOW_H
