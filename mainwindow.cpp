#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QAction>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    this->setWindowState(this->windowState()|Qt::WindowMaximized);
    ui->setupUi(this);
    problem=new Problem();
    problem->AddBoundaryTable(ui->BoundaryTable);
    problem->AddInnerTable(ui->InnerTable);
    connect(ui->BoundaryPointInput,SIGNAL(PointAdd(float,float)),this,SLOT(on_boundary_point(float,float)));
    connect(ui->InnerPointInput,SIGNAL(PointAdd(float,float)),this,SLOT(on_inner_point(float,float)));
    connect(problem,SIGNAL(UpdateBoundaryPoints(std::vector<Point>)),ui->GL,SLOT(SetPoints(std::vector<Point>)));
    connect(problem,SIGNAL(UpdateInnerPoints(std::vector<Point>)),ui->GL,SLOT(SetInnerPoints(std::vector<Point>)));
    connect(problem,SIGNAL(Error(QString)),this,SLOT(on_Problem_Error(QString)));
    connect(ui->GL,SIGNAL(PointClicked(float,float)),this,SLOT(on_GL_clicked(float,float)));
    ui->BoundaryTable->setContextMenuPolicy(Qt::CustomContextMenu);
    ui->InnerTable->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(ui->BoundaryTable,SIGNAL(customContextMenuRequested(QPoint)),this,SLOT(on_Context_menu(QPoint)));
    connect(ui->InnerTable,SIGNAL(customContextMenuRequested(QPoint)),this,SLOT(on_Context_menu(QPoint)));
    connect(ui->BoundaryTable,SIGNAL(cellChanged(int,int)),this,SLOT(on_Cell_Changed(int,int)));
    connect(ui->InnerTable,SIGNAL(cellChanged(int,int)),this,SLOT(on_Cell_Changed(int,int)));
    QAction *act=ui->menuBar->addAction("");
    act->setShortcut(QKeySequence("Ctrl+O"));
    connect(act,SIGNAL(triggered(bool)),ui->actionLoad,SLOT(trigger()));
    ui->pushButton_3->setShortcut(QKeySequence("Ctrl+R"));
    ui->pushButton->setShortcut(QKeySequence("Space"));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    problem->Solve();
}


void MainWindow::on_actionLoad_triggered()
{
    QString filename=QFileDialog::getOpenFileName(this,tr("Load Data"), "/home", tr("*.txt"));
    if(filename!=NULL){
        problem->LoadFromFile(filename);
        loaded=true;
    }
}

void MainWindow::on_actionSave_triggered()
{
    if(loaded){
        QString filename=QFileDialog::getSaveFileName(this,tr("Save Data"),"/home",tr("*txt"));
        if(!filename.endsWith(".txt")){
            filename+=".txt";
        }
        if(filename!=NULL)
            problem->SaveToFile(filename);
    }
}

void MainWindow::on_actionClose_triggered()
{
    ui->GL->Clear();
    problem->Clear();
    loaded=false;
}

void MainWindow::on_boundary_point(float x,float y){
    problem->AddBoundaryPoint(x,y);
    loaded=true;
}

void MainWindow::on_inner_point(float x,float y){
    problem->AddInnerPoint(x,y);
    loaded=true;
}

void MainWindow::on_GL_clicked(float x,float y){
    ui->statusBar->showMessage(QString::number(x)+" , "+QString::number(y));
}

void MainWindow::on_pushButton_2_clicked()
{
    ui->GL->Clear();
    problem->Clear();
    loaded=false;
}

void MainWindow::on_Context_menu(QPoint p){
    QMenu contextMenu(tr("Context menu"), this);

       QAction action1("Remove Data Point", this);
       connect(&action1, SIGNAL(triggered()), this, SLOT(removeDataPoint()));
       contextMenu.addAction(&action1);
       contextMenu.exec(ui->BoundaryTable->mapToGlobal(p));
}

void MainWindow::removeDataPoint(){
    if(ui->tabWidget->currentIndex()==0){
        QList<QTableWidgetSelectionRange> list= ui->BoundaryTable->selectedRanges();

        foreach(QTableWidgetSelectionRange range,list){
            for(int i=0;i<range.rowCount();i++){
                problem->DeleteBoundaryPoint(range.topRow()+i);
            }
        }
    }
    else{
        QList<QTableWidgetSelectionRange> list= ui->InnerTable->selectedRanges();

        foreach(QTableWidgetSelectionRange range,list){
            for(int i=0;i<range.rowCount();i++){
                problem->DeleteInnerPoint(range.topRow()+i);
            }
        }
    }
    if(ui->BoundaryTable->rowCount()==0&&ui->InnerTable->rowCount()==0){
        loaded=false;
    }
}

void MainWindow::on_Cell_Changed(int row,int col){
    if(ui->tabWidget->currentIndex()==0){
        problem->UpdateBoundaryPoint(row,ui->BoundaryTable->item(row,0)->text().toFloat(),
                                     ui->BoundaryTable->item(row,1)->text().toFloat(),
                                     ui->BoundaryTable->item(row,2)->checkState()==Qt::Checked?0:1,
                                     ui->BoundaryTable->item(row,3)->text().toInt(),
                                     ui->BoundaryTable->item(row,4)->text().toFloat());
    }
    else{
        problem->UpdateInnerPoint(row,ui->InnerTable->item(row,0)->text().toFloat(),
                                  ui->InnerTable->item(row,1)->text().toFloat());
    }
}

void MainWindow::on_checkBox_clicked(bool checked)
{
    ui->GL->SetGridShow(checked);
}

void MainWindow::on_pushButton_3_clicked()
{
    ui->GL->Restore();
}

void MainWindow::on_Problem_Error(QString error){
    QMessageBox msgbox(this);
    msgbox.setText(error);
    msgbox.setWindowTitle("Error");
    msgbox.setStandardButtons(QMessageBox::Ok);
    msgbox.setFont(QFont("Arial",14));
    msgbox.exec();
}
