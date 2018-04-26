#include "pointinput.h"
#include "ui_pointinput.h"

PointInput::PointInput(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PointInput)
{
    ui->setupUi(this);
}

PointInput::~PointInput()
{
    delete ui;
}


void PointInput::on_ADD_clicked()
{
    emit PointAdd(ui->X->value(),ui->Y->value());
}
