#ifndef POINTINPUT_H
#define POINTINPUT_H

#include <QWidget>

namespace Ui {
class PointInput;
}

class PointInput : public QWidget
{
    Q_OBJECT

public:
    explicit PointInput(QWidget *parent = 0);
    ~PointInput();

private slots:

    void on_ADD_clicked();

signals:
    void PointAdd(float x,float y);
private:
    Ui::PointInput *ui;
};

#endif // POINTINPUT_H
