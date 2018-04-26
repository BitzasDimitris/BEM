/********************************************************************************
** Form generated from reading UI file 'pointinput.ui'
**
** Created by: Qt User Interface Compiler version 5.9.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_POINTINPUT_H
#define UI_POINTINPUT_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_PointInput
{
public:
    QHBoxLayout *horizontalLayout_2;
    QLabel *label;
    QDoubleSpinBox *X;
    QLabel *label_2;
    QDoubleSpinBox *Y;
    QPushButton *ADD;

    void setupUi(QWidget *PointInput)
    {
        if (PointInput->objectName().isEmpty())
            PointInput->setObjectName(QStringLiteral("PointInput"));
        PointInput->resize(400, 43);
        horizontalLayout_2 = new QHBoxLayout(PointInput);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label = new QLabel(PointInput);
        label->setObjectName(QStringLiteral("label"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy);

        horizontalLayout_2->addWidget(label);

        X = new QDoubleSpinBox(PointInput);
        X->setObjectName(QStringLiteral("X"));
        X->setDecimals(3);
        X->setMinimum(-1000);
        X->setMaximum(1000);
        X->setSingleStep(0.25);

        horizontalLayout_2->addWidget(X);

        label_2 = new QLabel(PointInput);
        label_2->setObjectName(QStringLiteral("label_2"));
        sizePolicy.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy);

        horizontalLayout_2->addWidget(label_2);

        Y = new QDoubleSpinBox(PointInput);
        Y->setObjectName(QStringLiteral("Y"));
        Y->setDecimals(3);
        Y->setMinimum(-1000);
        Y->setMaximum(1000);
        Y->setSingleStep(0.25);

        horizontalLayout_2->addWidget(Y);

        ADD = new QPushButton(PointInput);
        ADD->setObjectName(QStringLiteral("ADD"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(ADD->sizePolicy().hasHeightForWidth());
        ADD->setSizePolicy(sizePolicy1);

        horizontalLayout_2->addWidget(ADD);


        retranslateUi(PointInput);

        QMetaObject::connectSlotsByName(PointInput);
    } // setupUi

    void retranslateUi(QWidget *PointInput)
    {
        PointInput->setWindowTitle(QApplication::translate("PointInput", "Form", Q_NULLPTR));
        label->setText(QApplication::translate("PointInput", "X:", Q_NULLPTR));
        label_2->setText(QApplication::translate("PointInput", "Y:", Q_NULLPTR));
        ADD->setText(QApplication::translate("PointInput", "ADD", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class PointInput: public Ui_PointInput {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_POINTINPUT_H
