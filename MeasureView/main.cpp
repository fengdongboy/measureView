#include "MeasuredataView.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
	MeasureDataView w;
    w.show();
    return a.exec();
}
