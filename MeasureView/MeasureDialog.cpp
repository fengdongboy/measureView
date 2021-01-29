#include "MeasureDialog.h"

MeasureDialog::MeasureDialog(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);
}

MeasureDialog::~MeasureDialog()
{
}

void MeasureDialog::sltCapturePoints(const QVector3D& point)
{
	QVector3D p;
	if (!mLastPointList.isEmpty())
		p = mLastPointList.last();
	mLastPointList.push_back(point);
	QString str = QString("%1 %2 %3").arg(point.x()).arg(point.y()).arg(point.z());
	if (mLastPointList.size()%2 == 1)
	{
		ui.lineEdit_point1->setText(str);
	}
	else
	{
		QString str1 = QString("%1 %2 %3").arg(p.x()).arg(p.y()).arg(p.z());
		ui.lineEdit_point2->setText(str1);

		float d = p.distanceToPoint(point);
		ui.lineEdit_value->setText(QString::number(d));

		ManualMeasurePoint ps;
		ps.point1 = p;
		ps.point2 = point;
		ps.value = d;
		mPointList.push_back(ps);
		ui.listWidget->addItem(QString::number(mPointList.size() - 1));
	}
}

void MeasureDialog::on_listWidget_itemClicked(QListWidgetItem *item)
{
	int n = ui.listWidget->row(item);
	if (n < mPointList.size())
	{
		QString str1 = QString("%1 %2 %3").arg(mPointList[n].point1.x()).arg(mPointList[n].point1.y()).arg(mPointList[n].point1.z());
		ui.lineEdit_point1->setText(str1);
		QString str2 = QString("%1 %2 %3").arg(mPointList[n].point2.x()).arg(mPointList[n].point2.y()).arg(mPointList[n].point2.z());
		ui.lineEdit_point2->setText(str2);
		ui.lineEdit_value->setText(QString::number(mPointList[n].value));
	}
}

void MeasureDialog::closeEvent(QCloseEvent *event)
{
	emit sigVisible(false);
}

void MeasureDialog::showEvent(QShowEvent *event)
{
	sigVisible(true);
}
