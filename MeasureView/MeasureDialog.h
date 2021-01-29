#pragma once

#include <QWidget>
#include <QVector3D>
#include "ui_MeasureDialog.h"

struct ManualMeasurePoint
{
	QVector3D point1;
	QVector3D point2;
	float value;
};

class MeasureDialog : public QWidget
{
	Q_OBJECT

public:
	MeasureDialog(QWidget *parent = Q_NULLPTR);
	~MeasureDialog();

	public Q_SLOTS:
	void sltCapturePoints(const QVector3D& point);

	/// µã»÷
	void on_listWidget_itemClicked(QListWidgetItem *item);

	Q_SIGNALS:
	void sigVisible(bool);

protected:
	void closeEvent(QCloseEvent *event);

	void showEvent(QShowEvent *event);


private:
	Ui::MeasureDialog ui;

	QVector<QVector3D> mLastPointList;
	QVector<ManualMeasurePoint> mPointList;
};
