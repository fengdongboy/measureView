#pragma once
#ifndef __MEASUREDATAVIEW_H__
#define __MEASUREDATAVIEW_H__

#include "face/include/Pnt3.h"
#include "Face/include/TriMeshUtils.h"

#include <QMainWindow>
#include <QFile>
#include <QDataStream>
#include <QTextStream>

#include "ui_MeasureDataView.h"

#ifndef __EsLineData__
#define __EsLineData__
struct EsLineData
{
	std::string name;
	float value;
	bool CF;

	std::vector<omesh::Pnt3> Points;
	std::vector<omesh::TriVtx> Tris;

	EsLineData() :name(""), value(0.0f) {}
	EsLineData(const std::string& str, float lenght = 0.0f) :name(str), value(lenght)
	{
	}

	void setName(const std::string& name)
	{
		this->name = name;
	}
	const std::string getName(void) const
	{
		return name;
	}

	void setValue(float v)
	{
		value = v;
	}
	float getValue(void) const
	{
		return value;
	}
};
#endif



class MeasureDataView : public QMainWindow
{
	Q_OBJECT

public:
	MeasureDataView(QWidget *parent = Q_NULLPTR);
	~MeasureDataView();

	public Q_SLOTS:
	void on_action_open_triggered(void);
	void on_action_export_triggered(void);

	void on_action_openModel_triggered(void);

	void on_listWidget_measureDatas_itemSelectionChanged(void);

protected:
	bool readSTLFile(const char *filename, std::vector<omesh::Pnt3>& vtx, std::vector<omesh::Pnt3>& normals, std::vector<omesh::TriVtx>& tris, int format = 0);

	PtCloud& readPoissonStl(const QString& name);
	PtCloud& readObjModel(const QString& name);

private:
	Ui::MeasureDataView ui;

	std::vector<EsLineData> mLineData;
};
#endif