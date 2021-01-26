#include "MeasureDataView.h"

#include "ObjReadWrite.h"

#include <QFileDialog>
#include <fstream>
#include <QTableWidgetItem>

typedef float(*EsMeasure_Fun)(std::string filepath, bool gender, std::vector<EsLineData>  & data, float & score, bool ttt);

void readWriteEslineData(std::vector<EsLineData>&data, QDataStream& stream, float& score, bool read = true)
{
	int size = data.size();
	if (read)
	{
		stream >> size;
		data.resize(size);
		stream >> score;
	}
	else
	{
		stream << size;
		stream << score;
	}


	auto readLineData = [=](EsLineData& d, QDataStream& s)
	{
		int size = 0;
		s >> size;
		char buf[256];
		s.readRawData(buf, size);
		d.name = buf;
		s >> d.value;
		s >> d.CF;
		s >> size;
		d.Points.resize(size);
		for (int i = 0; i < size; i++)
		{
			s >> d.Points[i][0] >> d.Points[i][1] >> d.Points[i][2];
		}
		s >> size;
		d.Tris.resize(size);
		for (int i = 0; i < size; i++)
		{
			s >> d.Tris[i][0] >> d.Tris[i][1] >> d.Tris[i][2];
		}
	};

	auto writeLineData = [=](EsLineData& d, QDataStream& s)
	{
		int size = d.name.size() + 1;
		s << size;
		s.writeRawData(d.name.c_str(), size);
		s << d.getValue();
		s << d.CF;
		size = d.Points.size();
		s << size;
		for (int i = 0; i < size; i++)
		{
			s << d.Points[i][0] << d.Points[i][1] << d.Points[i][2];
		}
		size = d.Tris.size();
		s << size;
		for (int i = 0; i < size; i++)
		{
			s << d.Tris[i][0] << d.Tris[i][1] << d.Tris[i][2];
		}
	};

	for (int i = 0; i < size; i++)
	{
		if (read)
			readLineData(data[i], stream);
		else
			writeLineData(data[i], stream);
	}
}

bool writeEslineDatas(std::vector<EsLineData>& data, float& score, const QString& name)
{
	QFile f(name);
	if (f.open(QFile::WriteOnly))
	{
		QDataStream stream(&f);
		readWriteEslineData(data, stream, score, false);
		return true;
	}
	return false;
}

bool readEslineDatas(std::vector<EsLineData>& data, float& score, const QString& name)
{
	QFile f(name);
	if (f.open(QFile::ReadOnly))
	{
		QDataStream stream(&f);
		readWriteEslineData(data, stream, score, true);
		return true;
	}
	return false;
}

void readWriteEslineData(std::vector<EsLineData>&data, QTextStream& stream, float& score, bool read = true)
{
	int size = data.size();
	if (read)
	{
		stream >> size;
		data.resize(size);
		stream >> score;
	}
	else
	{
		stream << size << " ";
		stream << score << " ";
	}


	auto readLineData = [=](EsLineData& d, QTextStream& s)
	{
		int size = 0;
		s >> size;
		QString str;
		s >> str;
		//str = s.read(size);
		d.name = str.toLocal8Bit().data();
		s >> d.value;
		int b;
		s >> b;
		d.CF = b;
		//s >> (char)d.CF;
		s >> size;
		d.Points.resize(size);
		for (int i = 0; i < size; i++)
		{
			s >> d.Points[i][0] >> d.Points[i][1] >> d.Points[i][2];
		}
		s >> size;
		d.Tris.resize(size);
		for (int i = 0; i < size; i++)
		{
			s >> d.Tris[i][0] >> d.Tris[i][1] >> d.Tris[i][2];
		}
	};

	auto writeLineData = [=](EsLineData& d, QTextStream& s)
	{
		int size = d.name.size() + 1;
		QString str = QString::fromLocal8Bit(d.name.c_str());
		size = str.size() + 1;
		s << size << " ";
		//s.writeRawData(d.name.c_str(), size);
		
		s << str << " ";
		s << d.getValue() << " ";
		int b = d.CF;
		s << b << " ";
		//s << (char)d.CF;
		size = d.Points.size();
		s << size << " ";
		for (int i = 0; i < size; i++)
		{
			s << d.Points[i][0] <<" "<< d.Points[i][1]<<" " << d.Points[i][2]<<" ";
		}
		size = d.Tris.size();
		s << size << " ";
		for (int i = 0; i < size; i++)
		{
			s << d.Tris[i][0] << " "<< d.Tris[i][1]<<" " << d.Tris[i][2]<< " ";
		}
	};

	for (int i = 0; i < size; i++)
	{
		if (read)
			readLineData(data[i], stream);
		else
			writeLineData(data[i], stream);
	}
}

bool writeEslineDatasTxt(std::vector<EsLineData>& data, float& score, const QString& name)
{
	QFile f(name);
	if (f.open(QFile::WriteOnly))
	{
		QTextStream stream(&f);
		//stream.setCodec("gbk");
		readWriteEslineData(data, stream, score, false);
		return true;
	}
	return false;
}

bool readEslineDatasTxt(std::vector<EsLineData>& data, float& score, const QString& name)
{
	QFile f(name);
	if (f.open(QFile::ReadOnly))
	{
		QTextStream stream(&f);
		//stream.setCodec("gbk");
		readWriteEslineData(data, stream, score, true);
		return true;
	}
	return false;
}

MeasureDataView::MeasureDataView(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	QVector3D matral;
	matral.setX(242 / 255.0f);
	matral.setY(242 / 255.0f);
	matral.setZ(242 / 255.0f);
	ui.widget_scene->setColor(matral, matral, QVector3D());
}

MeasureDataView::~MeasureDataView()
{
}

void lineData2Ptcloud(const EsLineData&lines, PtCloud& cloud)
{
	//for (int i = 0; i < lines.Points.size(); i++)
	//{
	//	cloud.push_back(Vertex3Normal(lines.Points[i][0], lines.Points[i][1], lines.Points[i][2], 0, 0, 1));
	//}

	for (int i = 0; i < lines.Tris.size(); i ++)
	{
		const omesh::Pnt3& pnt1 = lines.Points[lines.Tris[i][0]];
		const omesh::Pnt3& pnt2 = lines.Points[lines.Tris[i][1]];
		const omesh::Pnt3& pnt3 = lines.Points[lines.Tris[i][2]];

		/// 一个三角形画三条线
		cloud.push_back(Vertex3Normal(pnt1[0], pnt1[1], pnt1[2], 0, 0, 1));
		cloud.push_back(Vertex3Normal(pnt2[0], pnt2[1], pnt2[2], 0, 0, 1));

		cloud.push_back(Vertex3Normal(pnt1[0], pnt1[1], pnt1[2], 0, 0, 1));
		cloud.push_back(Vertex3Normal(pnt3[0], pnt3[1], pnt3[2], 0, 0, 1));

		cloud.push_back(Vertex3Normal(pnt2[0], pnt2[1], pnt2[2], 0, 0, 1));
		cloud.push_back(Vertex3Normal(pnt3[0], pnt3[1], pnt3[2], 0, 0, 1));
	}
}

void MeasureDataView::on_action_open_triggered(void)
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
		"",
		tr("datas (*.data);; txt (*.txt )"));
	if (!fileName.isEmpty())
	{
		float score = 0.0;
		bool ret = false;
		if (fileName.right(3) == "txt")
			ret = readEslineDatasTxt(mLineData, score, fileName);
		else
			ret = readEslineDatas(mLineData, score, fileName);
		if (ret)
		{
			ui.widget_scene->clear();
			PtCloud cloud;
			int rowCount = 0;
			ui.listWidget_measureDatas->setRowCount(0);
			for (int i = 0; i < mLineData.size(); i++)
			{ 
				rowCount = ui.listWidget_measureDatas->rowCount();
				ui.listWidget_measureDatas->insertRow(rowCount);//增加一行

																//插入元素
				ui.listWidget_measureDatas->setItem(rowCount, 0, new QTableWidgetItem(QString::number(i)));
				ui.listWidget_measureDatas->setItem(rowCount, 1, new QTableWidgetItem(QString::fromLocal8Bit(mLineData[i].name.c_str())));
				ui.listWidget_measureDatas->setItem(rowCount, 2, new QTableWidgetItem(QString::number(mLineData[i].value / 10, 'f', 2)));

				if (i == 0)
				{
					ui.listWidget_measureDatas->setItem(rowCount, 2, new QTableWidgetItem(QString::number(mLineData[i].value / 10.0f, 'f', 2) + "kg"));
				}

				cloud.clear();
				lineData2Ptcloud(mLineData[i], cloud);

				//if (!cloud.empty())
				ui.widget_scene->addLineData(&cloud);
			}
		}
	}
}

void MeasureDataView::on_action_export_triggered(void)
{
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
		"",
		tr("datas (*.data);; txt(*.txt )"));
	if (!fileName.isEmpty())
	{
		bool ret = false;
		float score = 0.0f;
		if (fileName.right(3) == "txt")
		{
			writeEslineDatasTxt(mLineData, score, fileName);
		}
		else
		{
			writeEslineDatas(mLineData, score, fileName);
		}		
	}
}

void MeasureDataView::on_action_openModel_triggered(void)
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
		"",
		tr("datas (*.obj);; txt (*.stl)"));

	PtCloud cloud;
	if (!fileName.isEmpty())
	{
		QString str = fileName.right(3);
		str = str.toLower();
		if (str == "obj")
		{
			//PtCloud cloud = readObjModel(fileName);
			//ui.widget_scene->setPtCloud(&cloud, NULL);
			bool ret = readObjModel(fileName, mScanData);
			if (ret)
			{
				if (mScanData.modelData.hasImg)				
					cloud = mScanData.modelData.getTexcoordCloud();				
				else				
					cloud = mScanData.modelData.getCloud();				
			}
			ui.widget_scene->setModel(&cloud, mScanData.modelData.img, mScanData.modelData.hasImg);
		}
		else if (str == "stl")
		{
			//PtCloud cloud = readPoissonStl(fileName);
			//ui.widget_scene->setPtCloud(&cloud, NULL);

			bool ret = readStlModel(fileName, mScanData);
			if (ret)
			{
				cloud = mScanData.modelData.getCloud();
				ui.widget_scene->setModel(&cloud, mScanData.modelData.img, false);
			}
		}
	}
}

void MeasureDataView::on_listWidget_measureDatas_itemSelectionChanged(void)
{
	QList<QTableWidgetItem *> list = ui.listWidget_measureDatas->selectedItems();
	std::vector<int> indexs;
	foreach(QTableWidgetItem* item, list)
	{
		//if (item->column() == 1)
		indexs.push_back(item->row());
	}
	std::sort(indexs.begin(), indexs.end());
	indexs.erase(unique(indexs.begin(), indexs.end()), indexs.end());
	if (!indexs.empty() && indexs[0] == 0)	/// 选中体重，选中全部
	{
		indexs.clear();
		for (int i = 0; i < ui.listWidget_measureDatas->rowCount(); i++)
		{
			indexs.push_back(i);
		}
	}
	ui.widget_scene->setSelectLines(indexs);
}

void MeasureDataView::on_action_openPorj_triggered(void)
{
#if 0
	static QString last;
	QString dir = QFileDialog::getExistingDirectory(this, QStringLiteral("请选择要打开的文件夹"),
		last,
		QFileDialog::ShowDirsOnly
		| QFileDialog::DontResolveSymlinks);
	if (dir.isEmpty())
		return;
	last = dir;
	QStringList list = getFileList(dir);
	float score = 0;
	for (int i = 0; i < list.size(); i++)
	{
		if (modelFileType[i] == en_stl)
		{
			readStlModel(list[i], mScanData[i]);			
		}
		else 
		{
			readObjModel(list[i], mScanData[i]);
		}
		readEslineDatas(mScanData[i].mLineData, score, list[i] + ".data");
	}
#else
#endif
}

void MeasureDataView::on_action_exportModel_triggered(void)
{

}

bool MeasureDataView::readSTLFile(const char *filename, std::vector<omesh::Pnt3>& vtx, std::vector<omesh::Pnt3>& normals, std::vector<omesh::TriVtx>& tris, int format /*= 0*/)
{
	if (format == 0) {
		ifstream ifs(filename, ios::binary);
		if (!ifs) { ifs.close(); return false; }

		vtx.clear();
		normals.clear();
		tris.clear();

		int intSize = sizeof(int);
		int floatSize = sizeof(float);
		ifs.ignore(80);

		// 面的个数
		int num_tris;
		ifs.read((char*)(&num_tris), intSize);

		float tn0, tn1, tn2;
		float v0, v1, v2;

		omesh::Pnt3 nortemp;

		for (int i = 0; i < num_tris; i++) {
			ifs.read((char*)(&(tn0)), floatSize);
			ifs.read((char*)(&(tn1)), floatSize);
			ifs.read((char*)(&(tn2)), floatSize);

			ifs.read((char*)(&(v0)), floatSize);
			ifs.read((char*)(&(v1)), floatSize);
			ifs.read((char*)(&(v2)), floatSize);
			vtx.push_back(omesh::Pnt3(v0, v1, v2));
			ifs.read((char*)(&(v0)), floatSize);
			ifs.read((char*)(&(v1)), floatSize);
			ifs.read((char*)(&(v2)), floatSize);
			vtx.push_back(omesh::Pnt3(v0, v1, v2));
			ifs.read((char*)(&(v0)), floatSize);
			ifs.read((char*)(&(v1)), floatSize);
			ifs.read((char*)(&(v2)), floatSize);
			vtx.push_back(omesh::Pnt3(v0, v1, v2));

			nortemp = omesh::normal(vtx[i * 3], vtx[i * 3 + 1], vtx[i * 3 + 2]);
			nortemp = nortemp.normalize();
			//nortemp = omesh::normal(tn0, tn1, tn2);
			normals.push_back(nortemp);
			normals.push_back(nortemp);
			normals.push_back(nortemp);

			tris.push_back(omesh::TriVtx(i * 3, i * 3 + 1, i * 3 + 2));

			ifs.ignore(2);
		}

		ifs.close();
		return true;

	}

	return false;
}

PtCloud& MeasureDataView::readPoissonStl(const QString& name)
{
	std::vector<omesh::Pnt3> vtx;
	std::vector<omesh::Pnt3> normals;
	std::vector<omesh::TriVtx> tris;
	bool ret = readSTLFile(name.toLocal8Bit().data(), vtx, normals, tris);
	static PtCloud cloud;
	cloud.clear();
	if (ret)
	{
		Vertex3Normal point;
		for (int i = 0; i < vtx.size(); i++)
		{
			point.x = vtx[i][0];
			point.y = vtx[i][1];
			point.z = vtx[i][2];

			point.nx = normals[i][0];
			point.ny = normals[i][1];
			point.nz = normals[i][2];
			cloud.push_back(point);
		}
	}
	return cloud;
}

PtCloud& MeasureDataView::readObjModel(const QString& name)
{
	static PtCloud mModelCloud;
	static ObjReadWrite w;
	w.clear();
	bool ret = w.readFile(name.toLocal8Bit().data());
	if (ret)
	{
		mModelCloud.clear();
		w.mNormals.resize(w.mPoints.size());

		auto fGetPoint = [=](int index)
		{
			return omesh::Pnt3(w.mPoints[index][0], w.mPoints[index][1], w.mPoints[index][2]);
		};

		omesh::Pnt3 p1, p2, p3, n;

		for (int i = 0; i < w.mTriangles.size(); i++)
		{
			p1 = fGetPoint(w.mTriangles[i][0] - 1);
			p2 = fGetPoint(w.mTriangles[i][1] - 1);
			p3 = fGetPoint(w.mTriangles[i][2] - 1);
			n = normal(p1, p2, p3);

			for (int j = 0; j < 3; j++)
			{
				int index = w.mTriangles[i][j] - 1;
				mModelCloud.push_back(Vertex3Normal(w.mPoints[index][0], w.mPoints[index][1], w.mPoints[index][2],
					n[0], n[1], n[2]));
			}
		}
	}
	return mModelCloud;
}

bool MeasureDataView::readObjModel(const QString& name, ScanData& model)
{
	static ObjReadWrite w;
	w.clear();
	bool ret = w.readFile(name.toLocal8Bit().data());
	if (ret)
	{
		model.modelData.trisVtx.resize(w.mTriangles.size());
		model.modelData.points.resize(w.mPoints.size());
		model.modelData.normals.resize(w.mNormals.size());
		model.modelData.texcoord.resize(w.mTextCoords.size());
		for (int i = 0; i < w.mTriangles.size(); i++)
		{
			for (int j = 0; j < 3; j++)
				model.modelData.trisVtx[i][j] = w.mTriangles[i][j] - 1;
		}

		auto conver = [](auto& from, auto& to)
		{
			for (int i = 0; i < from.size(); i++)
			{
				for (int j = 0; j < 3; j++)
					to[i][j] = from[i][j];
			}
		};
		conver(w.mPoints, model.modelData.points);
		conver(w.mNormals, model.modelData.normals);
		conver(w.mTextCoords, model.modelData.texcoord);
		model.modelData.mTextCoordIndexs = w.mTextCoordIndexs;

		/// 读取图片
		bool r = readModelTex(name, model.modelData.img);
		model.modelData.hasImg = r;
	}
	model.load = ret;
	return ret;
}

bool MeasureDataView::readModelTex(const QString& str, QImage& img)
{
	if (str.indexOf(".obj") == -1)
	{
		return false;
	}

	{
		QString fstr = str;
		fstr = fstr.replace(".obj", ".mtl");

		QFile f(fstr);
		QStringList list;
		if (f.open(QFile::ReadOnly | QFile::Text))
		{
			QTextStream stream(&f);
			QString strtemp;
			while (!stream.atEnd())
			{
				strtemp = stream.readLine();
				if (strtemp.indexOf("map_Kd") != -1)
				{
					list = strtemp.split(" ");
					break;
				}
			}
			f.close();
		}

		if (!list.isEmpty())
		{
			fstr = list[1];
			QFileInfo info(str);
			QString path = info.absolutePath();
			QString temp = path + "/" + fstr;

			bool b = img.load(temp);
			img = img.convertToFormat(QImage::Format_RGB888);
			img = img.mirrored();

			return true;
		}
	}
	return false;
}

QStringList MeasureDataView::getFileList(const QString& dir)
{
	QStringList list;
	QString str;
	QStringList l;
	l = dir.split("/");
	for (int i = 0; i < en_modelScanNull; i ++)
	{
		str = dir + "/" + l.last();
		str += EnScanTypeStr[i];
		str += EnScanTypeExtStr[i];
		list << str;
	}
	return list;
}

bool MeasureDataView::readStlModel(const QString& name, ScanData& model)
{
	bool ret = readSTLFile(name.toLocal8Bit().data(), model.modelData.points, model.modelData.normals, model.modelData.trisVtx);
	model.load = ret;
	return ret;
}
