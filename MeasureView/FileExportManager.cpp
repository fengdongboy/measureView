#include "FileExportManager.h"

#include "ObjReadWrite.h"

#include <QFile>
#include <QTextStream>


FileExportManager::FileExportManager()
{
}


FileExportManager::~FileExportManager()
{
}

FileExportManager* FileExportManager::getFileExportManager(void)
{
	static FileExportManager* p = new FileExportManager;
	return p;
}

bool FileExportManager::savePoints(const QString& name, const MeasureModel& modelData)
{
	QFile file(name);
	if (!file.open(QFile::WriteOnly | QFile::Text))
		return false;
	QTextStream stream(&file);
	for (int i = 0; i < modelData.points.size(); i ++)
	{
		const omesh::Pnt3& point = modelData.points[i];
		stream << point[0] << " " << point[1] << " " << point[2];
		if (!modelData.normals.empty())
		{
			stream<<" " << modelData.normals[i][0] << " " << modelData.normals[i][1] << " " << modelData.normals[i][2];
		}
		stream << endl;
	}
	return true;
}

bool FileExportManager::saveAsc(const QString& name, const MeasureModel& modelData)
{
	return savePoints(name, modelData);
}

bool FileExportManager::saveStl(const QString& name, const MeasureModel& modelData)
{
	std::vector<omesh::Pnt3> points = modelData.points;
	std::vector<omesh::Pnt3> normals = modelData.normals;
	std::vector<omesh::TriVtx> tris = modelData.trisVtx;
	if (normals.empty())
	{
		normals.resize(points.size());
		omesh::Pnt3 n;
		for (int i = 0; i < tris.size(); i ++)
		{
			n = omesh::normal(points[tris[i][0]], points[tris[i][1]], points[tris[i][2]]);
			normals[tris[i][0]] = normals[tris[i][1]] = normals[tris[i][2]] = n;
		}
	}
	
	return saveStlAssiiFile(name.toLocal8Bit().data(), points, normals, tris);
}

bool FileExportManager::saveXyz(const QString& name, const MeasureModel& modelData)
{
	return savePoints(name, modelData);
}

bool FileExportManager::saveObj(const QString& name, const MeasureModel& modelData)
{
	ObjReadWrite w;
	w.mTriangles.resize(modelData.trisVtx.size());
	w.mPoints.resize(modelData.points.size());
	w.mNormals.resize(modelData.normals.size());
	w.mTextCoords.resize(modelData.texcoord.size());
	for (int i = 0; i < w.mTriangles.size(); i++)
	{
		for (int j = 0; j < 3; j++)
			w.mTriangles[i][j] = modelData.trisVtx[i][j] + 1;
	}

	auto conver = [](auto& from, auto& to)
	{
		for (int i = 0; i < from.size(); i++)
		{
			for (int j = 0; j < 3; j++)
				to[i][j] = from[i][j];
		}
	};
	conver(modelData.points, w.mPoints);
	conver(modelData.normals, w.mNormals);
	conver(modelData.texcoord, w.mTextCoords);
	w.mTextCoordIndexs = modelData.mTextCoordIndexs;
	
	bool ret = w.writeFile(name.toLocal8Bit().data());
	if (modelData.hasImg)
	{
		QString bmg = name;
		bmg	= bmg.replace(".obj", ".bmp");
		bool r = modelData.img.save(bmg);
	}
	return ret;
}

bool FileExportManager::saveModel(const QString& name, const MeasureModel& modelData)
{
	QString ext = name.right(4);
	ext = ext.toLower();
	if (ext == ".asc")
	{
		return saveAsc(name, modelData);
	}
	else if (ext == ".stl")
	{
		return saveStl(name, modelData);
	}
	else if (ext == ".xyz")
	{
		return saveXyz(name, modelData);
	}
	else if (ext == ".obj")
	{
		return saveObj(name, modelData);
	}
	else return false;
}

bool FileExportManager::saveSTLBinaryFile(const char *filename, std::vector<omesh::Pnt3>& vtx, 
	std::vector<omesh::Pnt3>& nrm, std::vector<omesh::TriVtx>& tris, int format/*=0*/)
{
	format = 0;
	if (format == 0) {
		ofstream fs(filename, ios::binary);
		if (!fs) { fs.close(); return false; }

		int i;
		int intSize = sizeof(int);
		int floatSize = sizeof(float);
		// 文件头
		char* fileHead = "ROC";
		fs.write(fileHead, sizeof(char) * 3);
		// 附加信息
		char fileInfo[77];
		for (i = 0; i < 77; i++) fileInfo[i] = ' ';
		fs.write(fileInfo, sizeof(char) * 77);
		// 面的个数
		int num_tris = int(tris.size());
		fs.write((char*)(&num_tris), intSize);
		// 点列表，面列表
		char a[2];
		streamsize a_size = sizeof(char) * 2;
		omesh::Pnt3 tn;
		for (i = 0; i < num_tris; i++) {
			tn = nrm[tris[i][0]] + nrm[tris[i][1]] + nrm[tris[i][2]];
			fs.write((char*)(&(tn[0])), floatSize);
			fs.write((char*)(&(tn[1])), floatSize);
			fs.write((char*)(&(tn[2])), floatSize);
			fs.write((char*)(&(vtx[tris[i][0]][0])), floatSize);
			fs.write((char*)(&(vtx[tris[i][0]][1])), floatSize);
			fs.write((char*)(&(vtx[tris[i][0]][2])), floatSize);
			fs.write((char*)(&(vtx[tris[i][1]][0])), floatSize);
			fs.write((char*)(&(vtx[tris[i][1]][1])), floatSize);
			fs.write((char*)(&(vtx[tris[i][1]][2])), floatSize);
			fs.write((char*)(&(vtx[tris[i][2]][0])), floatSize);
			fs.write((char*)(&(vtx[tris[i][2]][1])), floatSize);
			fs.write((char*)(&(vtx[tris[i][2]][2])), floatSize);
			fs.write(a, a_size);
		}
		fs.close();
		return true;

	}
	if (format == 1) {
		int n, i;
		ofstream fs(filename);
		if (!fs) {
			fs.close();
			return false;
		}
		fs << "solid ods\n";
		n = int(tris.size());
		omesh::Pnt3 tn;
		for (i = 0; i < n; i++) {
			tn = nrm[tris[i][0]] + nrm[tris[i][1]] + nrm[tris[i][2]];
			fs << "facet normal " << tn[0] << " " << tn[1] << " " << tn[2] << "\nouter loop\n"
				<< "vertex " << vtx[tris[i][0]][0] << " " << vtx[tris[i][0]][1] << " " << vtx[tris[i][0]][2] << "\n"
				<< "vertex " << vtx[tris[i][1]][0] << " " << vtx[tris[i][1]][1] << " " << vtx[tris[i][1]][2] << "\n"
				<< "vertex " << vtx[tris[i][2]][0] << " " << vtx[tris[i][2]][1] << " " << vtx[tris[i][2]][2] << "\n"
				<< "endloop\nendfacet\n";
		}

		fs << "endsolid";
		fs.close();
		return true;
	}
	return false;
}

bool FileExportManager::saveStlAssiiFile(const char *filename, std::vector<omesh::Pnt3>& vtx, std::vector<omesh::Pnt3>& nrm, std::vector<omesh::TriVtx>& tris)
{
	std::ofstream* mf = new std::ofstream(filename);
	if (mf == NULL)
	{
		return false;
	}

	*mf << "solid " << "feng" << std::endl;

	omesh::Pnt3 normaltemp;
	for (int i = 0; i < tris.size(); i++)
	{
		int& i1 = tris[i][0];
		int& i2 = tris[i][1];
		int& i3 = tris[i][2];

		normaltemp = nrm[i1] + nrm[i2] + nrm[i3];

		*mf << "facet normal " << normaltemp[0] << " " << normaltemp[1] << " " << normaltemp[2] << std::endl;
		*mf << "outer loop" << std::endl;

		omesh::Pnt3& v1 = vtx[i1];
		omesh::Pnt3& v2 = vtx[i2];
		omesh::Pnt3& v3 = vtx[i3];

		*mf << "vertex " << v1[0] << " " << v1[1] << " " << v1[2] << std::endl;
		*mf << "vertex " << v2[0] << " " << v2[1] << " " << v2[2] << std::endl;
		*mf << "vertex " << v3[0] << " " << v3[1] << " " << v3[2] << std::endl;
		*mf << "endloop" << std::endl;
		*mf << "endfacet " << std::endl;
	}

	*mf << "endsolid " << std::endl;
	mf->close();
	delete mf;
	return true;
}
