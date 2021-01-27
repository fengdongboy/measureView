#pragma once
#ifndef __FILEEXPORTMANAGER_H__
#define __FILEEXPORTMANAGER_H__
/*!
 * \class FileExportManager
 *
 * \brief 实现文件的导出，支持 *.asc、*.stl、*.xyz、*.obj
 *
 * \author FDL
 * \date 一月 2021
 */

#include "types.h"
class FileExportManager
{
public:
	FileExportManager();
	~FileExportManager();

	static FileExportManager* getFileExportManager(void);

	bool savePoints(const QString& name, const MeasureModel& modelData);

	bool saveAsc(const QString& name, const MeasureModel& modelData);
	bool saveStl(const QString& name, const MeasureModel& modelData);
	bool saveXyz(const QString& name, const MeasureModel& modelData);
	bool saveObj(const QString& name, const MeasureModel& modelData);
	bool saveModel(const QString& name, const MeasureModel& modelData);

	// format 0 binary 1 ascii
	bool saveSTLBinaryFile(const char *filename, std::vector<omesh::Pnt3>& vtx, 
		std::vector<omesh::Pnt3>& nrm, std::vector<omesh::TriVtx>& tris, int format/*=0*/);

	bool saveStlAssiiFile(const char *filename, std::vector<omesh::Pnt3>& vtx, 
		std::vector<omesh::Pnt3>& nrm, std::vector<omesh::TriVtx>& tris);

};

#endif