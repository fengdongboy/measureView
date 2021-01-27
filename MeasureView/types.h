#pragma once
/*!
 * \class classname
 *
 * \brief 类型的头文件
 *
 * \author FDL
 * \date 一月 2021
 */
#ifndef __st_Vertex3Normal__
#define __st_Vertex3Normal__
#include <vector>
struct Vertex3Normal
{
	float x, y, z; // Position data
	float nx, ny, nz; //normal
	Vertex3Normal() : x(.0f), y(.0f), z(.0f), nx(0.0f), ny(1.0f), nz(0.0f)
	{	}

	Vertex3Normal(float x, float y, float z, float nx, float ny, float nz) : x(x), y(y), z(z), nx(nx), ny(ny), nz(nz)
	{	}
};

typedef std::vector<Vertex3Normal> PtCloud;
#endif

#include "Face/include/Bbox.h"
#include "Face/include/Pnt3.h"
#include "Face/include/TriMeshUtils.h"
#include "Face/include/tnt_matrix.h"
#include "ObjReadWrite.h"

#include <QImage>

struct MeasureModel
{
	std::vector<omesh::Pnt3> points;		/// 点
	std::vector<omesh::Pnt3> normals;		/// 法线
	std::vector<omesh::Pnt3> texcoord;		/// uv
	std::vector<omesh::TriVtx> trisVtx;		/// 三角形
	Vector3Is mTextCoordIndexs;

	/// obj的贴图
	QImage img;
	bool hasImg;

	MeasureModel():hasImg(false){}

	/// 三角形转换成点云
	PtCloud getCloud(void)
	{
		auto tranPnt3ToVer3 = [](const omesh::Pnt3& point, const omesh::Pnt3& normal, Vertex3Normal& verTex3)
		{
			verTex3.x = point[0];
			verTex3.y = point[1];
			verTex3.z = point[2];
			verTex3.nx = normal[0];
			verTex3.ny = normal[1];
			verTex3.nz = normal[2];
		};

		PtCloud cloud;
		if (trisVtx.empty())
		{
			cloud.resize(points.size());
			for (int i = 0; i < points.size(); i ++)
			{
				tranPnt3ToVer3(points[i], normals[i], cloud[i]);
			}
		}
		else
		{
			omesh::Pnt3 nor;
			cloud.resize(trisVtx.size() * 3);
			for (int i = 0; i < trisVtx.size(); i++)
			{
				if (normals.empty())
					nor = omesh::normal(points[trisVtx[i][0]], points[trisVtx[i][1]], points[trisVtx[i][2]]);
				for (int j = 0; j < 3; j++)
				{
					int index = trisVtx[i][j];
					if (normals.empty())
					tranPnt3ToVer3(points[index], nor, cloud[i * 3 + j]);
					else
					tranPnt3ToVer3(points[index], normals[index], cloud[i*3+j]);
				}				
			}
		}
		return cloud;
	}

	PtCloud getPointCloud(void)
	{
		auto tranPnt3ToVer3 = [](const omesh::Pnt3& point, const omesh::Pnt3& normal, Vertex3Normal& verTex3)
		{
			verTex3.x = point[0];
			verTex3.y = point[1];
			verTex3.z = point[2];
			verTex3.nx = normal[0];
			verTex3.ny = normal[1];
			verTex3.nz = normal[2];
		};

		PtCloud cloud;
		cloud.resize(points.size());
		for (int i = 0 ; i < points.size(); i ++)
		{
			tranPnt3ToVer3(points[i], normals[i], cloud[i]);
		}
		return cloud;
	}

	PtCloud getTexcoordCloud(void)
	{
		auto tranPnt3ToVer3 = [](const omesh::Pnt3& point, const omesh::Pnt3& normal, Vertex3Normal& verTex3)
		{
			verTex3.x = point[0];
			verTex3.y = point[1];
			verTex3.z = point[2];
			verTex3.nx = normal[0];
			verTex3.ny = normal[1];
			verTex3.nz = normal[2];
		};

		PtCloud cloud;
		if (trisVtx.empty())
		{
			cloud.resize(points.size());
			for (int i = 0; i < points.size(); i++)
			{
				tranPnt3ToVer3(points[i], normals[i], cloud[i]);
			}
		}
		else
		{
			cloud.resize(trisVtx.size() * 3);
			for (int i = 0; i < trisVtx.size(); i++)
			{
				for (int j = 0; j < 3; j++)
				{
					int index = trisVtx[i][j];
					int indext = mTextCoordIndexs[i][j]-1;
					tranPnt3ToVer3(points[index], texcoord[indext], cloud[i * 3 + j]);
				}
			}
		}
		return cloud;
	}
};

enum ModelScanType
{
	en_bodystand1,		/// 站1
	en_bodystand2,		/// 站2
	en_bodystand3,		/// 站3
	en_bodysit,			/// 坐
	en_headHead,		/// 头
	en_hand1,			/// 手1
	en_hand2,			/// 手2
	en_hand3,			/// 手3
	en_hand4,			/// 手4
	en_hand5,			/// 手5
	en_hand6,			/// 手6
	en_foot,				/// 脚
	en_modelScanNull
};

static const QString EnScanTypeStr[] = {"_a", "_b", "_t","_c", "_h", "_l1","_l2","_l3","_r1","_r2","_r3","_f" };

enum ModelType
{
	en_objTex,	/// 模型带贴图
	en_obj,		/// 单obj
	en_stl		/// stl 模型
};

static const ModelType modelFileType[] =
{
	en_stl,		/// 站1
	en_obj,		/// 站2
	en_stl,		/// 站3
	en_stl,			/// 坐
	en_objTex,		/// 头
	en_objTex,			/// 手1
	en_objTex,			/// 手2
	en_objTex,			/// 手3
	en_objTex,			/// 手4
	en_objTex,			/// 手5
	en_objTex,			/// 手6
	en_stl				/// 脚
};

static const QString EnScanTypeExtStr[] = 
{
	".stl",		/// 站1
	".obj",		/// 站2
	".stl",		/// 站3
	".stl",			/// 坐
	".obj",		/// 头
	".obj",			/// 手1
	".obj",			/// 手2
	".obj",			/// 手3
	".obj",			/// 手4
	".obj",			/// 手5
	".obj",			/// 手6
	".stl"				/// 脚
};


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


struct ScanData
{
	bool load;
	MeasureModel modelData;
	std::vector<EsLineData> mLineData;
	ScanData():load(false){}
};