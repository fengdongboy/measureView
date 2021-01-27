#pragma once
/*!
 * \class classname
 *
 * \brief ���͵�ͷ�ļ�
 *
 * \author FDL
 * \date һ�� 2021
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
	std::vector<omesh::Pnt3> points;		/// ��
	std::vector<omesh::Pnt3> normals;		/// ����
	std::vector<omesh::Pnt3> texcoord;		/// uv
	std::vector<omesh::TriVtx> trisVtx;		/// ������
	Vector3Is mTextCoordIndexs;

	/// obj����ͼ
	QImage img;
	bool hasImg;

	MeasureModel():hasImg(false){}

	/// ������ת���ɵ���
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
	en_bodystand1,		/// վ1
	en_bodystand2,		/// վ2
	en_bodystand3,		/// վ3
	en_bodysit,			/// ��
	en_headHead,		/// ͷ
	en_hand1,			/// ��1
	en_hand2,			/// ��2
	en_hand3,			/// ��3
	en_hand4,			/// ��4
	en_hand5,			/// ��5
	en_hand6,			/// ��6
	en_foot,				/// ��
	en_modelScanNull
};

static const QString EnScanTypeStr[] = {"_a", "_b", "_t","_c", "_h", "_l1","_l2","_l3","_r1","_r2","_r3","_f" };

enum ModelType
{
	en_objTex,	/// ģ�ʹ���ͼ
	en_obj,		/// ��obj
	en_stl		/// stl ģ��
};

static const ModelType modelFileType[] =
{
	en_stl,		/// վ1
	en_obj,		/// վ2
	en_stl,		/// վ3
	en_stl,			/// ��
	en_objTex,		/// ͷ
	en_objTex,			/// ��1
	en_objTex,			/// ��2
	en_objTex,			/// ��3
	en_objTex,			/// ��4
	en_objTex,			/// ��5
	en_objTex,			/// ��6
	en_stl				/// ��
};

static const QString EnScanTypeExtStr[] = 
{
	".stl",		/// վ1
	".obj",		/// վ2
	".stl",		/// վ3
	".stl",			/// ��
	".obj",		/// ͷ
	".obj",			/// ��1
	".obj",			/// ��2
	".obj",			/// ��3
	".obj",			/// ��4
	".obj",			/// ��5
	".obj",			/// ��6
	".stl"				/// ��
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