#pragma once
#ifndef __CLOUD_H__
#define __CLOUD_H__
/********************************************************************
	created:	2014/03/27
	created:	27:3:2014   11:09
	filename: 	d:\test program\Face\Face\Cloud.h
	file path:	d:\test program\Face\Face
	file base:	Cloud
	file ext:	h
	author:		feng
	
	purpose:	�������ݽṹ
*********************************************************************/

#include "stdafx.h"
#include "Xform.h"
#include "tnt_matrix.h"
#include "utils.h"
#include "EasyBMP.h"
#include "pif_format.h"
#include "RangeGrid.h"
#include "Mesh.h"
#include "PointProperty.h"

//#define File_Bmp	".bmp"  /// ��Ӧ��ͼƬ
//#define File_Txt	".txt"  /// ��Ӧ��txt
//#define File_xf	".xf"  /// ��Ӧ��txt
//#define File_infbin "_inf.bin" /// ��Ӧ��inf.bin
//#define File_roc	""		/// ��Ӧ��roc

#define File_Bmp	".bmp"  /// ��Ӧ��ͼƬ
#define File_Txt	".txt"  /// ��Ӧ��txt
#define File_xf	".xf"  /// ��Ӧ��txt
#define File_infbin "_inf.bin" /// ��Ӧ��inf.bin
#define File_roc	""		/// ��Ӧ��roc

#define wrongLenght 0.5		/// �������̾�����ж�

using omesh::Xform;
using TNT::Matrix;
using omesh::Mesh;
using omesh::RangeGrid;

typedef std::vector<ESTriangle> Triangles;
typedef std::vector<PointProperty> PointPropertys;

typedef std::list<int> Hole;
typedef std::vector<Hole> Holes;

const char wrlFileHeader[] = "#VRML V2.0 utf8\n\
#\n\
# Geomagic Studio\n\
#\n\
Viewpoint\n\
{\n\
	position -22.445450 -22.256793 -147.842509\n\
		orientation 0 0 1 0\n\
		fieldOfView 0.785398\n\
}\n\
DEF T Transform\n\
{\n\
	children\n\
		[\n\
			Shape\n\
			{\n\
				appearance Appearance\n\
				{\n\
					material Material\n\
					{\n\
						diffuseColor 0.7 0.7 1\n\
							ambientIntensity 0.9\n\
							specularColor 0.8 0.8 0.8\n\
					}\n\
				}\n\
				geometry IndexedFaceSet\n\
				{\n\
					solid FALSE\n\
";
const char wrlFileTail[] = "				}\n			}\n		]\n}\n";


struct TriInt
{
	int index1;
	int index2;
	int index3;
	int triIndex;
	/// �ͺ�����ԣ�����index ������Ϊ������index1��index2���ڱ� �п���������������һ��Ҳ�Ǹ���
	TriInt(int i1, int i2, int i3, int tri):index1(i1), index2(i2), index3(i3),triIndex(tri){}
};

typedef std::vector<TriInt> VecTriInt;

/// ������ͬ�㣬ɾ����ͬ��
struct OrderProducePoint
{
	bool h;
	bool r;
	int index;
	Ogre::Vector3 pos;
};
typedef std::vector<OrderProducePoint> VecOrderProducePoint;

/// ��¼�����У�ǰ���Ķ�Ӧ��ϵ��Ϊ�ӳٹ���׼��
struct LinkInt
{
	int hostIndex;
	int guestIndex;
	LinkInt(int h, int g):hostIndex(h), guestIndex(g) {}
};

//typedef std::vector<LinkInt> VecLinkInt;

typedef std::vector<std::pair<int, int> > VecLinkInt;



struct ProjectPoint
{
	bool mIsOverlappedIdx;		/// �Ƿ����غϵ�
	int mIndex;		/// ͶӰ�����������
	int mPointIndex;	/// ͶӰ�����
	double mUV[3];	/// ͶӰ��uv
	Ogre::Vector3 mPostion;		/// ͶӰ������	
	omesh::Pnt3 mColor;		/// ͶӰ����ɫ
	ProjectPoint()
		:mIndex(0)
		,mIsOverlappedIdx(false)
		,mPointIndex(0)
	{
		mUV[0] = 0.0f;
		mUV[1] = 0.0f;
		mUV[2] = 0.0f;
	}
};

//typedef std::map<int, ProjectPoint> MapPointProject;
typedef std::vector<std::pair<int, ProjectPoint> > MapPointProject;

struct ProducePoint
{
	bool replace;	/// ��������´����Ļ����滻�������´�����indexΪ��ǰ��������������������滻indexΪ������
	bool order;		/// ���������˳��
	int triIndex;	/// ��������������
	int srcPoint;	/// ԭ���ĵ�
	int line1;		/// �����ߵ�������1
	int line2;		/// �����ߵ�������2
	int targetIndex;		/// ��Ӧ�������
	int third;		/// �����߶��ŵĵ�
	float lenght;	/// ��ǰ��ĳ���
	Ogre::Vector3 pos;	/// ������λ��
	ProducePoint():replace(false),order(false),triIndex(0),srcPoint(0),line1(0),line2(0),targetIndex(0),lenght(0.0f){}

	void clear()
	{
		replace = false;
		order = false;
		triIndex = 0;
		srcPoint = 0;
		line1 = 0;
		line2 = 0;
		targetIndex = 0;
		lenght = 0.0f;
	}
};

struct HoleStatue
{
	int holeIndex;			/// �����ڵ�����
	bool statue;			/// ״̬
	int mapIndex;			/// ����ӳ���ϵ�λ��

	bool getStatue() { return statue; }
	void setStatue( bool s) { statue = s; }
	int getHoleIndex( void ) { return holeIndex; }
	void setHoleIndex( int i ) { holeIndex = i;}
	int getIndex( void ) { return mapIndex; };

	HoleStatue():holeIndex(-1), statue(false) {}
	HoleStatue(int i):holeIndex(-1), mapIndex(i), statue(false) {}
	HoleStatue(int i, bool s):holeIndex(i), statue(s) {}	
};

typedef std::map<int, HoleStatue> PointBorderLink;

typedef std::vector<ProducePoint> VecProducePoint;
typedef std::multimap<int, ProducePoint> MulProducePoint;
typedef std::pair<MulProducePoint::iterator, MulProducePoint::iterator> MulProPointInter;


class Cloud
{
public:
	Cloud(void)
		:mMesh(NULL)
	{

	}
	~Cloud(void)
	{
		if (mMesh)
		{
			delete mMesh;
			mMesh = NULL;
		}
	}

	void setNumber(int number );
	bool readData(const std::string& path, int subSamp=1);

	bool readData( const std::string& path, const std::string RocPaht, 
		const Matrix<float> &R, const Matrix<float> &G, const Matrix<float> &B, 
		const Xform<float>& XF, const Matrix<double>& KK, const Matrix<double>& KC, 
		const Matrix<double>& RM, const Matrix<double>& TM ,int subSamp/*=1*/ );

	bool readDataPCD(const std::string& path, int subSamp=1);

	bool transPointXF(const std::string& xf, std::vector<omesh::Pnt3>& pts_xyz, std::vector<omesh::Pnt3>& pts_normal);
	bool ReadTransformResult(Xform<float>& xf, string fileName);

	//bool transPointXF(const std::string& xf, std::vector<XYZ>& pts_xyz, std::vector<NORMAL>& pts_normal);

	Matrix<double> project2oulu(const vector <omesh::Pnt3> &vtx,const Matrix<double> &r,const Matrix<double> &t,double* fc,double* cc,double* kc);
	omesh::Pnt3 BilInter(float x, float y, const Matrix<float> &R, const Matrix<float> &G, const Matrix<float> &B);
	omesh::Pnt3 GetPointRGB(int x,int y,const Matrix<float> &R, const Matrix<float> &G, const Matrix<float> &B);
	bool GetKFC( double *kc,double *fc,double *cc,Matrix<double> KK_right_l2m,Matrix<double> kc_right_l2m);
	void getPreData(const std::string& path, Matrix<float> &R, Matrix<float> &G, Matrix<float> &B, 
		Xform<float>& XF, Matrix<double>& KK, Matrix<double>& KC, Matrix<double>& RM, Matrix<double>& TM);


	int getNumber( void )
	{
		return mNumber;
	}

	Mesh* getMesh( void )
	{
		return mMesh;
	}

	const Xform<float>& getXform( void )
	{
		 return mXform;
	}

	/// ��ʼ������
	void initData( void );

	/// ��ʼ��������
	void initPoints( void );

	/// ��ʼ������������
	void initTriangles( void );


	/// Ԥ�ȴ����������
	void makeTriangleArea( void );

	const PointPropertys& getPointPropertys( void ) const
	{
		return mPointPropertys;
	}

	PointPropertys& getPointPropertys( void )
	{
		return mPointPropertys;
	}

	const Triangles& getTriangles( void ) const
	{
		return mTriangles;
	}

	Triangles& getTriangles( void )
	{
		return mTriangles;
	}

	/// �õ��߽��ص�������
	std::vector<omesh::Pnt3>& getPntBorders( void )
	{
		return mPntBorders;
	}

	VectorIndex& getPntBorderIndexs( void )
	{
		return mPntBorderIndexs;
	}

	/*VecNewPoints& getNewPoints( void )
	{
		return mNewPoints;
	}*/

	/// �����wrl�ļ�
	/// fileName Ҫ������ļ���
	/// ���� true ����ɹ� false ����ʧ��
	bool saveWrl(const std::string& fileName);

	/// ��ȡwrl�ļ�
	/// fileName Ҫ������ļ���
	/// ���� true ��ȡ�ɹ� false ��ȡʧ��
	bool readWrl( const std::string& fileName );

	/// �ϲ�����һ������
	bool mergeOther(Cloud* pOther);

	/// ���ݵõ�����Ϣɾ��mesh�ĵ������Ϣ
	/// ���� true ,��ɾ������沢��ɾ���ɹ�
	/// ���� false û�п�ɾ���Ķ���
	bool deleteTriVtxFromMesh( void );

	/// �����ر����ر���ɫ
	/// type ,����
	void makeTestColor(const omesh::Pnt3& color);

	/// ɾ���ظ���
	void deleteSamePoint( void );

	/// ȡ��ĳ�������һ���߽��
	/// int curIndex ��ǰ��
	/// int lastIndex ǰһ����
	int getNextBorderPoint( int curIndex, int lastIndex );

	/// ����������ӵ�
	/// pos ������λ��
	/// return ���ӵ�������ڵ����е�����
	int addPoint(const Ogre::Vector3& pos);

	/// �����������������
	/// index1 ��������εĵ�һ��������
	/// index2 ��������εĵڶ���������
	/// index3 ��������εĵ�����������
	int addTriangle( int index1, int index2, int index3);

	/// ȡ���ͺ��������б�
	VecTriInt& getVecTriInt( void ) { return mVecTriInt; }
	const VecTriInt& getVecTriInt( void ) const { return mVecTriInt; }

	/// ����ؽ��ͺ������εĵ�������ι�ϵ
	/// size ����1�ĳߴ�
	/// vecTri �����ӵ������ι�ϵ
	void buildVecTriInt( int size, const VecTriInt& vecTri );

	/// �ؽ�ģ��
	void reBuildData( void );
	
	/// �ѵ������һ���㻻������һ����
	/// src ���滻�ĵ�
	/// target Ҫ�滻�ɵĵ�
	/// size ��Ҫ����ƫ��
	void replacePoint( int src, int target, int size);

	/// ɾ���ϲ�ʱ��������ظ���
	void deleteSameProducePoint( int size );

	VecLinkInt& getLinkInts( void )
	{
		return mLinkInts;
	}

	/// ���ҹ�ͬ���������
	void findSamePoints( void );

	/// �Ѻϲ��������θ��ƽ��б���
	/// size ��ͼ��ƫ��λ��
	void makeMergeTri( int size );

	/// �滻�ϲ���Ҫ�滻�ĵ�
	void makeMergePoint( int size );

	void addLinePoint(int line1, int line2)
	{
		mMegerKeyPoint.push_back(line1);
		mMegerKeyPoint.push_back(line2);
	}

	void addReplacePoints(int index )
	{
		mReplacePoints.push_back(index);
	}

	/// ��������
	void clear();

	//VectorIndex  mHostPointIndex;	/// �������ڿ͵���������������
	VectorIndex& getHostPointIndex( void )
	{
		return mHostPointIndex;
	}
	//VectorFloat mHostPointLenght;	/// �����Ƶ��͵��������ľ���
	VectorFloat& getHostPointLenght( void )
	{
		return mHostPointLenght;
	}

	/// ͶӰ��ͶӰ����Ϣ
	//MapPointProject mMapPointProperty;
	MapPointProject& getMapPointProperty( void )
	{
		return mMapPointProperty;
	}

	VectorIndex& getBordersFromDelete( void )
	{
		return mBordersFromDelete;
	}

	//KDTree *mKdtree;


	/// ȡ��������Ȩ����Ϣ��������ɫ���ɵ�ʱ���õ�
	float getPointWeight( int index ) const;

	bool compPoint(const int& a, const int&b);

	void reserve( int size, int triSize)
	{
		mTriangles.reserve(triSize);
		mPointPropertys.reserve(size);
		if (mMesh)
		{
			mMesh->vtx.reserve(size);
			mMesh->vtxRgb.reserve(size);
			mMesh->nrm.reserve(size);
			mMesh->tris.reserve(triSize);
		}
	}

	void saveCurDeletePointStatue( int number )
	{
		char buf[32] = {0};
		sprintf(buf, "%dpointDelete.txt", number);
		FILE* f = fopen(buf, "w");
		if ( f == NULL)
		{
			return;
		}
		for (int i = 0; i < mPointPropertys.size(); i ++)
		{
			fprintf(f, "%d id:%d %d\n", i, mPointPropertys.at(i).getIndex(), mPointPropertys.at(i).getNeedDelete());
		}
		fclose(f);
	}

	void getFromOtherMesh( Cloud* pCloud );

	/// ����mesh��Ϣ
	bool saveMesh( const std::string& fileName );

	/// ����ƽ������
	void makeAverageEdgeLenght( void );

	/// ȡ����������ڵ���е㣬�����������
	omesh::Pnt3 getMidBorderPos(int index);

	/// ȡ����������ڵ���е㣬�����������
	/// ƽ����Ե
	omesh::Pnt3 getSmoothBorderPos(int index, float f = 0.5);

	/// ��һ���ھӵ㲻��ͶӰ��
	bool getAnNeighborNotPropertyPoint( int index)
	{
		if (mPointPropertys.at(index).getPropertyPoint())
		{
			return false;
		}
		const VectorIndex& vecInts = mPointPropertys.at(index).getBorderList();
		for (int i = 0; i < vecInts.size(); i ++)
		{
			if (!mPointPropertys.at(i).getNeedDelete() && !mPointPropertys.at(i).getPropertyPoint())
			{
				return true;
			}
		}
		return false;
	}

	PointBorderLink& getPointBorderLink( void )
	{
		return mPointBorderLink;
	}

	Holes& getHoles( void )
	{
		return mHoles;
	}

	VectorIndex& getmWrongPoints( void )
	{
		return mWrongPoints;
	}

	void setMesh( omesh::Mesh* p)
	{
		if (p != NULL)
		{
			if (mMesh)
			{
				delete mMesh;
			}
			mMesh = p;;
		}
	}

	/// ȡ����һ������
	int getNextBorderIndex(int index, int last = -1);

protected:
	/// ��points��ĵ����Խ��и��£������е����������Լ�triSize�����е����Լ�pointSize
	/// points �����Լ���
	/// triSize ������ƫ�Ƴߴ�
	/// pointSize ��ƫ�Ƴߴ�
	void updatePoints(PointPropertys& points, int triSize, int pointSize) const;

	/// ����ͬ��ļ���vecɾ��ֻʣ1��
	/// vec ��ͬ�ĵ�ļ���
	/// size ��ߴ�
	void replaceVec( VectorIndex& vec, int size);

	/// ���������ƺϲ���Ĵ������кϲ�
	void megerKeyPoints( int size );

	/// ͨ�����ķ�ʽ����Ƿ��Ǳ߽�������
	bool isBuildBorderPoint( int index );

	/// �������ɵĵ��ؽ��߽��ϵ
	/// ���ʱ��߽�Ҳ���ǲ�׼ȷ�ģ����������ι�ϵ����׼ȷ��
	void reBuildKeyPoints( int size );

	/// ������޸����㷨�ǰ����ڵĶ��Ǳ߽��ĵ�ϳ�һ����
	void repairWrongPoints( void );

	/// ��n����ϳ�һ���㣬Ȼ��������յ�ƽ��������λ��
	void makeVecPointCenter( const VectorIndex& vec);

	/// ���½�����Щ��Ĺ�ϵ
	void rebuildVecPoints(VectorIndex& v, int size);

	void readLineData( int state, FILE* file );

	/// ɾ��mesh�����ͬ��
	/// ��Ҫ�����ݳ�ʼ��ǰ���������ݵ���ʱ�����ⲽ����Ȼ����
	void deleteMeshSamePoint( void );

	/// ɾ���ϳ�������ͼ���ʱ���������ͬ��
	void deleteSamePointMerge( void );
	public:
	/// ����ߺͶ��Ĺ�ϵ
	void initHoles( void );

private:
	/// ���ݵ����
	int mNumber;

	//float mAverageEdgeLenght;		/// ƽ���߳�����

	std::string mErrorString;		/// ��������

	/// ģ������
	Mesh* mMesh;

	//int mMaxLevel;		/// ������
	//int mMinLevel;		/// ��С����

	Xform<float> mXform;

	///// ��ɫ
	//Matrix<float> mImageR;
	//Matrix<float> mImageG;
	//Matrix<float> mImageB;

	/// �߽��ص�������
	std::vector<omesh::Pnt3> mPntBorders;

	/// �߽��ص�����ԭ����λ��
	/// ������������
	VectorIndex mPntBorderIndexs;

	PointPropertys mPointPropertys;		/// �������б�
	Triangles mTriangles;		/// ���б�

	/// ���ӵ��������б�
	VecTriInt mVecTriInt;

	VecLinkInt mLinkInts;

	/// �ϲ���ı߽��
	VectorIndex mMegerKeyPoint;

	VectorIndex mReplacePoints;
	
	/// ����㼯�ϣ����޸����Ǳ߽�����̵�
	VectorIndex mWrongPoints;

	/// ��Ϊɾ�������ı߽��
	VectorIndex mBordersFromDelete;

	PointBorderLink mPointBorderLink;
	Holes mHoles;

	///-------------------
	/// �Ӽ��������ù���������
	VectorIndex  mHostPointIndex;	/// �������ڿ͵���������������
	VectorFloat mHostPointLenght;	/// �����Ƶ��͵��������ľ���

	/// ͶӰ��ͶӰ����Ϣ
	MapPointProject mMapPointProperty;

	//KDTree *mKdtree;

};

#endif