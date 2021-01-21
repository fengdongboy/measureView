#ifndef __POINTPROPERTY_H__
#define __POINTPROPERTY_H__

#include "Tools.h"
//#include "OgreVector2.h"
//#include "OgreVector3.h"
//#include "OgreVector4.h"
#include < windows.h > 

class Cloud;
class ESTriangle;
/// ������
class PointProperty
{
	friend class Cloud;
	friend class ESTriangle;
public:
	PointProperty();
	~PointProperty();

	inline float dist(const PointProperty& other) const
	{
		/*float f1 = this->mPoints[0] - other.mPoints[0];
		float f2 = this->mPoints[1] - other.mPoints[1];
		float f3 = this->mPoints[2] - other.mPoints[2];
		return sqrtf( f1*f1 + f2*f2 + f3*f3);*/
		return 0;
	}

	/// ������ռ�����ε����Ȩ��
	/*inline void addAreaValue(float value)
	{
		mAreaValue += value;
	}*/

	inline void addAreaNormal(const omesh::Pnt3& normal, float value)
	{
		mNormalArea += value;
		mAreaNormal += (normal*value);
	}

	inline void endAreaNormal( void )
	{
		if (mNormalArea == 0.0f)
		{
			//printf("the %d of normalarea is zero\n", mIndex);
			mNormalArea = 1.0f;
		}
		if (mNormalArea > 0)
		{
			mAreaNormal /= mNormalArea;
			//mAreaNormal.normalise();
			mAreaNormal.normalize();
		}
	}

	inline void makeBelongTri( int index )
	{
		//EnterCriticalSection(&mBelongTriCS);
		mBelongTri.push_back( index );
		//LeaveCriticalSection(&mBelongTriCS);
	}

	/// �õ�������Ƿ��ǹ�����
	inline bool isAlonePoint( void )
	{
		return mBelongTri.empty();
	}

	inline void addBorderList(int index1, int index2)
	{
		/// ���������Ĵ�߽��
		//EnterCriticalSection(&mBorderListCS);
		mBorderList.push_back(index1);
		mBorderList.push_back(index2);
		//LeaveCriticalSection(&mBorderListCS);
	}

	void removeBorderlist( int index1, int index2 );

	/// �����Ƿ��Ǳ߽����Ϣ
	void buildBorder( void );

	bool isRingPoint( void );

	/// ��������Ϊ�߽��
	/// ��������˵�߽��ֻ���ڴ���ģ�͵�ʱ��õ���������ɾ�����ʱ���ھӵ�Ҳ���ܱ�Ϊ�߽��
	inline void setBorder(bool border)
	{
		mBorder = border;
	}

	/// ȡ��������Ƿ��Ǳ߽��
	inline bool isBorderPoint( void ) const
	{
		return mBorder;
	}
	
	inline omesh::Pnt3 getAreaNormal( void ) const
	{
		return mAreaNormal;
	}

	/*inline const Ogre::Vector3& getNormal( void ) const
	{
		return mNormal;
	}*/

	/// ȡ�ñ߽��ļ���
	inline const std::vector<int>& getBorderList( void ) const
	{
		return mBorderList;
	}

	/// ȡ�������������ε����
	inline const std::vector<int>& getBelongTri( void ) const
	{
		return mBelongTri;
	}

	inline bool getNeedDelete( void ) const
	{
		return mNeedDelete;
	}

	inline void setNeedDelete( bool need )
	{
		mNeedDelete = need;
	}

	inline void erasePoint( int index)
	{
		mBorderList.erase(std::remove(mBorderList.begin(), mBorderList.end(), index), mBorderList.end());
	}

	inline void setLevelLayer( int level = 0 )
	{
		if ( mLevelLayer == -1)
		{
			mLevelLayer = level;
		}
		else
		{
			mLevelLayer = (std::max)(level, mLevelLayer);
		}
	}
	inline int getLevelLayer( void ) const
	{
		return mLevelLayer;
	}

	/// �õ������ͽ����ĵ�������ɵ����Ƿ��Ǳ߽���
	int BorderLine( int value ) const;

	inline int getIndex( void ) const
	{
		return mIndex;
	}

	inline void addIndexBy( int index)
	{
		mIndex += index;

		/// �����������κ������㶼����
	}

	inline void merge(const PointProperty& p)
	{
		const VectorIndex& v = p.getBelongTri();
		mBelongTri.insert(mBelongTri.end(), v.begin(), v.end());
		Tools::deleteSameInt(mBelongTri);

		/*
		const VectorIndex& v1 = p.getBorderList();
		mBorderList.insert(mBorderList.end(), v1.begin(), v1.end());
		std::sort(mBorderList.begin(), mBorderList.end());*/
	}

	/// �����е���������������size
	void updateBelongTri( int size )
	{
		//VectorIndex::iterator it1, it2;
		//it1 = mBelongTri.begin();
		//it2 = mBelongTri.end();
		//for (; it1!= it2; ++it1)
		//{
		//	(*it1) += size;
		//}
		#pragma omp parallel for
		for ( int i = 0; i < mBelongTri.size(); i ++)
		{
			mBelongTri[i] += size;
		}
	}

	/// �����е����Ե�����size
	void updateBorderList( int size );

	void eraseTri( int triIndex );

	void replaceBorderPoint(int src, int target)
	{
		std::replace(mBorderList.begin(), mBorderList.end(), src, target);
	}

	void removeBorderList( int value )
	{
		mBorderList.erase(std::remove(mBorderList.begin(), mBorderList.end(), value), mBorderList.end());
	}

	/// ��յ��ϵ
	inline void clearBorderList( void )
	{
		mBorderList.clear();
		mBorder = false;
	}

	///�Ƿ���ͶӰ��
	inline void setPropertyPoint( bool pp)
	{
		mPropertyPoint = pp;
	}

	inline bool getPropertyPoint( void ) const
	{
		return mPropertyPoint;
	}

	void setBorderLevelLayer( float level );

	inline float getBorderLevelLayer( void ) const
	{
		return mBorderLevelLayer;
	}
	 
	void removeBelongTris(const VectorIndex& v);

	void setProjectTriIndex(int index)
	{
		mProjectTriIndex = index;
	}

	int getProjectTriIndex( void ) const
	{
		return mProjectTriIndex;
	}

private:
	unsigned int mIndex;	/// ��ǰ������
	int mProjectTriIndex;		/// ͶӰ����������

	bool mBorder;		/// �Ƿ�߽��
	bool mNeedDelete;	/// �Ƿ���Ҫɾ��
	bool mPropertyPoint;	/// �Ƿ���ͶӰ��

	omesh::Pnt3 mAreaNormal;	/// ���������
	float mNormalArea;			/// ���������
	//float mAreaValue;			/// ���������

	int mLevelLayer;				/// ��������
	float mBorderLevelLayer;			/// �߽��������

	std::vector<int> mBelongTri;		/// �����������ε����---------���������
	std::vector<int> mBorderList;		/// �߽�㼯��
};

/// ����������
class ESTriangle
{
public:
	ESTriangle()
		:mIndex(-1)
		,mArea(0.0f)
		,mNeedDelete(false)
		,mBorderOverlappTri(false)
	{

	}
	~ESTriangle()
	{

	}

	inline int setArea( float value )
	{
		mArea = value;
	}


	//���������
	inline float TriangleArea(const PointProperty &P1,const PointProperty &P2,const PointProperty &P3)
	{
		float L1=P1.dist(P2);
		float L2=P2.dist(P3);
		float L3=P3.dist(P1);
		float S=(L1+L2+L3)/2;
		return sqrt(S*(S-L1)*(S-L2)*(S-L3));
	}

	inline float TriangleArea( const omesh::Pnt3& p1, const omesh::Pnt3& p2, const omesh::Pnt3& p3)
	{
		float L1 = dist(p1, p2);
		float L2 = dist(p2, p3);
		float L3 = dist(p1, p3);
		float S = (L1 + L2 + L3)/2;
		return sqrt(S*(S-L1)*(S-L2)*(S-L3));
	}

	Ogre::Vector3 normal(const PointProperty &a, const PointProperty &b, const PointProperty &c);

	omesh::Pnt3 normal(const omesh::Pnt3& a, const omesh::Pnt3& b, const omesh::Pnt3&c );

	/// �������������
	/// ���������������к�
	void makeArea( std::vector<PointProperty>& points )
	{
		mArea = TriangleArea(points[ mPoints[0] ], points[ mPoints[1] ], points[ mPoints[2] ]);

		/// �õ�������
		//mNormal = normal(points[ mPoints[0] ], points[ mPoints[1] ], points[ mPoints[2] ]);
	}

	void makeArea( const std::vector<omesh::Pnt3>& points)
	{
		mArea = TriangleArea(points[ mPoints[0] ], points[ mPoints[1] ], points[ mPoints[2] ]);

		/// �õ�������
		mNormal = normal(points[ mPoints[0] ], points[ mPoints[1] ], points[ mPoints[2] ]);
	}

	void processArea( std::vector<PointProperty>& points)
	{
		for ( int i = 0; i < 3; i ++)
		{
			points[mPoints [i] ].addAreaNormal(mNormal, mArea);
		}
	}

	inline void setBorderOverlappTri(bool borOve =  true)
	{
		mBorderOverlappTri = borOve;
	}

	inline void setNeedDelete( bool idelete )
	{
		mNeedDelete = idelete;
	}

	inline bool getNeedDelete( void ) const
	{
		return mNeedDelete;
	}
	inline int getIndex( void ) const
	{
		return mIndex;
	}

	inline ESTriangle& operator=(const ESTriangle& other )
	{
		mNeedDelete = other.mNeedDelete;
		mIndex = other.mIndex;
		Tools::copyEveryTypePtr(mPoints, other.mPoints);
		return *this;
	}

//private:
	bool mBorderOverlappTri;			/// �Ƿ����ص��߽�������
	bool mNeedDelete;					/// �Ƿ���Ҫɾ������
	unsigned int mIndex;				/// ��ǰ����������
	unsigned int mPoints[3];			/// �����к�
	float mArea;			/// ���������
	omesh::Pnt3 mNormal;			/// ������ķ�����
};

#endif