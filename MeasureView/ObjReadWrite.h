#pragma once
#ifndef __OBJREADWRITE_H__
#define __OBJREADWRITE_H__
/********************************************************************
	created:	2017/06/27
	created:	27:6:2017   18:06
	filename: 	d:\work program\Face_3DScannerV2.0\Interface_module2\ObjReadWrite.h
	file path:	d:\work program\Face_3DScannerV2.0\Interface_module2
	file base:	ObjReadWrite
	file ext:	h
	author:		feng
	
	purpose:	事件简单版的obj读和写（不涉及图片和群组）
*********************************************************************/
#include <string>
#include <vector>

template <class T>
struct Vector3
{
	T x;
	T y;
	T z;

	T& operator[](int i)
	{
		if (i == 0)
		{
			return x;
		}
		if (i == 1)
		{
			return y;
		}
		if (i == 2)
		{
			return z;
		}
	}

	const T& operator[](int i) const 
	{
		if (i == 0)
		{
			return x;
		}
		if (i == 1)
		{
			return y;
		}
		if (i == 2)
		{
			return z;
		}
	}

	Vector3& operator=(const Vector3& other)
	{
		this->x = other.x;
		this->y = other.y;
		this->z = other.z;
		return *this;
	}
};


typedef Vector3<float> Vector3F;
typedef Vector3<int> Vector3I;

typedef std::vector<Vector3F> Vector3fs;
typedef std::vector<Vector3I> Vector3Is;

class ObjReadWrite
{
public:
	ObjReadWrite(void);
	~ObjReadWrite(void);

	/// 所有顶点集合
	Vector3fs mPoints;
	/// 所有uv坐标集合
	Vector3fs mTextCoords;
	/// 所有法线集合
	Vector3fs mNormals;

	/// 面索引集合
	Vector3Is mTriangles;
	/// 法线索引集合
	Vector3Is mNormalIndexs;
	/// uv 索引集合
	Vector3Is mTextCoordIndexs;

	/// obj 文件的名字
	std::string mObjName;
	std::string mTextName;

	std::string mFromFile;

	bool readFile(const std::string& fileName );

	bool writeFile( const std::string& fileName );

	bool writeMtl( const std::string& fileName );

	void clear();

	bool fileCopy(const std::string& from, const std::string& to);
};

#endif