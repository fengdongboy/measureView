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
	
	purpose:	�¼��򵥰��obj����д�����漰ͼƬ��Ⱥ�飩
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

	/// ���ж��㼯��
	Vector3fs mPoints;
	/// ����uv���꼯��
	Vector3fs mTextCoords;
	/// ���з��߼���
	Vector3fs mNormals;

	/// ����������
	Vector3Is mTriangles;
	/// ������������
	Vector3Is mNormalIndexs;
	/// uv ��������
	Vector3Is mTextCoordIndexs;

	/// obj �ļ�������
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