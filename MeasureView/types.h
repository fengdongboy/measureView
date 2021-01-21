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