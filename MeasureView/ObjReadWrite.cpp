#include "ObjReadWrite.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>


#define  _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>

ObjReadWrite::ObjReadWrite(void)
{
}


ObjReadWrite::~ObjReadWrite(void)
{
}

bool ObjReadWrite::readFile( const std::string& fileName )
{
	using namespace std;

	mFromFile = fileName;
	ifstream in;
	string line,word,ptname,ntname,fname;
	std::string txtPath;
	in.open(fileName.c_str());
	if(!in)
	{
		cout<<"Read obj error !"<<endl;
		return false;
	}
	while(getline(in,line))
	{
		if(line.size()==0||line[0]=='#') continue;
		istringstream is(line);
		is>>word;
		if(word=="v")
		{
			Vector3F p;
			is>>p.x>>p.y>>p.z;
			mPoints.push_back(p);
		}
		else if(word=="vt")
		{
			Vector3F p;
			is>>p.x>>p.y;
			mTextCoords.push_back(p);
		}
		else if(word=="vn")
		{
			Vector3F p;
			is>>p.x>>p.y>>p.z;
			mNormals.push_back(p);
		}
		else if(word=="o"||word=="g")
		{
			/// 不做那么复杂的考虑
			//if(!goname.empty()&&!faces.empty())
			//{
			//	Object obj(vertexs.begin(),vertexs.end(),texcoords.begin(),texcoords.end(),normals.begin(),normals.end(),faces.begin(),faces.end(),row,col,mtlname);
			//	while(n.count(goname)!=0)
			//	{
			//		goname.push_back('0');
			//	}
			//	m.insert(make_pair(goname,obj));
			//	n.insert(goname);
			//	faces.clear();
			//}
			//is>>goname;
		}
		else if(word=="f")
		{
			int r = 0,c = 0;
			Vector3I trangle;
			Vector3I normal;
			Vector3I textoors;
			while(is>>word)
			{
				c = count(word.begin(),word.end(),'/');
				if(c==0)
				{
					//faces.push_back(atoi(word.c_str()));
					trangle[r] = atoi(word.c_str());
				}
				else if(c==1)
				{
					//faces.push_back(atoi(string(word.begin(),word.begin()+word.find("/")).c_str()));
					//faces.push_back(atoi(string(word.begin()+word.find("/")+1,word.end()).c_str()));
					trangle[r] = atoi(string(word.begin(),word.begin()+word.find("/")).c_str());
					textoors[r] = atoi(string(word.begin()+word.find("/")+1,word.end()).c_str());
				}
				else if(c==2)
				{
					int a = word.find("/");
					int b = word.find("/",a+1);
					//faces.push_back(atoi(string(word.begin(),word.begin()+a).c_str()));
					//faces.push_back(atoi(string(word.begin()+a+1,word.begin()+b).c_str()));
					//faces.push_back(atoi(string(word.begin()+b+1,word.end()).c_str()));
					trangle[r] = atoi(string(word.begin(),word.begin()+a).c_str());
					textoors[r] = atoi(string(word.begin()+a+1,word.begin()+b).c_str());
					normal[r] = atoi(string(word.begin()+b+1,word.end()).c_str());
				}
				++r;

			}
			mTriangles.push_back(trangle);
			if (c > 0)
				mTextCoordIndexs.push_back(textoors);
			if (c > 1)
				mNormalIndexs.push_back(normal);
			//row = r;
			//col = c+1;
		}
		else if(word=="mtllib")
		{
			//is>>word;
			//ReadMtl(cd,word,matname);
		}
		else if(word=="usemtl")
		{
			//is>>mtlname;
		}
	}
}

bool ObjReadWrite::writeFile( const std::string& fileName )
{
	std::ofstream output;
	output.open ( fileName.c_str ( ) );

	if ( !output )
	{
		std::cerr << "\n";
		std::cerr << "OBJ_WRITE - Fatal error!\n";
		std::cerr << "  Could not open the output file \"" << fileName << "\".\n";
		return false;
	}

	output << "# " << fileName << "\n";
	output << "# created by feng\n";
	output << "\n";

	int pos = fileName.find_last_of("\\");
	if (pos == std::string::npos)
	{
		pos = fileName.find_last_of("/");
	}
	std::string outFile = fileName.substr(pos+1);
	std::string mtlName = outFile.substr(0, outFile.size() - 4);
	mtlName.append(".mtl");

	output << "mtllib " << mtlName << "\n";
	output << "usemtl material_0\n";

	std::string path = fileName.substr(0, pos);
	path += "\\";
	std::string txtName = outFile.substr(0, outFile.size() - 4);
	txtName.append(".bmp");
	mTextName = txtName;
	writeMtl(path + mtlName);

	/// 拷贝图像文件，图省事
	std::string img = mFromFile.substr(0, mFromFile.size() - 4);
	img = img + ".bmp";
	fileCopy(img, path + txtName);
	//saveTexture(txtName, _RMatrix, _GMatrix, _BMatrix);

	//int width = _RMatrix.num_cols();
	//int height = _RMatrix.num_rows();

	/// 写点坐标
	for (int i = 0; i < mPoints.size(); i ++)
	{
		output << "v " << mPoints[i].x << " " << mPoints[i].y << " " << mPoints[i].z << "\n";		/// 写一个点
	}
	/// 写法线
	for (int i = 0; i < mNormals.size(); i ++)
	{
		output << "vn "<< mNormals[i].x << " "<<mNormals[i].y<<" "<<mNormals[i].z<<"\n";
	}

	/// 写uv坐标
	for ( int i = 0; i < mTextCoords.size(); i ++)
	{
		output << "vt " << mTextCoords[i].x << " " << mTextCoords[i].y << "\n";
	}

	int mode = 0;
	if (!mTriangles.empty())
	{
		mode = (mode|1);
	}
	if (!mTextCoordIndexs.empty())
	{
		mode = mode | (1 << 1);
	}
	if (!mNormalIndexs.empty())
	{
		mode = mode | (1<<2);
	}
	//f Vertex1/Texture1/Normal1 Vertex2/Texture2/Normal2 Vertex3/Texture3/Normal3
	/// 写面信息
	for (int i = 0; i < mTriangles.size(); i ++)
	{
		/// 只有面
		if (mode == 1)
		{
			output << "f " << mTriangles[i].x << " "<<mTriangles[i].y<<" "<<mTriangles[i].z<<"\n"; 
		}
		/// 没有法线
		if (mode == 3)
		{
			output << "f "<< mTriangles[i].x<< "/" << mTextCoordIndexs[i].x <<"/" ;	/// 没有法线列
			output << " "<< mTriangles[i].y<< "/" << mTextCoordIndexs[i].y <<"/" ;	/// 没有法线列
			output << " "<< mTriangles[i].z<< "/" << mTextCoordIndexs[i].z <<"/ \n" ;	/// 没有法线列
		}
		/// 没有uv, 有法线
		if (mode == 5)
		{
			output << "f "<< mTriangles[i].x<< "//" << mNormalIndexs[i].x ;	/// 没有法线列
			output << " "<< mTriangles[i].y<< "//" << mNormalIndexs[i].y  ;	/// 没有法线列
			output << " "<< mTriangles[i].z<< "//" << mNormalIndexs[i].z <<"\n" ;	/// 没有法线列
		}
		/// 全部都有
		if (mode == 7)
		{
			output << "f "<< mTriangles[i].x<< "/"<< mTextCoordIndexs[i].x <<"/"<< mNormalIndexs[i].x ;	/// 没有法线列
			output << " "<< mTriangles[i].y<< "/" << mTextCoordIndexs[i].y <<"/"<< mNormalIndexs[i].y  ;	/// 没有法线列
			output << " "<< mTriangles[i].z<< "/" << mTextCoordIndexs[i].z <<"/"<< mNormalIndexs[i].z <<"\n" ;	/// 没有法线列
		}
	}

	output.close();
	return true;
}

bool ObjReadWrite::writeMtl( const std::string& fileName )
{
	std::ofstream fo(fileName.c_str());
	if (!fo.good()) {
		std::cout << "cannot open file: " << fileName.c_str() << std::endl;
		return false;
	}

	fo << "# Save by feng \n\n";

	fo << "newmtl material_0\n"; 
	fo << "Ka 0.200000 0.200000 0.200000\n";
	fo << "Kd 0.000000 0.000000 0.000000\n";
	fo << "Ks 1.000000 1.000000 1.000000\n";
	fo << "Tr 1.000000\n";
	fo << "illum 2\n";
	fo << "Ns 0.000000\n";
	fo << "map_Kd "<< mTextName;
	fo.close();
	return true;
}

void ObjReadWrite::clear()
{
	mPoints.clear();
	mTextCoords.clear();
	mNormals.clear();

	mTriangles.clear();
	mNormalIndexs.clear();
	mTextCoordIndexs.clear();
}

bool ObjReadWrite::fileCopy(const std::string& from, const std::string& to)
{
	FILE* ffrom = NULL;// fopen(from.c_str(), "rb");
	fopen_s(&ffrom, from.c_str(), "rb");
	FILE* fto = NULL;//fopen(to.c_str(), "wb");
	fopen_s(&fto, to.c_str(), "wb");
	if (ffrom == NULL || fto == NULL)
	{
		fclose(ffrom == NULL ? fto : ffrom);
		return false;
	}

	char buf[1024] = { 0 };

	int lenght = fread(buf, 1, 1024, ffrom);
	while (lenght > 0)
	{
		fwrite(buf, lenght, 1, fto);
		lenght = fread(buf, 1, 1024, ffrom);
	}
	fclose(ffrom);
	fclose(fto);
	return true;
}
