#ifndef _UTILS_H_
#define _UTILS_H_

#include "mesh.h"
#include "xform.h"
#include "tnt.h"

using namespace TNT;

namespace omesh
{

/*******************************************
����Xform
1.��r��t����Xformʱ��Xform.m = [r t]ת��
2.����Xformʱ��"*.xf"�ļ��б������[r t]
��Ϊ����saveXform(Xform<T>& xf, const char *name)�е�r����ת����r��ת�ã�
�������������[rת��ת�� t]��ʵ�ʵȼ���[r t]
3.template <class T> inline 
T& Xform<T>::operator()(int i, int j)
{
return m[j][i];
}
��� r[i][j] = m[j][i] = xf(i,j)
*******************************************/
template <class T>
	inline bool saveXformBin(Xform<T>& xf, const char *name)
{
	ofstream fs(name, ios::binary);
	if(!fs){
		fs.close();			
		return false;
	}

	T rt[4][4];
	for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
			rt[i][j] = xf(i,j);  //rt=[r t]

	int TSize = sizeof(T);

	fs.write((char*)(&(rt[0][0])), TSize);
	fs.write((char*)(&(rt[0][1])), TSize);
	fs.write((char*)(&(rt[0][2])), TSize);
	fs.write((char*)(&(rt[0][3])), TSize);

	fs.write((char*)(&(rt[1][0])), TSize);
	fs.write((char*)(&(rt[1][1])), TSize);
	fs.write((char*)(&(rt[1][2])), TSize);
	fs.write((char*)(&(rt[1][3])), TSize);

	fs.write((char*)(&(rt[2][0])), TSize);
	fs.write((char*)(&(rt[2][1])), TSize);
	fs.write((char*)(&(rt[2][2])), TSize);
	fs.write((char*)(&(rt[2][3])), TSize);

	fs.write((char*)(&(rt[3][0])), TSize);
	fs.write((char*)(&(rt[3][1])), TSize);
	fs.write((char*)(&(rt[3][2])), TSize);
	fs.write((char*)(&(rt[3][3])), TSize);

	fs.close();
	return true;
}

template <class T>
	inline bool readXformBin(Xform<T>& xf, const char *name)
{
	ifstream ifs(name, ios::binary);
	if(!ifs){
		ifs.close();
		return false;
	}

	T r[3][3], t[3];
	int TSize = sizeof(T);
	ifs.read((char*)(&r[0][0]), TSize);
	ifs.read((char*)(&r[0][1]), TSize);
	ifs.read((char*)(&r[0][2]), TSize);
	ifs.read((char*)(&t[0]), TSize);

	ifs.read((char*)(&r[1][0]), TSize);
	ifs.read((char*)(&r[1][1]), TSize);
	ifs.read((char*)(&r[1][2]), TSize);
	ifs.read((char*)(&t[1]), TSize);

	ifs.read((char*)(&r[2][0]), TSize);
	ifs.read((char*)(&r[2][1]), TSize);
	ifs.read((char*)(&r[2][2]), TSize);
	ifs.read((char*)(&t[2]), TSize);

	ifs.close();

	xf = Xform<T>(r,t);

	return true;
}

// ��Ļ�����xf.m
template <class T>
	inline void dispXf(const Xform<T>& xf)
{
	cout<<"Xform"<<endl;
	printf("%f ", xf(0,0));
	printf("%f ", xf(1,0));
	printf("%f ", xf(2,0));
	printf("%f\n", xf(3,0));

	printf("%f ", xf(0,1));
	printf("%f ", xf(1,1));
	printf("%f ", xf(2,1));
	printf("%f\n", xf(3,1));

	printf("%f ", xf(0,2));
	printf("%f ", xf(1,2));
	printf("%f ", xf(2,2));
	printf("%f\n", xf(3,2));

	printf("%f ", xf(0,3));
	printf("%f ", xf(1,3));
	printf("%f ", xf(2,3));
	printf("%f\n", xf(3,3));
}


template <class T>
inline bool readXform(Xform<T>& xf, const char *name)
{
	ifstream ifs(name);
	if(!ifs){
		ifs.close();
		return false;
	}

	T r[3][3], t[3];
	string bufline;	// read line buffer
	int lineIndex = 0;
	while(ifs.good()){
		getline(ifs, bufline, '\n');
		if(lineIndex<3 && bufline.length()!=0){
			int idx = -1;
			int lidx = 0;
			lidx = idx;	idx = int(bufline.find(' ', lidx+1));
			r[lineIndex][0] = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			lidx = idx;	idx = int(bufline.find(' ', lidx+1));
			r[lineIndex][1] = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			lidx = idx;	idx = int(bufline.find(' ', lidx+1));
			r[lineIndex][2] = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			lidx = idx;	idx = int(bufline.find(' ', lidx+1));
			t[lineIndex] = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
		}
		lineIndex ++;
	}
	ifs.close();

	xf = Xform<T>(r,t);

	return true;
}

template <class T>
inline bool readTxtfrom(Xform<T>& xf, const char *name)
{
	ifstream ifs(name);
	if(!ifs){
		ifs.close();
		return false;
	}
	int exponent = 0;
	T r[3][3], t[3];
	string bufline;	// read line buffer
	int lineIndex = 0;
    
	//r��t��index�������е�index
	int rtlineIndex = 0;
	while(ifs.good()){
		getline(ifs, bufline, '\n');
		if(lineIndex>1 && rtlineIndex<3 && bufline.length()!=0){
			int idx = -1;
			int lidx = 0;
			lidx = idx;	idx = int(bufline.find('e', lidx+1));
			r[rtlineIndex][0] = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			//lidx = idx;	idx++;
			lidx = idx;	idx = int(bufline.find(' ', lidx+1));
			exponent = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			if (exponent<0)
			{
				exponent = -exponent;
				for (float i=1;i<exponent+1;i++)
				{
					r[rtlineIndex][0] = r[rtlineIndex][0]/10;
				}
			}
			else if (exponent>0)
			{
				for (float i=1;i<exponent+1;i++)
				{
					r[rtlineIndex][0] = r[rtlineIndex][0]*10;
				}
			}
			/*else
			{
				
			}
			*/
			lidx = idx;	idx++;
			/// ÿ��������������� feng
			//cout<<r[rtlineIndex][0]<<" ";

			lidx = idx;	idx = int(bufline.find('e', lidx+1));
			r[rtlineIndex][1] = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			//lidx = idx;	idx++;
			lidx = idx;	idx = int(bufline.find(' ', lidx+1));
			exponent = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			if (exponent<0)
			{
				exponent = -exponent;
				for (float i=1;i<exponent+1;i++)
				{
					r[rtlineIndex][1] = r[rtlineIndex][1]/10;
				}
			}
			else if (exponent>0)
			{
				for (float i=1;i<exponent+1;i++)
				{
					r[rtlineIndex][1] = r[rtlineIndex][1]*10;
				}
			}
			/*else
			{
				
			}
			*/
			lidx = idx;	idx++;
			cout<<r[rtlineIndex][1]<<" ";


			lidx = idx;	idx = int(bufline.find('e', lidx+1));
			r[rtlineIndex][2] = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			//lidx = idx;	idx++;
			lidx = idx;	idx = int(bufline.find(' ', lidx+1));
			exponent = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			if (exponent<0)
			{
				exponent = -exponent;
				for (float i=1;i<exponent+1;i++)
				{
					r[rtlineIndex][2] = r[rtlineIndex][2]/10;
				}
			}
			else if (exponent>0)
			{
				for (float i=1;i<exponent+1;i++)
				{
					r[rtlineIndex][2] = r[rtlineIndex][2]*10;
				}
			}
			/*else
			{
				
			}
			*/
			lidx = idx;	idx++;
			cout<<r[rtlineIndex][2]<<" ";

			lidx = idx;	idx = int(bufline.find('e', lidx+1));
			t[rtlineIndex] = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			//lidx = idx;	idx++;
			lidx = idx;	idx = int(bufline.find(' ', lidx+1));
			exponent = T(atof((bufline.substr(lidx+1,idx-lidx-1)).c_str()));
			if (exponent<0)
			{
				exponent = -exponent;
				for (float i=1;i<exponent+1;i++)
				{
					t[rtlineIndex] = t[rtlineIndex]/10;
				}
			}
			else if (exponent>0)
			{
				for (float i=1;i<exponent+1;i++)
				{
					t[rtlineIndex] = t[rtlineIndex]*10;
				}
			}
			/*else
			{
				
			}
			*/
			lidx = idx;	idx++;
			cout<<t[rtlineIndex];
			cout<<endl;
            rtlineIndex ++;
		}
		lineIndex ++;
		
	}
	ifs.close();

	xf = Xform<T>(r,t);

	return true;
}

template <class T>
inline bool saveXform(Xform<T>& xf, const char *name)
{
	ofstream fs(name);
	if(!fs){
		fs.close();			
		return false;
	}

	T r[3][3], t[3];
	xf.get_e(0,r[0]); xf.get_e(1,r[1]); xf.get_e(2,r[2]);
	xf.getTranslation(t);

	fs << setprecision(16) << r[0][0] << " " << r[1][0] << " " << r[2][0] << " " << t[0] << endl;
	fs << setprecision(16) << r[0][1] << " " << r[1][1] << " " << r[2][1] << " " << t[1] << endl;
	fs << setprecision(16) << r[0][2] << " " << r[1][2] << " " << r[2][2] << " " << t[2] << endl;
	fs << setprecision(16) << 0 << " " << 0 << " " << 0 << " " << 1 << endl;

	fs.close();
	return true;
}

// format 0 binary 1 ascii
inline bool saveSTLFile(const char *filename, vector<Pnt3>& vtx, vector<Pnt3>& nrm, vector<TriVtx>& tris, int format=0)
{
	if(format == 0){
		ofstream fs(filename, ios::binary);
		if(!fs){ fs.close(); return false; }

		int i;
		int intSize = sizeof(int);
		int floatSize = sizeof(float);
		// �ļ�ͷ
		char* fileHead = "ROC";
		fs.write(fileHead, sizeof(char)*3);
		// ������Ϣ
		char fileInfo[77];
		for(i=0 ; i<77 ; i++) fileInfo[i]=' ';
		fs.write(fileInfo, sizeof(char)*77);
		// ��ĸ���
		int num_tris = int(tris.size());
		fs.write((char*)(&num_tris), intSize);
		// ���б����б�
		char a[2];
		streamsize a_size = sizeof(char)*2;
		Pnt3 tn;
		for(i=0 ; i<num_tris ; i++){
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
	if(format == 1){
		int n,i;
		ofstream fs(filename);
		if(!fs){
			fs.close();			
			return false;
		}
		fs << "solid ods\n";
		n = int(tris.size());
		Pnt3 tn;
		for(i=0 ; i<n ; i++){
			tn = nrm[tris[i][0]] + nrm[tris[i][1]] + nrm[tris[i][2]];
			fs << "facet normal " << tn[0] <<" "<< tn[1] <<" "<< tn[2] <<"\nouter loop\n"
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
};

inline bool readSTLFile(const char *filename, vector<Pnt3>& vtx, vector<TriVtx>& tris, int format=0)
{
	if(format == 0){
		ifstream ifs(filename, ios::binary);
		if(!ifs){ ifs.close(); return false; }

		vtx.clear();
		tris.clear();

		int intSize = sizeof(int);
		int floatSize = sizeof(float);
		ifs.ignore(80);

		// ��ĸ���
		int num_tris;
		ifs.read((char*)(&num_tris), intSize);
		
		float tn0, tn1, tn2;
		float v0, v1, v2;

		for(int i=0 ; i<num_tris ; i++){
			ifs.read((char*)(&(tn0)), floatSize);
			ifs.read((char*)(&(tn1)), floatSize);
			ifs.read((char*)(&(tn2)), floatSize);

			ifs.read((char*)(&(v0)), floatSize);
			ifs.read((char*)(&(v1)), floatSize);
			ifs.read((char*)(&(v2)), floatSize);
			vtx.push_back(Pnt3(v0,v1,v2));
			ifs.read((char*)(&(v0)), floatSize);
			ifs.read((char*)(&(v1)), floatSize);
			ifs.read((char*)(&(v2)), floatSize);
			vtx.push_back(Pnt3(v0,v1,v2));
			ifs.read((char*)(&(v0)), floatSize);
			ifs.read((char*)(&(v1)), floatSize);
			ifs.read((char*)(&(v2)), floatSize);
			vtx.push_back(Pnt3(v0,v1,v2));
			
			tris.push_back(TriVtx(i*3,i*3+1,i*3+2));

			ifs.ignore(2);
		}

		ifs.close();
		return true;

	}

	return false;
};

// �����ýṹ
struct VtxIdxSortItem {
	int i;
	Pnt3 value;
};
// ��������
inline bool VtxIdxSortItem_LessSecond_X(const VtxIdxSortItem & m1, const VtxIdxSortItem & m2) {
	return m1.value[0] < m2.value[0];
};
inline bool VtxIdxSortItem_LessSecond_Y(const VtxIdxSortItem & m1, const VtxIdxSortItem & m2) {
	return m1.value[1] < m2.value[1];
};
inline bool VtxIdxSortItem_LessSecond_Z(const VtxIdxSortItem & m1, const VtxIdxSortItem & m2) {
	return m1.value[2] < m2.value[2];
};
// ����XYZ����λ��������
inline bool VtxIdxSortItem_LessSecond(const VtxIdxSortItem & m1, const VtxIdxSortItem & m2) {
	bool b = false;

	if(m1.value[0] > m2.value[0])
		b = true;
	else if(m1.value[0] == m2.value[0])
		if(m1.value[1] > m2.value[1])
			b = true;
		else if(m1.value[1] == m2.value[1])
			if(m1.value[2] > m2.value[2])
				b = true;
			else
				b = false;
	
	return b;
};

// ����ֿ����
inline void UnifyDuplicatedVertices(vector<Pnt3> &vtx, vector<TriVtx> &tris)
{
	if(vtx.size() == 0)
		return;

	vector<Pnt3> vtxNew;
	vector<int> vtxIdxMap(vtx.size(),-1);

	vector<int> vtxIdx;
	for(int i=0 ; i<int(vtx.size()) ; i++){
		vtxIdx.push_back(i);
	}

	// ȥ��������
	vector<VtxIdxSortItem> vtxIdxSortItemList;
	VtxIdxSortItem item;
	for(int i=0; i<int(vtx.size()) ; i++){
		item.i = vtxIdx[i];
		item.value = vtx[i];
		vtxIdxSortItemList.push_back(item);
	}

	sort(vtxIdxSortItemList.begin(), vtxIdxSortItemList.end(), VtxIdxSortItem_LessSecond);
	
	vtxNew.push_back(vtxIdxSortItemList[0].value);
	vtxIdxMap[vtxIdxSortItemList[0].i] = int(vtxNew.size()-1);
	for(int i=1; i<int(vtxIdxSortItemList.size()) ; i++){
		if(vtxNew[vtxNew.size()-1] != vtxIdxSortItemList[i].value){
			vtxNew.push_back(vtxIdxSortItemList[i].value);
		}

		vtxIdxMap[vtxIdxSortItemList[i].i] = int(vtxNew.size()-1);
	}

	vtx = vtxNew;
	// �����µ��б�͸��������ε����
	for(int i=0 ; i<int(tris.size()) ; i++){
		tris[i][0] = vtxIdxMap[tris[i][0]];
		tris[i][1] = vtxIdxMap[tris[i][1]];
		tris[i][2] = vtxIdxMap[tris[i][2]];
	}

};

// ֻȥ�����
inline void UnifyDuplicatedVertices(vector<Pnt3> &vtx)
{
	if(vtx.size() == 0)
		return;

	vector<Pnt3> vtxNew;

	vector<int> vtxIdx;
	for(int i=0 ; i<int(vtx.size()) ; i++){
		vtxIdx.push_back(i);
	}

	// ȥ��������
	vector<VtxIdxSortItem> vtxIdxSortItemList;
	VtxIdxSortItem item;
	for(int i=0; i<int(vtx.size()) ; i++){
		item.i = vtxIdx[i];
		item.value = vtx[i];
		vtxIdxSortItemList.push_back(item);
	}

	sort(vtxIdxSortItemList.begin(), vtxIdxSortItemList.end(), VtxIdxSortItem_LessSecond);

	vtxNew.push_back(vtxIdxSortItemList[0].value);
	for(int i=1; i<int(vtxIdxSortItemList.size()) ; i++){
		if(vtxNew[vtxNew.size()-1] != vtxIdxSortItemList[i].value){
			vtxNew.push_back(vtxIdxSortItemList[i].value);
		}
	}

	vtx = vtxNew;
};

inline bool saveASCFile(const char *filename, vector<Pnt3>& vtx)
{
	int n,i;
	ofstream fs(filename);
	if(!fs){
		fs.close();			
		return false;
	}
	n = int(vtx.size());
	for(i=0 ; i<n ; i++){
		fs << vtx[i][0] <<" "<< vtx[i][1] <<" "<< vtx[i][2] << endl;
	}

	fs.close();
	return true;
};

// ɾ����
// del_list��vͬ��������1����ɾ����0��������
inline void delVtx(vector<Pnt3>& vtx, vector<TriVtx>& tris, vector<int>& del_list){

	// ɾ��������
	vector<TriVtx> old_tris = tris;
	tris.clear();
	for(int i= 0 ; i<int(old_tris.size()) ; i++){
		if(del_list[old_tris[i][0]]==0 && del_list[old_tris[i][1]]==0 && del_list[old_tris[i][2]]==0){
			tris.push_back(old_tris[i]);

		}
	}

	// ɾ����
	vector<Pnt3> old_vtx = vtx;
	vtx.clear();
	vector<int> vtx_index(old_vtx.size(), -1);
	int ind = 0;
	for(int i=0 ; i<int(old_vtx.size()) ; i++){
		if(del_list[i]==0){
			vtx.push_back(old_vtx[i]);
			vtx_index[i] = ind;
			ind ++;
		}
	}

	// �����������б�
	for(int i=0 ; i<int(tris.size()) ; i++){
		tris[i][0] = vtx_index[tris[i][0]];
		tris[i][1] = vtx_index[tris[i][1]];
		tris[i][2] = vtx_index[tris[i][2]];
	}

};

// �õ��ٽ����б�
inline vector< vector<int> > getNeighborVtx(vector<Pnt3>& vtx, vector<TriVtx>& tris)
{
	vector< vector<int> > neighborVtx;

	for(int i=0 ; i<int(vtx.size()) ; i++)
		neighborVtx.push_back(vector<int>());

	for(int i=0 ; i<int(tris.size()) ; i++){
		neighborVtx[tris[i][0]].push_back(tris[i][1]);
		neighborVtx[tris[i][0]].push_back(tris[i][2]);
		neighborVtx[tris[i][1]].push_back(tris[i][0]);
		neighborVtx[tris[i][1]].push_back(tris[i][2]);
		neighborVtx[tris[i][2]].push_back(tris[i][0]);
		neighborVtx[tris[i][2]].push_back(tris[i][1]);
	}

	// ȥ�������
	for(int i=0 ; i<int(neighborVtx.size()) ; i++){
		vector<int> t = neighborVtx[i];
		neighborVtx[i].clear();
		for(int j=0 ; j<int(t.size()) ; j++){
			bool a = true;
			for(int k=0 ; k<int(neighborVtx[i].size()) ; k++){
				if(t[j]==neighborVtx[i][k]){
					a = false;
					break;
				}
			}
			if(a==true)
				neighborVtx[i].push_back(t[j]);
		}		
	}

	return neighborVtx;
};

// �������ٽ��������б���Ϊ�㵽�����ε�����
inline vector< vector<int> > getVtxTris(vector<Pnt3>& vtx, vector<TriVtx>& tris)
{
	vector< vector<int> > vtxTris;

	for(int i=0 ; i<int(vtx.size()) ; i++)
		vtxTris.push_back(vector<int>());

	for(int i=0 ; i<int(tris.size()) ; i++){
		vtxTris[tris[i][0]].push_back(i);
		vtxTris[tris[i][1]].push_back(i);
		vtxTris[tris[i][2]].push_back(i);
	}

	return vtxTris;
}


// �õ��߽��б�
inline vector< vector<int> > getBoundaryBorder(vector<Pnt3>& vtx, vector<TriVtx>& tris)
{
	vector< vector<int> > boundaryList;

	vector< vector<int> > neighborVtx = getNeighborVtx(vtx, tris);
	vector< vector<int> > vtxTris = getVtxTris(vtx, tris);

	// ������ٽ��ĵ���ٽ��ĵ�����Ĳ�ֵ������0Ϊ�ڲ�����ǹ����㣬����0Ϊ�߽�㣬����1Ϊ���������Ľӵ�
	vector<int> vtxUsed;
	for(int i=0 ; i<int(vtx.size()) ; i++){
		vtxUsed.push_back( int( neighborVtx[i].size() - vtxTris[i].size() ) );
	}

	// ����߽�
	for(int N=0 ; N<100000 ; N++){
		if(N==100000){
			cout << "N�ﵽ���ѭ����" << endl;
			return boundaryList;
		}

		vector<int> boundary;
		// �ҳ���һ���߽��
		for(int i=0 ; i<int(vtxUsed.size()) ; i++){
			if(vtxUsed[i] == 1){ // ʹ�ñ߽����Ϊ��ʼ
				boundary.push_back(i);
				vtxUsed[i] --;
				break;
			}
		}
		if(boundary.empty())
			break;

		// Ѱ�ҵ�һ���߽�
		vector<int> t1 = vtxTris[boundary[0]];
		for(int i=0 ; i<int(t1.size()) ; i++){
			int nb;
			if(tris[t1[i]][0]==boundary[0]) nb=tris[t1[i]][1];
			if(tris[t1[i]][1]==boundary[0]) nb=tris[t1[i]][2];
			if(tris[t1[i]][2]==boundary[0]) nb=tris[t1[i]][0];

			// �ж��������Ƿ�ֻ�ǹ���һ��������
			if(vtxUsed[nb]>0){
				vector<int> t2 = vtxTris[nb];
				int s = 0;
				for(int m=0 ; m<int(t1.size()) ; m++)
					for(int n=0 ; n<int(t2.size()) ; n++)
						if(t1[m]==t2[n])
							s ++;
				if(s==1){
					boundary.push_back(nb);
					vtxUsed[nb] --;
					break;
				}
				if(s>2){
					cout << "���������ϵ������ι���" << endl;
					return boundaryList;
				}
			}
		}
		if(boundary.size()<2) break;

		for(int M=0 ; M<100000 ; M++){
			if(M==100000){
				cout << "M�ﵽ���ѭ����" << endl;
				return boundaryList;
			}

			int keeploop = 1;

			int nb;
			int nt;
			vector<int> u; // �Ѿ�ʹ�õ��ڽ����б�
			for(int K=0 ; K<100000 ; K++){
				if(K==100000){
					cout << "K�ﵽ���ѭ����" << endl;
					return boundaryList;
				}

				if(K == 0){
					nb = boundary[boundary.size()-2]; // �¼���ĵ�
					u.clear();
				}

				vector<int> t1 = vtxTris[boundary[boundary.size()-1]];
				vector<int> t2 = vtxTris[nb];
				int s=0;
				for(int m=0 ; m<int(t1.size()) ; m++){
					for(int n=0 ; n<int(t2.size()) ; n++){
						if(t1[m]==t2[n]){
							// �ж�t1[m]�Ƿ���u�б���
							vector<int>::iterator tmp;
							tmp = find(u.begin(), u.end(), t1[m]);
							if(tmp != u.end()) continue;

							nt = t1[m]; // �õ���������
							u.push_back(nt);
							s = s+1;
						}
					}
				}
				if(s>=2){
					cout << "���������ϵĹ��б߽��" << endl;
					return boundaryList;
				}
				// ���ֻ��һ���µ��ڽ�������
				if( s==1 ){
					// �õ��µ��ڽ���
					if( tris[nt][0]!=boundary[boundary.size()-1] && tris[nt][0]!=nb )
						nb = tris[nt][0];
					else if( tris[nt][1]!=boundary[boundary.size()-1] && tris[nt][1]!=nb )
						nb = tris[nt][1];
					else if( tris[nt][2]!=boundary[boundary.size()-1] && tris[nt][2]!=nb )
						nb = tris[nt][2];
					continue;
				}
				if( s==0 ){
					if( vtxUsed[nb] > 0 ){
						boundary.push_back(nb);
						vtxUsed[nb] --;
						break;
					}else{
						if( nb == boundary[0] ){
							keeploop = 0;
							break;
						}else{
							cout << "�߽粻���" << endl;
							return boundaryList;
						}
					}
				}
			}
			if(keeploop==0) 
				break;
		}

		boundaryList.push_back(boundary);

	}

	return boundaryList;
};

// ����������(Non-manifolds)����
inline void modifyNonmanifolds(vector<Pnt3>& vtx, vector<TriVtx>& tris)
{
	vector< vector<int> > neighborVtx = getNeighborVtx(vtx, tris);
	vector< vector<int> > vtxTris = getVtxTris(vtx, tris);

	// �õ����б�
	vector< vector<int> > pairList;
	for(int i=0 ; i<int(vtx.size()) ; i++){
		for(int j=0 ; j<int(neighborVtx[i].size()) ; j++){
			if(neighborVtx[i][j] > i){
				vector<int> t;
				t.push_back(i);
				t.push_back(neighborVtx[i][j]);
				pairList.push_back(t);
			}
		}
	}

	// �жϱ��ڽ��������θ���������������ɾ��
	vector<int> delList(tris.size(), 0);
	for(int i=0 ; i<int(pairList.size()) ; i++){
		vector<int> t1 = vtxTris[pairList[i][0]];
		vector<int> t2 = vtxTris[pairList[i][1]];
		int s = 0;
		for(int m=0 ; m<int(t1.size()) ; m++){
			for(int n=0 ; n<int(t2.size()) ; n++){
				if( t1[m] == t2[n] ){
					if( delList[t1[m]]==0 ){
						s = s+1;
						if( s>2 ) {
							delList[t1[m]] = 1;
							//	printf("del #%i tri\n", t1[m]);
						}
					}
				}
			}
		}
	}

	vector<TriVtx> tris_old = tris;
	tris.clear();
	for(int i=0 ; i<int(tris_old.size()) ; i++){
		if( delList[i] == 0 ){
			tris.push_back(tris_old[i]);
		}
	}

};

// ������˹ƽ��
inline void laplacianSmooth(vector<Pnt3>& vtx, vector<TriVtx>& tris)
{
	cout << "smooth ..." << endl;
	vector< vector<int> > neighborVtx = getNeighborVtx(vtx, tris);

	vector<Pnt3> vtx_old = vtx;
	for( int i=0 ; i<int(vtx_old.size()) ; i++){
		vector<int> nv = neighborVtx[i];
		if(int(nv.size())==0)
			continue;

		Pnt3 t(0.0f, 0.0f, 0.0f);
		for( int j=0 ; j<int(nv.size()) ; j++ ) {
			t = t + vtx_old[nv[j]];
		}
		vtx[i] = t / float(nv.size());
	}

};

// ���벻ͬ����ͨ����
inline vector< vector<int> > divideIntoPiece(vector<Pnt3> &vtx, vector<TriVtx> &tris)
{
	vector< vector<int> > neighborVtx = getNeighborVtx(vtx, tris);

	// ��������Ƭ�����
	vector<int> vtxPieceIndex(vtx.size(),-1);

	int pieceIndex = 0;
	vector< vector<int> > Piece;
	while(1){
		Piece.push_back(vector<int>());

		// Ѱ����ͬ������ʼ��
		vector<int> cadidateList;
		for(int i=0 ; i<int(vtx.size()) ; i++){
			if(vtxPieceIndex[i] == -1){
				cadidateList.push_back(i);
				vtxPieceIndex[i] = pieceIndex;
				Piece[pieceIndex].push_back(i);
				break;
			}
		}
		if(cadidateList.empty())
			break;

		while(1){
			vector<int> neighbor = neighborVtx[cadidateList[0]];
			cadidateList.erase(cadidateList.begin());

			for(int j=0 ; j<int(neighbor.size()) ; j++){
				if(vtxPieceIndex[neighbor[j]] == -1){
					cadidateList.push_back(neighbor[j]);
					vtxPieceIndex[neighbor[j]] = pieceIndex;
					Piece[pieceIndex].push_back(neighbor[j]);
				}
			}

			if(cadidateList.empty())
				break;
		}

		pieceIndex ++;
	}

	return Piece;
};

// ɾ��С��Ƭ
inline void removeSmallPiece(vector<Pnt3> &vtx, vector<TriVtx> &tris, int numThreshold = 10)
{
	vector< vector<int> > Piece = divideIntoPiece(vtx, tris);

	vector<int> del_list(vtx.size(), 0);
	for(int i=0 ; i<int(Piece.size()) ; i++){
		if(int(Piece[i].size()) <= numThreshold){
			for(int j=0 ; j<int(Piece[i].size()) ; j++){
				del_list[Piece[i][j]] = 1;
			}
		}
	}
	
	delVtx(vtx, tris, del_list);
};


}

#endif