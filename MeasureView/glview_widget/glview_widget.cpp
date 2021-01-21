
#include "glview_widget.h"
//#include "newMesh/include/TriMesh.h"
#include "data_mem.h"

#include <QPushButton>
#include "QMouseEvent"
using namespace std;
//using namespace trimesh;
using namespace omesh;


glview_widget::glview_widget(QWidget *parent) :
     QOpenGLWidget(parent),ptCloudMem(NULL),lineDataMem(NULL),\
     textureOpeMem(NULL),img_texture(NULL),mLineDataMem(NULL),mCoodinateVisible(false)
	, m_selectBall(0)
{

        //光照设置
        uLightPosition = QVector3D(0, 0, 100);
        uLightAmbient  = QVector3D(0.1, 0.1, 0.1);
        uLightDiffuse  = QVector3D(0.8, 0.8, 0.8);
        uLightSpecular = QVector3D(0.2, 0.2, 0.2);

        //模型材质参数 绘制点云//254 224 203
        uMaterialAmbient  = QVector3D(0.2, 0.2, 0.2);
        uMaterialDiffuse  = QVector3D(0.282*3.0, 0.224*3.0, 0.173*3.0);
        uMaterialSpecular = QVector3D(0.6,0.6, 0.6);
        uShininess = 68;

        //模型材质参数 绘制线条for line//254 224 203
        uMaterialAmbient_line  = QVector3D(0.7, 0.2, 0.2);
        uMaterialDiffuse_line  = QVector3D(0.282*3.0, 0.224*1.0, 0.173*1.0);
        uMaterialSpecular_line = QVector3D(0.9,0.4, 0.6);

        mPressedTag=1;

        scal_g=5.3;
        BOX_SCALL=1.2;
        PTVIEW_SCALL= 2;

        groupId=0;

		QSurfaceFormat surfaceFormat;
		surfaceFormat.setSamples(4);//多重采样
		setFormat(surfaceFormat);
}

glview_widget::~glview_widget(){

  delete ptCloudMem;
  delete lineDataMem;
    delete textureOpeMem;

	if (mLineDataMem)
		delete mLineDataMem;

	for (int i = 0; i < mLineDatas.size(); i ++)
	{
		delete mLineDatas[i];
	}
	mLineDatas.clear();

	if (mCoordinateMem)
	{
		delete mCoordinateMem[0];
		delete mCoordinateMem[1];
		delete mCoordinateMem[2];
	}
}

void glview_widget::setPtCloud(PtCloud* data,omesh::Bbox* box){

    points.clear();

    if (data->size()>0) {

#if 0
        while (!std::getline(inFile, s).eof())
        {
            point3d pt;
            char dat[100];
            strcpy(dat, s.c_str());
            sscanf(dat, "%f %f %f %f %f %f", &pt.x, &pt.y, &pt.z, &pt.nx, &pt.ny, &pt.nz);
            static int colorId = 0;
            pt.r = 0;
            pt.g = 100;
            pt.b = 40;
            points.push_back(pt);
            mBox.add(pt.x, pt.y, pt.z);

        }
#else

        if(data->size()){
            points.resize(data->size());
            memcpy(points.data(),data->data(),data->size()*sizeof(Vertex3Normal));
        }

#endif

        if(box){
			setBox(*box);
        }else{

            omesh::Bbox mBox;
            for(auto it=points.begin();it!=points.end();it++){

                mBox.add(it->x,it->y,it->z);
            }
			setBox(mBox);
        }

		rotation = QQuaternion();
		viewShift = QVector2D(0, 0);

#if 1

        if(ptCloudMem)
        {
			makeCurrent();
            ptCloudMem->ptCloudDateUpdate(&points);
            update();
			mCoodinateVisible = true;
        }
#endif
    }
}

void glview_widget::setPtCloud(PtCloud* data, omesh::Bbox* box, int pos)
{
	if (pos >= render_object_num || pos < 0)
		return;

	data_mem*& pCloudMem = mRenderObject[pos].ptCloudMem;
	if (pCloudMem == NULL)
	{
		pCloudMem = new data_mem(ePtData);
	}	

	if (data->size() > 0) {

		if (box) {
			setBox(*box);
		}
		else {

			omesh::Bbox mBox;
			for (auto it = data->begin(); it != data->end(); it++) {
				mBox.add(it->x, it->y, it->z);
			}
			setBox(mBox);
		}

		rotation = QQuaternion();
		viewShift = QVector2D(0, 0);

		if (pCloudMem)
		{
			makeCurrent();
			pCloudMem->ptCloudDateUpdate(data);
			update();
			mCoodinateVisible = true;
		}

	}
}

void glview_widget::setPtLines(PtCloud* data, omesh::Bbox* box)
{
	if (data->size() > 0) 
	{		
		omesh::Bbox Box;
		if (box)
		{
		}
		else 
		{			
			for (auto it = data->begin(); it != data->end(); it++)
			{
				Box.add(it->x, it->y, it->z);
			}
			box = &Box;
		}		
		setBox(*box);
		if (mLineDataMem)
		{
			makeCurrent();
			mLineDataMem->ptCloudDateUpdate(data);
			update();
			mCoodinateVisible = true;
		}
	}
}

void glview_widget::addLineData(PtCloud* lines)
{
	//if (!lines->empty())
	{		
		for (auto it = lines->begin(); it != lines->end(); it++)
		{
			mBox.add(it->x, it->y, it->z);
		}
		
		setBox(mBox);

		data_mem *l = new data_mem(ePtLinesData);
		makeCurrent();
		l->ptCloudDateUpdate(lines);
		mMutexLines.lock();
		mLineDatas.push_back(l);
		mMutexLines.unlock();
		update();
		mCoodinateVisible = true;
	}
}

void glview_widget::addLineData(PtCloud* lines, int pos)
{
	if (pos >= render_object_num || pos < 0)
		return;
	//if (!lines->empty())
	{
		std::vector<data_mem*>& total= mRenderObject[pos].mLineDatas;

		for (auto it = lines->begin(); it != lines->end(); it++)
		{
			mBox.add(it->x, it->y, it->z);
		}

		setBox(mBox);

		data_mem *l = new data_mem(ePtLinesData);
		makeCurrent();
		l->ptCloudDateUpdate(lines);
		total.push_back(l);
		update();
		mCoodinateVisible = true;
	}

}

void glview_widget::setPtCloudVisible(bool visible, int pos)
{
	if (pos >= render_object_num || pos < 0)
		return;
	mRenderObject[pos].visible = visible;
	update();
}

void glview_widget::setBox(omesh::Bbox& box)
{
	auto cent = box.center();
	auto dis = box.diag();
	pt3_g = cent;

	scal_g = 1.0 / dis;
	scal_g *= 2;
	setProperty("scal_g", scal_g);
	float sc = PTVIEW_SCALL / dis;// (15 / dis);
	box_forlinePts = box;
}

void glview_widget::addPtCloud(PtCloud* data, omesh::Bbox* box)
{
	if (data->size() > 0) 
	{
		points.insert(points.end(), data->begin(), data->end());
		if (box) 
		{
			box_forlinePts.add(*box);
		}
		else {

			omesh::Bbox mBox;
			for (auto it = points.begin(); it != points.end(); it++) {

				mBox.add(it->x, it->y, it->z);
			}

			auto cent = mBox.center();
			auto dis = mBox.diag();
			pt3_g = cent;

			scal_g = 1.0 / dis;
			scal_g *= 1.5;
			float sc = PTVIEW_SCALL / dis;// (15 / dis);
			setProperty("scal_g", scal_g);

										  //box of line same as points'box
			box_forlinePts = mBox;
		}
		if (ptCloudMem)
		{
			ptCloudMem->ptCloudDateUpdate(&points);
			update();
			mCoodinateVisible = true;
		}
	}
}


void glview_widget::clear(void)
{
	points.resize(1);
	if (ptCloudMem)
	{
		ptCloudMem->ptCloudDateUpdate(&points);
		update();
	}

	for (int i = 0; i < mLineDatas.size(); i ++)
	{
		delete mLineDatas[i];
	}
	mLineDatas.clear();
	mCoodinateVisible = false;
}

//! [4]

void glview_widget::setLinePts(PtCloud* data){

    line_pts.clear();

    if (data->size()>0) {

        line_pts=*data;
    }
}

void glview_widget::loadPly(QString xyzFileName){

    omesh::Bbox mBox;
    points.clear();
    string s;
    ifstream inFile(xyzFileName.toStdString().c_str());

    if (inFile.is_open()) {

        while (!std::getline(inFile, s).eof())
        {
             Vertex3Normal  pt;
            char dat[100];
            strcpy(dat, s.c_str());
            sscanf(dat, "%f %f %f %f %f %f", &pt.x, &pt.y, &pt.z, &pt.nx, &pt.ny, &pt.nz);
           // static int colorId = 0;

            //float va=rand()%10;
//            pt.x=va*0.1;
//            pt.y=va*0.1;
//            pt.z=va*0.1;
            pt.nx=0.8;
            pt.ny=0.3;
            pt.nz=0.3;
            points.push_back(pt);
            mBox.add(pt.x, pt.y, pt.z);

        }

        auto cent = mBox.center();
        auto dis =  mBox.diag();
        pt3_g=cent;
        mBox.makeCube(BOX_SCALL);
        //g_Box=mBox;

        float sc = PTVIEW_SCALL / dis;// (15 / dis);
        //相对原点放缩
        //initMatrix.scale(sc, sc, sc);

        //initMatrix.translate(-(cent[0] * sc),\
         //   - (cent[1] * sc), -(cent[2] * sc));

        //mGlType = eGlPoint;

       // transMatrixBk = transMatrix;
       // initMatrixBk = initMatrix;
        //sc_g=sc;
       // ce_x= -(cent[0] * sc);
       // ce_y= -(cent[1] * sc);
       // ce_z= -(cent[2] * sc);

        //int sz=points.size();
        if(ptCloudMem)
        {
            ptCloudMem->ptCloudDateUpdate(&points);
            update();
        }
    }

}

void glview_widget::loadTextures()
{
    // Load cube.png image
    if(!img_texture)
        img_texture = new QOpenGLTexture(QImage(":/cube.png").mirrored());
    else{

        img_texture->setData(QImage(":/cube.png"));
    }
    // Set nearest filtering mode for texture minification
    img_texture->setMinificationFilter(QOpenGLTexture::Nearest);

    // Set bilinear filtering mode for texture magnification
    img_texture->setMagnificationFilter(QOpenGLTexture::Linear);

    // Wrap texture coordinates by repeating
    // f.ex. texture coordinate (1.1, 1.2) is same as (0.1, 0.2)
    img_texture->setWrapMode(QOpenGLTexture::Repeat);

    //trig fbo,
    //programTexture.imgAxiUpdate();
    textureOpeMem->imgAxiUpdate();
    update();
}

void glview_widget::setGroupId(int id){

    groupId=id;
}

void glview_widget::setSelectLines(std::vector<int>& select)
{
	mSelectLineIndexs.swap(select);
	update();
}

void glview_widget::setSelectLines(std::vector<int>& select, int pos)
{
	if (pos < 0 && pos > render_object_num)
		return;
	mRenderObject[pos].mSelectLineIndexs = select;
}

void glview_widget::setColor(const QVector3D& MaterialAmbient, const QVector3D& MaterialDiffuse, const QVector3D& MaterialSpecular)
{
	uMaterialAmbient = MaterialAmbient;
	uMaterialDiffuse = MaterialDiffuse;
	uMaterialSpecular = uLightSpecular;
}


void glview_widget::addCoordinateData(void)
{
	PtCloud lx, lx1, lx2;
	PtCloud* lines = &lx;

	int size = 100;
	float scal_arrow = 1.2f;

	/// x
	lx.push_back(Vertex3Normal(-4.0*size, 0.0f, 0.0f,0, 0, 1));
	lx.push_back(Vertex3Normal(4.0*size, 0.0f, 0.0f, 0, 0, 1));

	// arrow
	lx.push_back(Vertex3Normal(4.0*size, 0.0f, 0.0f, 0, 0, 1));
	lx.push_back(Vertex3Normal(3.0*size*scal_arrow, 1.0f*size/ scal_arrow/2, 0.0f, 0, 0, 1));

	lx.push_back(Vertex3Normal(4.0*size, 0.0f, 0.0f, 0, 0, 1));
	lx.push_back(Vertex3Normal(3.0*size*scal_arrow, -1.0f*size / scal_arrow/2, 0.0f, 0, 0, 1));


	// y 
	lx1.push_back(Vertex3Normal(0.0, -4.0f*size, 0.0f, 0, 0, 1));
	lx1.push_back(Vertex3Normal(0.0, 4.0f*size, 0.0f, 0, 0, 1));

	// arrow
	lx1.push_back(Vertex3Normal(0.0, 4.0f*size, 0.0f, 0, 0, 1));
	lx1.push_back(Vertex3Normal(1.0*size / scal_arrow/2, 3.0f*size*scal_arrow, 0.0f, 0, 0, 1));

	lx1.push_back(Vertex3Normal(0.0, 4.0f*size, 0.0f, 0, 0, 1));
	lx1.push_back(Vertex3Normal(-1.0*size / scal_arrow/2, 3.0f*size*scal_arrow, 0.0f, 0, 0, 1));

	// z 
	//lx.push_back(Vertex3Normal(0.0, 0.0f, -4.0f*size, 0, 0, 1));
	lx2.push_back(Vertex3Normal(0.0, 0.0f, 0, 0, 0, 1));
	lx2.push_back(Vertex3Normal(0.0, 0.0f, 4.0f*size, 0, 0, 1));

	// arrow
	lx2.push_back(Vertex3Normal(0.0, 0.0f, 4.0f*size, 0, 0, 1));
	lx2.push_back(Vertex3Normal(0.0, 1.0f*size / scal_arrow/2, 3.0f*size*scal_arrow, 0, 0, 1));

	lx2.push_back(Vertex3Normal(0.0, 0.0f, 4.0f*size, 0, 0, 1));
	lx2.push_back(Vertex3Normal(0.0, -1.0f*size / scal_arrow/2, 3.0f*size*scal_arrow, 0, 0, 1));

	mCoordinateMem[0] = new data_mem(ePtLinesData);
	mCoordinateMem[0]->ptCloudDateUpdate(lines);

	mCoordinateMem[1] = new data_mem(ePtLinesData);
	mCoordinateMem[1]->ptCloudDateUpdate(&lx1);

	mCoordinateMem[2] = new data_mem(ePtLinesData);
	mCoordinateMem[2]->ptCloudDateUpdate(&lx2);
	update();
}

void glview_widget::clearPos(int pos /*= -1*/)
{
	memData tempData;
	tempData.resize(1);

	auto f = [=](int index, memData* p)
	{	
		if (mRenderObject[index].ptCloudMem)
		mRenderObject[index].ptCloudMem->ptCloudDateUpdate(p);
		for (int i = 0; i < mRenderObject[index].mLineDatas.size(); i++)
		{
			delete mRenderObject[index].mLineDatas[i];
		}
		mRenderObject[index].mLineDatas.clear();
		mCoodinateVisible = false;

	};
	if (pos >= 0 && pos < render_object_num)
	{
		f(pos, &tempData);
		auto iter = m_ballMap.find(pos);
		if (iter != m_ballMap.end())
		{
			vector<data_mem*> &mm = iter->second;
			for (size_t i = 0; i < mm.size(); i++)
			{
				delete mm[i];
			}
			m_ballMap.erase(iter);
		}
		update();
	}
	else if (pos == -1)
	{
		for (int i = 0; i < render_object_num; i ++)
		{
			f(i, &tempData);
		}
		for (auto iter = m_ballMap.begin(); iter !=m_ballMap.end(); iter++)
		{
			vector<data_mem*> &mm = iter->second;
			for (size_t i = 0; i < mm.size(); i++)
			{
				delete mm[i];
			}
		}
		m_ballMap.clear();
		update();
	}

	
		
}

void glview_widget::setBalls(vector<Ball> &balls, int pos)
{
// 	if (balls.size() == m_ballMems.size())
// 	{
// 		for (int i = 0; i < balls.size(); i++)
// 		{
// 			m_ballMems[i]->ptCloudDateUpdate(&balls[i].points);
// 		}
// 	}
// 	else 
// 	{
// 		for (size_t i = 0; i < m_ballMems.size(); i++)
// 		{
// 			delete m_ballMems[i];
// 		}
// 		m_ballMems.clear();
// 
// 		for (size_t i = 0; i < balls.size(); i++)
// 		{
// 			data_mem *dm = new data_mem(ePoints);
// 			dm->ptCloudDateUpdate(&balls[i].points);
// 			m_ballMems.push_back(dm);
// 		}
// 	}

	vector<data_mem*> ballMem;
	for (size_t i = 0; i < balls.size(); i++)
	{
		data_mem *dm = new data_mem(ePoints);
		dm->ptCloudDateUpdate(&balls[i].points);
		ballMem.push_back(dm);
	}
	m_ballMap[pos] = ballMem;
}

void glview_widget::setSelectBall(int index)
{
	m_selectBall = index;
	update();
}

void glview_widget::viewFront() {

    rotation=QQuaternion();
	diffAll = QVector2D(0, 0);
    update();
}

void glview_widget::viewLeft()
{
	rotation = QQuaternion::fromAxisAndAngle(0.0, 1.0, 0.0, 90);
	diffAll = QVector2D(0, 0);
	update();
}

void glview_widget::viewLeft60() {

   rotation=QQuaternion::fromAxisAndAngle(0.0,1.0,0.0, -40);
}

void glview_widget::viewRight60(){

   rotation=QQuaternion::fromAxisAndAngle(0.0,1.0,0.0, 40);
   update();
}

void glview_widget::viewUp60(){

    rotation=QQuaternion::fromAxisAndAngle(1.0,0.0,0.0, 40);
    update();
}

void glview_widget::viewDown60(){

     rotation=QQuaternion::fromAxisAndAngle(1.0,0.0,0.0,-40);
     update();
}


void glview_widget::viewReset()
{
	rotation = QQuaternion();
	viewShift = QVector2D(0, 0);
	diffAll = QVector2D(0, 0);
	scal_g = property("scal_g").toFloat();
	update();
}

void glview_widget::setCoordinateVisible(bool visible /*= true*/)
{
	mCoodinateVisible = visible;
	update();
}

void glview_widget::viewByQQuaternion(const QQuaternion& qqua)
{
	rotation = qqua;
	diffAll = QVector2D(0, 0);
	update();
}

void glview_widget::mousePressEvent(QMouseEvent *e)
{
    QOpenGLWidget::mousePressEvent(e);

    if(e->button()==Qt::LeftButton)
    {
        mPressedTag=1;

    }else if(e->button()==Qt::RightButton)
    {
        mPressedTag=0;

    }else if(e->button()==Qt::MidButton){
        mPressedTag=-1;
    }
    mousePressPosition=QVector2D(e->localPos());
    return;
    // Save mouse press position
    mousePressPosition = QVector2D(e->localPos());
}

void glview_widget::mouseReleaseEvent(QMouseEvent *e)
{
    QOpenGLWidget::mouseReleaseEvent(e);
    return;

}

void glview_widget::mouseMoveEvent(QMouseEvent *e)
{

    if(mPressedTag==1)
    {

#if 0

        auto NowMousePos=QVector2D(e->localPos());

        QVector2D diff = NowMousePos - mousePressPosition;

        // Rotation axis is perpendicular to the mouse position difference
        // vector
        QVector3D n = QVector3D(-diff.y(), diff.x(), 0.0).normalized();

        // Accelerate angular speed relative to the length of the mouse sweep
        qreal acc = diff.length() ;//diff.length() / 100.0;

        // Calculate new rotation axis as weighted sum
        rotationAxis = (rotationAxis * angularSpeed + n * acc).normalized();
        angularSpeed = acc;

        rotation = QQuaternion::fromAxisAndAngle(rotationAxis, angularSpeed) * rotation;

        // Request an update
        update();

        mousePressPosition=NowMousePos;
#else
        auto ps=QCursor::pos();
        QVector2D NowMousePos(ps);

        QVector2D diff =(NowMousePos - mousePressPosition);
        float xall=diffAll.x()+diff.x();
        float yall=diffAll.y()+diff.y();

        //QVector2D diff =NowMousePos - mousePressPosition;
        //float xall=diffAll.x()+diff.x();
        //float yall=diffAll.y()+diff.y();

        mousePressPosition=NowMousePos;

        if(diff.x()>100) return;
        if(diff.y()>100) return;
        diffAll.setX(xall);
        diffAll.setY(yall);
        auto rotationTemp0 = QQuaternion::fromAxisAndAngle(0.0,1.0,0.0,xall);
        auto rotationTemp1 = QQuaternion::fromAxisAndAngle(1.0,0.0,0.0,-yall);
        rotation=rotationTemp1*rotationTemp0;

		//mMatriViewCoodinate.rotate(rotation);

        update();
#endif

    }
    else if(mPressedTag==0)
    {
        auto NowMousePos=QVector2D(e->localPos());

        QVector2D diff = NowMousePos - mousePressPosition;

        viewShift+=diff;
        // Request an update
        update();

        mousePressPosition=NowMousePos;
    }
}

 void glview_widget::wheelEvent(QWheelEvent *event)
 {
    int ang=event->delta()/120;
   // ang++;
    scal_g*=(1+ang*0.05);
    update();
 }

 void glview_widget::mouseDoubleClickEvent(QMouseEvent *e){

     viewFront();
 }
//! [0]

//! [1]
//void es_gl_core::timerEvent(QTimerEvent *)
//{
//    // Decrease angular speed (friction)
//    angularSpeed *= 0.59;

//    // Stop rotation when speed goes below threshold
//    if (angularSpeed < 0.01) {
//        angularSpeed = 0.0;
//    } else {
//        // Update rotation
//        rotation = QQuaternion::fromAxisAndAngle(rotationAxis, angularSpeed) * rotation;

//        // Request an update
//        update();
//    }
//}
//! [1]


 void glview_widget::initializeGL()
 {
     #if 1
     initializeOpenGLFunctions();

     GLenum err = glGetError();
     printf(" initializeOpenGLFunctions err code:%d\n",err);

     //glClearColor(0.5, 0.5, 0.6, 1);

	 /// 浅蓝背景
	 glClearColor(0.749f, 0.835f, 0.96, 1);
     err = glGetError();
     printf("glClearColor err code:%d\n",err);

     initShaders();
     err = glGetError();
     printf("initShaders err code:%d\n",err);
     //initTextures();

 //! [2]
     // Enable depth buffer
     glEnable(GL_DEPTH_TEST);
     glFrontFace(GL_CCW);
     glCullFace(GL_BACK);

     // Enable back face culling
     glEnable(GL_CULL_FACE);
 //! [2]

     //glPolygonMode(GL_FRONT, GL_FILL);

    // glPointSize(3.0);
 #if  1
     err = glGetError();
     printf("GeometryEngPt begin err code:%d\n",err);

     //create point cloud or line mem,so can update some time
     ptCloudMem = new data_mem(ePtData);
     lineDataMem = new data_mem(eLineData);
     textureOpeMem=new data_mem(eTexture);
	 mLineDataMem = new data_mem(ePtLinesData);
     err = glGetError();
    printf("GeometryEngPt end err code:%d\n",err);

 #endif
     //axiDraw=new AxiDraw();

     // Use QBasicTimer because its faster than QTimer
     //timer.start(12, this);
 //loadPly();
     //
    #endif

	addCoordinateData();
 }

//! [5]
void glview_widget::resizeGL(int w, int h)
{
    #if 1
    // Calculate aspect ratio
    qreal aspect = qreal(w) / qreal(h ? h : 1);

    // Set near plane to 3.0, far plane to 7.0, field of view 45 degrees
    const qreal zNear = 3.0, zFar = 170.0, fov = 45.0;

    // Reset projection
    projection.setToIdentity();

    // Set perspective projection
    //projection.perspective(fov, aspect, zNear, zFar);
    projection.ortho(-aspect,aspect,-1.0,1.0,-100,100.0);

    viewShiftMutPara[0]=2.0*((float)w/(float)h)/((float)(w?w:1));
    viewShiftMutPara[1]=2.0*(1.0)/((float)(h?h:1));

    #else


    // Calculate aspect ratio
    qreal aspect = qreal(w) / qreal(h ? h : 1);

    // Set near plane to 3.0, far plane to 7.0, field of view 45 degrees
    const qreal zNear = 3.0, zFar = 7.0, fov = 45.0;

    // Reset projection
    projection.setToIdentity();

    // Set perspective projection
    projection.perspective(fov, aspect, zNear, zFar);
    #endif
}
//! [5]

void glview_widget::paintGL()
{
    //mingwe will broken
    //if(!points.size()) return;

    makeCurrent();
    #if 1
    //GLenum err = glGetError();
    //printf("err code:%d\n",err);
    // Clear color and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_MULTISAMPLE_ARB);

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	//glShadeModel(GL_SMOOTH);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //texture->bind();

    //! [6]
    // Calculate model view transformation
    QMatrix4x4 matrix;

    // initMatrix.x;
    // matrix.translate(0.0, 0.0, -0.0);

    //matrix.translate( ce_x, ce_y, ce_z);
    matrix.scale( scal_g );
    matrix.rotate( rotation );//pt3_g

    QMatrix4x4 matrixInit, matrixCoord;
	if (groupId == 0) {
		matrixInit.rotate(-90, 1.0, 0.0, 0.0);
		matrixInit.rotate(180.0, 0.0, 0.0, 1.0);
	}
	//else if (groupId == 1) {
	//	matrixInit.rotate(-180.0, 1.0, 0.0, 0.0);
	//	matrixInit.rotate(90.0, 0.0, 0.0, 1.0);
	//}
	matrixCoord = matrixInit;
    matrixInit.translate( -pt3_g[0], -pt3_g[1], -pt3_g[2]);

    QMatrix4x4 view;
    view.lookAt(QVector3D(0, 0, -10), QVector3D(/*viewShift[0]*0.01*/0,/*viewShift[1]*0.01*/0,0), QVector3D(0,1,0));

    QMatrix4x4 shitfMat;

    //float wit=this->width();
    //float hit=this->height();

    //float hor=2.0*(wit/hit)*(((float)viewShift[0])/wit);
    //float ver=2.0*(1)*(((float)viewShift[1])/hit);
    shitfMat.translate( -viewShift[0]*viewShiftMutPara[0], -viewShift[1]*viewShiftMutPara[1], 0);

    QMatrix4x4 modelview = view * shitfMat*matrix * matrixInit;

    QMatrix3x3 normalRot = modelview.normalMatrix();


	if (mCoordinateMem && mCoodinateVisible)
	{
		programAxes.bind();

		programAxes.setUniformValue("uModelView", modelview);
		programAxes.setUniformValue("uNormalRotation", normalRot);
		programAxes.setUniformValue("uModelViewProjection", projection*modelview /*projection * modelview*/);

		QMatrix4x4 matrix;
		matrix.rotate(-90, 1.0, 0.0, 0.0);
		matrix.rotate(180.0, 0.0, 0.0, 1.0);
		matrix.rotate(rotation);

		matrix.scale(100);
		//trant.translate(-0.9, -0.9, 0);
		//programPtcloud.setUniformValue("uModelViewProjection", trant* matrix /*projection * modelview*/);

		programAxes.setUniformValue("renderType", 3);
		programAxes.setUniformValue("inColor", QVector4D(0, 0, 1, 1));
		glLineWidth(2.5);
		mCoordinateMem[0]->drawWithProgram(&programAxes);

		programAxes.setUniformValue("renderType", 3);
		programAxes.setUniformValue("inColor", QVector4D(1, 0, 0, 1));
		mCoordinateMem[1]->drawWithProgram(&programAxes);

		programAxes.setUniformValue("renderType", 3);
		programAxes.setUniformValue("inColor", QVector4D(1, 1, 0, 1));
		mCoordinateMem[2]->drawWithProgram(&programAxes);
	}



    programPtcloud.bind();

    // Set modelview-projection matrix
    programPtcloud.setUniformValue("uModelView", modelview);
    programPtcloud.setUniformValue("uNormalRotation", normalRot);
    programPtcloud.setUniformValue("uModelViewProjection",projection*modelview /*projection * modelview*/);

    programPtcloud.setUniformValue("uLightPosition", uLightPosition);
    programPtcloud.setUniformValue("uLightAmbient", uLightAmbient);
    programPtcloud.setUniformValue("uLightDiffuse", uLightDiffuse);
    programPtcloud.setUniformValue("uLightSpecular", uLightSpecular);

    programPtcloud.setUniformValue("uMaterialAmbient", uMaterialAmbient);
    programPtcloud.setUniformValue("uMaterialDiffuse", uMaterialDiffuse);
    programPtcloud.setUniformValue("uMaterialSpecular", uMaterialSpecular);
    programPtcloud.setUniformValue("uShininess", uShininess);
    //! [6]

    // Use texture unit 0 which contains cube.png
    //program.setUniformValue("texture", 0);

    // Draw cube geometry


	int index = 0;
	for (int i = 0; i < mLineDatas.size(); i ++)
	{		
		if (index<mSelectLineIndexs.size()&& mSelectLineIndexs[index] == i)
		{			
			/// 设置选中的颜色
			programPtcloud.setUniformValue("renderType", 2);
			index++;
			glLineWidth(4);
			mLineDatas[i]->drawWithProgram(&programPtcloud);
			
		}
		else
		{
			continue;
			glLineWidth(2);
			programPtcloud.setUniformValue("renderType", 1);
			mLineDatas[i]->drawWithProgram(&programPtcloud);
		}		
	}

	if (ptCloudMem)
	{
		int type = mLineDatas.empty() ? 0 : 3;
		programPtcloud.setUniformValue("renderType", type);
		if (type == 3)
			ptCloudMem->drawWithProgram(&programPtcloud);
		else
			ptCloudMem->drawWithProgramPoints(&programPtcloud);
	}

	if (mLineDataMem)
		mLineDataMem->drawWithProgram(&programPtcloud);

//	renderBalls(TODO);

	renderSubSence();


    #endif

    //img texure
    #if 1

    if(img_texture)
    {

        //for img test
        // Clear color and depth buffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        img_texture->bind();

        programTexture.bind();
        //! [6]
        // Calculate model view transformation
        QMatrix4x4 matrix;
        matrix.translate(0.0, 0.0, -5.0);
        matrix.rotate(rotation);

        // Set modelview-projection matrix
        programTexture.setUniformValue("mvp_matrix", projection /** matrix*/);
        //! [6]

        // Use texture unit 0 which contains cube.png
        programTexture.setUniformValue("texture", 0);

        // Draw cube geometry
        //programTexture->drawWithProgram(&program);
        if(textureOpeMem){
            textureOpeMem->drawWithProgram(&programTexture);
         }
    }
    #endif		
}


void glview_widget::renderBalls(int pos)
{
// 	for (size_t i = 0; i < m_ballMems.size(); i++)
// 	{
// 		if (i != m_selectBall)
// 		{
// 			programPtcloud.setUniformValue("renderType", 4);
// 		}
// 		else
// 			programPtcloud.setUniformValue("renderType", 5);
// 		
// 		m_ballMems[i]->drawWithProgram(&programPtcloud);
// 	}
	auto iter = m_ballMap.find(pos);
	if (iter != m_ballMap.end())
	{
		vector<data_mem*> &mems = iter->second;
		for (size_t i = 0; i < mems.size(); i++)
		{
			if (i != m_selectBall)
			{
				programPtcloud.setUniformValue("renderType", 4);
			}
			else
				programPtcloud.setUniformValue("renderType", 5);

			mems[i]->drawWithProgram(&programPtcloud);
		}
	}
}

//! [3]
void glview_widget::initShaders()
{
    // Compile vertex shader



    //img texture
    {

            if (!programTexture.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/shader/vshader.glsl"))
            //if (!program.addShaderFromSourceCode(QOpenGLShader::Vertex, shaderVP))
            {
                QString errs = programTexture.log();

                GLenum err = glGetError();
                printf(" vp.glsl err code:%d,%s\n",err,errs.toStdString().c_str());
                exit(-1);
            };

            if (!programTexture.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/shader/fshader.glsl"))
            //if (!program.addShaderFromSourceCode(QOpenGLShader::Vertex, shaderVP))
            {
                QString errs = programTexture.log();

                GLenum err = glGetError();
                printf(" vp.glsl err code:%d,%s\n",err,errs.toStdString().c_str());
                exit(-1);
            };

            // Link shader pipeline
            if (!programTexture.link())
            {
                QString errs = programPtcloud.log();
                GLenum err = glGetError();
                printf(" glsl link err code:%d ,%s\n",err,errs.toStdString().c_str());
                  exit(-1);
            };

             //errs = programPtcloud.log();
    }
    //point cloud
    {
        if (!programPtcloud.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/shader/vp_pt.glsl"))
        //if (!program.addShaderFromSourceCode(QOpenGLShader::Vertex, shaderVP))
        {
            QString errs = programPtcloud.log();

            GLenum err = glGetError();
            printf(" vp.glsl err code:%d,%s\n",err,errs.toStdString().c_str());
            exit(-1);
        };

        QString errs = programPtcloud.log();
        // Compile fragment shader
        if (!programPtcloud.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/shader/vf_pt.glsl"))
        //if (!program.addShaderFromSourceCode(QOpenGLShader::Fragment, shaderFP))
        {
            QString errs = programPtcloud.log();
            GLenum err = glGetError();
            printf(" vf.glsl err code:%d ,%s\n",err,errs.toStdString().c_str());
            exit(-1);
        };

         errs = programPtcloud.log();


        // Link shader pipeline
        if (!programPtcloud.link())
        {
            QString errs = programPtcloud.log();
            GLenum err = glGetError();
            printf(" glsl link err code:%d ,%s\n",err,errs.toStdString().c_str());
              exit(-1);
        };

         errs = programPtcloud.log();

        // Bind shader pipeline for use
        if (!programPtcloud.bind())
        {
            QString errs = programPtcloud.log();
            GLenum err = glGetError();
            printf(" glsl bind err code:%d ,%s\n",err,errs.toStdString().c_str());
              exit(-1);
        };

         errs = programPtcloud.log();
    }

	/// axes
	{
		if (!programAxes.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/shader/vp_pt_axes.glsl"))
		{
			QString errs = programAxes.log();

			GLenum err = glGetError();
			printf(" vp.glsl err code:%d,%s\n", err, errs.toStdString().c_str());
			exit(-1);
		};

		QString errs = programAxes.log();
		// Compile fragment shader
		if (!programAxes.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/shader/vf_pt_axes.glsl"))
		{
			QString errs = programAxes.log();
			GLenum err = glGetError();
			printf(" vf.glsl err code:%d ,%s\n", err, errs.toStdString().c_str());
			exit(-1);
		};
		errs = programAxes.log();

		// Link shader pipeline
		if (!programAxes.link())
		{
			QString errs = programAxes.log();
			GLenum err = glGetError();
			printf(" glsl link err code:%d ,%s\n", err, errs.toStdString().c_str());
			exit(-1);
		};

		errs = programAxes.log();

		// Bind shader pipeline for use
		if (!programAxes.bind())
		{
			QString errs = programAxes.log();
			GLenum err = glGetError();
			printf(" glsl bind err code:%d ,%s\n", err, errs.toStdString().c_str());
			exit(-1);
		};
		errs = programAxes.log();
	}


     //line shader
     {
         if (!programLine.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/shader/vp_line.glsl"))
         //if (!program.addShaderFromSourceCode(QOpenGLShader::Vertex, shaderVP))
         {
             QString errs = programLine.log();

             GLenum err = glGetError();
             printf(" vp.glsl err code:%d,%s\n",err,errs.toStdString().c_str());
             exit(-1);
         };

         QString errs = programLine.log();

         // Compile fragment shader
         if (!programLine.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/shader/vf_line.glsl"))
         //if (!program.addShaderFromSourceCode(QOpenGLShader::Fragment, shaderFP))
         {
             QString errs = programLine.log();
             GLenum err = glGetError();
             printf(" vf.glsl err code:%d ,%s\n",err,errs.toStdString().c_str());
             exit(-1);
         };

          errs = programLine.log();


         // Link shader pipeline
         if (!programLine.link())
         {
             QString errs = programLine.log();
             GLenum err = glGetError();
             printf(" glsl link err code:%d ,%s\n",err,errs.toStdString().c_str());
               exit(-1);
         };

          errs = programLine.log();
     }
}

void glview_widget::renderSubSence(void)
{
	auto f = [=](int pos) 
	{	
		if (mRenderObject[pos].visible == false)		
			return;
		
		std::vector<data_mem*>& LineDatas = mRenderObject[pos].mLineDatas;
		data_mem *pCloudMem = mRenderObject[pos].ptCloudMem;
		int index = 0;
	for (int i = 0; i < LineDatas.size(); i++)
	{
		if (index < mRenderObject[pos].mSelectLineIndexs.size() && mRenderObject[pos].mSelectLineIndexs[index] == i)
		{
			/// 设置选中的颜色
			programPtcloud.setUniformValue("renderType", 2);
			index++;
			glLineWidth(4);
			LineDatas[i]->drawWithProgram(&programPtcloud);

		}
		else
		{
			continue;
			glLineWidth(2);
			programPtcloud.setUniformValue("renderType", 1);
			LineDatas[i]->drawWithProgram(&programPtcloud);
		}
	}

	if (pCloudMem)
	{
		int type = LineDatas.empty() ? 0 : 3;
		programPtcloud.setUniformValue("renderType", type);
		if (type == 3)
			pCloudMem->drawWithProgram(&programPtcloud);
		else
			pCloudMem->drawWithProgramPoints(&programPtcloud);
	}
	renderBalls(pos);
	//if (mLineDataMem)
	//	mLineDataMem->drawWithProgram(&programPtcloud);
	};

	for (int i = 0; i < render_object_num; i ++)
	{
		f(i);
	}

}

//! [3]

//! [4]


