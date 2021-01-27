
#include "glview_widget.h"
//#include "newMesh/include/TriMesh.h"
#include "data_mem.h"

#include <QPushButton>
#include "QMouseEvent"
using namespace std;
//using namespace trimesh;
using namespace omesh;


glview_widget::glview_widget(QWidget *parent) :
     QOpenGLWidget(parent),mCoodinateVisible(false)
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


	if (mCoordinateMem)
	{
		delete mCoordinateMem[0];
		delete mCoordinateMem[1];
		delete mCoordinateMem[2];
	}
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


void glview_widget::addLineData( PtCloud *lines)
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
	mRenderModel.mLineDatas.push_back(l);
	mMutexLines.unlock();
	update();
	mCoodinateVisible = true;
}

void glview_widget::clear(void)
{
	delete mRenderModel.ptCloudMem;
	mRenderModel.ptCloudMem = NULL;
	for (int i = 0; i < mRenderModel.mLineDatas.size(); i++)
	{
		delete mRenderModel.mLineDatas[i];
	}
	mRenderModel.mLineDatas.clear();
}

//! [4]


void glview_widget::loadTextures()
{
#if 0
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

    update();
#endif
}

void glview_widget::setGroupId(int id){

    groupId=id;
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


void glview_widget::setModel( PtCloud* cloud, const QImage& img, bool hasImg /*= false*/)
{
	for (auto it = cloud->begin(); it != cloud->end(); it++)
	{
		mBox.add(it->x, it->y, it->z);
	}
	setBox(mBox);

	data_mem *&textureOpeMem = mRenderModel.ptCloudMem;
	if (textureOpeMem == NULL)
	{
		textureOpeMem = new data_mem(ePtData);
	}

	textureOpeMem->ptCloudDateUpdate(cloud);
	mRenderModel.type = en_renderFace;

	if (hasImg)
	{
		textureOpeMem->setDrawType(eTexture);
		mRenderModel.type = en_renderTexture;
		QOpenGLTexture *&img_texture = mRenderModel.img_texture;
		if (img_texture)
		{
			delete img_texture;
			img_texture = NULL;
		}

		if (!img_texture)
		{
			img_texture = new QOpenGLTexture(img);
			img_texture->setMinificationFilter(QOpenGLTexture::Nearest);
			img_texture->setMagnificationFilter(QOpenGLTexture::Linear);
			img_texture->setWrapMode(QOpenGLTexture::Repeat);
		}
		else
		{
			img_texture->setData(img);
		}
	}
	update();
}

void glview_widget::setSelectLines(const std::vector<int>& indexs)
{
	mRenderModel.mSelectLineIndexs = indexs;
	update();
}

void glview_widget::renderModel(void)
{
	/// 画线
	programPtcloud.bind();
	std::vector<data_mem*>& LineDatas = mRenderModel.mLineDatas;
	std::vector<int>& mSelectLineIndexs = mRenderModel.mSelectLineIndexs;
	int index = 0;
	for (int i = 0; i < LineDatas.size(); i++)
	{
		if (index < mSelectLineIndexs.size() && mSelectLineIndexs[index] == i)
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

	/// 画模型
	if (mRenderModel.type == en_renderFace)
	{
		programPtcloud.bind();
		programPtcloud.setUniformValue("renderType", 3);
		if (mRenderModel.ptCloudMem)
			mRenderModel.ptCloudMem->drawWithProgram(&programPtcloud);
	}
	else if (mRenderModel.type == en_renderTexture)
	{
		if (mRenderModel.img_texture)
			mRenderModel.img_texture->bind();
		programTexture.bind();
		programTexture.setUniformValue("texture", 0);
		if (mRenderModel.ptCloudMem)
			mRenderModel.ptCloudMem->drawWithProgram(&programTexture);
	}


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

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_MULTISAMPLE_ARB);

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	//glShadeModel(GL_SMOOTH);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    QMatrix4x4 matrix;


    matrix.scale( scal_g );
    matrix.rotate( rotation );//pt3_g

    QMatrix4x4 matrixInit, matrixCoord;
	if (groupId == 0) {
		matrixInit.rotate(-90, 1.0, 0.0, 0.0);
		matrixInit.rotate(180.0, 0.0, 0.0, 1.0);
	}

	matrixCoord = matrixInit;
    matrixInit.translate( -pt3_g[0], -pt3_g[1], -pt3_g[2]);

    QMatrix4x4 view;
    view.lookAt(QVector3D(0, 0, -10), QVector3D(/*viewShift[0]*0.01*/0,/*viewShift[1]*0.01*/0,0), QVector3D(0,1,0));

    QMatrix4x4 shitfMat;

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




	


	programTexture.bind();
	QMatrix4x4 mvMat;
	mvMat.translate(0.0f, 0.0f, -3.0f);

	//QMatrix4x4 matrixInit1;
	//matrixInit1.rotate(180.0, 0.0, 0.0, 1.0);
	//matrixInit1.translate(-pt3_g[0], -pt3_g[1], -pt3_g[2]);

	//QMatrix4x4 modelview1 = view * shitfMat*matrix*matrixInit1;

	//programTexture.setUniformValue("uModelViewProjection", projection * modelview1);
	programTexture.setUniformValue("uModelViewProjection", projection * modelview);

	programTexture.setUniformValue("texture", 0);
	renderModel();
#else
#endif
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


//! [3]

//! [4]


