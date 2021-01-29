#ifndef GLVIEW_WIDGET_H
#define GLVIEW_WIDGET_H

#include "glview_widget_global.h"
#include "../types.h"

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector2D>
#include <QBasicTimer>
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QMutex>


class data_mem;

#define render_object_num 6
#define PI 3.141592653589793

/// 实现渲染3次

enum renderFaceType
{
	en_renderFace,
	en_renderTexture
};

struct RenderModel
{
	int type;		/// 模型的类型，是否带贴图
	data_mem *ptCloudMem;
	QOpenGLTexture *img_texture;

	std::vector<int> mSelectLineIndexs;
	std::vector<data_mem*> mLineDatas;
	RenderModel():type(0),ptCloudMem(NULL),img_texture(NULL){}
};

struct RenderPoint
{
	data_mem *points;
	data_mem *line;

	PtCloud clouds;
	RenderPoint():points(NULL),line(NULL){}

	void addPoint(const QPoint& point, const QMatrix4x4& projection,int w, int h);

	void draw(QOpenGLShaderProgram *program);

	void setTranst(const QMatrix4x4& m)
	{
		transtMatrix = m;
	}

	QMatrix4x4 transtMatrix;
};


//本模块最好能承担显示二维图，三维图，线条，面，可重入释放，可以多个界面并存的特点。
//class PtCloud;


class GLVIEW_WIDGETSHARED_EXPORT glview_widget:public QOpenGLWidget, protected QOpenGLFunctions
{

    Q_OBJECT
public:
    explicit glview_widget(QWidget *parent = 0);

    ~glview_widget();


#if 1
	void addLineData( PtCloud *cloud);

	void clear(void);

    void loadTextures();
    void setGroupId(int id);

	void setBox(omesh::Bbox& box);

	void setColor(const QVector3D& MaterialAmbient, const QVector3D& MaterialDiffuse, const QVector3D& MaterialSpecular);

	void addCoordinateData(void);

	void setModel( PtCloud* points, const QImage& img, bool hasImg = false);
	void setSelectLines(const std::vector<int>& indexs);

	void renderModel(void);
public slots:
    //below no scal work
    void viewFront();
	void viewLeft();
    void viewLeft60();
    void viewRight60();
    void viewUp60();
    void viewDown60();
	void viewReset();
	void setCoordinateVisible(bool visible = true);
	void viewByQQuaternion(const QQuaternion& qqua);

protected:
    void mousePressEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    //void timerEvent(QTimerEvent *e) override;
    void mouseMoveEvent(QMouseEvent *e) override;
    void wheelEvent(QWheelEvent *event) override;
    void mouseDoubleClickEvent(QMouseEvent *e) override;
    //init gl eviroment,load shader
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void initShaders();
private:
    QBasicTimer timer;
    QOpenGLShaderProgram programPtcloud;
    QOpenGLShaderProgram programLine;
    QOpenGLShaderProgram programTexture;
	QOpenGLShaderProgram programAxes;

	data_mem *mCoordinateMem[3];		/// 坐标轴的点

	omesh::Bbox mBox;
	QMutex mMutexLines;					/// 线条的数据锁


    QMatrix4x4 projection;

    QVector2D mousePressPosition;
    QVector3D rotationAxis;
    qreal angularSpeed;
    QQuaternion rotation;

private:
    //光照设置
    QVector3D uLightPosition;
    QVector3D uLightAmbient;
    QVector3D uLightDiffuse;
    QVector3D uLightSpecular;
    //模型材质参数
    QVector3D uMaterialAmbient;
    QVector3D uMaterialDiffuse;
    QVector3D uMaterialSpecular;
    //模型材质参数
    QVector3D uMaterialAmbient_line;
    QVector3D uMaterialDiffuse_line;
    QVector3D uMaterialSpecular_line;
    float uShininess;

    int mPressedTag;//0:right ,1:left,-1:midle
    //may be useless
    QVector2D viewShift;
    QVector2D viewShiftMutPara;

    omesh::Bbox box_forlinePts;
    float BOX_SCALL;//=1.2
    float PTVIEW_SCALL;// 2

    float scal_g;
    omesh::Pnt3 pt3_g;
    QVector2D diffAll;
    int groupId;

	/// 坐标轴是否显示
	bool mCoodinateVisible;

	QMatrix4x4 mMatriViewCoodinate;

	/// 渲染的模型
	RenderModel mRenderModel;

	RenderPoint mRenderPoint;

#endif
};

#endif // GLVIEW_WIDGET_H
