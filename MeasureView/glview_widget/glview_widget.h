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
struct RenderObject
{
	bool visible;
	data_mem *ptCloudMem;
	std::vector<int> mSelectLineIndexs;
	std::vector<data_mem*> mLineDatas;

	RenderObject():visible(false),ptCloudMem(NULL){}
};

struct Ball
{
	bool visible;
	float radius;
	PtCloud points;
	std::vector<QVector3D> indices;
	
	Ball(float centerx, float centery, float centerz,int radius = 5) :visible(false)
	{
		int stacks = 100;
		int Slices = 50;
		Vertex3Normal point;
// 		point.r = 255;
// 		point.g = point.b = 0;
		for (int i = 0; i <= stacks; ++i) {

			float V = i / (float)stacks;
			float phi = V * PI;

			// Loop Through Slices
			for (int j = 0; j <= Slices; ++j) {

				float U = j / (float)Slices;
				float theta = U * (PI * 2);

				// Calc The Vertex Positions
				float x = cosf(theta) * sinf(phi);
				float y = cosf(phi);
				float z = sinf(theta) * sinf(phi);

				// Push Back Vertex Data
				point.x = x*radius + centerx;
				point.y = y*radius + centery;
				point.z = z*radius + centerz;
				points.push_back(point);
			}
		}

		// Calc The Index Positions
		QVector3D ind;
		for (int i = 0; i < Slices * stacks + Slices; ++i)
		{
			indices.push_back(QVector3D(i, i + Slices + 1, i + Slices));
			indices.push_back(QVector3D(i + Slices + 1, i, i + 1));
		}
	}
};


//本模块最好能承担显示二维图，三维图，线条，面，可重入释放，可以多个界面并存的特点。
//class PtCloud;


class GLVIEW_WIDGETSHARED_EXPORT glview_widget:public QOpenGLWidget, protected QOpenGLFunctions
{

    Q_OBJECT
public:
    explicit glview_widget(QWidget *parent = 0);

    ~glview_widget();
//    Q_OBJECT
//public:
//    explicit es_gl_core(QWidget *parent = 0);
//    ~es_gl_core();

#if 1

    void setPtCloud(PtCloud* data,omesh::Bbox* box);
	void setPtLines(PtCloud* data, omesh::Bbox* box);
	void addLineData(PtCloud* lines);					/// 画线
	void addPtCloud(PtCloud* data, omesh::Bbox* box);
	void clear(void);
    void setLinePts(PtCloud* data);
    void loadPly(QString xyzFileName ="/home/hl/smb_share/3.xyz");
    void loadTextures();
    void setGroupId(int id);
	void setSelectLines(std::vector<int>& select);

	/// 多组姿势
	void setPtCloud(PtCloud* data, omesh::Bbox* box, int pos);
	void addLineData(PtCloud* lines, int pos);
	void setPtCloudVisible(bool visible, int pos);
	void setSelectLines(std::vector<int>& select, int pos);

	void setBox(omesh::Bbox& box);

	void setColor(const QVector3D& MaterialAmbient, const QVector3D& MaterialDiffuse, const QVector3D& MaterialSpecular);

	void addCoordinateData(void);

	void clearPos(int pos = -1);

	void setBalls(vector<Ball> &balls, int  pos);

	void setSelectBall(int index);

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

	void renderBalls(int pos);

    void initShaders();

	/// 封装渲染子场景的函数
	void renderSubSence(void);

private:
    QBasicTimer timer;
    QOpenGLShaderProgram programPtcloud;
    QOpenGLShaderProgram programLine;
    QOpenGLShaderProgram programTexture;
	QOpenGLShaderProgram programAxes;
    data_mem *ptCloudMem;
    data_mem *lineDataMem;
	data_mem *mLineDataMem;			/// 用点的shader画线
    data_mem *textureOpeMem;
	data_mem *mCoordinateMem[3];		/// 坐标轴的点

	vector<data_mem*> m_ballMems;  //球显示
	map<int, vector<data_mem*>> m_ballMap;

	omesh::Bbox mBox;
	QMutex mMutexLines;					/// 线条的数据锁
	std::vector<data_mem*> mLineDatas;

	RenderObject mRenderObject[render_object_num];
    QOpenGLTexture *img_texture;
    //AxiDraw* axiDraw;

	std::vector<int> mSelectLineIndexs;

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


    //支持同时显示点云和线条
    PtCloud points;
    PtCloud line_pts;

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

	int m_selectBall;
#endif
};

#endif // GLVIEW_WIDGET_H
