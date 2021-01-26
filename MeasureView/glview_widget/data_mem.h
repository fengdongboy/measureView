#ifndef GEOMETRYENGPT_H
#define GEOMETRYENGPT_H

#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>

#include "../types.h"

typedef PtCloud memData;

enum eDrawData{
    eNullType,ePtData,eLineData,eTexture,ePtLinesData,ePoints
};

class data_mem : protected QOpenGLFunctions
{
public:
    data_mem(eDrawData tp);
    virtual ~data_mem();

     void drawWithProgram(QOpenGLShaderProgram *program);
	 void drawWithProgramPoints(QOpenGLShaderProgram *program);

     void ptCloudDateUpdate(memData* ptData);
     void lineDataUpdate(memData* lineData);
     void imgAxiUpdate();

	 int getPointSize(void) const;

	 void setDrawType(eDrawData tp);
private:

    int meshSize;
    int pointSize;
    QOpenGLBuffer arrayBuf;
    QOpenGLBuffer indexBuf;

    eDrawData drawType;
    //QOpenGLBuffer arrayBuf;
    //QOpenGLBuffer indexBuf;
};

#endif // GEOMETRYENGPT_H
