#include "data_mem.h"

#include <QVector2D>
#include <QVector3D>

//#include "newMesh/include/TriMesh.h"
//#include "newMesh/include/XForm.h"
//#include "newMesh/include/Box.h"
#include "../types.h"

//#include "./Shader/GeometryObject.h"

using namespace std;
//using namespace trimesh;
using namespace omesh;

struct VertexData
{
    QVector3D position;
    QVector3D normal;
};

//! [0]
data_mem::data_mem(eDrawData tp)
    : indexBuf(QOpenGLBuffer::IndexBuffer),drawType(eNullType),pointSize(0)
{
    initializeOpenGLFunctions();

    // Generate 2 VBOs
    arrayBuf.create();
    indexBuf.create();

    // Initializes cube geometry and transfers it to VBOs
   // initCubeGeometry();

    drawType=tp;
    //why this must need??????????????????????
//    if(ePtData==tp){

//        memData dat;
//        ptCloudDateUpdate(&dat);
//    }

//    if(eLineData==tp){

//        memData dat;
//        lineDataUpdate(&dat);
//    }

//    if(eTexture==tp){

//    }

}

data_mem::~data_mem()
{
    arrayBuf.destroy();
    indexBuf.destroy();
}
//! [0]

//extern vector< yxCoda:: > points;
//extern omesh::XForm<float> transMatrix;
//extern omesh::XForm<float> initMatrix;

void data_mem::ptCloudDateUpdate(memData* ptData)
{
    // For cube we would need only 8 vertices but we have to
    // duplicate vertex for each face because texture coordinate
    // is different.

    //drawType=ePtData;
    //isPtShow=1;
#if 0
    yxCoda::GeometryGenerator gg;

    std::vector< yxCoda::Vertex3NormalUV2 > sphere;
    std::vector< yxCoda::Triangle > sphereMesh;
    gg.Sphere(sphere, sphereMesh);

    meshSize = sphereMesh.size()*3;
    pointSize = sphere.size();

    // Transfer vertex data to VBO 0
    arrayBuf.bind();
    arrayBuf.allocate(sphere.data(), sphere.size() * sizeof(yxCoda::Vertex3NormalUV2));

    // Transfer index data to VBO 1
    indexBuf.bind();
    indexBuf.allocate(sphereMesh.data(), sphereMesh.size() * sizeof(yxCoda::Triangle));
#else
   // extern vector< yxCoda::Vertex3Normal > points;
    memData sphere;
    sphere=*ptData;
    pointSize = ptData->size();

    // Transfer vertex data to VBO 0
    bool ret = arrayBuf.bind();
	if (ret == false)	
		printf("gl can not bind %d\n", pointSize);
	
    arrayBuf.allocate(sphere.data(), sphere.size() * sizeof(Vertex3Normal));

    //arrayBuf.write(0,sphere.data(),sphere.size() * sizeof(yxCoda::Vertex3Normal));
    //indexBuf.bind();
    // Transfer index data to VBO 1
    //indexBuf.bind();
    //indexBuf.allocate(sphereMesh.data(), sphereMesh.size() * sizeof(yxCoda::Triangle));

#endif


    #if 1

//    // Transfer vertex data to VBO 0
//    arrayBuf.bind();
//    arrayBuf.allocate(vertices, 24 * sizeof(VertexData));

//    // Transfer index data to VBO 1
//    indexBuf.bind();
//    indexBuf.allocate(cubeMesh.data(), cubeMesh.size() * sizeof(int));

    #else

    int siz=100000;//points.size()/100;
    VertexData verticesBig[siz];
    GLushort indicesBig[siz];

    int pointSiz=points.size();
    for(int ii=0;ii<siz && ii< pointSiz;ii++)
    {
        verticesBig[ii].position.setX(points[ii].x);
        verticesBig[ii].position.setY(points[ii].y);
        verticesBig[ii].position.setZ(points[ii].z);
        verticesBig[ii].texCoord.setX(0.2);
        verticesBig[ii].texCoord.setY(0.2);

        indicesBig[ii]=ii;
    }

    arrayBuf.bind();
    arrayBuf.allocate(verticesBig, siz * sizeof(VertexData));

    // Transfer index data to VBO 1
    indexBuf.bind();
    indexBuf.allocate(indicesBig, siz * sizeof(GLushort));
    #endif
}

void data_mem::lineDataUpdate(memData* lineData){

    //isPtShow=0;
    drawType=eLineData;
    //extern vector< yxCoda::Vertex3Normal > line_pts;
    memData sphere;
    sphere=*lineData/*line_pts*/;
    pointSize = sphere.size();

    // Transfer vertex data to VBO 0
    bool ret = arrayBuf.bind();
	if (ret == false)
		printf("%s can not bind\n", __FUNCTION__);
    arrayBuf.allocate(sphere.data(), sphere.size() * sizeof(Vertex3Normal));
}

void data_mem::imgAxiUpdate()
{
    // For cube we would need only 8 vertices but we have to
    // duplicate vertex for each face because texture coordinate
    // is different.
    VertexData vertices[] = {
        // Vertex data for face 0
//        {QVector3D(-1.0f, -1.0f,  1.0f), QVector2D(0.0f, 0.0f)},  // v0
//        {QVector3D( 1.0f, -1.0f,  1.0f), QVector2D(0.33f, 0.0f)}, // v1
//        {QVector3D(-1.0f,  1.0f,  1.0f), QVector2D(0.0f, 0.5f)},  // v2
//        {QVector3D( 1.0f,  1.0f,  1.0f), QVector2D(0.33f, 0.5f)}, // v3

        {QVector3D(-1.0f, -1.0f,  1.0f), QVector2D(0.0f, 0.0f)},  // v0
        {QVector3D( 1.0f, -1.0f,  1.0f), QVector2D(1.0f, 0.0f)}, // v1
        {QVector3D(-1.0f,  1.0f,  1.0f), QVector2D(0.0f, 1.0f)},  // v2
        {QVector3D( 1.0f,  1.0f,  1.0f), QVector2D(1.0f, 1.0f)}, // v3

        // Vertex data for face 1
        {QVector3D( 1.0f, -1.0f,  1.0f), QVector2D( 0.0f, 0.5f)}, // v4
        {QVector3D( 1.0f, -1.0f, -1.0f), QVector2D(0.33f, 0.5f)}, // v5
        {QVector3D( 1.0f,  1.0f,  1.0f), QVector2D(0.0f, 1.0f)},  // v6
        {QVector3D( 1.0f,  1.0f, -1.0f), QVector2D(0.33f, 1.0f)}, // v7

        // Vertex data for face 2
        {QVector3D( 1.0f, -1.0f, -1.0f), QVector2D(0.66f, 0.5f)}, // v8
        {QVector3D(-1.0f, -1.0f, -1.0f), QVector2D(1.0f, 0.5f)},  // v9
        {QVector3D( 1.0f,  1.0f, -1.0f), QVector2D(0.66f, 1.0f)}, // v10
        {QVector3D(-1.0f,  1.0f, -1.0f), QVector2D(1.0f, 1.0f)},  // v11

        // Vertex data for face 3
        {QVector3D(-1.0f, -1.0f, -1.0f), QVector2D(0.66f, 0.0f)}, // v12
        {QVector3D(-1.0f, -1.0f,  1.0f), QVector2D(1.0f, 0.0f)},  // v13
        {QVector3D(-1.0f,  1.0f, -1.0f), QVector2D(0.66f, 0.5f)}, // v14
        {QVector3D(-1.0f,  1.0f,  1.0f), QVector2D(1.0f, 0.5f)},  // v15

        // Vertex data for face 4
        {QVector3D(-1.0f, -1.0f, -1.0f), QVector2D(0.33f, 0.0f)}, // v16
        {QVector3D( 1.0f, -1.0f, -1.0f), QVector2D(0.66f, 0.0f)}, // v17
        {QVector3D(-1.0f, -1.0f,  1.0f), QVector2D(0.33f, 0.5f)}, // v18
        {QVector3D( 1.0f, -1.0f,  1.0f), QVector2D(0.66f, 0.5f)}, // v19

        // Vertex data for face 5
        {QVector3D(-1.0f,  1.0f,  1.0f), QVector2D(0.33f, 0.5f)}, // v20
        {QVector3D( 1.0f,  1.0f,  1.0f), QVector2D(0.66f, 0.5f)}, // v21
        {QVector3D(-1.0f,  1.0f, -1.0f), QVector2D(0.33f, 1.0f)}, // v22
        {QVector3D( 1.0f,  1.0f, -1.0f), QVector2D(0.66f, 1.0f)}  // v23
    };

    // Indices for drawing cube faces using triangle strips.
    // Triangle strips can be connected by duplicating indices
    // between the strips. If connecting strips have opposite
    // vertex order then last index of the first strip and first
    // index of the second strip needs to be duplicated. If
    // connecting strips have same vertex order then only last
    // index of the first strip needs to be duplicated.
    GLushort indices[] = {
         0,  1,  2,  3,  3,     // Face 0 - triangle strip ( v0,  v1,  v2,  v3)
         4,  4,  5,  6,  7,  7, // Face 1 - triangle strip ( v4,  v5,  v6,  v7)
         8,  8,  9, 10, 11, 11, // Face 2 - triangle strip ( v8,  v9, v10, v11)
        12, 12, 13, 14, 15, 15, // Face 3 - triangle strip (v12, v13, v14, v15)
        16, 16, 17, 18, 19, 19, // Face 4 - triangle strip (v16, v17, v18, v19)
        20, 20, 21, 22, 23      // Face 5 - triangle strip (v20, v21, v22, v23)
    };

//! [1]
    // Transfer vertex data to VBO 0
    arrayBuf.bind();
    arrayBuf.allocate(vertices, 24 * sizeof(VertexData));

    // Transfer index data to VBO 1
    indexBuf.bind();
    indexBuf.allocate(indices, 34 * sizeof(GLushort));

    drawType=eTexture;
//! [1]
}

int data_mem::getPointSize(void) const
{
	return pointSize;
}

void data_mem::setDrawType(eDrawData tp)
{
	drawType = tp;
}

//! [2]
void data_mem::drawWithProgram(QOpenGLShaderProgram *program)
{
    if(drawType==eNullType || pointSize == 0) {

        return;
    }

    //here is also ok
    //imgAxiUpdate();
    if(drawType==ePtData || drawType == ePtLinesData)
    {
        // Tell OpenGL which VBOs to use
        bool ret = arrayBuf.bind();
		if (ret == false)
			printf("%s can not bind\n", __FUNCTION__);
        indexBuf.bind();

        // Offset for position
        int offset = 0;

        // Tell OpenGL programmable pipeline how to locate vertex position data
        //int vertexLocationTs =sizeof(yxCoda::Vertex3Normal);
        int vertexLocation = program->attributeLocation("vaPosition");
        program->enableAttributeArray(vertexLocation);
        program->setAttributeBuffer(vertexLocation, GL_FLOAT, offset, 3, sizeof(Vertex3Normal));

        // Offset for texture coordinate
        offset += sizeof(float)*3;

        // Tell OpenGL programmable pipeline how to locate vertex texture coordinate data
        int normalLocation = program->attributeLocation("vaNormal");
        program->enableAttributeArray(normalLocation);
        program->setAttributeBuffer(normalLocation, GL_FLOAT, offset, 3, sizeof(Vertex3Normal));
		
    //    offset += sizeof(float)*3;

    //    int texcoordLocation = program->attributeLocation("vaTexcoord");
    //    program->enableAttributeArray(texcoordLocation);
    //    program->setAttributeBuffer(texcoordLocation, GL_FLOAT, offset, 2, sizeof(yxCoda::Vertex3NormalUV2));

        // Draw cube geometry using indices from VBO 1
        //glDrawElements(GL_POINTS, 100000, GL_UNSIGNED_SHORT, 0);
         //glDrawElements(GL_TRIANGLES, 24, GL_UNSIGNED_SHORT, 0);
        //glDrawElements(GL_TRIANGLES, meshSize, GL_UNSIGNED_INT, 0);
		if (drawType == ePtData)
		{
			glDrawArrays(GL_TRIANGLES, 0, pointSize);
		}
		else
		{			
			glDrawArrays(GL_LINES, 0, pointSize);
		}
        
    }
    else if(drawType==eLineData)
    {
        // Tell OpenGL which VBOs to use
        arrayBuf.bind();
        indexBuf.bind();

        // Offset for position
        int offset = 0;

        // Tell OpenGL programmable pipeline how to locate vertex position data
        int vertexLocation = program->attributeLocation("vaPosition");
        program->enableAttributeArray(vertexLocation);
        program->setAttributeBuffer(vertexLocation, GL_FLOAT, offset, 3, sizeof(Vertex3Normal));

        // Offset for texture coordinate
        offset += sizeof(float)*3;

        // Tell OpenGL programmable pipeline how to locate vertex texture coordinate data
        int normalLocation = program->attributeLocation("vaNormal");
        program->enableAttributeArray(normalLocation);
        program->setAttributeBuffer(normalLocation, GL_FLOAT, offset, 3, sizeof(Vertex3Normal));

    //    offset += sizeof(float)*3;

    //    int texcoordLocation = program->attributeLocation("vaTexcoord");
    //    program->enableAttributeArray(texcoordLocation);
    //    program->setAttributeBuffer(texcoordLocation, GL_FLOAT, offset, 2, sizeof(yxCoda::Vertex3NormalUV2));

        // Draw cube geometry using indices from VBO 1
        //glDrawElements(GL_POINTS, 100000, GL_UNSIGNED_SHORT, 0);
         //glDrawElements(GL_TRIANGLES, 24, GL_UNSIGNED_SHORT, 0);
        //glDrawElements(GL_TRIANGLES, meshSize, GL_UNSIGNED_INT, 0);

        glDrawArrays(GL_LINES, 0, pointSize);
    }
    else if(drawType==eTexture)
    {
        //glDrawArrays(GL_LINES, 0, pointSize);
        // Tell OpenGL which VBOs to use
        arrayBuf.bind();
        indexBuf.bind();

        // Offset for position
        quintptr offset = 0;

        // Tell OpenGL programmable pipeline how to locate vertex position data
        int vertexLocation = program->attributeLocation("a_position");
        program->enableAttributeArray(vertexLocation);
        program->setAttributeBuffer(vertexLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

        // Offset for texture coordinate
        offset += sizeof(QVector3D);

        // Tell OpenGL programmable pipeline how to locate vertex texture coordinate data
        int texcoordLocation = program->attributeLocation("a_texcoord");
        program->enableAttributeArray(texcoordLocation);
        program->setAttributeBuffer(texcoordLocation, GL_FLOAT, offset, 2, sizeof(VertexData));

        // Draw cube geometry using indices from VBO 1
        //glDrawElements(GL_TRIANGLE_STRIP,/* 34*/4, GL_UNSIGNED_SHORT, 0);
		glDrawArrays(GL_TRIANGLES, 0, pointSize);
		glDrawArrays(GL_POINTS, 0, pointSize);
    }
	else if (drawType == ePoints)
	{
		bool ret = arrayBuf.bind();
		if (ret == false)
			printf("%s can not bind\n", __FUNCTION__);
		indexBuf.bind();

		// Offset for position
		int offset = 0;
		int vertexLocation = program->attributeLocation("vaPosition");
		program->enableAttributeArray(vertexLocation);
		program->setAttributeBuffer(vertexLocation, GL_FLOAT, offset, 3, sizeof(Vertex3Normal));

		// Offset for texture coordinate
		offset += sizeof(float) * 3;

		// Tell OpenGL programmable pipeline how to locate vertex texture coordinate data
		int normalLocation = program->attributeLocation("vaNormal");
		program->enableAttributeArray(normalLocation);
		program->setAttributeBuffer(normalLocation, GL_FLOAT, offset, 3, sizeof(Vertex3Normal));
		
		glDrawArrays(GL_POINTS, 0, pointSize);
	}
}

void data_mem::drawWithProgramPoints(QOpenGLShaderProgram *program)
{
	arrayBuf.bind();
	indexBuf.bind();

	int offset = 0;

	// Tell OpenGL programmable pipeline how to locate vertex position data
	//int vertexLocationTs =sizeof(yxCoda::Vertex3Normal);
	int vertexLocation = program->attributeLocation("vaPosition");
	program->enableAttributeArray(vertexLocation);
	program->setAttributeBuffer(vertexLocation, GL_FLOAT, offset, 3, sizeof(Vertex3Normal));

	// Offset for texture coordinate
	offset += sizeof(float) * 3;

	// Tell OpenGL programmable pipeline how to locate vertex texture coordinate data
	int normalLocation = program->attributeLocation("vaNormal");
	program->enableAttributeArray(normalLocation);
	program->setAttributeBuffer(normalLocation, GL_FLOAT, offset, 3, sizeof(Vertex3Normal));

	glDrawArrays(GL_POINTS, 0, pointSize);
}

//! [2]

