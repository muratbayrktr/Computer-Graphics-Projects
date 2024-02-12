#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	// Get the required matrices for each camera
	Matrix4 finalTransformationMatrix = getIdentityMatrix();
	Matrix4 modelingTransformationMatrix = getIdentityMatrix();
	Matrix4 cameraTransformationMatrix = getCameraTransformation(camera);
	Matrix4 projectionTransformationMatrix = getProjectionTransformation(camera);
	Matrix4 viewportTransformationMatrix = getViewportTransformation(camera);
	
	// For each mesh consist of triangles do the forward rendering pipeline
	for (int i=0; i<(this->meshes).size(); i++)
	{
		// Get the model transformation matrix for each mesh -> 
		modelingTransformationMatrix = getModelingTransformation(this->meshes[i], this->scalings, this->rotations, this->translations);

		// Final transformation matrix should be -> M_vp * M_per * M_proj * M_cam * M_model
		finalTransformationMatrix = multiplyMatrixWithMatrix(projectionTransformationMatrix, multiplyMatrixWithMatrix(cameraTransformationMatrix, modelingTransformationMatrix));

		// For each triangle in the mesh do the pipeline:
		// 1. Project the vertices of the triangle
		// 2. Do the perspective divide 
		// 3. Do the culling
		// 4. Do the clipping if wireFrame mode is enabled
			// 4.1.1 If we do clipping, we should construct lines, do clipping with lines
			// 4.1.2 Do the viewport transformation
			// 4.1.3 Do the rasterization with lines
		// 4.2.1 If we do not do clipping, do viewport transformation directly
		// 4.2.2 Do the rasterization with triangles
		for(int j=0; j<this->meshes[i]->triangles.size(); j++)
		{
			// Project the vertices of the triangle
			std::vector<Vec4> projectedVertices = getProjectedVertices((this->meshes[i]->triangles[j]), this->vertices, finalTransformationMatrix);
			
			// Do the perspective divide
			if(camera->projectionType == PERSPECTIVE_PROJECTION)
			{
				projectedVertices = perspectiveDivision(projectedVertices);
			}
			
			// Do the culling
			if(Culling(this->cullingEnabled, projectedVertices, camera->position)) // ?? Check this!
			{
				//std::cout << "Culling is done!" << std::endl;
				continue;
			}

			// Do not clip if the mesh is solid, do viewport transformation and rasterization directly
			if(this->meshes[i]->type == SOLID_MESH)
			{
				// Do the viewport transformation
				for(int k=0; k<projectedVertices.size(); k++)
				{
					projectedVertices[k] = multiplyMatrixWithVec4(viewportTransformationMatrix, projectedVertices[k]);
				}
				
				//std::cout << "Viewport transformation is done for SOLIDMESH!" << std::endl;
				
				triangleRasterization(this->image, this->depth, this->colorsOfVertices, projectedVertices, camera->horRes, camera->verRes);
			}
			else if(this->meshes[i]->type == WIREFRAME_MESH)
			{
				// When clipping, we should also interpolate the colors of the vertices
				Vec4 v0 = projectedVertices[0];
				Vec4 v1 = projectedVertices[1];
				Color c0 = *(this->colorsOfVertices[v0.colorId-1]);
				Color c1 = *(this->colorsOfVertices[v1.colorId-1]);

				// Do the clipping
				bool v0v1Line = clippingLiangBarsky(v0, v1, c0, c1);
				//std::cout << v0v1Line << std::endl;
				
				// Do the viewport transformation
				v0 = multiplyMatrixWithVec4(viewportTransformationMatrix, v0);
				v1 = multiplyMatrixWithVec4(viewportTransformationMatrix, v1);
				
				// Do the rasterization with lines, including depth buffer algorithm
				if(v0v1Line)
				{
					//std::cout << "v0v1Line is true!" << std::endl;
					lineRasterization(this->image, this->depth, v0, v1, c0, c1);
					
				}
				v0 = projectedVertices[1];
				v1 = projectedVertices[2];
				c0 = *(this->colorsOfVertices[v0.colorId-1]);
				c1 = *(this->colorsOfVertices[v1.colorId-1]);
				bool v1v2Line = clippingLiangBarsky(v0, v1, c0, c1);
				//std::cout << v1v2Line << std::endl;
				v0 = multiplyMatrixWithVec4(viewportTransformationMatrix, v0);
				v1 = multiplyMatrixWithVec4(viewportTransformationMatrix, v1);
				if(v1v2Line)
				{
					//std::cout << "v1v2Line is true!" << std::endl;
					lineRasterization(this->image, this->depth, v0, v1, c0, c1);
					
				}
				v0 = projectedVertices[2];
				v1 = projectedVertices[0];
				c0 = *(this->colorsOfVertices[v0.colorId-1]);
				c1 = *(this->colorsOfVertices[v1.colorId-1]);
				bool v2v0Line = clippingLiangBarsky(v0, v1, c0, c1);
				//std::cout << v2v0Line << std::endl;
				v0 = multiplyMatrixWithVec4(viewportTransformationMatrix, v0);
				v1 = multiplyMatrixWithVec4(viewportTransformationMatrix, v1);
				if(v2v0Line)
				{
					//std::cout << "v2v0Line is true!" << std::endl;
					lineRasterization(this->image, this->depth, v0, v1, c0, c1);
					
				}
			}
		}
	}
}	

/*
	Forward Rendering Pipeline step by step (from lectures and slides):

	1. Model Transformation
	2. Camera Transformation
	3. Projection Transformation (Orthographic or Perspective)
	4. Perspective Divide
	5. Culling
	6. Clipping if wireFrame mode is enabled
	7. Viewport Transformation
	8. Rasterization

*/

