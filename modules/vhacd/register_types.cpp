/*************************************************************************/
/*  register_types.cpp                                                   */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2021 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2021 Godot Engine contributors (cf. AUTHORS.md).   */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#include "register_types.h"
#include "core/math/geometry.h"
#include "scene/resources/mesh.h"
#include "thirdparty/vhacd/inc/btConvexHullComputer.h"
#include "thirdparty/vhacd/public/VHACD.h"

static Vector<Vector<Face3> > convex_decompose(const Vector<Face3> &p_faces) {

	Vector<float> vertices;
	vertices.resize(p_faces.size() * 9);
	Vector<uint32_t> indices;
	indices.resize(p_faces.size() * 3);

	for (int i = 0; i < p_faces.size(); i++) {
		for (int j = 0; j < 3; j++) {
			vertices.write[i * 9 + j * 3 + 0] = p_faces[i].vertex[j].x;
			vertices.write[i * 9 + j * 3 + 1] = p_faces[i].vertex[j].y;
			vertices.write[i * 9 + j * 3 + 2] = p_faces[i].vertex[j].z;
			indices.write[i * 3 + j] = i * 3 + j;
		}
	}

	VHACD::IVHACD *decomposer = VHACD::CreateVHACD();
	VHACD::IVHACD::Parameters params;
	decomposer->Compute(vertices.ptr(), vertices.size() / 3, indices.ptr(), indices.size() / 3, params);

	int hull_count = decomposer->GetNConvexHulls();

	Vector<Vector<Face3> > ret;

	for (int i = 0; i < hull_count; i++) {
		Vector<Face3> triangles;
		VHACD::IVHACD::ConvexHull hull;
		decomposer->GetConvexHull(i, hull);
		triangles.resize(hull.m_nTriangles);
		for (uint32_t j = 0; j < hull.m_nTriangles; j++) {
			Face3 f;
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					f.vertex[k][l] = hull.m_points[hull.m_triangles[j * 3 + k] * 3 + l];
				}
			}
			triangles.write[j] = f;
		}
		ret.push_back(triangles);
	}

	decomposer->Clean();
	decomposer->Release();

	return ret;
}

static Geometry::MeshData convex_hull(const Vector<Vector3> &p_points) {
	// build the convex hull using the convex hull computer from bullet.
	// simply copies the data over to Godot's types
	if (p_points.size() == 0)
		return Geometry::MeshData();

	VHACD::btConvexHullComputer ch;
	ch.compute(&p_points.ptr()[0][0], sizeof(p_points.ptr()[0]), p_points.size(), -1.0, -1.0);

	Geometry::MeshData ret;
	ret.vertices.resize(ch.vertices.size());
	for (int i = 0; i < ch.vertices.size(); i++) {
		ret.vertices.write[i].x = ch.vertices[i].getX();
		ret.vertices.write[i].y = ch.vertices[i].getY();
		ret.vertices.write[i].z = ch.vertices[i].getZ();
	}

	ret.edges.resize(ch.edges.size());
	for (int i = 0; i < ch.edges.size(); i++) {
		ret.edges.write[i].a = (&ch.edges[i])->getSourceVertex();
		ret.edges.write[i].b = (&ch.edges[i])->getTargetVertex();
	}

	ret.faces.resize(ch.faces.size());
	for (int i = 0; i < ch.faces.size(); i++) {
		const VHACD::btConvexHullComputer::Edge *e_start = &ch.edges[ch.faces[i]];
		const VHACD::btConvexHullComputer::Edge *e = e_start;
		Geometry::MeshData::Face &face = ret.faces.write[i];

		do {
			face.indices.push_back(e->getTargetVertex());

			e = e->getNextEdgeOfFace();
		} while (e != e_start);

		// compute normal
		if (face.indices.size() >= 3) {
			face.plane = Plane(ret.vertices[face.indices[0]], ret.vertices[face.indices[2]], ret.vertices[face.indices[1]]);
		} else {
			WARN_PRINT("Too few vertices per face.");
		}
	}

	return ret;
}

void register_vhacd_types() {
	Mesh::convex_decomposition_function = convex_decompose;
	Geometry::convex_hull_function = convex_hull;
}

void unregister_vhacd_types() {
	Mesh::convex_decomposition_function = NULL;
	Geometry::convex_hull_function = NULL;
}
