#include "WalkMesh.hpp"

#include "read_write_chunk.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <glm/gtx/hash.hpp> //allows the use of 'uvec2' as an unordered_map key

#define ERROR_VAL 0.0333f

WalkMesh::WalkMesh(std::vector< glm::vec3 > const& vertices_, std::vector< glm::vec3 > const& normals_, std::vector< glm::uvec3 > const& triangles_)
	: vertices(vertices_), normals(normals_), triangles(triangles_) {

	//construct next_vertex map (maps each edge to the next vertex in the triangle):
	next_vertex.reserve(triangles.size() * 3);
	auto do_next = [this](uint32_t a, uint32_t b, uint32_t c) {
		auto ret = next_vertex.insert(std::make_pair(glm::uvec2(a, b), c));
		assert(ret.second);
	};
	for (auto const& tri : triangles) {
		do_next(tri.x, tri.y, tri.z);
		do_next(tri.y, tri.z, tri.x);
		do_next(tri.z, tri.x, tri.y);
	}

	//DEBUG: are vertex normals consistent with geometric normals?
	for (auto const& tri : triangles) {
		glm::vec3 const& a = vertices[tri.x];
		glm::vec3 const& b = vertices[tri.y];
		glm::vec3 const& c = vertices[tri.z];
		glm::vec3 out = glm::normalize(glm::cross(b - a, c - a));

		float da = glm::dot(out, normals[tri.x]);
		float db = glm::dot(out, normals[tri.y]);
		float dc = glm::dot(out, normals[tri.z]);

		assert(da > 0.1f && db > 0.1f && dc > 0.1f);
	}
}


//project pt to the plane of triangle a,b,c and return the barycentric weights of the projected point:
glm::vec3 barycentric_weights(glm::vec3 const& a, glm::vec3 const& b, glm::vec3 const& c, glm::vec3 const& pt) {
	//Source: My own submission for lesson 8
	//Reference https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates
	auto norm = [](glm::vec3 v) {
		return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	};
	//https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d
	glm::vec3 pointVec = pt - a;
	glm::vec3 normal = -glm::normalize(glm::cross(b - a, c - a));
	float distance = glm::dot(pointVec, normal);
	glm::vec3 projPt = pt - glm::vec3(distance) * normal;

	float ABC = norm(glm::cross(b - a, c - a)) / 2.f;
	float A = norm(glm::cross(b - projPt, c - projPt)) / 2.f;
	float B = norm(glm::cross(projPt - a, c - a)) / 2.f;
	float C = norm(glm::cross(b - a, projPt - a)) / 2.f;
	float u = B / ABC; float v = C / ABC; float w = A / ABC;
	if (abs(-u - v + w - 1.0f) < 0.0001) return glm::vec3(w, -u, -v);
	if (abs(-u + v - w - 1.0f) < 0.0001) return glm::vec3(-w, -u, v);
	if (abs(u - v - w - 1.0f) < 0.0001) return glm::vec3(-w, u, -v);
	if (abs(-u + v + w - 1.0f) < 0.0001) return glm::vec3(w, -u, v);
	if (abs(u - v + w - 1.0f) < 0.0001) return glm::vec3(w, u, -v);
	if (abs(u + v - w - 1.0f) < 0.0001) return glm::vec3(-w, u, v);
	return glm::vec3(w, u, v);
}

WalkPoint WalkMesh::nearest_walk_point(glm::vec3 const& world_point) const {
	assert(!triangles.empty() && "Cannot start on an empty walkmesh");

	WalkPoint closest;
	float closest_dis2 = std::numeric_limits< float >::infinity();

	for (auto const& tri : triangles) {
		//find closest point on triangle:

		glm::vec3 const& a = vertices[tri.x];
		glm::vec3 const& b = vertices[tri.y];
		glm::vec3 const& c = vertices[tri.z];

		//get barycentric coordinates of closest point in the plane of (a,b,c):
		glm::vec3 coords = barycentric_weights(a, b, c, world_point);

		//is that point inside the triangle?
		if (coords.x >= 0.0f && coords.y >= 0.0f && coords.z >= 0.0f) {
			//yes, point is inside triangle.
			float dis2 = glm::length2(world_point - to_world_point(WalkPoint(tri, coords)));
			if (dis2 < closest_dis2) {
				closest_dis2 = dis2;
				closest.indices = tri;
				closest.weights = coords;
			}
		}
		else {
			//check triangle vertices and edges:
			auto check_edge = [&world_point, &closest, &closest_dis2, this](uint32_t ai, uint32_t bi, uint32_t ci) {
				glm::vec3 const& a = vertices[ai];
				glm::vec3 const& b = vertices[bi];

				//find closest point on line segment ab:
				float along = glm::dot(world_point - a, b - a);
				float max = glm::dot(b - a, b - a);
				glm::vec3 pt;
				glm::vec3 coords;
				if (along < 0.0f) {
					pt = a;
					coords = glm::vec3(1.0f, 0.0f, 0.0f);
				}
				else if (along > max) {
					pt = b;
					coords = glm::vec3(0.0f, 1.0f, 0.0f);
				}
				else {
					float amt = along / max;
					pt = glm::mix(a, b, amt);
					coords = glm::vec3(1.0f - amt, amt, 0.0f);
				}

				float dis2 = glm::length2(world_point - pt);
				if (dis2 < closest_dis2) {
					closest_dis2 = dis2;
					closest.indices = glm::uvec3(ai, bi, ci);
					closest.weights = coords;
				}
			};
			check_edge(tri.x, tri.y, tri.z);
			check_edge(tri.y, tri.z, tri.x);
			check_edge(tri.z, tri.x, tri.y);
		}
	}
	assert(closest.indices.x < vertices.size());
	assert(closest.indices.y < vertices.size());
	assert(closest.indices.z < vertices.size());
	return closest;
}


void WalkMesh::walk_in_triangle(WalkPoint const& start, glm::vec3 const& step, WalkPoint* end_, float* time_) const {

	//Source: Code shared for lesson 8 in course discord

	assert(end_);
	auto& end = *end_;
	assert(time_);
	auto& time = *time_;

	assert(vertices.size() > start.indices.x);
	assert(vertices.size() > start.indices.y);
	assert(vertices.size() > start.indices.z);
	glm::vec3 const& a = vertices[start.indices.x];
	glm::vec3 const& b = vertices[start.indices.y];
	glm::vec3 const& c = vertices[start.indices.z];
	assert(step != glm::vec3(0.f));

	//TODO: transform 'step' into a barycentric velocity on (a,b,c)
	glm::vec3 unitNorm = glm::normalize(glm::cross((b - a), (c - a)));
	//project step on triangle plane, assume step and normal at the same side of triangle
	glm::vec3 step_proj = step - glm::vec3(dot(step, unitNorm)) * unitNorm;


	//can use cross product to rotate 90 degree
	assert(a != b && b != c && c != a);
	assert(step_proj != glm::vec3(0.0f));
	glm::vec3 scaleA = (a - b) - glm::dot(a - b, glm::normalize(c - b)) * glm::normalize(c - b);
	glm::vec3 scaleB = (b - a) - glm::dot(b - a, glm::normalize(c - a)) * glm::normalize(c - a);
	glm::vec3 scaleC = (c - a) - glm::dot(c - a, glm::normalize(b - a)) * glm::normalize(b - a);

	float stepA = glm::dot(step_proj, glm::normalize(scaleA)) / glm::length(scaleA);
	float stepB = glm::dot(step_proj, glm::normalize(scaleB)) / glm::length(scaleB);
	float stepC = glm::dot(step_proj, glm::normalize(scaleC)) / glm::length(scaleC);

	glm::vec3 step_bary = glm::vec3(stepA, stepB, stepC);

	//TODO: check when/if this velocity pushes start.weights into an edge
	//weight already in barycentric coordinate
	glm::vec3 move = start.weights + step_bary * 1.0f; //Calculate max timestep as 1.0f;
	assert(step_bary != glm::vec3(0.f) && move != start.weights);

	float minTime = 1.0f; glm::uvec2 edge = glm::uvec2(start.indices.y, start.indices.x);
	if (move.x < 0.0f) {
		auto next_edge = next_vertex.find(edge);
		assert(step_bary.x != 0);
		minTime = -start.weights.x / step_bary.x;
		move = start.weights + step_bary * -start.weights.x / step_bary.x;
		if (minTime != 1.0f)
			assert(!(move.x > ERROR_VAL && move.y > ERROR_VAL && move.z > ERROR_VAL)); //here
		assert(move != glm::vec3(0.0f));
	}
	if (move.y < 0.0f) {
		auto next_edge = next_vertex.find(edge);
		assert(step_bary.y != 0);
		minTime = -start.weights.y / step_bary.y;
		move = start.weights + step_bary * -start.weights.y / step_bary.y;
		if (minTime != 1.0f)
			assert(!(move.x > ERROR_VAL && move.y > ERROR_VAL && move.z > ERROR_VAL));
		assert(move != glm::vec3(0.0f));
	}
	if (move.z < 0.0f) {
		auto next_edge = next_vertex.find(edge);
		assert(step_bary.z != 0);
		minTime = -start.weights.z / step_bary.z;
		move = start.weights + step_bary * -start.weights.z / step_bary.z;
		if (minTime != 1.0f)
			assert(!(move.x > ERROR_VAL && move.y > ERROR_VAL && move.z > ERROR_VAL));
		assert(move != glm::vec3(0.0f));
	}


	auto rotateEnd = [this](WalkPoint* end_, float time) {
		if (time == 1.0f) {
			return;
		}
		while (end_->weights.z > ERROR_VAL) {
			std::swap(end_->weights.x, end_->weights.y);
			std::swap(end_->weights.y, end_->weights.z);
			std::swap(end_->indices.x, end_->indices.y);
			std::swap(end_->indices.y, end_->indices.z);
		}
		assert(end_->weights.z <= ERROR_VAL);
	};

	*time_ = minTime;
	end_->weights = move;
	end_->indices = start.indices;
	rotateEnd(end_, *time_);
}




std::vector<uint32_t> WalkMesh::getRegion(glm::uvec2 startEdge, std::vector<bool>* claimedVerts) const{
	std::vector<uint32_t> retVerts; //Newly found vertices
	retVerts.reserve(11);

	//Then add to a counter once function is done, and update colors, back in PlayMode

	auto outerLoop = [this](glm::uvec2 nextEdge, glm::uvec2 startEdge, std::vector<uint32_t> retVerts, glm::uvec3 triangle,bool* done, std::vector<bool> *claimedVerts) {

		if ((*claimedVerts)[nextEdge.x] == false) {
			(*claimedVerts)[nextEdge.x] = true;
			retVerts.emplace_back(nextEdge.x);
		}
		if ((*claimedVerts)[nextEdge.y] == false) {
			(*claimedVerts)[nextEdge.y] = true;
			retVerts.emplace_back(nextEdge.y);
		}
		auto nextEdgeIter = next_vertex.find(nextEdge);
		if (nextEdgeIter == next_vertex.end()) { //Not a triangle
			*done = true;
			return std::make_pair( retVerts, nextEdge);
		}
		nextEdge = glm::uvec2(nextEdge.y, next_vertex.at(nextEdge)); //Next edge in triangle
		if (((*claimedVerts))[nextEdge.y] == false) {
			(*claimedVerts)[nextEdge.y] = true;
			retVerts.emplace_back(nextEdge.y);
		}
		auto inTriangle = [this](glm::uvec2 edge, glm::uvec3 triangle) { //Checks if flip edge is the triangle, if so, need the next edge after current, 
			if (edge.x == triangle.y && edge.y == triangle.x) return true; //Then can get the flip edge
			if (edge.x == triangle.z && edge.y == triangle.y) return true;
			if (edge.x == triangle.x && edge.y == triangle.z) return true;
			return false;
		};
		if(inTriangle(nextEdge,triangle))
			nextEdge = glm::uvec2(nextEdge.y, next_vertex.at(nextEdge));

		glm::uvec2 edge = glm::uvec2(startEdge.y, startEdge.x); //Gets flip edge
		auto next_edgeIter = next_vertex.find(edge);
		if (next_edgeIter == next_vertex.end()) { //If not a triangle, just return what has already been found
			*done = true;
			return std::make_pair(retVerts, edge);
		}
		nextEdge = glm::uvec2(edge.y, next_vertex.at(edge)); //What is next edge here for?
		if (nextEdge == startEdge) {
			*done = true;
			return std::make_pair(retVerts, edge);
		}
		*done = false;
		return std::make_pair(retVerts, edge); 
	};

	if ((*claimedVerts)[startEdge.x] == false) {
		(*claimedVerts)[startEdge.x] = true;
		retVerts.emplace_back(startEdge.x);
	}
	if ((*claimedVerts)[startEdge.y] == false) {
		(*claimedVerts)[startEdge.y] = true;
		retVerts.emplace_back(startEdge.y);
	}

	glm::uvec2 edge = glm::uvec2(startEdge.y, startEdge.x); //Get flipped edge
	auto next_edgeIter = next_vertex.find(edge);
	if (next_edgeIter == next_vertex.end()) //If no triangle, return
		return retVerts;


	auto triZ = next_vertex.at(startEdge);
	if ((*claimedVerts)[triZ] == false) {
		(*claimedVerts)[triZ] = true;
		retVerts.emplace_back(triZ);
	}

	glm::uvec3 triangle;
	triangle.x = startEdge.x; triangle.y = startEdge.y; triangle.z = triZ;

	glm::uvec2 next_edge = glm::uvec2(edge.y, next_vertex.at(edge)); //Get edge that will be first found by outerloop
	next_edge = glm::uvec2(next_edge.y,next_vertex.at(next_edge)); //of this triangle to signify when it should stop
	bool done = false;
	while (!done) {
		std::pair<std::vector<uint32_t>,glm::uvec2> retPair = outerLoop(edge, next_edge, retVerts, triangle, &done, claimedVerts);
		edge = retPair.second;
		retVerts = retPair.first;
	}
	return retVerts;
}

bool WalkMesh::cross_edge(WalkPoint const& start, WalkPoint* end_, glm::quat* rotation_) const {


	assert(end_);

	assert(rotation_);

	assert(start.weights.z <= ERROR_VAL); //*must* be on an edge.

	// check if edge (start.indices.x, start.indices.y) has a triangle on the other side:
	glm::uvec2 edge = glm::uvec2(start.indices.y, start.indices.x);
	auto next_edge = next_vertex.find(edge);
	if (next_edge == next_vertex.end())
		return false;

	// if there is another triangle:
	// set end's weights and indicies on that triangle:
	end_->indices = glm::vec3(start.indices.y, start.indices.x, next_edge->second);
	end_->weights = glm::vec3(start.weights.y, start.weights.x, 0);

	// compute rotation that takes starting triangle's normal to ending triangle's normal:
	assert(vertices.size() > start.indices.x);
	assert(vertices.size() > start.indices.y);
	assert(vertices.size() > start.indices.z);
	assert(vertices.size() > next_edge->second);
	glm::vec3 const& b = vertices[start.indices.x];
	glm::vec3 const& c = vertices[start.indices.y];
	glm::vec3 const& a = vertices[start.indices.z];
	glm::vec3 const& d = vertices[next_edge->second];

	glm::vec3 start_norm = normalize(cross(c - b, a - b));
	glm::vec3 end_norm = normalize(cross(d - b, c - b));

	*rotation_ = glm::rotation(start_norm, end_norm);

	//return 'true' if there was another triangle, 'false' otherwise:
	return true;
}


WalkMeshes::WalkMeshes(std::string const& filename) {
	std::ifstream file(filename, std::ios::binary);

	std::vector< glm::vec3 > vertices;
	read_chunk(file, "p...", &vertices);

	std::vector< glm::vec3 > normals;
	read_chunk(file, "n...", &normals);

	std::vector< glm::uvec3 > triangles;
	read_chunk(file, "tri0", &triangles);

	std::vector< char > names;
	read_chunk(file, "str0", &names);

	struct IndexEntry {
		uint32_t name_begin, name_end;
		uint32_t vertex_begin, vertex_end;
		uint32_t triangle_begin, triangle_end;
	};

	std::vector< IndexEntry > index;
	read_chunk(file, "idxA", &index);

	if (file.peek() != EOF) {
		std::cerr << "WARNING: trailing data in walkmesh file '" << filename << "'" << std::endl;
	}

	//-----------------

	if (vertices.size() != normals.size()) {
		throw std::runtime_error("Mis-matched position and normal sizes in '" + filename + "'");
	}

	for (auto const& e : index) {
		if (!(e.name_begin <= e.name_end && e.name_end <= names.size())) {
			throw std::runtime_error("Invalid name indices in index of '" + filename + "'");
		}
		if (!(e.vertex_begin <= e.vertex_end && e.vertex_end <= vertices.size())) {
			throw std::runtime_error("Invalid vertex indices in index of '" + filename + "'");
		}
		if (!(e.triangle_begin <= e.triangle_end && e.triangle_end <= triangles.size())) {
			throw std::runtime_error("Invalid triangle indices in index of '" + filename + "'");
		}

		//copy vertices/normals:
		std::vector< glm::vec3 > wm_vertices(vertices.begin() + e.vertex_begin, vertices.begin() + e.vertex_end);
		std::vector< glm::vec3 > wm_normals(normals.begin() + e.vertex_begin, normals.begin() + e.vertex_end);

		//remap triangles:
		std::vector< glm::uvec3 > wm_triangles; wm_triangles.reserve(e.triangle_end - e.triangle_begin);
		for (uint32_t ti = e.triangle_begin; ti != e.triangle_end; ++ti) {
			if (!((e.vertex_begin <= triangles[ti].x && triangles[ti].x < e.vertex_end)
				&& (e.vertex_begin <= triangles[ti].y && triangles[ti].y < e.vertex_end)
				&& (e.vertex_begin <= triangles[ti].z && triangles[ti].z < e.vertex_end))) {
				throw std::runtime_error("Invalid triangle in '" + filename + "'");
			}
			wm_triangles.emplace_back(
				triangles[ti].x - e.vertex_begin,
				triangles[ti].y - e.vertex_begin,
				triangles[ti].z - e.vertex_begin
			);
		}

		std::string name(names.begin() + e.name_begin, names.begin() + e.name_end);

		auto ret = meshes.emplace(name, WalkMesh(wm_vertices, wm_normals, wm_triangles));
		if (!ret.second) {
			throw std::runtime_error("WalkMesh with duplicated name '" + name + "' in '" + filename + "'");
		}

	}
}

WalkMesh const& WalkMeshes::lookup(std::string const& name) const {
	auto f = meshes.find(name);
	if (f == meshes.end()) {
		throw std::runtime_error("WalkMesh with name '" + name + "' not found.");
	}
	return f->second;
}