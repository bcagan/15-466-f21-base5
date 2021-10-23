#include "PlayMode.hpp"

#include "LitColorTextureProgram.hpp"

#include "DrawLines.hpp"
#include "Mesh.hpp"
#include "Load.hpp"
#include "gl_errors.hpp"
#include "data_path.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>

#include <random>
#include <assert.h>
#include <stdio.h>
#include <glm/gtx/hash.hpp>
#define ERROR_VAL 0.00003f
#define ORANGE_COLOR glm::u8vec4(255, 132, 0,255)
#define BLUE_COLOR glm::u8vec4(0, 50, 255,255)

//To do: 
/*
Move VBO for scene into playmoed and make it only for visualwalkmesh
Create new vao for visualwalkmesh. Handle all this in playmode (ie, don't rely on mesh.cpp)
Make it so tjhere is only one copy of vertices stored in playmode, the one for visual walk mesh
Add code to copy walkmesh into a new drawable manually instead of making second mesh*/


GLuint scene1_meshes_for_lit_color_texture_program = 0;
Load< MeshBuffer > scene1_meshes(LoadTagDefault, []() -> MeshBuffer const * {
	MeshBuffer const *ret = new MeshBuffer(data_path("scene1.pnct"));
	scene1_meshes_for_lit_color_texture_program = ret->make_vao_for_program(lit_color_texture_program->program);
	return ret;
});

Load< Scene > scene1_scene(LoadTagDefault, []() -> Scene const * {
	return new Scene(data_path("scene1.scene"), [&](Scene &scene, Scene::Transform *transform, std::string const &mesh_name){
		Mesh const &mesh = scene1_meshes->lookup(mesh_name);

		scene.drawables.emplace_back(transform);
		Scene::Drawable &drawable = scene.drawables.back();
		
		drawable.pipeline = lit_color_texture_program_pipeline;

		drawable.pipeline.vao = scene1_meshes_for_lit_color_texture_program;
		drawable.pipeline.type = mesh.type;
		drawable.pipeline.start = mesh.start;
		drawable.pipeline.count = mesh.count;
	});
});

WalkMesh const *walkmesh = nullptr;
Load< WalkMeshes > scene1_walkmeshes(LoadTagDefault, []() -> WalkMeshes const * {
	WalkMeshes *ret = new WalkMeshes(data_path("scene1.w"));
	walkmesh = &ret->lookup("WalkMesh");
	return ret;
});

PlayMode::PlayMode() : scene(*scene1_scene) {
	//create a player transform:
	
	auto createVAO = [this](GLuint program) {

		glGenBuffers(1, &walkmeshBuffer);
		GL_ERRORS();

		//Create VBO in the same way as mesh
		glBindBuffer(GL_ARRAY_BUFFER, walkmeshBuffer);
		GL_ERRORS();
		glBufferData(GL_ARRAY_BUFFER, walkmeshData.size() * sizeof(Mesh::Vertex), walkmeshData.data(), GL_DYNAMIC_DRAW);
		GL_ERRORS();
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		GL_ERRORS();
		MeshBuffer::Attrib Position = MeshBuffer::Attrib(3, GL_FLOAT, GL_FALSE, sizeof(Mesh::Vertex), offsetof(Mesh::Vertex, Position));
		MeshBuffer::Attrib Normal = MeshBuffer::Attrib(3, GL_FLOAT, GL_FALSE, sizeof(Mesh::Vertex), offsetof(Mesh::Vertex, Normal));
		MeshBuffer::Attrib Color = MeshBuffer::Attrib(4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Mesh::Vertex), offsetof(Mesh::Vertex, Color));
		MeshBuffer::Attrib TexCoord = MeshBuffer::Attrib(2, GL_FLOAT, GL_FALSE, sizeof(Mesh::Vertex), offsetof(Mesh::Vertex, TexCoord));

		//Create  VAO in the same way as mesh and return
		GLuint vao = 0;
		glGenVertexArrays(1, &vao);
		GL_ERRORS();
		glBindVertexArray(vao);
		GL_ERRORS();

		//Try to bind all attributes in this buffer:
		glBindBuffer(GL_ARRAY_BUFFER, walkmeshBuffer);
		GL_ERRORS();
		auto bind_attribute = [&](char const* name, MeshBuffer::Attrib const& attrib, GLint program) {
			if (attrib.size == 0) return; //don't bind empty attribs
			GLint location = glGetAttribLocation(program, name);
			GL_ERRORS();
			if (location == -1) return; //can't bind missing attribs
			glVertexAttribPointer(location, attrib.size, attrib.type, attrib.normalized, attrib.stride, (GLbyte*)0 + attrib.offset);
			GL_ERRORS();
			glEnableVertexAttribArray(location);
			GL_ERRORS();
		};
		bind_attribute("Position", Position, program);
		bind_attribute("Normal", Normal, program);
		bind_attribute("Color", Color, program);
		bind_attribute("TexCoord", TexCoord, program);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		GL_ERRORS();
		glBindVertexArray(0);
		GL_ERRORS();

		return vao;

	};
	//Creating transform causes issues
	scene.transforms.emplace_back();
	Scene::Transform* transformPoint = &scene.transforms.back();
	Scene::Drawable newDrawable = Scene::Drawable(transformPoint);
	scene.drawables.push_back(newDrawable);
	walkmeshDrawable= &scene.drawables.back();

	scene.transforms.emplace_back();
	player.transform = &scene.transforms.back();

	//Create vertices
	walkmeshDataVerts.reserve(walkmesh->vertices.size());
	for (size_t ind = 0; ind < walkmesh->vertices.size(); ind++) {
		Mesh::Vertex newVert;
		newVert.Position = walkmesh->vertices[ind];
		newVert.Normal = walkmesh->normals[ind];
		newVert.Color = ORANGE_COLOR;
		newVert.TexCoord = glm::uvec2(1.0f, 1.0f); //TexCord isn't used
		walkmeshDataVerts.push_back( newVert);
	}
	walkmeshData.reserve(walkmesh->triangles.size() * 3);
	for (size_t ind = 0; ind < walkmesh->triangles.size(); ind++) {
		Mesh::Vertex vert0 = walkmeshDataVerts[walkmesh->triangles[ind].x];
		Mesh::Vertex vert1 = walkmeshDataVerts[walkmesh->triangles[ind].y];
		Mesh::Vertex vert2 = walkmeshDataVerts[walkmesh->triangles[ind].z];
		walkmeshData.push_back(vert0);
		walkmeshData.push_back(vert1);
		walkmeshData.push_back(vert2);
	}
	
	walkmeshDrawable->pipeline = lit_color_texture_program_pipeline;
	walkmeshDrawable->pipeline.program = lit_color_texture_program->program;
	walkmeshDrawable->pipeline.vao = createVAO(lit_color_texture_program->program);
	walkmeshDrawable->pipeline.type = scene.drawables.front().pipeline.type;
	walkmeshDrawable->pipeline.count = (GLuint)(walkmesh->triangles.size() * 3);
	walkmeshDrawable->pipeline.start = 0;
	walkmeshDrawable->transform->name = std::string("walkmeshDrawable");
	totalVertSize = walkmesh->vertices.size();

	//create a player camera attached to a child of the player transform:
	scene.transforms.emplace_back();
	scene.cameras.emplace_back(&scene.transforms.back());
	player.camera = &scene.cameras.back();
	player.camera->fovy = glm::radians(60.0f);
	player.camera->near = 0.01f;
	player.camera->transform->parent = player.transform;

	//player's eyes are 1.8 units above the ground:
	player.camera->transform->position = glm::vec3(0.0f, 0.0f, 1.8f);

	//rotate camera facing direction (-z) to player facing direction (+y):
	player.camera->transform->rotation = glm::angleAxis(glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));

	//start player walking at nearest walk point:
	player.at = walkmesh->nearest_walk_point(player.transform->position);

	claimedVerts.reserve(walkmesh->vertices.size());
	for (size_t c = 0; c < walkmesh->vertices.size(); c++)
		claimedVerts[c] = false;

	stageWinRatio = 0.8f;

}

PlayMode::~PlayMode() {
}

bool PlayMode::handle_event(SDL_Event const &evt, glm::uvec2 const &window_size) {

	if (evt.type == SDL_KEYDOWN) {
		if (evt.key.keysym.sym == SDLK_ESCAPE) {
			SDL_SetRelativeMouseMode(SDL_FALSE);
			return true;
		} else if (evt.key.keysym.sym == SDLK_a) {
			left.downs += 1;
			left.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.downs += 1;
			right.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			up.downs += 1;
			up.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			down.downs += 1;
			down.pressed = true;
			return true;
		}
	} else if (evt.type == SDL_KEYUP) {
		if (evt.key.keysym.sym == SDLK_a) {
			left.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			up.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			down.pressed = false;
			return true;
		}
	} else if (evt.type == SDL_MOUSEBUTTONDOWN) {
		if (SDL_GetRelativeMouseMode() == SDL_FALSE) {
			SDL_SetRelativeMouseMode(SDL_TRUE);
			return true;
		}
	} else if (evt.type == SDL_MOUSEMOTION) {
		if (SDL_GetRelativeMouseMode() == SDL_TRUE) {
			glm::vec2 motion = glm::vec2(
				evt.motion.xrel / float(window_size.y),
				-evt.motion.yrel / float(window_size.y)
			);
			glm::vec3 up = walkmesh->to_world_smooth_normal(player.at);
			player.transform->rotation = glm::angleAxis(-motion.x * player.camera->fovy, up) * player.transform->rotation;

			float pitch = glm::pitch(player.camera->transform->rotation);
			pitch += motion.y * player.camera->fovy;
			//camera looks down -z (basically at the player's feet) when pitch is at zero.
			pitch = std::min(pitch, 0.95f * 3.1415926f);
			pitch = std::max(pitch, 0.05f * 3.1415926f);
			player.camera->transform->rotation = glm::angleAxis(pitch, glm::vec3(1.0f, 0.0f, 0.0f));

			return true;
		}
	}

	return false;
}

void PlayMode::update(float elapsed) {
	if(!ifWon) timer += elapsed;
	if (timer >= timeLimit) ifWon = true;
	//player walking:
	{
		//combine inputs into a move:
		constexpr float PlayerSpeed = 3.0f;
		glm::vec2 move = glm::vec2(0.0f);
		if (left.pressed && !right.pressed) move.x =-1.0f;
		if (!left.pressed && right.pressed) move.x = 1.0f;
		if (down.pressed && !up.pressed) move.y =-1.0f;
		if (!down.pressed && up.pressed) move.y = 1.0f;

		glm::vec3 startPos = player.transform->position;
		//make it so that moving diagonally doesn't go faster:
		if (move != glm::vec2(0.0f)) move = glm::normalize(move) * PlayerSpeed * elapsed;

		//get move in world coordinate system:
		glm::vec3 remain = player.transform->make_local_to_world() * glm::vec4(move.x, move.y, 0.0f, 0.0f);

		//using a for() instead of a while() here so that if walkpoint gets stuck in
		// some awkward case, code will not infinite loop:
		for (uint32_t iter = 0; iter < 100; ++iter) {
			if (remain == glm::vec3(0.0f)) break;
			WalkPoint end;
			float time = 0.f;
			walkmesh->walk_in_triangle(player.at, remain, &end, &time);
			player.at = end; //crash is here
			if (time == 1.0f) {
				//finished within triangle:
				remain = glm::vec3(0.0f);
				break;
			}
			assert(player.at.weights.z <= ERROR_VAL);
			//some step remains:
			remain *= glm::vec3(1.0f - time);
			//try to step over edge:
			glm::quat rotation;
			if (walkmesh->cross_edge(player.at, &end, &rotation)) {
				//stepped to a new triangle:
				player.at = end;
				//rotate step to follow surface:
				remain = rotation * remain;
			} else {
				//ran into a wall, bounce / slide along it:
				glm::vec3 const &a = walkmesh->vertices[player.at.indices.x];
				glm::vec3 const &b = walkmesh->vertices[player.at.indices.y];
				glm::vec3 const &c = walkmesh->vertices[player.at.indices.z];
				glm::vec3 along = glm::normalize(b-a);
				glm::vec3 normal = glm::normalize(glm::cross(b-a, c-a));
				glm::vec3 in = glm::cross(normal, along);

				//check how much 'remain' is pointing out of the triangle:
				float d = glm::dot(remain, in);
				if (d < 0.0f) {
					//bounce off of the wall:
					remain += (-1.25f * d) * in;
				} else {
					//if it's just pointing along the edge, bend slightly away from wall:
					remain += 0.01f * d * in;
				}
			}
		}


		if (remain != glm::vec3(0.0f)) {
			std::cout << "NOTE: code used full iteration budget for walking." << std::endl;
		}

		//update player's position to respect walking:
		player.transform->position = walkmesh->to_world_point(player.at);
		//pos not being updates
		//assert(move == glm::vec2(0.f) ||  player.transform->position != startPos);

		{ //update player's rotation to respect local (smooth) up-vector:
			
			glm::quat adjust = glm::rotation(
				player.transform->rotation * glm::vec3(0.0f, 0.0f, 1.0f), //current up vector
				walkmesh->to_world_smooth_normal(player.at) //smoothed up vector at walk location
			);
			player.transform->rotation = glm::normalize(adjust * player.transform->rotation);
		}

		//update vertex colors/score
		std::vector<uint32_t> claimed = walkmesh->getRegion(glm::uvec2(player.at.indices.x, player.at.indices.y), &claimedVerts);
		numClaimed += (claimed.size());
		for (size_t c = 0; c < claimed.size(); c++) {
			walkmeshDataVerts[claimed[c]].Color = BLUE_COLOR;
		}
	}

	//reset button press counters:
	left.downs = 0;
	right.downs = 0;
	up.downs = 0;
	down.downs = 0;

	if ((float)numClaimed / (float)totalVertSize >= stageWinRatio) ifWon = true;

}

void PlayMode::draw(glm::uvec2 const &drawable_size) {
	//update camera aspect ratio for drawable:
	player.camera->aspect = float(drawable_size.x) / float(drawable_size.y);

	//set up light type and position for lit_color_texture_program:
	// TODO: consider using the Light(s) in the scene to do this
	glUseProgram(lit_color_texture_program->program);
	glUniform1i(lit_color_texture_program->LIGHT_TYPE_int, 1);
	glUniform3fv(lit_color_texture_program->LIGHT_DIRECTION_vec3, 1, glm::value_ptr(glm::vec3(0.0f, 0.0f,-1.0f)));
	glUniform3fv(lit_color_texture_program->LIGHT_ENERGY_vec3, 1, glm::value_ptr(glm::vec3(1.0f, 1.0f, 0.95f)));
	glUseProgram(0);

	glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
	glClearDepth(1.0f); //1.0 is actually the default value to clear the depth buffer to, but FYI you can change it.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS); //this is the default depth comparison function, but FYI you can change it.

	//Update walkmesh triangle data
	for (size_t ind = 0; ind < walkmesh->triangles.size(); ind++) {
		Mesh::Vertex vert0 = walkmeshDataVerts[walkmesh->triangles[ind].x];
		Mesh::Vertex vert1 = walkmeshDataVerts[walkmesh->triangles[ind].y];
		Mesh::Vertex vert2 = walkmeshDataVerts[walkmesh->triangles[ind].z];
		walkmeshData[3 * ind] = vert0;
		walkmeshData[3 * ind + 1] = vert1;
		walkmeshData[3 * ind + 2] = vert2;
	}

	//Update vertex colors for walkmesh
	glBindBuffer(GL_ARRAY_BUFFER, walkmeshBuffer);
	GL_ERRORS();
	void* bufferLoc = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
	GL_ERRORS();
	assert(bufferLoc != nullptr);
	memcpy(bufferLoc, walkmeshData.data(), walkmeshData.size() * sizeof(Mesh::Vertex));
	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	scene.draw(*player.camera);
	

	{ //use DrawLines to overlay some text:
		glDisable(GL_DEPTH_TEST);
		float aspect = float(drawable_size.x) / float(drawable_size.y);
		DrawLines lines(glm::mat4(
			1.0f / aspect, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		));

		constexpr float H = 0.09f;
		float currentPercentage = (float)numClaimed / (float)totalVertSize; 
		std::string textString = std::string("");
		if (!ifWon) textString = std::string("Mouse looks; WASD moves; escape ungrabs mouse; Percentage: ").append(std::to_string(currentPercentage))
			.append(std::string("; To win: ")).append(std::to_string(stageWinRatio)).append("Time: ").append(std::to_string((int)(timeLimit - timer)));
		else if (timer < timeLimit) textString = std::string("You won! Percentage: ").append(std::to_string(currentPercentage));
		else textString = std::string("You Lost.");
		lines.draw_text(textString.c_str(),
			glm::vec3(-aspect + 0.1f * H, -1.0 + 0.1f * H, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0x00, 0x00, 0x00, 0x00));
		float ofs = 2.0f / drawable_size.y;
		lines.draw_text(textString.c_str(),
			glm::vec3(-aspect + 0.1f * H + ofs, -1.0 + + 0.1f * H + ofs, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0xff, 0xff, 0xff, 0x00));
	}
	GL_ERRORS();
}
