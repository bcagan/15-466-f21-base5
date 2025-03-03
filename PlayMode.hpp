#include "Mode.hpp"

#include "Scene.hpp"
#include "WalkMesh.hpp"

#include <glm/glm.hpp>

#include <vector>
#include <deque>
#include <glm/gtx/hash.hpp> 

struct PlayMode : Mode {
	PlayMode();
	virtual ~PlayMode();

	//functions called by main loop:
	virtual bool handle_event(SDL_Event const &, glm::uvec2 const &window_size) override;
	virtual void update(float elapsed) override;
	virtual void draw(glm::uvec2 const &drawable_size) override;

	//----- game state -----

	//input tracking:
	struct Button {
		uint8_t downs = 0;
		uint8_t pressed = 0;
	} left, right, down, up;

	//local copy of the game scene (so code can change it during gameplay):
	Scene scene;

	//player info:
	struct Player {
		WalkPoint at;
		//transform is at player's feet and will be yawed by mouse left/right motion:
		Scene::Transform *transform = nullptr;
		//camera is at player's head and will be pitched by mouse up/down motion:
		Scene::Camera *camera = nullptr;
	} player;

	//Vertex information
	size_t numClaimed = 0;
	std::vector<bool> claimedVerts;
	size_t totalVertSize = 0;
	float stageWinRatio = 1.0f;
	bool ifWon = false;
	float timer = 0.f;
	float timeLimit = 15.f;
	Scene::Transform transform;
	Scene::Drawable *walkmeshDrawable;
	GLuint walkmeshBuffer = 0;
	std::vector< Mesh::Vertex > walkmeshDataVerts; //Verts are updated quickly, then they are passed to the triangle data
	std::vector< Mesh::Vertex > walkmeshData;
};
