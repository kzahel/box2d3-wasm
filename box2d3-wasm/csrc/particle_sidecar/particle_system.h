#pragma once

#include "particle_bridge.h"

#include <box2cpp/box2cpp.h>

#include <vector>

namespace physbox3::particle_sidecar {

struct ParticleSystemDef
{
	float radius = 0.1f;
	float density = 1.0f;
	float gravityScale = 1.0f;
	int initialCapacity = 256;
	int maxParticles = 0;
};

class ParticleSystemSidecar
{
public:
	ParticleSystemSidecar( b2::World& world, const ParticleSystemDef& def );

	void Destroy();
	bool IsValid() const;

	int GetParticleCount() const;
	float GetParticleRadius() const;
	const b2Vec2* GetPositionData() const;

	int SpawnParticlesInCircle( b2Vec2 center, float radius, float spacing, b2Vec2 initialVelocity );
	void Step( float timeStep );
	void Clear();

private:
	bool CanCreateParticleAt( b2Vec2 position );
	void ResolveBodyContacts( b2Vec2* position, b2Vec2* velocity );
	void QueryNearbyShapes( b2Vec2 position, float radius );

	bool valid_ = true;
	Box2DParticleWorldBridge bridge_;
	float particleRadius_ = 0.1f;
	float particleMass_ = 0.0f;
	float gravityScale_ = 1.0f;
	int maxParticles_ = 0;
	std::vector<b2Vec2> positions_;
	std::vector<b2Vec2> velocities_;
	std::vector<b2ShapeId> queryScratch_;
};

ParticleSystemSidecar* CreateParticleSystem( b2::World& world, const ParticleSystemDef& def );
void DestroyParticleSystem( ParticleSystemSidecar* system );

} // namespace physbox3::particle_sidecar
