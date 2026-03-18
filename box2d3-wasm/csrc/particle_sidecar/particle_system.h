#pragma once

#include "particle_bridge.h"

#include <box2cpp/box2cpp.h>

#include <cstdint>
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

struct ParticleProxy
{
	std::uint64_t tag = 0;
	int index = 0;
};

struct ParticleCellRange
{
	std::uint64_t tag = 0;
	int begin = 0;
	int end = 0;
};

struct ParticleContact
{
	int indexA = 0;
	int indexB = 0;
	b2Vec2 normal = { 1.0f, 0.0f };
	float distance = 0.0f;
	float weight = 0.0f;
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
	void StepWithWorld( float timeStep, int subStepCount );
	void Clear();

private:
	void ResizeSolverScratch();
	void UpdateParticleContacts();
	void SolveParticlePressure( float timeStep );
	void SolveParticleDamping( float timeStep );
	void UpdateVelocitiesFromPositions( float invTimeStep );
	bool RefreshContactGeometry( ParticleContact* contact ) const;
	static int ComputeCellCoordinate( float value, float cellSize );
	static std::uint64_t ComputeCellTag( int cellX, int cellY );
	const ParticleCellRange* FindCellRange( std::uint64_t tag ) const;
	bool CanCreateParticleAt( b2Vec2 position );
	void ResolveBodyContacts( b2Vec2* position, b2Vec2* velocity );
	void QueryNearbyShapes( b2Vec2 position, float radius );

	bool valid_ = true;
	Box2DParticleWorldBridge bridge_;
	float particleRadius_ = 0.1f;
	float particleDiameter_ = 0.2f;
	float particleMass_ = 0.0f;
	float contactRadius_ = 0.3f;
	float gravityScale_ = 1.0f;
	int maxParticles_ = 0;
	std::vector<b2Vec2> positions_;
	std::vector<b2Vec2> velocities_;
	std::vector<b2Vec2> previousPositions_;
	std::vector<float> densityScratch_;
	std::vector<float> nearDensityScratch_;
	std::vector<ParticleProxy> proxies_;
	std::vector<ParticleCellRange> cellRanges_;
	std::vector<ParticleContact> contacts_;
	std::vector<b2ShapeId> queryScratch_;
};

ParticleSystemSidecar* CreateParticleSystem( b2::World& world, const ParticleSystemDef& def );
void DestroyParticleSystem( ParticleSystemSidecar* system );

} // namespace physbox3::particle_sidecar
