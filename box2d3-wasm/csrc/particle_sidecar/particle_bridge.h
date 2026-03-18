#pragma once

#include <box2d/box2d.h>

namespace physbox3::particle_sidecar {

using ParticleShapeQueryCallback = bool ( * )( b2ShapeId shapeId, void* context );

struct ParticleShapeContact
{
	b2ShapeId shapeId = b2_nullShapeId;
	b2BodyId bodyId = b2_nullBodyId;
	b2Vec2 point = { 0.0f, 0.0f };
	b2Vec2 normal = { 0.0f, 1.0f };
	float signedDistance = 0.0f;
};

// The reduced solver talks only to this interface. Box2D v3 flat-C adaptation
// stays on the bridge side so future upstream updates are localized here.
class ParticleWorldBridge
{
public:
	virtual ~ParticleWorldBridge() = default;

	virtual b2Vec2 GetGravity() const = 0;
	virtual void QueryShapesInAABB( const b2AABB& aabb, ParticleShapeQueryCallback callback, void* context ) const = 0;
	virtual bool ComputeParticleShapeContact( b2ShapeId shapeId, b2Vec2 particlePosition, float particleRadius,
											  ParticleShapeContact* outContact ) const = 0;
	virtual b2Vec2 GetBodyPointVelocity( b2BodyId bodyId, b2Vec2 worldPoint ) const = 0;
	virtual void ApplyBodyLinearImpulse( b2BodyId bodyId, b2Vec2 impulse, b2Vec2 worldPoint, bool wake ) const = 0;
};

class Box2DParticleWorldBridge final : public ParticleWorldBridge
{
public:
	explicit Box2DParticleWorldBridge( b2WorldId worldId );

	b2WorldId GetWorldId() const;

	b2Vec2 GetGravity() const override;
	void QueryShapesInAABB( const b2AABB& aabb, ParticleShapeQueryCallback callback, void* context ) const override;
	bool ComputeParticleShapeContact( b2ShapeId shapeId, b2Vec2 particlePosition, float particleRadius,
									  ParticleShapeContact* outContact ) const override;
	b2Vec2 GetBodyPointVelocity( b2BodyId bodyId, b2Vec2 worldPoint ) const override;
	void ApplyBodyLinearImpulse( b2BodyId bodyId, b2Vec2 impulse, b2Vec2 worldPoint, bool wake ) const override;

private:
	b2WorldId world_id_;
};

} // namespace physbox3::particle_sidecar
