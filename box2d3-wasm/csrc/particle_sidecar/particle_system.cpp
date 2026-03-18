#include "particle_system.h"

#include <box2d/math_functions.h>

#include <numbers>

namespace physbox3::particle_sidecar {

namespace {

static_assert( sizeof( b2Vec2 ) == sizeof( float ) * 2, "Particle position buffer assumes packed b2Vec2 values." );

struct ShapeCollector
{
	std::vector<b2ShapeId>* shapes = nullptr;
};

bool CollectShapeId( b2ShapeId shapeId, void* context )
{
	ShapeCollector* collector = static_cast<ShapeCollector*>( context );
	collector->shapes->push_back( shapeId );
	return true;
}

} // namespace

ParticleSystemSidecar::ParticleSystemSidecar( b2::World& world, const ParticleSystemDef& def ) : bridge_( world.Handle() )
{
	particleRadius_ = def.radius > 0.0f ? def.radius : 0.1f;
	particleMass_ = def.density * std::numbers::pi_v<float> * particleRadius_ * particleRadius_;
	gravityScale_ = def.gravityScale;
	maxParticles_ = def.maxParticles;

	if ( def.initialCapacity > 0 )
	{
		positions_.reserve( def.initialCapacity );
		velocities_.reserve( def.initialCapacity );
	}
	queryScratch_.reserve( 16 );
}

void ParticleSystemSidecar::Destroy()
{
	Clear();
	valid_ = false;
}

bool ParticleSystemSidecar::IsValid() const
{
	return valid_ && b2World_IsValid( bridge_.GetWorldId() );
}

int ParticleSystemSidecar::GetParticleCount() const
{
	return static_cast<int>( positions_.size() );
}

float ParticleSystemSidecar::GetParticleRadius() const
{
	return particleRadius_;
}

const b2Vec2* ParticleSystemSidecar::GetPositionData() const
{
	return positions_.data();
}

int ParticleSystemSidecar::SpawnParticlesInCircle( b2Vec2 center, float radius, float spacing, b2Vec2 initialVelocity )
{
	if ( !IsValid() || radius <= 0.0f )
	{
		return 0;
	}

	const float stride = spacing > 0.0f ? spacing : particleRadius_ * 2.0f;
	int created = 0;

	for ( float y = -radius; y <= radius; y += stride )
	{
		for ( float x = -radius; x <= radius; x += stride )
		{
			if ( maxParticles_ > 0 && static_cast<int>( positions_.size() ) >= maxParticles_ )
			{
				return created;
			}

			const b2Vec2 local = { x, y };
			if ( b2Dot( local, local ) > radius * radius )
			{
				continue;
			}

			const b2Vec2 position = { center.x + x, center.y + y };
			if ( !CanCreateParticleAt( position ) )
			{
				continue;
			}

			positions_.push_back( position );
			velocities_.push_back( initialVelocity );
			++created;
		}
	}

	return created;
}

void ParticleSystemSidecar::Step( float timeStep )
{
	if ( !IsValid() || timeStep <= 0.0f )
	{
		return;
	}

	const b2Vec2 gravity = b2MulSV( gravityScale_, bridge_.GetGravity() );

	for ( std::size_t i = 0; i < positions_.size(); ++i )
	{
		b2Vec2 velocity = velocities_[i];
		velocity = b2MulAdd( velocity, timeStep, gravity );

		b2Vec2 position = positions_[i];
		position = b2MulAdd( position, timeStep, velocity );

		ResolveBodyContacts( &position, &velocity );

		positions_[i] = position;
		velocities_[i] = velocity;
	}
}

void ParticleSystemSidecar::Clear()
{
	positions_.clear();
	velocities_.clear();
	queryScratch_.clear();
}

bool ParticleSystemSidecar::CanCreateParticleAt( b2Vec2 position )
{
	QueryNearbyShapes( position, particleRadius_ );

	for ( const b2ShapeId shapeId : queryScratch_ )
	{
		ParticleShapeContact contact;
		if ( bridge_.ComputeParticleShapeContact( shapeId, position, particleRadius_, &contact ) )
		{
			return false;
		}
	}

	return true;
}

void ParticleSystemSidecar::ResolveBodyContacts( b2Vec2* position, b2Vec2* velocity )
{
	QueryNearbyShapes( *position, particleRadius_ );

	for ( const b2ShapeId shapeId : queryScratch_ )
	{
		ParticleShapeContact contact;
		if ( !bridge_.ComputeParticleShapeContact( shapeId, *position, particleRadius_, &contact ) )
		{
			continue;
		}

		const float correction = particleRadius_ - contact.signedDistance;
		if ( correction <= 0.0f )
		{
			continue;
		}

		*position = b2MulAdd( *position, correction, contact.normal );

		const b2Vec2 bodyVelocity = bridge_.GetBodyPointVelocity( contact.bodyId, contact.point );
		const b2Vec2 relativeVelocity = b2Sub( *velocity, bodyVelocity );
		const float normalSpeed = b2Dot( relativeVelocity, contact.normal );

		if ( normalSpeed < 0.0f )
		{
			*velocity = b2MulAdd( *velocity, -normalSpeed, contact.normal );
			const b2Vec2 bodyImpulse = b2MulSV( particleMass_ * normalSpeed, contact.normal );
			bridge_.ApplyBodyLinearImpulse( contact.bodyId, bodyImpulse, contact.point, true );
		}
	}
}

void ParticleSystemSidecar::QueryNearbyShapes( b2Vec2 position, float radius )
{
	queryScratch_.clear();

	b2AABB aabb;
	aabb.lowerBound = { position.x - radius, position.y - radius };
	aabb.upperBound = { position.x + radius, position.y + radius };

	ShapeCollector collector = {
		.shapes = &queryScratch_,
	};
	bridge_.QueryShapesInAABB( aabb, CollectShapeId, &collector );
}

ParticleSystemSidecar* CreateParticleSystem( b2::World& world, const ParticleSystemDef& def )
{
	return new ParticleSystemSidecar( world, def );
}

void DestroyParticleSystem( ParticleSystemSidecar* system )
{
	delete system;
}

} // namespace physbox3::particle_sidecar
