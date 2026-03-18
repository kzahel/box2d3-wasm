#include "particle_system.h"

#include <box2d/math_functions.h>

#include <algorithm>
#include <cmath>
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

constexpr float kContactRadiusMultiplier = 3.0f;
constexpr int kPressureIterations = 2;
constexpr float kRestDensity = 0.45f;
constexpr float kPressureStrength = 140.0f;
constexpr float kNearPressureStrength = 260.0f;
constexpr float kDampingStrength = 0.18f;
constexpr float kMaxPressureDisplacementMultiplier = 0.75f;
constexpr float kContactEpsilon = 1.0e-5f;

} // namespace

ParticleSystemSidecar::ParticleSystemSidecar( b2::World& world, const ParticleSystemDef& def ) : bridge_( world.Handle() )
{
	particleRadius_ = def.radius > 0.0f ? def.radius : 0.1f;
	particleDiameter_ = particleRadius_ * 2.0f;
	particleMass_ = def.density * std::numbers::pi_v<float> * particleRadius_ * particleRadius_;
	contactRadius_ = particleRadius_ * kContactRadiusMultiplier;
	gravityScale_ = def.gravityScale;
	maxParticles_ = def.maxParticles;

	if ( def.initialCapacity > 0 )
	{
		positions_.reserve( def.initialCapacity );
		velocities_.reserve( def.initialCapacity );
		previousPositions_.reserve( def.initialCapacity );
		densityScratch_.reserve( def.initialCapacity );
		nearDensityScratch_.reserve( def.initialCapacity );
		proxies_.reserve( def.initialCapacity );
		contacts_.reserve( def.initialCapacity * 8 );
	}
	queryScratch_.reserve( 16 );
	cellRanges_.reserve( 64 );
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

int ParticleSystemSidecar::GetMaxParticles() const
{
	return maxParticles_;
}

void ParticleSystemSidecar::SetMaxParticles( int maxParticles )
{
	maxParticles_ = maxParticles > 0 ? maxParticles : 0;
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

int ParticleSystemSidecar::DestroyParticlesInCircle( b2Vec2 center, float radius )
{
	if ( !IsValid() || radius <= 0.0f || positions_.empty() )
	{
		return 0;
	}

	const float radiusSq = radius * radius;
	std::size_t writeIndex = 0;
	int removed = 0;

	for ( std::size_t readIndex = 0; readIndex < positions_.size(); ++readIndex )
	{
		const b2Vec2 delta = b2Sub( positions_[readIndex], center );
		if ( b2Dot( delta, delta ) <= radiusSq )
		{
			++removed;
			continue;
		}

		if ( writeIndex != readIndex )
		{
			positions_[writeIndex] = positions_[readIndex];
			velocities_[writeIndex] = velocities_[readIndex];
		}
		++writeIndex;
	}

	if ( removed == 0 )
	{
		return 0;
	}

	positions_.resize( writeIndex );
	velocities_.resize( writeIndex );
	previousPositions_.clear();
	densityScratch_.clear();
	nearDensityScratch_.clear();
	proxies_.clear();
	cellRanges_.clear();
	contacts_.clear();

	return removed;
}

void ParticleSystemSidecar::Step( float timeStep )
{
	if ( !IsValid() || timeStep <= 0.0f )
	{
		return;
	}

	ResizeSolverScratch();
	const b2Vec2 gravity = b2MulSV( gravityScale_, bridge_.GetGravity() );

	for ( std::size_t i = 0; i < positions_.size(); ++i )
	{
		previousPositions_[i] = positions_[i];
		velocities_[i] = b2MulAdd( velocities_[i], timeStep, gravity );
		positions_[i] = b2MulAdd( positions_[i], timeStep, velocities_[i] );
	}

	UpdateParticleContacts();
	SolveParticlePressure( timeStep );
	UpdateVelocitiesFromPositions( 1.0f / timeStep );
	SolveParticleDamping( timeStep );

	for ( std::size_t i = 0; i < positions_.size(); ++i )
	{
		ResolveBodyContacts( &positions_[i], &velocities_[i] );
	}
}

void ParticleSystemSidecar::StepWithWorld( float timeStep, int subStepCount )
{
	if ( !IsValid() || timeStep <= 0.0f )
	{
		return;
	}

	const int clampedSubStepCount = subStepCount > 0 ? subStepCount : 1;
	const float subTimeStep = timeStep / static_cast<float>( clampedSubStepCount );

	for ( int i = 0; i < clampedSubStepCount; ++i )
	{
		Step( subTimeStep );
		bridge_.StepWorld( subTimeStep, 1 );
	}
}

void ParticleSystemSidecar::Clear()
{
	positions_.clear();
	velocities_.clear();
	previousPositions_.clear();
	densityScratch_.clear();
	nearDensityScratch_.clear();
	proxies_.clear();
	cellRanges_.clear();
	contacts_.clear();
	queryScratch_.clear();
}

void ParticleSystemSidecar::ResizeSolverScratch()
{
	const std::size_t particleCount = positions_.size();
	previousPositions_.resize( particleCount );
	densityScratch_.resize( particleCount );
	nearDensityScratch_.resize( particleCount );
	proxies_.resize( particleCount );
}

void ParticleSystemSidecar::UpdateParticleContacts()
{
	contacts_.clear();
	cellRanges_.clear();

	const int particleCount = static_cast<int>( positions_.size() );
	if ( particleCount < 2 )
	{
		return;
	}

	for ( int i = 0; i < particleCount; ++i )
	{
		const b2Vec2 position = positions_[i];
		const int cellX = ComputeCellCoordinate( position.x, contactRadius_ );
		const int cellY = ComputeCellCoordinate( position.y, contactRadius_ );
		proxies_[i] = {
			.tag = ComputeCellTag( cellX, cellY ),
			.index = i,
		};
	}

	std::sort( proxies_.begin(), proxies_.end(), []( const ParticleProxy& lhs, const ParticleProxy& rhs ) {
		if ( lhs.tag != rhs.tag )
		{
			return lhs.tag < rhs.tag;
		}

		return lhs.index < rhs.index;
	} );

	for ( int begin = 0; begin < particleCount; )
	{
		const std::uint64_t tag = proxies_[begin].tag;
		int end = begin + 1;
		while ( end < particleCount && proxies_[end].tag == tag )
		{
			++end;
		}

		cellRanges_.push_back( {
			.tag = tag,
			.begin = begin,
			.end = end,
		} );

		begin = end;
	}

	const float maxDistanceSq = contactRadius_ * contactRadius_;
	contacts_.reserve( positions_.size() * 8 );

	for ( int i = 0; i < particleCount; ++i )
	{
		const b2Vec2 position = positions_[i];
		const int cellX = ComputeCellCoordinate( position.x, contactRadius_ );
		const int cellY = ComputeCellCoordinate( position.y, contactRadius_ );

		for ( int offsetY = -1; offsetY <= 1; ++offsetY )
		{
			for ( int offsetX = -1; offsetX <= 1; ++offsetX )
			{
				const ParticleCellRange* cell = FindCellRange( ComputeCellTag( cellX + offsetX, cellY + offsetY ) );
				if ( cell == nullptr )
				{
					continue;
				}

				for ( int proxyIndex = cell->begin; proxyIndex < cell->end; ++proxyIndex )
				{
					const int otherIndex = proxies_[proxyIndex].index;
					if ( otherIndex <= i )
					{
						continue;
					}

					const b2Vec2 delta = b2Sub( positions_[otherIndex], position );
					const float distanceSq = b2Dot( delta, delta );
					if ( distanceSq > maxDistanceSq )
					{
						continue;
					}

					ParticleContact contact = {
						.indexA = i,
						.indexB = otherIndex,
					};
					if ( RefreshContactGeometry( &contact ) )
					{
						contacts_.push_back( contact );
					}
				}
			}
		}
	}
}

void ParticleSystemSidecar::SolveParticlePressure( float timeStep )
{
	if ( contacts_.empty() )
	{
		return;
	}

	const float maxDisplacement = particleRadius_ * kMaxPressureDisplacementMultiplier;
	const float pressureStep = timeStep * timeStep;

	for ( int iteration = 0; iteration < kPressureIterations; ++iteration )
	{
		std::fill( densityScratch_.begin(), densityScratch_.end(), 0.0f );
		std::fill( nearDensityScratch_.begin(), nearDensityScratch_.end(), 0.0f );

		for ( ParticleContact& contact : contacts_ )
		{
			if ( !RefreshContactGeometry( &contact ) )
			{
				continue;
			}

			const float q = contact.weight;
			const float q2 = q * q;
			densityScratch_[contact.indexA] += q2;
			densityScratch_[contact.indexB] += q2;
			nearDensityScratch_[contact.indexA] += q2 * q;
			nearDensityScratch_[contact.indexB] += q2 * q;
		}

		for ( ParticleContact& contact : contacts_ )
		{
			if ( contact.weight <= 0.0f )
			{
				continue;
			}

			const float q = contact.weight;
			const float pressureA = kPressureStrength * std::max( 0.0f, densityScratch_[contact.indexA] - kRestDensity );
			const float pressureB = kPressureStrength * std::max( 0.0f, densityScratch_[contact.indexB] - kRestDensity );
			const float nearPressureA = kNearPressureStrength * nearDensityScratch_[contact.indexA];
			const float nearPressureB = kNearPressureStrength * nearDensityScratch_[contact.indexB];
			float displacementMagnitude =
				pressureStep * ( ( pressureA + pressureB ) * q + ( nearPressureA + nearPressureB ) * q * q );
			displacementMagnitude = std::min( displacementMagnitude, maxDisplacement );
			if ( displacementMagnitude <= 0.0f )
			{
				continue;
			}

			const b2Vec2 offset = b2MulSV( 0.5f * displacementMagnitude, contact.normal );
			positions_[contact.indexA] = b2Sub( positions_[contact.indexA], offset );
			positions_[contact.indexB] = b2Add( positions_[contact.indexB], offset );
		}
	}
}

void ParticleSystemSidecar::SolveParticleDamping( float timeStep )
{
	if ( contacts_.empty() )
	{
		return;
	}

	const float velocityScale = particleDiameter_ > 0.0f ? timeStep / particleDiameter_ : 0.0f;

	for ( ParticleContact& contact : contacts_ )
	{
		if ( !RefreshContactGeometry( &contact ) )
		{
			continue;
		}

		const b2Vec2 relativeVelocity = b2Sub( velocities_[contact.indexB], velocities_[contact.indexA] );
		const float normalSpeed = b2Dot( relativeVelocity, contact.normal );
		if ( normalSpeed >= 0.0f )
		{
			continue;
		}

		const float damping = b2ClampFloat(
			kDampingStrength * contact.weight + ( -normalSpeed * velocityScale ) * 0.02f,
			0.0f,
			0.5f
		);
		const b2Vec2 impulse = b2MulSV( 0.5f * damping * normalSpeed, contact.normal );
		velocities_[contact.indexA] = b2Add( velocities_[contact.indexA], impulse );
		velocities_[contact.indexB] = b2Sub( velocities_[contact.indexB], impulse );
	}
}

void ParticleSystemSidecar::UpdateVelocitiesFromPositions( float invTimeStep )
{
	for ( std::size_t i = 0; i < positions_.size(); ++i )
	{
		velocities_[i] = b2MulSV( invTimeStep, b2Sub( positions_[i], previousPositions_[i] ) );
	}
}

bool ParticleSystemSidecar::RefreshContactGeometry( ParticleContact* contact ) const
{
	if ( contact == nullptr )
	{
		return false;
	}

	const b2Vec2 delta = b2Sub( positions_[contact->indexB], positions_[contact->indexA] );
	const float distanceSq = b2Dot( delta, delta );
	const float maxDistanceSq = contactRadius_ * contactRadius_;
	if ( distanceSq > maxDistanceSq )
	{
		contact->distance = std::sqrt( distanceSq );
		contact->weight = 0.0f;
		return false;
	}

	if ( distanceSq > kContactEpsilon * kContactEpsilon )
	{
		contact->distance = std::sqrt( distanceSq );
		contact->normal = b2MulSV( 1.0f / contact->distance, delta );
	}
	else
	{
		contact->distance = 0.0f;
		contact->normal = { 1.0f, 0.0f };
	}

	contact->weight = std::max( 0.0f, 1.0f - contact->distance / contactRadius_ );
	return contact->weight > 0.0f;
}

int ParticleSystemSidecar::ComputeCellCoordinate( float value, float cellSize )
{
	return static_cast<int>( std::floor( value / cellSize ) );
}

std::uint64_t ParticleSystemSidecar::ComputeCellTag( int cellX, int cellY )
{
	return ( static_cast<std::uint64_t>( static_cast<std::uint32_t>( cellX ) ) << 32 ) |
		   static_cast<std::uint32_t>( cellY );
}

const ParticleCellRange* ParticleSystemSidecar::FindCellRange( std::uint64_t tag ) const
{
	const auto it = std::lower_bound(
		cellRanges_.begin(),
		cellRanges_.end(),
		tag,
		[]( const ParticleCellRange& cell, std::uint64_t value ) { return cell.tag < value; }
	);
	if ( it == cellRanges_.end() || it->tag != tag )
	{
		return nullptr;
	}

	return &*it;
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
