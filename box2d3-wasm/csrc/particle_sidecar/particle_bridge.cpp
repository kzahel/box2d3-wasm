#include "particle_bridge.h"

#include <box2d/collision.h>
#include <box2d/id.h>
#include <box2d/math_functions.h>

#include <cfloat>
#include <cmath>

namespace physbox3::particle_sidecar {

namespace {

struct QueryForwarder
{
	ParticleShapeQueryCallback callback = nullptr;
	void* context = nullptr;
};

bool ForwardOverlapResult( b2ShapeId shapeId, void* context )
{
	QueryForwarder* forwarder = static_cast<QueryForwarder*>( context );
	return forwarder->callback == nullptr ? true : forwarder->callback( shapeId, forwarder->context );
}

constexpr float kEpsilon = 1.0e-5f;

float Length( b2Vec2 value )
{
	return std::sqrt( value.x * value.x + value.y * value.y );
}

b2Vec2 NormalizeOr( b2Vec2 value, b2Vec2 fallback )
{
	const float length = Length( value );
	if ( length > kEpsilon )
	{
		const float invLength = 1.0f / length;
		return { value.x * invLength, value.y * invLength };
	}

	return fallback;
}

b2Vec2 ClosestPointOnSegment( b2Vec2 a, b2Vec2 b, b2Vec2 point )
{
	const b2Vec2 edge = b2Sub( b, a );
	const float edgeLengthSq = edge.x * edge.x + edge.y * edge.y;
	if ( edgeLengthSq <= kEpsilon * kEpsilon )
	{
		return a;
	}

	float t = b2Dot( b2Sub( point, a ), edge ) / edgeLengthSq;
	t = b2ClampFloat( t, 0.0f, 1.0f );
	return { a.x + t * edge.x, a.y + t * edge.y };
}

b2Vec2 RightNormal( b2Vec2 edge )
{
	return NormalizeOr( { edge.y, -edge.x }, { 0.0f, 1.0f } );
}

ParticleShapeContact MakeOutsideContact( b2ShapeId shapeId, b2BodyId bodyId, b2Vec2 point, b2Vec2 closestPoint )
{
	const b2Vec2 delta = b2Sub( point, closestPoint );
	const b2Vec2 normal = NormalizeOr( delta, { 0.0f, 1.0f } );
	return {
		.shapeId = shapeId,
		.bodyId = bodyId,
		.point = closestPoint,
		.normal = normal,
		.signedDistance = Length( delta ),
	};
}

ParticleShapeContact MakeCircleInsideContact( b2ShapeId shapeId, b2BodyId bodyId, const b2Circle& circle, b2Vec2 point )
{
	const b2Vec2 center = b2Body_GetWorldPoint( bodyId, circle.center );
	const b2Vec2 delta = b2Sub( point, center );
	const b2Vec2 normal = NormalizeOr( delta, { 0.0f, 1.0f } );
	const float signedDistance = Length( delta ) - circle.radius;
	return {
		.shapeId = shapeId,
		.bodyId = bodyId,
		.point = b2Sub( point, b2MulSV( signedDistance, normal ) ),
		.normal = normal,
		.signedDistance = signedDistance,
	};
}

ParticleShapeContact MakeCapsuleInsideContact( b2ShapeId shapeId, b2BodyId bodyId, const b2Capsule& capsule, b2Vec2 point )
{
	const b2Vec2 center1 = b2Body_GetWorldPoint( bodyId, capsule.center1 );
	const b2Vec2 center2 = b2Body_GetWorldPoint( bodyId, capsule.center2 );
	const b2Vec2 center = ClosestPointOnSegment( center1, center2, point );
	const b2Vec2 fallbackNormal = RightNormal( b2Sub( center2, center1 ) );
	const b2Vec2 delta = b2Sub( point, center );
	const b2Vec2 normal = NormalizeOr( delta, fallbackNormal );
	const float signedDistance = Length( delta ) - capsule.radius;
	return {
		.shapeId = shapeId,
		.bodyId = bodyId,
		.point = b2Sub( point, b2MulSV( signedDistance, normal ) ),
		.normal = normal,
		.signedDistance = signedDistance,
	};
}

ParticleShapeContact MakePolygonInsideContact( b2ShapeId shapeId, b2BodyId bodyId, const b2Polygon& polygon, b2Vec2 point )
{
	float maxSeparation = -FLT_MAX;
	b2Vec2 bestNormal = { 0.0f, 1.0f };

	for ( int i = 0; i < polygon.count; ++i )
	{
		const b2Vec2 worldVertex = b2Body_GetWorldPoint( bodyId, polygon.vertices[i] );
		const b2Vec2 worldNormal = NormalizeOr( b2Body_GetWorldVector( bodyId, polygon.normals[i] ), bestNormal );
		const float separation = b2Dot( worldNormal, b2Sub( point, worldVertex ) );
		if ( separation > maxSeparation )
		{
			maxSeparation = separation;
			bestNormal = worldNormal;
		}
	}

	const float signedDistance = maxSeparation - polygon.radius;
	return {
		.shapeId = shapeId,
		.bodyId = bodyId,
		.point = b2Sub( point, b2MulSV( signedDistance, bestNormal ) ),
		.normal = bestNormal,
		.signedDistance = signedDistance,
	};
}

ParticleShapeContact MakeSegmentContact( b2ShapeId shapeId, b2BodyId bodyId, b2Vec2 point1, b2Vec2 point2, b2Vec2 point )
{
	const b2Vec2 closestPoint = ClosestPointOnSegment( point1, point2, point );
	const b2Vec2 fallbackNormal = RightNormal( b2Sub( point2, point1 ) );
	const b2Vec2 delta = b2Sub( point, closestPoint );
	return {
		.shapeId = shapeId,
		.bodyId = bodyId,
		.point = closestPoint,
		.normal = NormalizeOr( delta, fallbackNormal ),
		.signedDistance = Length( delta ),
	};
}

bool MakeChainSegmentContact( b2ShapeId shapeId, b2BodyId bodyId, const b2ChainSegment& chainSegment, b2Vec2 point,
							  ParticleShapeContact* outContact )
{
	const b2Vec2 point1 = b2Body_GetWorldPoint( bodyId, chainSegment.segment.point1 );
	const b2Vec2 point2 = b2Body_GetWorldPoint( bodyId, chainSegment.segment.point2 );
	const b2Vec2 segment = b2Sub( point2, point1 );
	const b2Vec2 surfaceNormal = RightNormal( segment );

	if ( b2Dot( b2Sub( point, point1 ), surfaceNormal ) <= 0.0f )
	{
		return false;
	}

	ParticleShapeContact contact = MakeSegmentContact( shapeId, bodyId, point1, point2, point );
	if ( b2Dot( contact.normal, surfaceNormal ) < 0.0f )
	{
		contact.normal = surfaceNormal;
	}

	*outContact = contact;
	return true;
}

} // namespace

Box2DParticleWorldBridge::Box2DParticleWorldBridge( b2WorldId worldId ) : world_id_( worldId )
{
}

b2WorldId Box2DParticleWorldBridge::GetWorldId() const
{
	return world_id_;
}

b2Vec2 Box2DParticleWorldBridge::GetGravity() const
{
	return b2World_GetGravity( world_id_ );
}

void Box2DParticleWorldBridge::StepWorld( float timeStep, int subStepCount ) const
{
	b2World_Step( world_id_, timeStep, subStepCount );
}

void Box2DParticleWorldBridge::QueryShapesInAABB( const b2AABB& aabb, ParticleShapeQueryCallback callback, void* context ) const
{
	QueryForwarder forwarder = {
		.callback = callback,
		.context = context,
	};

	b2World_OverlapAABB( world_id_, aabb, b2DefaultQueryFilter(), ForwardOverlapResult, &forwarder );
}

bool Box2DParticleWorldBridge::ComputeParticleShapeContact( b2ShapeId shapeId, b2Vec2 particlePosition, float particleRadius,
															ParticleShapeContact* outContact ) const
{
	if ( outContact == nullptr || !B2_IS_NON_NULL( shapeId ) || b2Shape_IsSensor( shapeId ) )
	{
		return false;
	}

	const b2BodyId bodyId = b2Shape_GetBody( shapeId );
	const b2ShapeType shapeType = b2Shape_GetType( shapeId );
	ParticleShapeContact contact;

	switch ( shapeType )
	{
		case b2_circleShape:
		{
			if ( b2Shape_TestPoint( shapeId, particlePosition ) )
			{
				contact = MakeCircleInsideContact( shapeId, bodyId, b2Shape_GetCircle( shapeId ), particlePosition );
			}
			else
			{
				contact = MakeOutsideContact( shapeId, bodyId, particlePosition, b2Shape_GetClosestPoint( shapeId, particlePosition ) );
			}
			break;
		}

		case b2_capsuleShape:
		{
			if ( b2Shape_TestPoint( shapeId, particlePosition ) )
			{
				contact = MakeCapsuleInsideContact( shapeId, bodyId, b2Shape_GetCapsule( shapeId ), particlePosition );
			}
			else
			{
				contact = MakeOutsideContact( shapeId, bodyId, particlePosition, b2Shape_GetClosestPoint( shapeId, particlePosition ) );
			}
			break;
		}

		case b2_polygonShape:
		{
			if ( b2Shape_TestPoint( shapeId, particlePosition ) )
			{
				contact = MakePolygonInsideContact( shapeId, bodyId, b2Shape_GetPolygon( shapeId ), particlePosition );
			}
			else
			{
				contact = MakeOutsideContact( shapeId, bodyId, particlePosition, b2Shape_GetClosestPoint( shapeId, particlePosition ) );
			}
			break;
		}

		case b2_segmentShape:
		{
			const b2Segment segment = b2Shape_GetSegment( shapeId );
			const b2Vec2 point1 = b2Body_GetWorldPoint( bodyId, segment.point1 );
			const b2Vec2 point2 = b2Body_GetWorldPoint( bodyId, segment.point2 );
			contact = MakeSegmentContact( shapeId, bodyId, point1, point2, particlePosition );
			break;
		}

		case b2_chainSegmentShape:
		{
			if ( !MakeChainSegmentContact( shapeId, bodyId, b2Shape_GetChainSegment( shapeId ), particlePosition, &contact ) )
			{
				return false;
			}
			break;
		}

		default:
			return false;
	}

	if ( contact.signedDistance > particleRadius )
	{
		return false;
	}

	*outContact = contact;
	return true;
}

b2Vec2 Box2DParticleWorldBridge::GetBodyPointVelocity( b2BodyId bodyId, b2Vec2 worldPoint ) const
{
	if ( !B2_IS_NON_NULL( bodyId ) )
	{
		return { 0.0f, 0.0f };
	}

	return b2Body_GetWorldPointVelocity( bodyId, worldPoint );
}

void Box2DParticleWorldBridge::ApplyBodyLinearImpulse( b2BodyId bodyId, b2Vec2 impulse, b2Vec2 worldPoint, bool wake ) const
{
	if ( !B2_IS_NON_NULL( bodyId ) )
	{
		return;
	}

	b2Body_ApplyLinearImpulse( bodyId, impulse, worldPoint, wake );
}

} // namespace physbox3::particle_sidecar
