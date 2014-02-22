

#ifndef CORE_INTERSECTION
#define CORE_INTERSECTION

#include "Object.h"
#include "Ray.h"
#include <stdlib.h>
#include <stdio.h>

namespace Core
{
	/***************************************************************
	Intersection class which stores information about a 
	ray-object intersection.
	***************************************************************/
	class Intersection
	{
	public:
		Intersection( Object obj = NULL, Ray ray=NULL, float dist = 0.0f )
		{
			this->mObject = obj;
			this->mRay = ray;
			this->mDist = dist;
		}


		Object object() const {return mObject;}
		Ray ray() const {return mRay;}
		float distance const {return mDist;}

	protected:
		Object			mObject;			// The object the ray hit.
		Ray				mRay;				// The ray that hit an object.
		float			mDist;				// The distance from the ray origin to the hit point.
	};


	/***************************************************************
	IntersectionDetector class.
	The class is used to detect ray-scene intersections.
	***************************************************************/
	class IntersectionDetector
	{
		public detectIntersect( Ray ray, /*Scene secen*/ )
		{

		}
	};
}

#endif