

#ifndef CORE_RAY
#define CORE_RAY
#include "Math/Vector3.h"


namespace Core
{
	/***************************************************************
	Ray class.
		A ray is defined by two features: the origin point and the 
	direction of the ray.
	***************************************************************/
	class Ray
	{
	public:
		Ray(){}

		Ray(Math::Vector3 o, Math::Vector dir)
		{
			this->mOrigin = o;
			this->mDirect = dir;
		}

		Math::Vector3 origin() const {return this->mOrigin;}
		Math::Vector3 direct() const {return this->mDirect;}

	private:
		Math::Vector3 mOrigin;
		Math::Vector3 mDirect;
	};
}


#endif