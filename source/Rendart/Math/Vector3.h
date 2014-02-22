/* Vector3 class
**
** Author: Yinan Zhang
** Time: Feb. 22, 2014
**
*/

#ifndef MATH_VECTOR3
#define MATH_VECTOR3

#include "glm/glm.hpp"

namespace Math
{
	class Vector3
	{
	public:
		Vector3( float x, float y, float z )
		{
			this->mVector = glm::vec3( x, y, z );
		}

		Vector3 operator+( const Vector3 &vect ) const
		{
			return Vector3( this->mVector + vect.mVector );
		}

		Vector3 operator-( const Vector3 &vect ) const
		{
			return Vector3( this->mVector - vect.mVector );
		}

		Vector3 operator*( const float val ) const
		{
			return Vector3( this->mVector * val );
		}

		Vector3 operator*( const Vector3 &vect ) const
		{
			return Vector3( this->mVector * vect.mVector);
		}

		float dot( const Vector3 &vect ) const
		{
			return glm::dot( this->mVector, vect.mVector );
		}

		Vector3 cross( const Vector3 &vect ) const
		{
			return Vector3( glm::cross( this->mVector, vect.mVector) );
		}

		/*
		Return the normalized value. The function will not change the value of the instalce
		*/
		Vector3 normalize( void ) const
		{
			return Vector3( glm::normalize( this->mVector ) );
		}

		float X() const{ return this->mVector.x; }
		float Y() const{ return this->mVector.y; }
		float Z() const{ return this->mVector.z; }

	private:
		Vector3( glm::vec3 vect )
		{
			this->mVector = vect;
		}

		glm::vec3 mVector;
	};

}

#endif