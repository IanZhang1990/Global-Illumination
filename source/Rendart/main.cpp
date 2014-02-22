

#include <stdlib.h>
#include "Math/Vector3.h"
#include "Math/glm/glm.hpp"


int main()
{
	Math::Vector3 vect( 1,1,1 );

	glm::vec3 vect1(1,2,3);
	glm::vec3 vect2(4,5,6);

	float a = glm::dot( vect1, vect2 );

	return 0;
}