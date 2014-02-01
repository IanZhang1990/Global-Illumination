/****************************************************************
** smallpt
** explain: https://docs.google.com/file/d/0B8g97JkuSSBwUENiWTJXeGtTOHFmSm51UC01YWtCZw/edit
** code : http://kevinbeason.com/smallpt/explicit.cpp
** created time: Jan 31th, 2014.
** note:
**		For viewing the image: http://www.sciweavers.org/free-online-image-converter
*****************************************************************/



#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double M_PI = 3.14159265354;
double M_1_PI = 1.0 / M_PI;

double erand48( unsigned short xsubi[3] )
{
	return (double)rand() / double(RAND_MAX);
}

struct Vec3
{
	double x,y,z;
	Vec3 (double x_ = 0, double y_ = 0, double z_ = 0) { this->x = x_; this->y=y_; this->z = z_; }
	Vec3 operator+(const Vec3 &b) const { return Vec3(b.x+this->x, b.y+this->y, b.z+this->z); }
	Vec3 operator-(const Vec3 &b) const {return Vec3(this->x-b.x, this->y-b.y,this->z-b.z);}
	Vec3 operator*(double b) const {return Vec3(this->x*b,  this->y*b, this->z*b);}
	Vec3 mult( const Vec3 &b ) const {return Vec3(x*b.x, y*b.y, z*b.z );}
	Vec3& norm() { return *this=*this*( 1/sqrt( x*x + y*y + z*z ) ); }
	double dot(const Vec3 &b) const { return x*b.x+y*b.y+z*b.z; }
	Vec3 cross( Vec3 &b ){return Vec3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);} // Cross;
};

struct Ray
{
	Vec3 o, dir;
	Ray( Vec3 o_, Vec3 dir_ ):o(o_),dir(dir_){}
};

enum Refl_t {DIFF, SPEC, REFR}; // Material types, used in radiance

struct Sphere 
{
	double r;		// radius
	Vec3 p, e, c;	// position, emission, color;
	Refl_t refl;	// reflection type (DIFFuse, SPECular, REFRactive)

	Sphere( double r_, Vec3 p_, Vec3 e_, Vec3 c_, Refl_t refl_ ):r(r_), p(p_),e(e_),c(c_), refl(refl_){}

	// Returns distance , 0 if not hit
	double intersect( const Ray &ray ) const
	{
		Vec3 op = p-ray.o; 
		// Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		double t, 
				   eps	=	1e-4, 
				   b	=	op.dot(ray.dir), 
				   det	=	b*b-op.dot(op)+r*r;

		if(det<0) return 0; else det=sqrt(det);
		return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
	}
};

struct Camera
{
	Vec3 pos;	// Position of the camera 
	Vec3 dir;	// Direction of the camera
	Vec3 cx;	// x direction increment ( uses implicit 0 for y, z )
	Vec3 cy;   // y direction increment ( note cross product )

	Camera( Vec3 pos_, Vec3 dir_, Vec3 cx_ ):
		pos(pos_), dir(dir_), cx(cx_)
	{
		this->cy = (cx.cross(this->dir)).norm() * 0.5135;
	}
};

struct Image
{
	Vec3* pixels;
	int width, height;
	Image(int width_=0, int height_=0):pixels(NULL)
	{
		this->width = width_;
		this->height = height_;
		this->pixels = new Vec3[width, height];
	}
};


Sphere spheres[] = 
{//Scene: radius, position, emission, color, material
	//		radius,            position				emission,				color,                   material
	Sphere(1e5, Vec3( 1e5+1,40.8,81.6),		Vec3(),					Vec3(.75,.25,.25),	DIFF),		//Left
	Sphere(1e5, Vec3(-1e5+99,40.8,81.6),		Vec3(),					Vec3(.25,.25,.75),	DIFF),		//Rght
	Sphere(1e5, Vec3(50,40.8, 1e5),			Vec3(),					Vec3(.75,.75,.75),	DIFF),		//Back
	Sphere(1e5, Vec3(50,40.8,-1e5+170),		Vec3(),					Vec3(),				DIFF),		//Frnt
	Sphere(1e5, Vec3(50, 1e5, 81.6),			Vec3(),					Vec3(.75,.75,.75),	DIFF),		//Botm
	Sphere(1e5, Vec3(50,-1e5+81.6,81.6),		Vec3(),					Vec3(.75,.75,.75),	DIFF),		//Top
	Sphere(16.5,Vec3(27,16.5,47),			Vec3(),					Vec3(1,1,1)*.999,	SPEC),		//Mirr
	Sphere(16.5,Vec3(73,16.5,78),			Vec3(),					Vec3(1,1,1)*.999,	REFR),		//Glas
	Sphere(1.5,  Vec3(50,81.6-16.5,81.6),		Vec3(4,4,4)*100,		Vec3(),				DIFF),		//Lite
};

int numSpheres = sizeof(spheres)/sizeof(Sphere);

// CLAMP function
inline double clamp(double x){ return x<0 ? 0 : x>1? 1 : x; }

// The output of the "radiance" function is a set of unbounded colors.
// This has to be converted to between 0 and 255 for display.
// Converts floats to integers to be saved in PPM file
inline int toInt( double x )
{ 
	double gamma = 2.2;
	return int( pow( clamp(x), 1/ gamma ) * 255 + 0.5 ); 
}

// Test if ray intersect with sphere[id]
inline bool intersect( const Ray &ray, double &t, int &id )
{
	double n = sizeof(spheres)/sizeof(Sphere);
	double d = -1.0;
	double inf = t = 1e20;
	for (int i= int(n); i >= 0; i-- )
		if ( (d=spheres[i].intersect(ray)) && d < t )
		{
			t = d;	id = i;
		}

	return t < inf;
}

double tentFilter( double rnd  )
{
	double result = rnd < 1 ? sqrt(rnd)-1 : 1-sqrt(2-rnd);
	return result;
}

// Compute the radiance estimate along a ray
// @param ray:	the ray we are casting
// @param depth:	the ray depth. (termination condition)
// @param Xi:		random number seed
// @param E:		whether to include emissive color
// return value: Vec3 the radiance estimation
Vec3 radiance( const Ray &ray, int depth, unsigned short* Xi, int E=1	)
{
	double t;															// distance to intersection
	int	id = 0;															// id of intersected object
	if ( !intersect( ray, t, id )  ) return Vec3();						// if miss, return black
	const Sphere &obj = spheres[id];									// get the  hit object
	Vec3 x = ray.o + ray.dir*t;										// the hit point
	Vec3 n = ( x - obj.p).norm();										// normal vector of the point in the surface
	Vec3 nl = n.dot(ray.dir) < 0? n:n*(-1);							// properly oriented surface
																		// when a ray hits a glass surface, the ray tracer must 
																		// determine if it is entering or leaving the glass to compute refraction ray.
	Vec3 f = obj.c;													// object color (BRDF modulator)
	
	// Stop the recursion randomly based on the surface reflectivity.
	//		- Use the max component (r,g,b) of the surface color
	//		- Don't do Russian Roulette until after depth 5
	
	// Use maximum reflectivity amount for Russian roulette
	double p = f.x > f.y && f.x > f.z ? f.x : f.y>f.z? f.y: f.z;	// Get the max reflectivity
	if ( ++depth > 5 || !p )
		if ( erand48(Xi) < p )
			f = f * ( 1/p );
		else 
			return obj.e * E;

	if ( obj.refl == DIFF )											// Ideal DIFFUSE reflection
	{
		// Construct random ray:
		//		- Get random angle (r1)
		//		- Get random distance from center (r2s)
		//		- Use normal to create orthonormal coordinate fram (w,u,v)
		double r1 = 2 * M_PI * erand48(Xi);			// angle arbound
		double r2 = erand48(Xi);
		double r2s = sqrt(r2);								// distance from center
		// TODO: 
		// Following four lines of code for generating random reflection ray doesn't seem necessary to me.
		// In a hemisphere coordinate system, two angles and one position can already define the position of a new ray.
		Vec3 w = nl;
		Vec3 u = ( fabs( w.x ) > 0.1 ? Vec3(0,1) : Vec3(1)).cross( w ).norm(); // u is perpendicular to w
		Vec3 v = w.cross(u);														// v is perpendicular to u and w
		Vec3 d = (u*cos(r1)*r2s + v*sin(r1)*r2s +w*sqrt(1-r2)).norm();			// d is random reflection ray
		
		// Loop over any lights
		Vec3 e;
		for ( int i = 0; i < numSpheres; i++ )
		{
			const Sphere& s = spheres[i];
			if ( s.e.x <= 0 && s.e.y <= 0 && s.e.z <= 0 )
				continue;												// Skip non lights
		}
	}
	else if( obj.refl == SPEC )										// Ideal SPECULAR reflection
	{
	}
	else																// Ideal dielectric Refraction
	{

	}




	return Vec3();
}

int WIDTH = 512;
int HEIGHT = 384;

void main( int argc, char * argv[] )
{
	int samps = argc == 2? atoi( argv[1] )/4 : 1; // # Samples (default of 1)

	Camera cam( Vec3(50,52,295.6), Vec3(0, -0.042612, -1 ).norm(), Vec3( WIDTH*0.5135/HEIGHT ) );
	Vec3 r;	// Used for colors of samples.
	Image img(WIDTH, HEIGHT); // Image to display

	#pragma omp parallel for schedule (dynamic, 1) private (r) // OpenMP
	
	// Loop over all image pixels
	for ( int y = 0; y<HEIGHT; y++ )	// Loop rows
	{
		fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4, 100.*y/(HEIGHT-1));
		unsigned short Xi[3] = { 0,0, y*y*y };								
		for ( unsigned short x=0; x<WIDTH; x++ )
			for (int sy = 0; sy<2; sy++)												//|| Each Pixel is devided into 2x2 sub pixels
			{																			//|| Each sub-pixel will be sampled 'samps' times. 
				int i = ( HEIGHT-y-1 )*WIDTH + x; // pixel[i] ;
				for (int sx=0; sx<2; sx++)
				{
					r = Vec3();
					for ( int s=0; s<samps; s++ )
					{
						// Tent filter 
						double dx = tentFilter( 2*erand48(Xi) );		// | dx and dy are random filtered values
						double dy = tentFilter( 2*erand48(Xi) );		// | used to determine location of sample withing pixel

						Vec3 rayDir = cam.cx * ( ( ( sx+0.5+dx )/2 + x ) / WIDTH - 0.5 ) +
											  cam.cy * ( ( ( sy+0.5+dy )/2 + y ) / HEIGHT -0.5 ) + cam.dir;
						// Camera rays are pushed ^^^^ forward to start in interior
						r = r + radiance( Ray(cam.pos + rayDir * 140, rayDir.norm() ), 0, Xi ) * (1.0 / samps);
					}
				}
			}
	}

	return;
}
