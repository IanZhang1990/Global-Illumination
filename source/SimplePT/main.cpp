/****************************************************************
** smallpt
** explain: https://docs.google.com/file/d/0B8g97JkuSSBwUENiWTJXeGtTOHFmSm51UC01YWtCZw/edit
** code : http://kevinbeason.com/smallpt/explicit.cpp
** created time: Jan 31th, 2014.
** note:
**		To view the image: http://www.sciweavers.org/free-online-image-converter
*****************************************************************/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ctime>

#ifndef M_PI
#define M_PI  3.14159265354
#endif

#ifndef M_1_PI
#define M_1_PI (1.0/M_PI)
#endif

#define MAXDEPTH 5

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

enum Refl_t {DIFF, SPEC, REFR};						// Material types, used in radiance

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
	Vec3 pos;			// Position of the camera 
	Vec3 dir;			// Direction of the camera
	Vec3 cx;			// x direction increment ( uses implicit 0 for y, z )
	Vec3 cy;			// y direction increment ( note cross product )

	Camera( Vec3 pos_, Vec3 dir_, Vec3 cx_ ):
		pos(pos_), dir(dir_), cx(cx_)
	{
		this->cy = (cx.cross(this->dir)).norm() * 0.5135;
	}
};

inline int toInt( double x );

struct Image
{
	Vec3* pixels;
	int width, height;
	Image(int width_=0, int height_=0):pixels(NULL)
	{
		this->width = width_;
		this->height = height_;
		this->pixels = new Vec3[width * height];
	}

	Vec3 getPixel( int x, int y )
	{
		int i = ( this->height-y-1 )*this->width + x;
		return this->pixels[ i ];
	}

	void setPixel( int x, int y, Vec3 value )
	{
		int i = ( this->height-y-1 )*this->width + x;
		this->pixels[i] = value;
	}

	void saveImg(char* filename)
	{
		FILE *f = fopen( filename , "w");         // Write image to PPM file.
		fprintf(f, "P3\n%d %d\n%d\n", this->width, this->height, 255);

		for (int i=0; i<this->width*this->height; i++)
			fprintf(f,"%d %d %d ", toInt(pixels[i].x), toInt(pixels[i].y), toInt(pixels[i].z));

	}
};


Sphere spheres[] = 
{//Scene: radius, position, emission, color, material
	//	 radius,        position			   emission,				color,                 material
	Sphere(1e5, Vec3( 1e5+1,40.8,81.6),			Vec3(),					Vec3(.75,.25,.25),		DIFF),		//Left
	Sphere(1e5, Vec3(-1e5+99,40.8,81.6),		Vec3(),					Vec3(.25,.25,.75),		DIFF),		//Rght
	Sphere(1e5, Vec3(50,40.8, 1e5),				Vec3(),					Vec3(.75,.75,.75),		DIFF),		//Back
	Sphere(1e5, Vec3(50,40.8,-1e5+170),			Vec3(),					Vec3(),					DIFF),		//Frnt
	Sphere(1e5, Vec3(50, 1e5, 81.6),			Vec3(),					Vec3(.75,.75,.75),		DIFF),		//Botm
	Sphere(1e5, Vec3(50,-1e5+81.6,81.6),		Vec3(),					Vec3(.75,.75,.75),		DIFF),		//Top
	Sphere(16.5,Vec3(27,16.5,47),				Vec3(),					Vec3(1,1,1)*.999,		SPEC),		//Mirr
	Sphere(16.5,Vec3(73,25.0,78),				Vec3(),					Vec3(1,1,1)*.999,		REFR),		//Glas
	Sphere(1.5,  Vec3(50,81.6-16.5,81.6),		Vec3(4,4,4)*100,		Vec3(),					DIFF),		//Lite
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
	// TODO: does this exceed the range of sphere? what is n? 9? 8?
	// Test result is n=9;  so you can't use "int i = int(n)""
	//for (int i= int(n); i >= 0; i-- )
	for (int i= 0; i < int(n); i++ )
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
	Vec3 x = ray.o + ray.dir*t;											// the hit point
	Vec3 n = ( x - obj.p).norm();										// normal vector of the point in the surface
	Vec3 nl = n.dot(ray.dir) < 0? n:n*(-1);								// properly oriented surface
																		// when a ray hits a glass surface, the ray tracer must 
																		// determine if it is entering or leaving the glass to compute refraction ray.
	Vec3 f = obj.c;														// object color (BRDF modulator)
	depth++;	

	// Stop the recursion randomly based on the surface reflectivity.
	//		- Use the max component (r,g,b) of the surface color
	//		- Don't do Russian Roulette until after depth 5
	
	// Use maximum reflectivity amount for Russian roulette
	double p = f.x > f.y && f.x > f.z ? f.x : f.y>f.z? f.y: f.z;		// Get the max reflectivity
	if ( depth > MAXDEPTH || !p )
		if ( erand48(Xi) < p )
			f = f * ( 1/p );
		else 
			return obj.e * E;

	if ( obj.refl == DIFF )												// Ideal DIFFUSE reflection
	{
		// Construct random ray:
		//		- Get random angle (r1)
		//		- Get random distance from center (r2s)
		//		- Use normal to create orthonormal coordinate fram (w,u,v)
		double r1 = 2 * M_PI * erand48(Xi);					// angle arbound
		double r2 = erand48(Xi);
		double r2s = sqrt(r2);								// distance from center
		// TODO: 
		// Following four lines of code for generating random reflection ray doesn't seem necessary to me.
		// In a hemisphere coordinate system, two angles and one position can already define the position of a new ray.
		Vec3 w = nl;
		Vec3 u = ( fabs( w.x ) > 0.1 ? Vec3(0,1) : Vec3(1)).cross( w ).norm();		// u is perpendicular to w
		Vec3 v = w.cross(u);														// v is perpendicular to u and w
		Vec3 d = (u*cos(r1)*r2s + v*sin(r1)*r2s +w*sqrt(1-r2)).norm();				// d is random reflection ray
		
		// Loop over any lights. TODO: Create a light class to handle lights
		Vec3 e;
		for ( int i = 0; i < numSpheres; i++ )
		{
			const Sphere& s = spheres[i];
			if ( s.e.x <= 0 && s.e.y <= 0 && s.e.z <= 0 )
				continue;															// Skip non lights
			
			// Create random direction towards sphere using method from <<Realistic Ray Tracing>>.
			Vec3 sw=s.p-x, su=((fabs(sw.x)>0.1?Vec3(0,1):Vec3(1)).cross(sw)).norm(), sv=sw.cross(su);
			double cos_a_max = sqrt(1-s.r*s.r/(x-s.p).dot(x-s.p));
			double eps1 = erand48(Xi), eps2 = erand48(Xi);
			double cos_a = 1-eps1+eps1*cos_a_max;
			double sin_a = sqrt(1-cos_a*cos_a);
			double phi = 2*M_PI*eps2;
			Vec3 l = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
			l.norm();

			// Shoot shadow ray
			if ( intersect(Ray( x, l ), t, id) && id == i )						// Check for occlusion with shadow ray
			{
				double omega = 2 * M_PI * ( 1 - cos_a_max );					// Compute 1/probability with respect to solid angle
				e = e + f.mult( s.e * l.dot( nl ) * omega ) * M_1_PI;			// 1/pi for brdf. Calculate lighting and add to current value;	
			}

			// Make recursive call with random ray direction computed earlier
			return obj.e*E + e + f.mult( radiance( Ray(x,d), depth, Xi, 0 ) );	// The 0 parameter at the end turns off the emissive term at the next recursion level
		}
	}
	else if( obj.refl == SPEC )													// Ideal SPECULAR reflection
	{
		// Mirror reflection. 
		// TODO: Mirror reflection should not return 100% of lights, instead, energy should decrease
		return obj.e + f.mult( radiance( Ray( x, ray.dir-n*2*n.dot( ray.dir ) ), depth, Xi ) ) * 0.95;
	}
	else																		// Ideal dielectric Refraction
	{
		// We have a dielectric (glass) surface
		// Glass is both reflective and refractive
		Ray reflRay( x, ray.dir - n*2*n.dot(ray.dir) );							// Ideal dielectric reflection ray
		bool into = n.dot(nl) > 0;												// Is the ray from outside and going into the glass?
		double nc=1, nt = 1.5, nnt = into? nc/nt : nt/nc;						// IOR for glass is 1.5. nnt is either 1.5 or 1/1.5
		double ddn = ray.dir.dot(nl), cos2t;

		// If total internal reflection, reflect. (Fresnel effect)
		if ( ( cos2t = 1-nnt*nnt*(1-ddn*ddn) ) < 0 )							// total internal reflection
			return obj.e + f.mult( radiance( reflRay, depth, Xi ) );
		
		// Otherwise, choose reflection or refraction
		Vec3 tdir = (ray.dir*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
		double a=nt-nc, b=nt+nc;
		double R0=a*a/(b*b);													// reflectance at normal incidence based on IOR
		double c = 1-(into?-ddn:tdir.dot(n));									// 1- cos(theta)
		double Re=R0+(1-R0)*c*c*c*c*c;											// Fresnel reflectance
		double Tr=1-Re;					
		double P=.25+.5*Re;														// Probability of reflecting
		double RP=Re/P;												
		double TP=Tr/(1-P); 

		// Finally make 1 or 2 recursive calls
		//		- Make 2 if depth is <= 2;
		//		- Make 1 randomly if depth > 2
		return obj.e + f.mult( depth > 2 ?
			(erand48(Xi) < P ?  radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
			radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr );

	}
}

int WIDTH = 1366;
int HEIGHT = 950;


void showCurrentTime()
{
	// current date/time based on current system
   time_t now = time(0);
   
   // convert now to string form
   char* dt = ctime(&now);

   std::cout << "\nThe local date and time is: " << dt << std::endl;
}


int main( int argc, char * argv[] )
{
	int samps = argc == 2? atoi( argv[1] )/4 : 1; // # Samples (default of 1)

	Camera cam( Vec3(50,52,295.6), Vec3(0, -0.042612, -1 ).norm(), Vec3( WIDTH*0.5135/HEIGHT ) );
	Vec3 r;	// Used for colors of samples.
	Image img(WIDTH, HEIGHT); // Image to display

	showCurrentTime();

	// Loop over all image pixels
	#pragma omp parallel for schedule (dynamic, 1) private (r) // OpenMP
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
						double dx = tentFilter( 2*erand48(Xi) );						// | dx and dy are random filtered values
						double dy = tentFilter( 2*erand48(Xi) );						// | used to determine location of sample withing pixel

						Vec3 rayDir = cam.cx * ( ( ( sx+0.5+dx )/2 + x ) / WIDTH - 0.5 ) +
										 cam.cy * ( ( ( sy+0.5+dy )/2 + y ) / HEIGHT - 0.5 ) + cam.dir;
						// Camera rays are pushed ^^^^ forward to start in interior
						r = r + radiance( Ray(cam.pos + rayDir * 140, rayDir.norm() ), 0, Xi ) * (1.0 / samps);
					}
					Vec3 oldVal = img.getPixel( x, y );
					double gamma = 0.25; 
					Vec3 newVal = oldVal + Vec3( clamp( r.x ), clamp(r.y), clamp(r.z) ) * gamma;
					img.setPixel(x, y, newVal);
				}
			}
	}

	showCurrentTime();

	img.saveImg( "image.ppm" );

	return 0;
}
