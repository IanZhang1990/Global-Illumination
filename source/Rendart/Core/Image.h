

#ifndef CORE_IMAGE
#define CORE_IMAGE

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

namespace Core
{
	/***************************************************************
	Pixel class 
	***************************************************************/
	struct Pixel
	{
		Pixel( float _r = 0.0f, float _g = 00.0f, float _b = 0.0f )
		{
			this->r = _r;
			this->g = _g;
			this->b = _b;
		}

		float r;		// R value 
		float g;		// G value
		float b;		// B value
	};

	/***************************************************************
	Abstract Image class 
	***************************************************************/
	class Image
	{
	public:
		Image( int width = 0, int height = 0 )
			:mPixels( NULL ){};

		// Get pixel at position (x,y)
		virtual Pixel pixel( int x, int y ) const = 0;
		// Set pixel value at position (x,y)
		virtual void setPixel( int x, int y, Pixel val ) = 0;
		// Save the image to a file
		virtual void saveImage( char* filename ) = 0;

	protected:
		Pixel* mPixels;
		int mWidth;
		int mHeight;
	};

	// CLAMP function
	inline double clamp(float x){ return x<0 ? 0 : x>1? 1 : x; }

	// The output of the "radiance" function is a set of unbounded colors.
	// This has to be converted to between 0 and 255 for display.
	// Converts floats to integers to be saved in PPM file
	inline int toInt( float x, float gamma )
	{ 
		return int( pow( clamp(x), 1.0/gamma ) * 255 + 0.5 ); 
	}


	/***************************************************************
	PPM Image class which implements Image abstract class 
	***************************************************************/
	class PPMImage: public Image
	{
		/*
		Constructor of PPM image class.
		The class will write image to .ppm format
		*/
		PPMImage(int width=0, int height=0)
		{
			this->mWidth = width;
			this->mHeight = height;
			this->mPixels = new Pixel[width * height];
		}

		~PPMImage()
		{
			delete this->mPixels;
		}

		/*
		Get pixel at (x, y) position.
		*/
		Pixel pixel( int x, int y )
		{
			int i = ( this->mHeight-y-1 )*this->mWidth + x;
			return this->mPixels[ i ];
		}

		/*
		Set pixel value at (x, y) position.
		*/
		void setPixel( int x, int y, Pixel value )
		{
			int i = ( this->mHeight-y-1 )*this->mWidth + x;
			this->mPixels[i] = value;
		}

		/*
		Save the image to a .ppm file.
		*/
		void saveImage(char* filename)
		{
			FILE *f = fopen( filename , "w");         // Write image to PPM file.
			fprintf(f, "P3\n%d %d\n%d\n", this->mWidth, this->mHeight, 255);

			float gamma = 2.2;
			for (int i=0; i<this->mWidth*this->mHeight; i++)
				fprintf(f,"%d %d %d ", toInt(mPixels[i].r, gamma), toInt(mPixels[i].g, gamma), toInt(mPixels[i].b, gamma));
		}
	};
}


#endif