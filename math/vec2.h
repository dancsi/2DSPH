#pragma once
#include <math.h>
#include <ostream>
#include <memory>

namespace math
{
	class vec2
	{
	public:
		double x,y;
	public:
		vec2(void){}

		vec2(double length)
		{
			int xSign = rand()%10>5 ? 1 : -1;
			int ySign = rand()%10>5 ? 1 : -1;

			x = (rand()%100)*xSign, 
				y = (rand()%100)*ySign;

			Normalise();
			x*=length;
			y*=length;
		}

		vec2(double Ix,double Iy): x(Ix), y(Iy){}

		void Set(double x1, double y1) {x=x1; y=y1;};
		void Rescale(double length) { Normalise(); (*this)*= length; } 
		vec2 &operator /=(const double Scalar)  { x /= Scalar; y /= Scalar;		return *this; }

		vec2 &operator *=(const double Scalar)  { x *= Scalar; y *= Scalar;		return *this; }

		vec2 &operator +=(const vec2 &Other) { x += Other.x;	y += Other.y;	return *this; }

		vec2 &operator -=(const vec2 &Other) { x -= Other.x;	y -= Other.y;	return *this; }

		double operator ^ (const vec2 &V) const {	return (x * V.y) - (y * V.x); } // cross product

		double operator * (const vec2 &V) const {	return (x*V.x) + (y*V.y); } // dot product

		vec2 operator * (double  s)			const	{	return vec2(x*s, y*s); }

		vec2 operator / (double  s)			const	{	return vec2(x/s, y/s); }

		vec2 operator + (const vec2 &V)	const	{	return vec2(x+V.x, y+V.y); }

		vec2 operator - (const vec2 &V)	const	{	return vec2(x-V.x, y-V.y); }

		bool operator ==(const vec2 &V) const { return x==V.x && y==V.y; } 
		bool operator !=(const vec2 &V) const { return x!=V.x || y!=V.y; } 

		friend vec2 operator * (double k, const vec2& V) {	return vec2(V.x*k, V.y*k); }

		friend std::ostream &operator<<( std::ostream &out, const vec2 &v )
		{
			out << "x: " << v.x << " y: " << v.y << std::endl;
			return out;
		}

		vec2 operator -(void) const { return vec2(-x, -y); }

		double Length(void) const { return (double) sqrt(x*x + y*y); }
		double LengthSq() const { return x*x + y*y; }

		double Normalise(void) 
		{	
			double fLength = Length();	

			if (fLength == 0.0f) 
				return 0.0f; 

			(*this) /= fLength;

			return fLength;	
		}

		vec2 Direction(void) const
		{
			vec2 Temp(*this);

			Temp.Normalise();
			return Temp;
		}

		void Perp()
		{
			double temp = x;
			x = -y;
			y = temp;
		}

		void Clear()
		{
			Set(0,0);
		}

		vec2 GetPerp() const
		{
			vec2 a= *this;
			return vec2(-a.y, a.x);
		}

		vec2& rotate(double angle)
		{
			double tx = x;
			x =  x * cos(angle) - y * sin(angle);
			y = tx * sin(angle) + y * cos(angle);
			return *this;
		}

		vec2& transform(const vec2& trans, double rot)
		{
			vec2 D = *this;
			D.rotate(rot);
			*this = trans + D;
			return *this;
		}
		operator std::string()
		{
			char buf[100];
			sprintf(buf, "(%.2lf, %.2lf)", x, y);
			return std::string(buf);
		}
	};
}