#include <cmath>
#include "three_d_vector.h"
#include <iostream>
#include <stdio.h> 
#include <cfloat>
using namespace std;

ThreeDVector::ThreeDVector(long double _x, long double _y, long double _z) {
	x = _x;
	y = _y;
	z = _z;
}

long double ThreeDVector::magnitude(){
	long double x = this->x;
	long double y = this->y;
	long double z = this->z;
	return sqrt(x*x+y*y+z*z);
}

ThreeDVector* ThreeDVector::normalize(){
	//cout << "BEFORE: " << this->x << ", " << this->y << ", " << this->z << endl;
	long double mag = magnitude();
	if (mag == 0) {
		return new ThreeDVector();
	} else {
		return new ThreeDVector(this->x / mag, this->y / mag, this->z / mag);
	}
	//cout << "AFTER: " << this->x << ", " << this->y << ", " << this->z << endl;
}

void ThreeDVector::normalize_bang(){
	//cout << "BEFORE: " << this->x << ", " << this->y << ", " << this->z << endl;
	long double mag = magnitude();
	this->x /= mag;
	this->y /= mag;
	this->z /= mag;
	//cout << "AFTER: " << this->x << ", " << this->y << ", " << this->z << endl;
}



long double ThreeDVector::dot_product(ThreeDVector* v){
	return this->x * v->x + this->y * v->y + this->z * v->z;
}


ThreeDVector* ThreeDVector::scalar_multiply(long double k){	
	return new ThreeDVector(this->x * k, this->y * k, this->z * k);
}

void ThreeDVector::scalar_multiply_bang(long double k){	
	this->x *= k;
	this->y *= k;
	this->z *= k;
}

ThreeDVector* ThreeDVector::vector_add(ThreeDVector* v){
	return new ThreeDVector(this->x + v->x, this->y + v->y, this->z + v->z);
}

void ThreeDVector::vector_add_bang(ThreeDVector* v){
	this->x += v->x;
	this->y += v->y;
	this->z += v->z;
}

ThreeDVector* ThreeDVector::vector_subtract(ThreeDVector* v){
	return new ThreeDVector(this->x - v->x, this->y - v->y, this->z - v->z);
}

ThreeDVector* ThreeDVector::vector_multiply(ThreeDVector* v){
	return new ThreeDVector(v->x * this->x, v->y * this->y, v->z * this->z);
}

ThreeDVector* ThreeDVector::cross_product(ThreeDVector* v){
	return new ThreeDVector(this->y * v->z - this->z * v->y, this->z * v->x - this->x * v->z, this->x * v->y - this->y * v->x);                                                                         
}

long double ThreeDVector::distance(ThreeDVector* v) {
	long double delta_x = this->x - v->x;
	long double delta_y = this->y - v->y;
	long double delta_z = this->z - v->z;
	return sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
}

ThreeDVector* ThreeDVector::midpoint(ThreeDVector* v) {
	return new ThreeDVector((this->x + v->x)/2, (this->y + v->y)/2, (this->z + v->z)/2);
}

ThreeDVector* ThreeDVector::clone(){
	return new ThreeDVector(this->x, this->y, this->z);
}

void ThreeDVector::transform_bang(Eigen::Matrix4f transformation, bool point) {
	Eigen::Vector4f old;
	if (point) {
		old = Eigen::Vector4f(this->x, this->y, this->z, 1);
	} else {
		old = Eigen::Vector4f(this->x, this->y, this->z, 0);
	}
	Eigen::Vector4f transformed = transformation * old;

	this->x = transformed[0];
	this->y = transformed[1];
	this->z = transformed[2];
}

char* ThreeDVector::repr() {
	char* buffer = new char[1000];
	sprintf(buffer, "<ThreeDVector, x = %0.2Lf, y = %0.2Lf, z = %0.2Lf>", this->x, this->y, this->z);
	return buffer;
}

char* ThreeDVector::print() {
	char* buffer = new char[1000];
	sprintf(buffer, "%0.2Lf %0.2Lf %0.2Lf", this->x, this->y, this->z);
	return buffer;
}

/*
//Override Hash and Equality Operations
bool ThreeDVector::operator==(const ThreeDVector &other) const {
	return this->x == other.x && this->y == other.y & this->z == other.z;    
}

namespace std
{
    template<>
    struct hash<ThreeDVector>
    {
    public:
        std::size_t operator()(ThreeDVector const& v) const 
        {
            return ((hash<long double>()(v.x)
               ^ (hash<long double>()(v.y) << 1)) >> 1)
               ^ (hash<long double>()(v.z) << 1);
        }
    };
}*/
