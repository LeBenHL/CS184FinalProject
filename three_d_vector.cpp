#include <cmath>
#include "three_d_vector.h"
#include <iostream>
#include <stdio.h> 
#include <cfloat>
using namespace std;

ThreeDVector::ThreeDVector(float _x, float _y, float _z) {
	x = _x;
	y = _y;
	z = _z;
}

float ThreeDVector::magnitude(){
	float x = this->x;
	float y = this->y;
	float z = this->z;
	return sqrt(x*x+y*y+z*z);
}

ThreeDVector* ThreeDVector::normalize(){
	//cout << "BEFORE: " << this->x << ", " << this->y << ", " << this->z << endl;
	float mag = magnitude();
	if (mag == 0) {
		return new ThreeDVector();
	} else {
		return new ThreeDVector(this->x / mag, this->y / mag, this->z / mag);
	}
	//cout << "AFTER: " << this->x << ", " << this->y << ", " << this->z << endl;
}

void ThreeDVector::normalize_bang(){
	//cout << "BEFORE: " << this->x << ", " << this->y << ", " << this->z << endl;
	float mag = magnitude();
	if (mag != 0) {
		this->x /= mag;
		this->y /= mag;
		this->z /= mag;
	}
	//cout << "AFTER: " << this->x << ", " << this->y << ", " << this->z << endl;
}



float ThreeDVector::dot_product(ThreeDVector* v){
	return this->x * v->x + this->y * v->y + this->z * v->z;
}


ThreeDVector* ThreeDVector::scalar_multiply(float k){	
	return new ThreeDVector(this->x * k, this->y * k, this->z * k);
}

void ThreeDVector::scalar_multiply_bang(float k){	
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

float ThreeDVector::distance(ThreeDVector* v) {
	float delta_x = this->x - v->x;
	float delta_y = this->y - v->y;
	float delta_z = this->z - v->z;
	return sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
}

ThreeDVector* ThreeDVector::midpoint(ThreeDVector* v) {
	return new ThreeDVector((this->x + v->x)/2, (this->y + v->y)/2, (this->z + v->z)/2);
}

ThreeDVector* ThreeDVector::clone(){
	return new ThreeDVector(this->x, this->y, this->z);
}

char* ThreeDVector::repr() {
	char* buffer = new char[1000];
	sprintf(buffer, "<ThreeDVector, x = %0.2f, y = %0.2f, z = %0.2f>", this->x, this->y, this->z);
	return buffer;
}

char* ThreeDVector::print() {
	char* buffer = new char[1000];
	sprintf(buffer, "%0.2f %0.2f %0.2f", this->x, this->y, this->z);
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
            return ((hash<float>()(v.x)
               ^ (hash<float>()(v.y) << 1)) >> 1)
               ^ (hash<float>()(v.z) << 1);
        }
    };
}*/
