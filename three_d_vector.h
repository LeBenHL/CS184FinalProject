#ifndef THREEDVECTOR_H
#define THREEDVECTOR_H
#include <Eigen/Dense>

class ThreeDVector{
	public:
		long double x;
		long double y;
		long double z;
		ThreeDVector(long double=0, long double=0, long double=0);
		long double magnitude();
		ThreeDVector* normalize();
		void normalize_bang();
		long double dot_product(ThreeDVector*);
		ThreeDVector* scalar_multiply(long double k);
		void scalar_multiply_bang(long double k);
		ThreeDVector* vector_add(ThreeDVector*);
		void vector_add_bang(ThreeDVector*);
		ThreeDVector* vector_subtract(ThreeDVector*);
		ThreeDVector* vector_multiply(ThreeDVector*);
		ThreeDVector* cross_product(ThreeDVector*);
		long double distance(ThreeDVector*);
		ThreeDVector* midpoint(ThreeDVector*);
		ThreeDVector* clone();
		void transform_bang(Eigen::Matrix4f transformation, bool point);
		char* repr();
		char* print();
};

#endif