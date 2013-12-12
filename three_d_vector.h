#ifndef THREEDVECTOR_H
#define THREEDVECTOR_H
using namespace std;

class ThreeDVector{
	public:
		float x;
		float y;
		float z;
		ThreeDVector(float=0, float=0, float=0);
		float magnitude();
		ThreeDVector* normalize();
		void normalize_bang();
		float dot_product(ThreeDVector*);
		ThreeDVector* scalar_multiply(float k);
		void scalar_multiply_bang(float k);
		ThreeDVector* vector_add(ThreeDVector*);
		void vector_add_bang(ThreeDVector*);
		ThreeDVector* vector_subtract(ThreeDVector*);
		ThreeDVector* vector_multiply(ThreeDVector*);
		ThreeDVector* cross_product(ThreeDVector*);
		float distance(ThreeDVector*);
		float distance(float x, float y, float z);
		ThreeDVector* midpoint(ThreeDVector*);
		ThreeDVector* clone();
		char* repr();
		char* print();
		bool operator==(const ThreeDVector &other) const;
};

class comparator {
    public:
        bool operator()(const ThreeDVector* lhs, const ThreeDVector* rhs) const {
        	const ThreeDVector* left_v = (lhs);
        	const ThreeDVector* right_v = (rhs);
            return left_v->x < right_v->x
                || ( left_v->x == right_v->x && ( left_v->y < right_v->y
                || ( left_v->y == right_v->y && left_v->z < right_v->z)));
        }
};

#endif