/*************************************************************************/
/*  convex_hull.cpp                                                      */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2021 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2021 Godot Engine contributors (cf. AUTHORS.md).   */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

/*
Copyright (c) 2011 Ole Kniemeyer, MAXON, www.maxon.net

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#include <string.h>

#include "aabb.h"
#include "convex_hull.h"

#include "core/error_macros.h"
#include "core/math/math_defs.h"
#include "core/os/memory.h"

//#define DEBUG_CONVEX_HULL
//#define SHOW_ITERATIONS

// -- GODOT start --
// Assembly optimizations are not used at the moment.
//#define USE_X86_64_ASM
// -- GODOT end --

#ifdef DEBUG_ENABLED
#define CHULL_ASSERT(m_cond) do { if (unlikely(!(m_cond))) { ERR_PRINT("Assertion \"" _STR(m_cond) "\" failed.") } } while(0)
#else
#define CHULL_ASSERT(m_cond) do { } while(0)
#endif

#if defined(DEBUG_CONVEX_HULL) || defined(SHOW_ITERATIONS)
#include <stdio.h>
#endif

// Convex hull implementation based on Preparata and Hong
// Ole Kniemeyer, MAXON Computer GmbH
class ConvexHullInternal {
public:
	class Point64 {
	public:
		int64_t x;
		int64_t y;
		int64_t z;

		Point64(int64_t p_x, int64_t p_y, int64_t p_z) {
			x = p_x;
			y = p_y;
			z = p_z;
		}

		bool is_zero() {
			return (x == 0) && (y == 0) && (z == 0);
		}

		int64_t dot(const Point64 &b) const {
			return x * b.x + y * b.y + z * b.z;
		}
	};

	class Point32 {
	public:
		int32_t x = 0;
		int32_t y = 0;
		int32_t z = 0;
		int32_t index = -1;

		Point32() {
		}

		Point32(int32_t p_x, int32_t p_y, int32_t p_z) {
			x = p_x;
			y = p_y;
			z = p_z;
		}

		bool operator==(const Point32 &b) const {
			return (x == b.x) && (y == b.y) && (z == b.z);
		}

		bool operator!=(const Point32 &b) const {
			return (x != b.x) || (y != b.y) || (z != b.z);
		}

		bool is_zero() {
			return (x == 0) && (y == 0) && (z == 0);
		}

		Point64 cross(const Point32 &b) const {
			return Point64(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
		}

		Point64 cross(const Point64 &b) const {
			return Point64(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
		}

		int64_t dot(const Point32 &b) const {
			return x * b.x + y * b.y + z * b.z;
		}

		int64_t dot(const Point64 &b) const {
			return x * b.x + y * b.y + z * b.z;
		}

		Point32 operator+(const Point32 &b) const {
			return Point32(x + b.x, y + b.y, z + b.z);
		}

		Point32 operator-(const Point32 &b) const {
			return Point32(x - b.x, y - b.y, z - b.z);
		}
	};

	class Int128 {
	public:
		uint64_t low = 0;
		uint64_t high = 0;

		Int128() {
		}

		Int128(uint64_t p_low, uint64_t p_high) {
			low = p_low;
			high = p_high;
		}

		Int128(uint64_t p_low) {
			low = p_low;
			high = 0;
		}

		Int128(int64_t p_value) {
			low = p_value;
			if (p_value >= 0) {
				high = 0;
			} else {
				high = (uint64_t)-1LL;
			}
		}

		static Int128 mul(int64_t a, int64_t b);

		static Int128 mul(uint64_t a, uint64_t b);

		Int128 operator-() const {
			return Int128((uint64_t) - (int64_t)low, ~high + (low == 0));
		}

		Int128 operator+(const Int128 &b) const {
#ifdef USE_X86_64_ASM
			Int128 result;
			__asm__("addq %[bl], %[rl]\n\t"
					"adcq %[bh], %[rh]\n\t"
					: [rl] "=r"(result.low), [rh] "=r"(result.high)
					: "0"(low), "1"(high), [bl] "g"(b.low), [bh] "g"(b.high)
					: "cc");
			return result;
#else
			uint64_t lo = low + b.low;
			return Int128(lo, high + b.high + (lo < low));
#endif
		}

		Int128 operator-(const Int128 &b) const {
#ifdef USE_X86_64_ASM
			Int128 result;
			__asm__("subq %[bl], %[rl]\n\t"
					"sbbq %[bh], %[rh]\n\t"
					: [rl] "=r"(result.low), [rh] "=r"(result.high)
					: "0"(low), "1"(high), [bl] "g"(b.low), [bh] "g"(b.high)
					: "cc");
			return result;
#else
			return *this + -b;
#endif
		}

		Int128 &operator+=(const Int128 &b) {
#ifdef USE_X86_64_ASM
			__asm__("addq %[bl], %[rl]\n\t"
					"adcq %[bh], %[rh]\n\t"
					: [rl] "=r"(low), [rh] "=r"(high)
					: "0"(low), "1"(high), [bl] "g"(b.low), [bh] "g"(b.high)
					: "cc");
#else
			uint64_t lo = low + b.low;
			if (lo < low) {
				++high;
			}
			low = lo;
			high += b.high;
#endif
			return *this;
		}

		Int128 &operator++() {
			if (++low == 0) {
				++high;
			}
			return *this;
		}

		Int128 operator*(int64_t b) const;

		real_t to_scalar() const {
			return ((int64_t)high >= 0) ? real_t(high) * (real_t(0x100000000LL) * real_t(0x100000000LL)) + real_t(low) : -(-*this).to_scalar();
		}

		int32_t get_sign() const {
			return ((int64_t)high < 0) ? -1 : (high || low) ? 1 :
																0;
		}

		bool operator<(const Int128 &b) const {
			return (high < b.high) || ((high == b.high) && (low < b.low));
		}

		int32_t ucmp(const Int128 &b) const {
			if (high < b.high) {
				return -1;
			}
			if (high > b.high) {
				return 1;
			}
			if (low < b.low) {
				return -1;
			}
			if (low > b.low) {
				return 1;
			}
			return 0;
		}
	};

	class Rational64 {
	private:
		uint64_t numerator;
		uint64_t denominator;
		int32_t sign;

	public:
		Rational64(int64_t p_numerator, int64_t p_denominator) {
			if (p_numerator > 0) {
				sign = 1;
				numerator = (uint64_t)p_numerator;
			} else if (p_numerator < 0) {
				sign = -1;
				numerator = (uint64_t)-p_numerator;
			} else {
				sign = 0;
				numerator = 0;
			}
			if (p_denominator > 0) {
				denominator = (uint64_t)p_denominator;
			} else if (p_denominator < 0) {
				sign = -sign;
				denominator = (uint64_t)-p_denominator;
			} else {
				denominator = 0;
			}
		}

		bool is_negative_infinity() const {
			return (sign < 0) && (denominator == 0);
		}

		bool is_nan() const {
			return (sign == 0) && (denominator == 0);
		}

		int32_t compare(const Rational64 &b) const;

		real_t to_scalar() const {
			return sign * ((denominator == 0) ? FLT_MAX : (real_t)numerator / denominator);
		}
	};

	class Rational128 {
	private:
		Int128 numerator;
		Int128 denominator;
		int32_t sign;
		bool is_int_64;

	public:
		Rational128(int64_t p_value) {
			if (p_value > 0) {
				sign = 1;
				this->numerator = p_value;
			} else if (p_value < 0) {
				sign = -1;
				this->numerator = -p_value;
			} else {
				sign = 0;
				this->numerator = (uint64_t)0;
			}
			this->denominator = (uint64_t)1;
			is_int_64 = true;
		}

		Rational128(const Int128 &p_numerator, const Int128 &p_denominator) {
			sign = p_numerator.get_sign();
			if (sign >= 0) {
				this->numerator = p_numerator;
			} else {
				this->numerator = -p_numerator;
			}
			int32_t dsign = p_denominator.get_sign();
			if (dsign >= 0) {
				this->denominator = p_denominator;
			} else {
				sign = -sign;
				this->denominator = -p_denominator;
			}
			is_int_64 = false;
		}

		int32_t compare(const Rational128 &b) const;

		int32_t compare(int64_t b) const;

		real_t to_scalar() const {
			return sign * ((denominator.get_sign() == 0) ? FLT_MAX : numerator.to_scalar() / denominator.to_scalar());
		}
	};

	class PointR128 {
	public:
		Int128 x;
		Int128 y;
		Int128 z;
		Int128 denominator;

		PointR128() {
		}

		PointR128(Int128 p_x, Int128 p_y, Int128 p_z, Int128 p_denominator) {
			x = p_x;
			y = p_y;
			z = p_z;
			denominator = p_denominator;
		}

		real_t xvalue() const {
			return x.to_scalar() / denominator.to_scalar();
		}

		real_t yvalue() const {
			return y.to_scalar() / denominator.to_scalar();
		}

		real_t zvalue() const {
			return z.to_scalar() / denominator.to_scalar();
		}
	};

	class Edge;
	class Face;

	class Vertex {
	public:
		Vertex *next = nullptr;
		Vertex *prev = nullptr;
		Edge *edges = nullptr;
		Face *first_nearby_face = nullptr;
		Face *last_nearby_face = nullptr;
		PointR128 point128;
		Point32 point;
		int32_t copy = -1;

		Vertex() {
		}

#ifdef DEBUG_CONVEX_HULL
		void print() {
			printf("V%d (%d, %d, %d)", point.index, point.x, point.y, point.z);
		}

		void print_graph();
#endif

		Point32 operator-(const Vertex &b) const {
			return point - b.point;
		}

		Rational128 dot(const Point64 &b) const {
			return (point.index >= 0) ? Rational128(point.dot(b)) : Rational128(point128.x * b.x + point128.y * b.y + point128.z * b.z, point128.denominator);
		}

		real_t xvalue() const {
			return (point.index >= 0) ? real_t(point.x) : point128.xvalue();
		}

		real_t yvalue() const {
			return (point.index >= 0) ? real_t(point.y) : point128.yvalue();
		}

		real_t zvalue() const {
			return (point.index >= 0) ? real_t(point.z) : point128.zvalue();
		}

		void receive_nearby_faces(Vertex *p_src) {
			if (last_nearby_face) {
				last_nearby_face->next_with_same_nearby_vertex = p_src->first_nearby_face;
			} else {
				first_nearby_face = p_src->first_nearby_face;
			}
			if (p_src->last_nearby_face) {
				last_nearby_face = p_src->last_nearby_face;
			}
			for (Face *f = p_src->first_nearby_face; f; f = f->next_with_same_nearby_vertex) {
				CHULL_ASSERT(f->nearby_vertex == p_src);
				f->nearby_vertex = this;
			}
			p_src->first_nearby_face = nullptr;
			p_src->last_nearby_face = nullptr;
		}
	};

	class Edge {
	public:
		Edge *next = nullptr;
		Edge *prev = nullptr;
		Edge *reverse = nullptr;
		Vertex *target = nullptr;
		Face *face = nullptr;
		int32_t copy = -1;

		~Edge() {
			next = nullptr;
			prev = nullptr;
			reverse = nullptr;
			target = nullptr;
			face = nullptr;
		}

		void link(Edge *n) {
			CHULL_ASSERT(reverse->target == n->reverse->target);
			next = n;
			n->prev = this;
		}

#ifdef DEBUG_CONVEX_HULL
		void print() {
			printf("E%p : %d -> %d,  n=%p p=%p   (0 %d\t%d\t%d) -> (%d %d %d)", this, reverse->target->point.index, target->point.index, next, prev,
					reverse->target->point.x, reverse->target->point.y, reverse->target->point.z, target->point.x, target->point.y, target->point.z);
		}
#endif
	};

	class Face {
	public:
		Face *next = nullptr;
		Vertex *nearby_vertex = nullptr;
		Face *next_with_same_nearby_vertex = nullptr;
		Point32 origin;
		Point32 dir0;
		Point32 dir1;

		Face() {
		}

		void init(Vertex *p_a, Vertex *p_b, Vertex *p_c) {
			nearby_vertex = p_a;
			origin = p_a->point;
			dir0 = *p_b - *p_a;
			dir1 = *p_c - *p_a;
			if (p_a->last_nearby_face) {
				p_a->last_nearby_face->next_with_same_nearby_vertex = this;
			} else {
				p_a->first_nearby_face = this;
			}
			p_a->last_nearby_face = this;
		}

		Point64 getNormal() {
			return dir0.cross(dir1);
		}
	};

	template <typename UWord, typename UHWord>
	class DMul {
	private:
		static uint32_t high(uint64_t p_value) {
			return (uint32_t)(p_value >> 32);
		}

		static uint32_t low(uint64_t p_value) {
			return (uint32_t)p_value;
		}

		static uint64_t mul(uint32_t a, uint32_t b) {
			return (uint64_t)a * (uint64_t)b;
		}

		static void shl_half(uint64_t &p_value) {
			p_value <<= 32;
		}

		static uint64_t high(Int128 p_value) {
			return p_value.high;
		}

		static uint64_t low(Int128 p_value) {
			return p_value.low;
		}

		static Int128 mul(uint64_t a, uint64_t b) {
			return Int128::mul(a, b);
		}

		static void shl_half(Int128 &p_value) {
			p_value.high = p_value.low;
			p_value.low = 0;
		}

	public:
		static void mul(UWord p_a, UWord p_b, UWord &r_low, UWord &r_high) {
			UWord p00 = mul(low(p_a), low(p_b));
			UWord p01 = mul(low(p_a), high(p_b));
			UWord p10 = mul(high(p_a), low(p_b));
			UWord p11 = mul(high(p_a), high(p_b));
			UWord p0110 = UWord(low(p01)) + UWord(low(p10));
			p11 += high(p01);
			p11 += high(p10);
			p11 += high(p0110);
			shl_half(p0110);
			p00 += p0110;
			if (p00 < p0110) {
				++p11;
			}
			r_low = p00;
			r_high = p11;
		}
	};

private:
	class IntermediateHull {
	public:
		Vertex *min_xy = nullptr;
		Vertex *max_xy = nullptr;
		Vertex *min_yx = nullptr;
		Vertex *max_yx = nullptr;

		IntermediateHull() {
		}

		void print();
	};

	enum Orientation { NONE,
		CLOCKWISE,
		COUNTER_CLOCKWISE };

	template <typename T>
	class PoolArray {
	private:
		T *array;
		int32_t size;

	public:
		PoolArray<T> *next;

		PoolArray(int32_t p_size) {
			size = p_size;
			next = nullptr;
			array = (T *)memalloc(sizeof(T) * p_size);
		}

		~PoolArray() {
			memfree(array);
		}

		T *init() {
			T *o = array;
			for (int32_t i = 0; i < size; i++, o++) {
				o->next = (i + 1 < size) ? o + 1 : nullptr;
			}
			return array;
		}
	};

	template <typename T>
	class Pool {
	private:
		PoolArray<T> *arrays = nullptr;
		PoolArray<T> *next_array = nullptr;
		T *free_objects = nullptr;
		int32_t array_size = 256;

	public:
		Pool() {
		}

		~Pool() {
			while (arrays) {
				PoolArray<T> *p = arrays;
				arrays = p->next;
				memdelete(p);
			}
		}

		void reset() {
			next_array = arrays;
			free_objects = nullptr;
		}

		void set_array_size(int32_t p_array_size) {
			this->array_size = p_array_size;
		}

		T *new_object() {
			T *o = free_objects;
			if (!o) {
				PoolArray<T> *p = next_array;
				if (p) {
					next_array = p->next;
				} else {
					p = memnew(PoolArray<T>(array_size));
					p->next = arrays;
					arrays = p;
				}
				o = p->init();
			}
			free_objects = o->next;
			return memnew_placement(o, T);
		};

		void free_object(T *p_object) {
			p_object->~T();
			p_object->next = free_objects;
			free_objects = p_object;
		}
	};

	Vector3 scaling;
	Vector3 center;
	Pool<Vertex> vertex_pool;
	Pool<Edge> edge_pool;
	Pool<Face> face_pool;
	LocalVector<Vertex *> original_vertices;
	int32_t merge_stamp = 0;
	int32_t min_axis = 0;
	int32_t med_axis = 0;
	int32_t max_axis = 0;
	int32_t used_edge_pairs = 0;
	int32_t max_used_edge_pairs = 0;

	static Orientation get_orientation(const Edge *p_prev, const Edge *p_next, const Point32 &p_s, const Point32 &p_t);
	Edge *find_max_angle(bool p_ccw, const Vertex *p_start, const Point32 &p_s, const Point64 &p_rxs, const Point64 &p_ssxrxs, Rational64 &p_min_cot);
	void find_edge_for_coplanar_faces(Vertex *p_c0, Vertex *p_c1, Edge *&p_e0, Edge *&p_e1, Vertex *p_stop0, Vertex *p_stop1);

	Edge *new_edge_pair(Vertex *p_from, Vertex *p_to);

	void remove_edge_pair(Edge *p_edge) {
		Edge *n = p_edge->next;
		Edge *r = p_edge->reverse;

		CHULL_ASSERT(p_edge->target && r->target);

		if (n != p_edge) {
			n->prev = p_edge->prev;
			p_edge->prev->next = n;
			r->target->edges = n;
		} else {
			r->target->edges = nullptr;
		}

		n = r->next;

		if (n != r) {
			n->prev = r->prev;
			r->prev->next = n;
			p_edge->target->edges = n;
		} else {
			p_edge->target->edges = nullptr;
		}

		edge_pool.free_object(p_edge);
		edge_pool.free_object(r);
		used_edge_pairs--;
	}

	void compute_internal(int32_t p_start, int32_t p_end, IntermediateHull &r_result);

	bool merge_projection(IntermediateHull &p_h0, IntermediateHull &p_h1, Vertex *&r_c0, Vertex *&r_c1);

	void merge(IntermediateHull &p_h0, IntermediateHull &p_h1);

	Vector3 to_gd_vector(const Point32 &p_v);

	Vector3 get_gd_normal(Face *p_face);

	bool shift_face(Face *p_face, real_t p_amount, LocalVector<Vertex *> p_stack);

public:
	Vertex *vertex_list;

	void compute(const Vector3 *p_coords, int32_t p_count);

	Vector3 get_coordinates(const Vertex *p_v);

	real_t shrink(real_t amount, real_t p_clamp_amount);
};

ConvexHullInternal::Int128 ConvexHullInternal::Int128::operator*(int64_t b) const {
	bool negative = (int64_t)high < 0;
	Int128 a = negative ? -*this : *this;
	if (b < 0) {
		negative = !negative;
		b = -b;
	}
	Int128 result = mul(a.low, (uint64_t)b);
	result.high += a.high * (uint64_t)b;
	return negative ? -result : result;
}

ConvexHullInternal::Int128 ConvexHullInternal::Int128::mul(int64_t a, int64_t b) {
	Int128 result;

#ifdef USE_X86_64_ASM
	__asm__("imulq %[b]"
			: "=a"(result.low), "=d"(result.high)
			: "0"(a), [b] "r"(b)
			: "cc");
	return result;

#else
	bool negative = a < 0;
	if (negative) {
		a = -a;
	}
	if (b < 0) {
		negative = !negative;
		b = -b;
	}
	DMul<uint64_t, uint32_t>::mul((uint64_t)a, (uint64_t)b, result.low, result.high);
	return negative ? -result : result;
#endif
}

ConvexHullInternal::Int128 ConvexHullInternal::Int128::mul(uint64_t a, uint64_t b) {
	Int128 result;

#ifdef USE_X86_64_ASM
	__asm__("mulq %[b]"
			: "=a"(result.low), "=d"(result.high)
			: "0"(a), [b] "r"(b)
			: "cc");

#else
	DMul<uint64_t, uint32_t>::mul(a, b, result.low, result.high);
#endif

	return result;
}

int32_t ConvexHullInternal::Rational64::compare(const Rational64 &b) const {
	if (sign != b.sign) {
		return sign - b.sign;
	} else if (sign == 0) {
		return 0;
	}

	//	return (numerator * b.denominator > b.numerator * denominator) ? sign : (numerator * b.denominator < b.numerator * denominator) ? -sign : 0;

#ifdef USE_X86_64_ASM

	int32_t result;
	int64_t tmp;
	int64_t dummy;
	__asm__("mulq %[bn]\n\t"
			"movq %%rax, %[tmp]\n\t"
			"movq %%rdx, %%rbx\n\t"
			"movq %[tn], %%rax\n\t"
			"mulq %[bd]\n\t"
			"subq %[tmp], %%rax\n\t"
			"sbbq %%rbx, %%rdx\n\t" // rdx:rax contains 128-bit-difference "numerator*b.denominator - b.numerator*denominator"
			"setnsb %%bh\n\t" // bh=1 if difference is non-negative, bh=0 otherwise
			"orq %%rdx, %%rax\n\t"
			"setnzb %%bl\n\t" // bl=1 if difference if non-zero, bl=0 if it is zero
			"decb %%bh\n\t" // now bx=0x0000 if difference is zero, 0xff01 if it is negative, 0x0001 if it is positive (i.e., same sign as difference)
			"shll $16, %%ebx\n\t" // ebx has same sign as difference
			: "=&b"(result), [tmp] "=&r"(tmp), "=a"(dummy)
			: "a"(denominator), [bn] "g"(b.numerator), [tn] "g"(numerator), [bd] "g"(b.denominator)
			: "%rdx", "cc");
	return result ? result ^ sign // if sign is +1, only bit 0 of result is inverted, which does not change the sign of result (and cannot result in zero)
					// if sign is -1, all bits of result are inverted, which changes the sign of result (and again cannot result in zero)
					:
					  0;

#else

	return sign * Int128::mul(numerator, b.denominator).ucmp(Int128::mul(denominator, b.numerator));

#endif
}

int32_t ConvexHullInternal::Rational128::compare(const Rational128 &b) const {
	if (sign != b.sign) {
		return sign - b.sign;
	} else if (sign == 0) {
		return 0;
	}
	if (is_int_64) {
		return -b.compare(sign * (int64_t)numerator.low);
	}

	Int128 nbdLow, nbdHigh, dbnLow, dbnHigh;
	DMul<Int128, uint64_t>::mul(numerator, b.denominator, nbdLow, nbdHigh);
	DMul<Int128, uint64_t>::mul(denominator, b.numerator, dbnLow, dbnHigh);

	int32_t cmp = nbdHigh.ucmp(dbnHigh);
	if (cmp) {
		return cmp * sign;
	}
	return nbdLow.ucmp(dbnLow) * sign;
}

int32_t ConvexHullInternal::Rational128::compare(int64_t b) const {
	if (is_int_64) {
		int64_t a = sign * (int64_t)numerator.low;
		return (a > b) ? 1 : (a < b) ? -1 :
										 0;
	}
	if (b > 0) {
		if (sign <= 0) {
			return -1;
		}
	} else if (b < 0) {
		if (sign >= 0) {
			return 1;
		}
		b = -b;
	} else {
		return sign;
	}

	return numerator.ucmp(denominator * b) * sign;
}

ConvexHullInternal::Edge *ConvexHullInternal::new_edge_pair(Vertex *p_from, Vertex *p_to) {
	CHULL_ASSERT(p_from && p_to);
	Edge *e = edge_pool.new_object();
	Edge *r = edge_pool.new_object();
	e->reverse = r;
	r->reverse = e;
	e->copy = merge_stamp;
	r->copy = merge_stamp;
	e->target = p_to;
	r->target = p_from;
	e->face = nullptr;
	r->face = nullptr;
	used_edge_pairs++;
	if (used_edge_pairs > max_used_edge_pairs) {
		max_used_edge_pairs = used_edge_pairs;
	}
	return e;
}

bool ConvexHullInternal::merge_projection(IntermediateHull &r_h0, IntermediateHull &r_h1, Vertex *&r_c0, Vertex *&r_c1) {
	Vertex *v0 = r_h0.max_yx;
	Vertex *v1 = r_h1.min_yx;
	if ((v0->point.x == v1->point.x) && (v0->point.y == v1->point.y)) {
		CHULL_ASSERT(v0->point.z < v1->point.z);
		Vertex *v1p = v1->prev;
		if (v1p == v1) {
			r_c0 = v0;
			if (v1->edges) {
				CHULL_ASSERT(v1->edges->next == v1->edges);
				v1 = v1->edges->target;
				CHULL_ASSERT(v1->edges->next == v1->edges);
			}
			r_c1 = v1;
			return false;
		}
		Vertex *v1n = v1->next;
		v1p->next = v1n;
		v1n->prev = v1p;
		if (v1 == r_h1.min_xy) {
			if ((v1n->point.x < v1p->point.x) || ((v1n->point.x == v1p->point.x) && (v1n->point.y < v1p->point.y))) {
				r_h1.min_xy = v1n;
			} else {
				r_h1.min_xy = v1p;
			}
		}
		if (v1 == r_h1.max_xy) {
			if ((v1n->point.x > v1p->point.x) || ((v1n->point.x == v1p->point.x) && (v1n->point.y > v1p->point.y))) {
				r_h1.max_xy = v1n;
			} else {
				r_h1.max_xy = v1p;
			}
		}
	}

	v0 = r_h0.max_xy;
	v1 = r_h1.max_xy;
	Vertex *v00 = nullptr;
	Vertex *v10 = nullptr;
	int32_t sign = 1;

	for (int32_t side = 0; side <= 1; side++) {
		int32_t dx = (v1->point.x - v0->point.x) * sign;
		if (dx > 0) {
			while (true) {
				int32_t dy = v1->point.y - v0->point.y;

				Vertex *w0 = side ? v0->next : v0->prev;
				if (w0 != v0) {
					int32_t dx0 = (w0->point.x - v0->point.x) * sign;
					int32_t dy0 = w0->point.y - v0->point.y;
					if ((dy0 <= 0) && ((dx0 == 0) || ((dx0 < 0) && (dy0 * dx <= dy * dx0)))) {
						v0 = w0;
						dx = (v1->point.x - v0->point.x) * sign;
						continue;
					}
				}

				Vertex *w1 = side ? v1->next : v1->prev;
				if (w1 != v1) {
					int32_t dx1 = (w1->point.x - v1->point.x) * sign;
					int32_t dy1 = w1->point.y - v1->point.y;
					int32_t dxn = (w1->point.x - v0->point.x) * sign;
					if ((dxn > 0) && (dy1 < 0) && ((dx1 == 0) || ((dx1 < 0) && (dy1 * dx < dy * dx1)))) {
						v1 = w1;
						dx = dxn;
						continue;
					}
				}

				break;
			}
		} else if (dx < 0) {
			while (true) {
				int32_t dy = v1->point.y - v0->point.y;

				Vertex *w1 = side ? v1->prev : v1->next;
				if (w1 != v1) {
					int32_t dx1 = (w1->point.x - v1->point.x) * sign;
					int32_t dy1 = w1->point.y - v1->point.y;
					if ((dy1 >= 0) && ((dx1 == 0) || ((dx1 < 0) && (dy1 * dx <= dy * dx1)))) {
						v1 = w1;
						dx = (v1->point.x - v0->point.x) * sign;
						continue;
					}
				}

				Vertex *w0 = side ? v0->prev : v0->next;
				if (w0 != v0) {
					int32_t dx0 = (w0->point.x - v0->point.x) * sign;
					int32_t dy0 = w0->point.y - v0->point.y;
					int32_t dxn = (v1->point.x - w0->point.x) * sign;
					if ((dxn < 0) && (dy0 > 0) && ((dx0 == 0) || ((dx0 < 0) && (dy0 * dx < dy * dx0)))) {
						v0 = w0;
						dx = dxn;
						continue;
					}
				}

				break;
			}
		} else {
			int32_t x = v0->point.x;
			int32_t y0 = v0->point.y;
			Vertex *w0 = v0;
			Vertex *t;
			while (((t = side ? w0->next : w0->prev) != v0) && (t->point.x == x) && (t->point.y <= y0)) {
				w0 = t;
				y0 = t->point.y;
			}
			v0 = w0;

			int32_t y1 = v1->point.y;
			Vertex *w1 = v1;
			while (((t = side ? w1->prev : w1->next) != v1) && (t->point.x == x) && (t->point.y >= y1)) {
				w1 = t;
				y1 = t->point.y;
			}
			v1 = w1;
		}

		if (side == 0) {
			v00 = v0;
			v10 = v1;

			v0 = r_h0.min_xy;
			v1 = r_h1.min_xy;
			sign = -1;
		}
	}

	v0->prev = v1;
	v1->next = v0;

	v00->next = v10;
	v10->prev = v00;

	if (r_h1.min_xy->point.x < r_h0.min_xy->point.x) {
		r_h0.min_xy = r_h1.min_xy;
	}
	if (r_h1.max_xy->point.x >= r_h0.max_xy->point.x) {
		r_h0.max_xy = r_h1.max_xy;
	}

	r_h0.max_yx = r_h1.max_yx;

	r_c0 = v00;
	r_c1 = v10;

	return true;
}

void ConvexHullInternal::compute_internal(int32_t p_start, int32_t p_end, IntermediateHull &r_result) {
	int32_t n = p_end - p_start;
	switch (n) {
		case 0:
			r_result.min_xy = nullptr;
			r_result.max_xy = nullptr;
			r_result.min_yx = nullptr;
			r_result.max_yx = nullptr;
			return;
		case 2: {
			Vertex *v = original_vertices[p_start];
			Vertex *w = v + 1;
			if (v->point != w->point) {
				int32_t dx = v->point.x - w->point.x;
				int32_t dy = v->point.y - w->point.y;

				if ((dx == 0) && (dy == 0)) {
					if (v->point.z > w->point.z) {
						Vertex *t = w;
						w = v;
						v = t;
					}
					CHULL_ASSERT(v->point.z < w->point.z);
					v->next = v;
					v->prev = v;
					r_result.min_xy = v;
					r_result.max_xy = v;
					r_result.min_yx = v;
					r_result.max_yx = v;
				} else {
					v->next = w;
					v->prev = w;
					w->next = v;
					w->prev = v;

					if ((dx < 0) || ((dx == 0) && (dy < 0))) {
						r_result.min_xy = v;
						r_result.max_xy = w;
					} else {
						r_result.min_xy = w;
						r_result.max_xy = v;
					}

					if ((dy < 0) || ((dy == 0) && (dx < 0))) {
						r_result.min_yx = v;
						r_result.max_yx = w;
					} else {
						r_result.min_yx = w;
						r_result.max_yx = v;
					}
				}

				Edge *e = new_edge_pair(v, w);
				e->link(e);
				v->edges = e;

				e = e->reverse;
				e->link(e);
				w->edges = e;

				return;
			}
		}
		// lint -fallthrough
		case 1: {
			Vertex *v = original_vertices[p_start];
			v->edges = nullptr;
			v->next = v;
			v->prev = v;

			r_result.min_xy = v;
			r_result.max_xy = v;
			r_result.min_yx = v;
			r_result.max_yx = v;

			return;
		}
	}

	int32_t split0 = p_start + n / 2;
	Point32 p = original_vertices[split0 - 1]->point;
	int32_t split1 = split0;
	while ((split1 < p_end) && (original_vertices[split1]->point == p)) {
		split1++;
	}
	compute_internal(p_start, split0, r_result);
	IntermediateHull hull1;
	compute_internal(split1, p_end, hull1);
#ifdef DEBUG_CONVEX_HULL
	printf("\n\nMerge\n");
	r_result.print();
	hull1.print();
#endif
	merge(r_result, hull1);
#ifdef DEBUG_CONVEX_HULL
	printf("\n  Result\n");
	r_result.print();
#endif
}

#ifdef DEBUG_CONVEX_HULL
void ConvexHullInternal::IntermediateHull::print() {
	printf("    Hull\n");
	for (Vertex *v = min_xy; v;) {
		printf("      ");
		v->print();
		if (v == max_xy) {
			printf(" max_xy");
		}
		if (v == min_yx) {
			printf(" min_yx");
		}
		if (v == max_yx) {
			printf(" max_yx");
		}
		if (v->next->prev != v) {
			printf(" Inconsistency");
		}
		printf("\n");
		v = v->next;
		if (v == min_xy) {
			break;
		}
	}
	if (min_xy) {
		min_xy->copy = (min_xy->copy == -1) ? -2 : -1;
		min_xy->print_graph();
	}
}

void ConvexHullInternal::Vertex::print_graph() {
	print();
	printf("\nEdges\n");
	Edge *e = edges;
	if (e) {
		do {
			e->print();
			printf("\n");
			e = e->next;
		} while (e != edges);
		do {
			Vertex *v = e->target;
			if (v->copy != copy) {
				v->copy = copy;
				v->print_graph();
			}
			e = e->next;
		} while (e != edges);
	}
}
#endif

ConvexHullInternal::Orientation ConvexHullInternal::get_orientation(const Edge *p_prev, const Edge *p_next, const Point32 &p_s, const Point32 &p_t) {
	CHULL_ASSERT(p_prev->reverse->target == p_next->reverse->target);
	if (p_prev->next == p_next) {
		if (p_prev->prev == p_next) {
			Point64 n = p_t.cross(p_s);
			Point64 m = (*p_prev->target - *p_next->reverse->target).cross(*p_next->target - *p_next->reverse->target);
			CHULL_ASSERT(!m.is_zero());
			int64_t dot = n.dot(m);
			CHULL_ASSERT(dot != 0);
			return (dot > 0) ? COUNTER_CLOCKWISE : CLOCKWISE;
		}
		return COUNTER_CLOCKWISE;
	} else if (p_prev->prev == p_next) {
		return CLOCKWISE;
	} else {
		return NONE;
	}
}

ConvexHullInternal::Edge *ConvexHullInternal::find_max_angle(bool p_ccw, const Vertex *p_start, const Point32 &p_s, const Point64 &p_rxs, const Point64 &p_sxrxs, Rational64 &p_min_cot) {
	Edge *min_edge = nullptr;

#ifdef DEBUG_CONVEX_HULL
	printf("find max edge for %d\n", p_start->point.index);
#endif
	Edge *e = p_start->edges;
	if (e) {
		do {
			if (e->copy > merge_stamp) {
				Point32 t = *e->target - *p_start;
				Rational64 cot(t.dot(p_sxrxs), t.dot(p_rxs));
#ifdef DEBUG_CONVEX_HULL
				printf("      Angle is %f (%d) for ", Math::atan(cot.to_scalar()), (int32_t)cot.is_nan());
				e->print();
#endif
				if (cot.is_nan()) {
					CHULL_ASSERT(p_ccw ? (t.dot(p_s) < 0) : (t.dot(p_s) > 0));
				} else {
					int32_t cmp;
					if (min_edge == nullptr) {
						p_min_cot = cot;
						min_edge = e;
					} else if ((cmp = cot.compare(p_min_cot)) < 0) {
						p_min_cot = cot;
						min_edge = e;
					} else if ((cmp == 0) && (p_ccw == (get_orientation(min_edge, e, p_s, t) == COUNTER_CLOCKWISE))) {
						min_edge = e;
					}
				}
#ifdef DEBUG_CONVEX_HULL
				printf("\n");
#endif
			}
			e = e->next;
		} while (e != p_start->edges);
	}
	return min_edge;
}

void ConvexHullInternal::find_edge_for_coplanar_faces(Vertex *p_c0, Vertex *p_c1, Edge *&p_e0, Edge *&p_e1, Vertex *p_stop0, Vertex *p_stop1) {
	Edge *start0 = p_e0;
	Edge *start1 = p_e1;
	Point32 et0 = start0 ? start0->target->point : p_c0->point;
	Point32 et1 = start1 ? start1->target->point : p_c1->point;
	Point32 s = p_c1->point - p_c0->point;
	Point64 normal = ((start0 ? start0 : start1)->target->point - p_c0->point).cross(s);
	int64_t dist = p_c0->point.dot(normal);
	CHULL_ASSERT(!start1 || (start1->target->point.dot(normal) == dist));
	Point64 perp = s.cross(normal);
	CHULL_ASSERT(!perp.is_zero());

#ifdef DEBUG_CONVEX_HULL
	printf("   Advancing %d %d  (%p %p, %d %d)\n", p_c0->point.index, p_c1->point.index, start0, start1, start0 ? start0->target->point.index : -1, start1 ? start1->target->point.index : -1);
#endif

	int64_t maxDot0 = et0.dot(perp);
	if (p_e0) {
		while (p_e0->target != p_stop0) {
			Edge *e = p_e0->reverse->prev;
			if (e->target->point.dot(normal) < dist) {
				break;
			}
			CHULL_ASSERT(e->target->point.dot(normal) == dist);
			if (e->copy == merge_stamp) {
				break;
			}
			int64_t dot = e->target->point.dot(perp);
			if (dot <= maxDot0) {
				break;
			}
			maxDot0 = dot;
			p_e0 = e;
			et0 = e->target->point;
		}
	}

	int64_t maxDot1 = et1.dot(perp);
	if (p_e1) {
		while (p_e1->target != p_stop1) {
			Edge *e = p_e1->reverse->next;
			if (e->target->point.dot(normal) < dist) {
				break;
			}
			CHULL_ASSERT(e->target->point.dot(normal) == dist);
			if (e->copy == merge_stamp) {
				break;
			}
			int64_t dot = e->target->point.dot(perp);
			if (dot <= maxDot1) {
				break;
			}
			maxDot1 = dot;
			p_e1 = e;
			et1 = e->target->point;
		}
	}

#ifdef DEBUG_CONVEX_HULL
	printf("   Starting at %d %d\n", et0.index, et1.index);
#endif

	int64_t dx = maxDot1 - maxDot0;
	if (dx > 0) {
		while (true) {
			int64_t dy = (et1 - et0).dot(s);

			if (p_e0 && (p_e0->target != p_stop0)) {
				Edge *f0 = p_e0->next->reverse;
				if (f0->copy > merge_stamp) {
					int64_t dx0 = (f0->target->point - et0).dot(perp);
					int64_t dy0 = (f0->target->point - et0).dot(s);
					if ((dx0 == 0) ? (dy0 < 0) : ((dx0 < 0) && (Rational64(dy0, dx0).compare(Rational64(dy, dx)) >= 0))) {
						et0 = f0->target->point;
						dx = (et1 - et0).dot(perp);
						p_e0 = (p_e0 == start0) ? nullptr : f0;
						continue;
					}
				}
			}

			if (p_e1 && (p_e1->target != p_stop1)) {
				Edge *f1 = p_e1->reverse->next;
				if (f1->copy > merge_stamp) {
					Point32 d1 = f1->target->point - et1;
					if (d1.dot(normal) == 0) {
						int64_t dx1 = d1.dot(perp);
						int64_t dy1 = d1.dot(s);
						int64_t dxn = (f1->target->point - et0).dot(perp);
						if ((dxn > 0) && ((dx1 == 0) ? (dy1 < 0) : ((dx1 < 0) && (Rational64(dy1, dx1).compare(Rational64(dy, dx)) > 0)))) {
							p_e1 = f1;
							et1 = p_e1->target->point;
							dx = dxn;
							continue;
						}
					} else {
						CHULL_ASSERT((p_e1 == start1) && (d1.dot(normal) < 0));
					}
				}
			}

			break;
		}
	} else if (dx < 0) {
		while (true) {
			int64_t dy = (et1 - et0).dot(s);

			if (p_e1 && (p_e1->target != p_stop1)) {
				Edge *f1 = p_e1->prev->reverse;
				if (f1->copy > merge_stamp) {
					int64_t dx1 = (f1->target->point - et1).dot(perp);
					int64_t dy1 = (f1->target->point - et1).dot(s);
					if ((dx1 == 0) ? (dy1 > 0) : ((dx1 < 0) && (Rational64(dy1, dx1).compare(Rational64(dy, dx)) <= 0))) {
						et1 = f1->target->point;
						dx = (et1 - et0).dot(perp);
						p_e1 = (p_e1 == start1) ? nullptr : f1;
						continue;
					}
				}
			}

			if (p_e0 && (p_e0->target != p_stop0)) {
				Edge *f0 = p_e0->reverse->prev;
				if (f0->copy > merge_stamp) {
					Point32 d0 = f0->target->point - et0;
					if (d0.dot(normal) == 0) {
						int64_t dx0 = d0.dot(perp);
						int64_t dy0 = d0.dot(s);
						int64_t dxn = (et1 - f0->target->point).dot(perp);
						if ((dxn < 0) && ((dx0 == 0) ? (dy0 > 0) : ((dx0 < 0) && (Rational64(dy0, dx0).compare(Rational64(dy, dx)) < 0)))) {
							p_e0 = f0;
							et0 = p_e0->target->point;
							dx = dxn;
							continue;
						}
					} else {
						CHULL_ASSERT((p_e0 == start0) && (d0.dot(normal) < 0));
					}
				}
			}

			break;
		}
	}
#ifdef DEBUG_CONVEX_HULL
	printf("   Advanced edges to %d %d\n", et0.index, et1.index);
#endif
}

void ConvexHullInternal::merge(IntermediateHull &p_h0, IntermediateHull &p_h1) {
	if (!p_h1.max_xy) {
		return;
	}
	if (!p_h0.max_xy) {
		p_h0 = p_h1;
		return;
	}

	merge_stamp--;

	Vertex *c0 = nullptr;
	Edge *toPrev0 = nullptr;
	Edge *firstNew0 = nullptr;
	Edge *pendingHead0 = nullptr;
	Edge *pendingTail0 = nullptr;
	Vertex *c1 = nullptr;
	Edge *toPrev1 = nullptr;
	Edge *firstNew1 = nullptr;
	Edge *pendingHead1 = nullptr;
	Edge *pendingTail1 = nullptr;
	Point32 prevPoint;

	if (merge_projection(p_h0, p_h1, c0, c1)) {
		Point32 s = *c1 - *c0;
		Point64 normal = Point32(0, 0, -1).cross(s);
		Point64 t = s.cross(normal);
		CHULL_ASSERT(!t.is_zero());

		Edge *e = c0->edges;
		Edge *start0 = nullptr;
		if (e) {
			do {
				int64_t dot = (*e->target - *c0).dot(normal);
				CHULL_ASSERT(dot <= 0);
				if ((dot == 0) && ((*e->target - *c0).dot(t) > 0)) {
					if (!start0 || (get_orientation(start0, e, s, Point32(0, 0, -1)) == CLOCKWISE)) {
						start0 = e;
					}
				}
				e = e->next;
			} while (e != c0->edges);
		}

		e = c1->edges;
		Edge *start1 = nullptr;
		if (e) {
			do {
				int64_t dot = (*e->target - *c1).dot(normal);
				CHULL_ASSERT(dot <= 0);
				if ((dot == 0) && ((*e->target - *c1).dot(t) > 0)) {
					if (!start1 || (get_orientation(start1, e, s, Point32(0, 0, -1)) == COUNTER_CLOCKWISE)) {
						start1 = e;
					}
				}
				e = e->next;
			} while (e != c1->edges);
		}

		if (start0 || start1) {
			find_edge_for_coplanar_faces(c0, c1, start0, start1, nullptr, nullptr);
			if (start0) {
				c0 = start0->target;
			}
			if (start1) {
				c1 = start1->target;
			}
		}

		prevPoint = c1->point;
		prevPoint.z++;
	} else {
		prevPoint = c1->point;
		prevPoint.x++;
	}

	Vertex *first0 = c0;
	Vertex *first1 = c1;
	bool firstRun = true;

	while (true) {
		Point32 s = *c1 - *c0;
		Point32 r = prevPoint - c0->point;
		Point64 rxs = r.cross(s);
		Point64 sxrxs = s.cross(rxs);

#ifdef DEBUG_CONVEX_HULL
		printf("\n  Checking %d %d\n", c0->point.index, c1->point.index);
#endif
		Rational64 minCot0(0, 0);
		Edge *min0 = find_max_angle(false, c0, s, rxs, sxrxs, minCot0);
		Rational64 minCot1(0, 0);
		Edge *min1 = find_max_angle(true, c1, s, rxs, sxrxs, minCot1);
		if (!min0 && !min1) {
			Edge *e = new_edge_pair(c0, c1);
			e->link(e);
			c0->edges = e;

			e = e->reverse;
			e->link(e);
			c1->edges = e;
			return;
		} else {
			int32_t cmp = !min0 ? 1 : !min1 ? -1 :
												minCot0.compare(minCot1);
#ifdef DEBUG_CONVEX_HULL
			printf("    -> Result %d\n", cmp);
#endif
			if (firstRun || ((cmp >= 0) ? !minCot1.is_negative_infinity() : !minCot0.is_negative_infinity())) {
				Edge *e = new_edge_pair(c0, c1);
				if (pendingTail0) {
					pendingTail0->prev = e;
				} else {
					pendingHead0 = e;
				}
				e->next = pendingTail0;
				pendingTail0 = e;

				e = e->reverse;
				if (pendingTail1) {
					pendingTail1->next = e;
				} else {
					pendingHead1 = e;
				}
				e->prev = pendingTail1;
				pendingTail1 = e;
			}

			Edge *e0 = min0;
			Edge *e1 = min1;

#ifdef DEBUG_CONVEX_HULL
			printf("   Found min edges to %d %d\n", e0 ? e0->target->point.index : -1, e1 ? e1->target->point.index : -1);
#endif

			if (cmp == 0) {
				find_edge_for_coplanar_faces(c0, c1, e0, e1, nullptr, nullptr);
			}

			if ((cmp >= 0) && e1) {
				if (toPrev1) {
					for (Edge *e = toPrev1->next, *n = nullptr; e != min1; e = n) {
						n = e->next;
						remove_edge_pair(e);
					}
				}

				if (pendingTail1) {
					if (toPrev1) {
						toPrev1->link(pendingHead1);
					} else {
						min1->prev->link(pendingHead1);
						firstNew1 = pendingHead1;
					}
					pendingTail1->link(min1);
					pendingHead1 = nullptr;
					pendingTail1 = nullptr;
				} else if (!toPrev1) {
					firstNew1 = min1;
				}

				prevPoint = c1->point;
				c1 = e1->target;
				toPrev1 = e1->reverse;
			}

			if ((cmp <= 0) && e0) {
				if (toPrev0) {
					for (Edge *e = toPrev0->prev, *n = nullptr; e != min0; e = n) {
						n = e->prev;
						remove_edge_pair(e);
					}
				}

				if (pendingTail0) {
					if (toPrev0) {
						pendingHead0->link(toPrev0);
					} else {
						pendingHead0->link(min0->next);
						firstNew0 = pendingHead0;
					}
					min0->link(pendingTail0);
					pendingHead0 = nullptr;
					pendingTail0 = nullptr;
				} else if (!toPrev0) {
					firstNew0 = min0;
				}

				prevPoint = c0->point;
				c0 = e0->target;
				toPrev0 = e0->reverse;
			}
		}

		if ((c0 == first0) && (c1 == first1)) {
			if (toPrev0 == nullptr) {
				pendingHead0->link(pendingTail0);
				c0->edges = pendingTail0;
			} else {
				for (Edge *e = toPrev0->prev, *n = nullptr; e != firstNew0; e = n) {
					n = e->prev;
					remove_edge_pair(e);
				}
				if (pendingTail0) {
					pendingHead0->link(toPrev0);
					firstNew0->link(pendingTail0);
				}
			}

			if (toPrev1 == nullptr) {
				pendingTail1->link(pendingHead1);
				c1->edges = pendingTail1;
			} else {
				for (Edge *e = toPrev1->next, *n = nullptr; e != firstNew1; e = n) {
					n = e->next;
					remove_edge_pair(e);
				}
				if (pendingTail1) {
					toPrev1->link(pendingHead1);
					pendingTail1->link(firstNew1);
				}
			}

			return;
		}

		firstRun = false;
	}
}

struct PointComparator {
	_FORCE_INLINE_ bool operator()(const ConvexHullInternal::Point32 &p, const ConvexHullInternal::Point32 &q) const {
		return (p.y < q.y) || ((p.y == q.y) && ((p.x < q.x) || ((p.x == q.x) && (p.z < q.z))));
	}
};

void ConvexHullInternal::compute(const Vector3 *p_coords, int32_t p_count) {
	AABB aabb;
	for (int32_t i = 0; i < p_count; i++) {
		Vector3 p = p_coords[i];
		if (i == 0) {
			aabb.position = p;
		} else {
			aabb.expand_to(p);
		}
	}

	Vector3 s = aabb.size;
	max_axis = s.max_axis();
	min_axis = s.min_axis();
	if (min_axis == max_axis) {
		min_axis = (max_axis + 1) % 3;
	}
	med_axis = 3 - max_axis - min_axis;

	s /= real_t(10216);
	if (((med_axis + 1) % 3) != max_axis) {
		s *= -1;
	}
	scaling = s;

	if (s[0] != 0) {
		s[0] = real_t(1) / s[0];
	}
	if (s[1] != 0) {
		s[1] = real_t(1) / s[1];
	}
	if (s[2] != 0) {
		s[2] = real_t(1) / s[2];
	}

	center = aabb.position;

	LocalVector<Point32> points;
	points.resize(p_count);
	for (int32_t i = 0; i < p_count; i++) {
		Vector3 p = p_coords[i];
		p = (p - center) * s;
		points[i].x = (int32_t)p[med_axis];
		points[i].y = (int32_t)p[max_axis];
		points[i].z = (int32_t)p[min_axis];
		points[i].index = i;
	}

	points.sort_custom<PointComparator>();

	vertex_pool.reset();
	vertex_pool.set_array_size(p_count);
	original_vertices.resize(p_count);
	for (int32_t i = 0; i < p_count; i++) {
		Vertex *v = vertex_pool.new_object();
		v->edges = nullptr;
		v->point = points[i];
		v->copy = -1;
		original_vertices[i] = v;
	}

	points.clear();

	edge_pool.reset();
	edge_pool.set_array_size(6 * p_count);

	used_edge_pairs = 0;
	max_used_edge_pairs = 0;

	merge_stamp = -3;

	IntermediateHull hull;
	compute_internal(0, p_count, hull);
	vertex_list = hull.min_xy;
#ifdef DEBUG_CONVEX_HULL
	printf("max. edges %d (3v = %d)", max_used_edge_pairs, 3 * p_count);
#endif
}

Vector3 ConvexHullInternal::to_gd_vector(const Point32 &p_v) {
	Vector3 p;
	p[med_axis] = real_t(p_v.x);
	p[max_axis] = real_t(p_v.y);
	p[min_axis] = real_t(p_v.z);
	return p * scaling;
}

Vector3 ConvexHullInternal::get_gd_normal(Face *p_face) {
	return to_gd_vector(p_face->dir0).cross(to_gd_vector(p_face->dir1)).normalized();
}

Vector3 ConvexHullInternal::get_coordinates(const Vertex *p_v) {
	Vector3 p;
	p[med_axis] = p_v->xvalue();
	p[max_axis] = p_v->yvalue();
	p[min_axis] = p_v->zvalue();
	return p * scaling + center;
}

real_t ConvexHullInternal::shrink(real_t p_amount, real_t p_clamp_amount) {
	if (!vertex_list) {
		return 0;
	}
	int32_t stamp = --merge_stamp;
	LocalVector<Vertex *> stack;
	vertex_list->copy = stamp;
	stack.push_back(vertex_list);
	LocalVector<Face *> faces;

	Point32 ref = vertex_list->point;
	Int128 hullCenterX(0, 0);
	Int128 hullCenterY(0, 0);
	Int128 hullCenterZ(0, 0);
	Int128 volume(0, 0);

	while (stack.size() > 0) {
		Vertex *v = stack[stack.size() - 1];
		stack.remove(stack.size() - 1);
		Edge *e = v->edges;
		if (e) {
			do {
				if (e->target->copy != stamp) {
					e->target->copy = stamp;
					stack.push_back(e->target);
				}
				if (e->copy != stamp) {
					Face *face = face_pool.new_object();
					face->init(e->target, e->reverse->prev->target, v);
					faces.push_back(face);
					Edge *f = e;

					Vertex *a = nullptr;
					Vertex *b = nullptr;
					do {
						if (a && b) {
							int64_t vol = (v->point - ref).dot((a->point - ref).cross(b->point - ref));
							CHULL_ASSERT(vol >= 0);
							Point32 c = v->point + a->point + b->point + ref;
							hullCenterX += vol * c.x;
							hullCenterY += vol * c.y;
							hullCenterZ += vol * c.z;
							volume += vol;
						}

						CHULL_ASSERT(f->copy != stamp);
						f->copy = stamp;
						f->face = face;

						a = b;
						b = f->target;

						f = f->reverse->prev;
					} while (f != e);
				}
				e = e->next;
			} while (e != v->edges);
		}
	}

	if (volume.get_sign() <= 0) {
		return 0;
	}

	Vector3 hullCenter;
	hullCenter[med_axis] = hullCenterX.to_scalar();
	hullCenter[max_axis] = hullCenterY.to_scalar();
	hullCenter[min_axis] = hullCenterZ.to_scalar();
	hullCenter /= 4 * volume.to_scalar();
	hullCenter *= scaling;

	int32_t faceCount = faces.size();

	if (p_clamp_amount > 0) {
		real_t minDist = FLT_MAX;
		for (int32_t i = 0; i < faceCount; i++) {
			Vector3 normal = get_gd_normal(faces[i]);
			real_t dist = normal.dot(to_gd_vector(faces[i]->origin) - hullCenter);
			if (dist < minDist) {
				minDist = dist;
			}
		}

		if (minDist <= 0) {
			return 0;
		}

		p_amount = MIN(p_amount, minDist * p_clamp_amount);
	}

	uint32_t seed = 243703;
	for (int32_t i = 0; i < faceCount; i++, seed = 1664525 * seed + 1013904223) {
		SWAP(faces[i], faces[seed % faceCount]);
	}

	for (int32_t i = 0; i < faceCount; i++) {
		if (!shift_face(faces[i], p_amount, stack)) {
			return -p_amount;
		}
	}

	return p_amount;
}

bool ConvexHullInternal::shift_face(Face *p_face, real_t p_amount, LocalVector<Vertex *> p_stack) {
	Vector3 origShift = get_gd_normal(p_face) * -p_amount;
	if (scaling[0] != 0) {
		origShift[0] /= scaling[0];
	}
	if (scaling[1] != 0) {
		origShift[1] /= scaling[1];
	}
	if (scaling[2] != 0) {
		origShift[2] /= scaling[2];
	}
	Point32 shift((int32_t)origShift[med_axis], (int32_t)origShift[max_axis], (int32_t)origShift[min_axis]);
	if (shift.is_zero()) {
		return true;
	}
	Point64 normal = p_face->getNormal();
#ifdef DEBUG_CONVEX_HULL
	printf("\nShrinking p_face (%d %d %d) (%d %d %d) (%d %d %d) by (%d %d %d)\n",
			p_face->origin.x, p_face->origin.y, p_face->origin.z, p_face->dir0.x, p_face->dir0.y, p_face->dir0.z, p_face->dir1.x, p_face->dir1.y, p_face->dir1.z, shift.x, shift.y, shift.z);
#endif
	int64_t origDot = p_face->origin.dot(normal);
	Point32 shiftedOrigin = p_face->origin + shift;
	int64_t shiftedDot = shiftedOrigin.dot(normal);
	CHULL_ASSERT(shiftedDot <= origDot);
	if (shiftedDot >= origDot) {
		return false;
	}

	Edge *intersection = nullptr;

	Edge *startEdge = p_face->nearby_vertex->edges;
#ifdef DEBUG_CONVEX_HULL
	printf("Start edge is ");
	startEdge->print();
	printf(", normal is (%lld %lld %lld), shifted dot is %lld\n", normal.x, normal.y, normal.z, shiftedDot);
#endif
	Rational128 optDot = p_face->nearby_vertex->dot(normal);
	int32_t cmp = optDot.compare(shiftedDot);
#ifdef SHOW_ITERATIONS
	int32_t n = 0;
#endif
	if (cmp >= 0) {
		Edge *e = startEdge;
		do {
#ifdef SHOW_ITERATIONS
			n++;
#endif
			Rational128 dot = e->target->dot(normal);
			CHULL_ASSERT(dot.compare(origDot) <= 0);
#ifdef DEBUG_CONVEX_HULL
			printf("Moving downwards, edge is ");
			e->print();
			printf(", dot is %f (%f %lld)\n", (float)dot.to_scalar(), (float)optDot.to_scalar(), shiftedDot);
#endif
			if (dot.compare(optDot) < 0) {
				int32_t c = dot.compare(shiftedDot);
				optDot = dot;
				e = e->reverse;
				startEdge = e;
				if (c < 0) {
					intersection = e;
					break;
				}
				cmp = c;
			}
			e = e->prev;
		} while (e != startEdge);

		if (!intersection) {
			return false;
		}
	} else {
		Edge *e = startEdge;
		do {
#ifdef SHOW_ITERATIONS
			n++;
#endif
			Rational128 dot = e->target->dot(normal);
			CHULL_ASSERT(dot.compare(origDot) <= 0);
#ifdef DEBUG_CONVEX_HULL
			printf("Moving upwards, edge is ");
			e->print();
			printf(", dot is %f (%f %lld)\n", (float)dot.to_scalar(), (float)optDot.to_scalar(), shiftedDot);
#endif
			if (dot.compare(optDot) > 0) {
				cmp = dot.compare(shiftedDot);
				if (cmp >= 0) {
					intersection = e;
					break;
				}
				optDot = dot;
				e = e->reverse;
				startEdge = e;
			}
			e = e->prev;
		} while (e != startEdge);

		if (!intersection) {
			return true;
		}
	}

#ifdef SHOW_ITERATIONS
	printf("Needed %d iterations to find initial intersection\n", n);
#endif

	if (cmp == 0) {
		Edge *e = intersection->reverse->next;
#ifdef SHOW_ITERATIONS
		n = 0;
#endif
		while (e->target->dot(normal).compare(shiftedDot) <= 0) {
#ifdef SHOW_ITERATIONS
			n++;
#endif
			e = e->next;
			if (e == intersection->reverse) {
				return true;
			}
#ifdef DEBUG_CONVEX_HULL
			printf("Checking for outwards edge, current edge is ");
			e->print();
			printf("\n");
#endif
		}
#ifdef SHOW_ITERATIONS
		printf("Needed %d iterations to check for complete containment\n", n);
#endif
	}

	Edge *firstIntersection = nullptr;
	Edge *faceEdge = nullptr;
	Edge *firstFaceEdge = nullptr;

#ifdef SHOW_ITERATIONS
	int32_t m = 0;
#endif
	while (true) {
#ifdef SHOW_ITERATIONS
		m++;
#endif
#ifdef DEBUG_CONVEX_HULL
		printf("Intersecting edge is ");
		intersection->print();
		printf("\n");
#endif
		if (cmp == 0) {
			Edge *e = intersection->reverse->next;
			startEdge = e;
#ifdef SHOW_ITERATIONS
			n = 0;
#endif
			while (true) {
#ifdef SHOW_ITERATIONS
				n++;
#endif
				if (e->target->dot(normal).compare(shiftedDot) >= 0) {
					break;
				}
				intersection = e->reverse;
				e = e->next;
				if (e == startEdge) {
					return true;
				}
			}
#ifdef SHOW_ITERATIONS
			printf("Needed %d iterations to advance intersection\n", n);
#endif
		}

#ifdef DEBUG_CONVEX_HULL
		printf("Advanced intersecting edge to ");
		intersection->print();
		printf(", cmp = %d\n", cmp);
#endif

		if (!firstIntersection) {
			firstIntersection = intersection;
		} else if (intersection == firstIntersection) {
			break;
		}

		int32_t prevCmp = cmp;
		Edge *prevIntersection = intersection;
		Edge *prevFaceEdge = faceEdge;

		Edge *e = intersection->reverse;
#ifdef SHOW_ITERATIONS
		n = 0;
#endif
		while (true) {
#ifdef SHOW_ITERATIONS
			n++;
#endif
			e = e->reverse->prev;
			CHULL_ASSERT(e != intersection->reverse);
			cmp = e->target->dot(normal).compare(shiftedDot);
#ifdef DEBUG_CONVEX_HULL
			printf("Testing edge ");
			e->print();
			printf(" -> cmp = %d\n", cmp);
#endif
			if (cmp >= 0) {
				intersection = e;
				break;
			}
		}
#ifdef SHOW_ITERATIONS
		printf("Needed %d iterations to find other intersection of p_face\n", n);
#endif

		if (cmp > 0) {
			Vertex *removed = intersection->target;
			e = intersection->reverse;
			if (e->prev == e) {
				removed->edges = nullptr;
			} else {
				removed->edges = e->prev;
				e->prev->link(e->next);
				e->link(e);
			}
#ifdef DEBUG_CONVEX_HULL
			printf("1: Removed part contains (%d %d %d)\n", removed->point.x, removed->point.y, removed->point.z);
#endif

			Point64 n0 = intersection->face->getNormal();
			Point64 n1 = intersection->reverse->face->getNormal();
			int64_t m00 = p_face->dir0.dot(n0);
			int64_t m01 = p_face->dir1.dot(n0);
			int64_t m10 = p_face->dir0.dot(n1);
			int64_t m11 = p_face->dir1.dot(n1);
			int64_t r0 = (intersection->face->origin - shiftedOrigin).dot(n0);
			int64_t r1 = (intersection->reverse->face->origin - shiftedOrigin).dot(n1);
			Int128 det = Int128::mul(m00, m11) - Int128::mul(m01, m10);
			CHULL_ASSERT(det.get_sign() != 0);
			Vertex *v = vertex_pool.new_object();
			v->point.index = -1;
			v->copy = -1;
			v->point128 = PointR128(Int128::mul(p_face->dir0.x * r0, m11) - Int128::mul(p_face->dir0.x * r1, m01) + Int128::mul(p_face->dir1.x * r1, m00) - Int128::mul(p_face->dir1.x * r0, m10) + det * shiftedOrigin.x,
					Int128::mul(p_face->dir0.y * r0, m11) - Int128::mul(p_face->dir0.y * r1, m01) + Int128::mul(p_face->dir1.y * r1, m00) - Int128::mul(p_face->dir1.y * r0, m10) + det * shiftedOrigin.y,
					Int128::mul(p_face->dir0.z * r0, m11) - Int128::mul(p_face->dir0.z * r1, m01) + Int128::mul(p_face->dir1.z * r1, m00) - Int128::mul(p_face->dir1.z * r0, m10) + det * shiftedOrigin.z,
					det);
			v->point.x = (int32_t)v->point128.xvalue();
			v->point.y = (int32_t)v->point128.yvalue();
			v->point.z = (int32_t)v->point128.zvalue();
			intersection->target = v;
			v->edges = e;

			p_stack.push_back(v);
			p_stack.push_back(removed);
			p_stack.push_back(nullptr);
		}

		if (cmp || prevCmp || (prevIntersection->reverse->next->target != intersection->target)) {
			faceEdge = new_edge_pair(prevIntersection->target, intersection->target);
			if (prevCmp == 0) {
				faceEdge->link(prevIntersection->reverse->next);
			}
			if ((prevCmp == 0) || prevFaceEdge) {
				prevIntersection->reverse->link(faceEdge);
			}
			if (cmp == 0) {
				intersection->reverse->prev->link(faceEdge->reverse);
			}
			faceEdge->reverse->link(intersection->reverse);
		} else {
			faceEdge = prevIntersection->reverse->next;
		}

		if (prevFaceEdge) {
			if (prevCmp > 0) {
				faceEdge->link(prevFaceEdge->reverse);
			} else if (faceEdge != prevFaceEdge->reverse) {
				p_stack.push_back(prevFaceEdge->target);
				while (faceEdge->next != prevFaceEdge->reverse) {
					Vertex *removed = faceEdge->next->target;
					remove_edge_pair(faceEdge->next);
					p_stack.push_back(removed);
#ifdef DEBUG_CONVEX_HULL
					printf("2: Removed part contains (%d %d %d)\n", removed->point.x, removed->point.y, removed->point.z);
#endif
				}
				p_stack.push_back(nullptr);
			}
		}
		faceEdge->face = p_face;
		faceEdge->reverse->face = intersection->face;

		if (!firstFaceEdge) {
			firstFaceEdge = faceEdge;
		}
	}
#ifdef SHOW_ITERATIONS
	printf("Needed %d iterations to process all intersections\n", m);
#endif

	if (cmp > 0) {
		firstFaceEdge->reverse->target = faceEdge->target;
		firstIntersection->reverse->link(firstFaceEdge);
		firstFaceEdge->link(faceEdge->reverse);
	} else if (firstFaceEdge != faceEdge->reverse) {
		p_stack.push_back(faceEdge->target);
		while (firstFaceEdge->next != faceEdge->reverse) {
			Vertex *removed = firstFaceEdge->next->target;
			remove_edge_pair(firstFaceEdge->next);
			p_stack.push_back(removed);
#ifdef DEBUG_CONVEX_HULL
			printf("3: Removed part contains (%d %d %d)\n", removed->point.x, removed->point.y, removed->point.z);
#endif
		}
		p_stack.push_back(nullptr);
	}

	CHULL_ASSERT(p_stack.size() > 0);
	vertex_list = p_stack[0];

#ifdef DEBUG_CONVEX_HULL
	printf("Removing part\n");
#endif
#ifdef SHOW_ITERATIONS
	n = 0;
#endif
	uint32_t pos = 0;
	while (pos < p_stack.size()) {
		uint32_t end = p_stack.size();
		while (pos < end) {
			Vertex *kept = p_stack[pos++];
#ifdef DEBUG_CONVEX_HULL
			kept->print();
#endif
			bool deeper = false;
			Vertex *removed;
			while ((removed = p_stack[pos++]) != nullptr) {
#ifdef SHOW_ITERATIONS
				n++;
#endif
				kept->receive_nearby_faces(removed);
				while (removed->edges) {
					if (!deeper) {
						deeper = true;
						p_stack.push_back(kept);
					}
					p_stack.push_back(removed->edges->target);
					remove_edge_pair(removed->edges);
				}
			}
			if (deeper) {
				p_stack.push_back(nullptr);
			}
		}
	}
#ifdef SHOW_ITERATIONS
	printf("Needed %d iterations to remove part\n", n);
#endif

	p_stack.resize(0);
	p_face->origin = shiftedOrigin;

	return true;
}

static int32_t get_vertex_copy(ConvexHullInternal::Vertex *p_vertex, LocalVector<ConvexHullInternal::Vertex *> &p_vertices) {
	int32_t index = p_vertex->copy;
	if (index < 0) {
		index = p_vertices.size();
		p_vertex->copy = index;
		p_vertices.push_back(p_vertex);
#ifdef DEBUG_CONVEX_HULL
		printf("Vertex %d gets index *%d\n", p_vertex->point.index, index);
#endif
	}
	return index;
}

real_t ConvexHullComputer::compute(const Vector3 *p_coords, int32_t p_count, real_t p_shrink, real_t p_shrink_clamp) {
	if (p_count <= 0) {
		vertices.clear();
		edges.clear();
		faces.clear();
		return 0;
	}

	ConvexHullInternal hull;
	hull.compute(p_coords, p_count);

	real_t shift = 0;
	if ((p_shrink > 0) && ((shift = hull.shrink(p_shrink, p_shrink_clamp)) < 0)) {
		vertices.clear();
		edges.clear();
		faces.clear();
		return shift;
	}

	vertices.resize(0);
	edges.resize(0);
	faces.resize(0);

	LocalVector<ConvexHullInternal::Vertex *> oldVertices;
	get_vertex_copy(hull.vertex_list, oldVertices);
	int32_t copied = 0;
	while (copied < (int32_t)oldVertices.size()) {
		ConvexHullInternal::Vertex *v = oldVertices[copied];
		vertices.push_back(hull.get_coordinates(v));
		ConvexHullInternal::Edge *firstEdge = v->edges;
		if (firstEdge) {
			int32_t firstCopy = -1;
			int32_t prevCopy = -1;
			ConvexHullInternal::Edge *e = firstEdge;
			do {
				if (e->copy < 0) {
					int32_t s = edges.size();
					edges.push_back(Edge());
					edges.push_back(Edge());
					Edge *c = &edges[s];
					Edge *r = &edges[s + 1];
					e->copy = s;
					e->reverse->copy = s + 1;
					c->reverse = 1;
					r->reverse = -1;
					c->target_vertex = get_vertex_copy(e->target, oldVertices);
					r->target_vertex = copied;
#ifdef DEBUG_CONVEX_HULL
					printf("      CREATE: Vertex *%d has edge to *%d\n", copied, c->get_target_vertex());
#endif
				}
				if (prevCopy >= 0) {
					edges[e->copy].next = prevCopy - e->copy;
				} else {
					firstCopy = e->copy;
				}
				prevCopy = e->copy;
				e = e->next;
			} while (e != firstEdge);
			edges[firstCopy].next = prevCopy - firstCopy;
		}
		copied++;
	}

	for (int32_t i = 0; i < copied; i++) {
		ConvexHullInternal::Vertex *v = oldVertices[i];
		ConvexHullInternal::Edge *firstEdge = v->edges;
		if (firstEdge) {
			ConvexHullInternal::Edge *e = firstEdge;
			do {
				if (e->copy >= 0) {
#ifdef DEBUG_CONVEX_HULL
					printf("Vertex *%d has edge to *%d\n", i, edges[e->copy].get_target_vertex());
#endif
					faces.push_back(e->copy);
					ConvexHullInternal::Edge *f = e;
					do {
#ifdef DEBUG_CONVEX_HULL
						printf("   Face *%d\n", edges[f->copy].get_target_vertex());
#endif
						f->copy = -1;
						f = f->reverse->prev;
					} while (f != e);
				}
				e = e->next;
			} while (e != firstEdge);
		}
	}

	return shift;
}

Error ConvexHullComputer::convex_hull(const Vector<Vector3> &p_points, Geometry::MeshData &r_mesh) {
	r_mesh = Geometry::MeshData(); // clear

	if (p_points.size() == 0) {
		return FAILED; // matches QuickHull
	}

	ConvexHullComputer ch;
	ch.compute(p_points.ptr(), p_points.size(), -1.0, -1.0);

	r_mesh.vertices = ch.vertices;

	r_mesh.edges.resize(ch.edges.size());
	for (uint32_t i = 0; i < ch.edges.size(); i++) {
		r_mesh.edges.write[i].a = (&ch.edges[i])->get_source_vertex();
		r_mesh.edges.write[i].b = (&ch.edges[i])->get_target_vertex();
	}

	r_mesh.faces.resize(ch.faces.size());
	for (uint32_t i = 0; i < ch.faces.size(); i++) {
		const Edge *e_start = &ch.edges[ch.faces[i]];
		const Edge *e = e_start;
		Geometry::MeshData::Face &face = r_mesh.faces.write[i];

		do {
			face.indices.push_back(e->get_target_vertex());

			e = e->get_next_edge_of_face();
		} while (e != e_start);

		// compute normal
		if (face.indices.size() >= 3) {
			face.plane = Plane(r_mesh.vertices[face.indices[0]], r_mesh.vertices[face.indices[2]], r_mesh.vertices[face.indices[1]]);
		} else {
			WARN_PRINT("Too few vertices per face.");
		}
	}

	return OK;
}
