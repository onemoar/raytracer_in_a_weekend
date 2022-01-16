#if !defined(M_PI)
#define M_PI 3.14159265359f
#endif
#define DEG_TO_RAD (M_PI / 180)

union quaternion;

union vec3 {
	struct {
		float x, y, z;
	};

	struct {
		float r, g, b;
	};

	float E[3];

	vec3() {}
	vec3(float _x, float _y, float _z): x(_x), y(_y), z(_z) {}
	vec3(const vec3 &v) {x = v.x; y = v.y; z = v.z;}

	float& operator[](int i) {
		return E[i];
	}

};

#define v3_zero vec3(0,0,0)
#define v3_forward vec3(0,0,1)
#define v3_up vec3(0,1,0)
#define v3_right vec3(1,0,0)
#define v3_one vec3(1,1,1)

union vec4 {
	struct {
		float x, y, z, w;
	};

	struct {
		float r, g, b, a;
	};
};

union quaternion {
	struct {
		float x;
		float y;
		float z;
		float w;
	};

	struct {
		vec3 xyz;
		float w;
	};

	quaternion() {}
	quaternion(float _x, float _y, float _z, float _w){
		x = _x;
		y = _y;
		z = _z;
		w = _w;
	}
	quaternion(const quaternion& q) : x(q.x), y(q.y), z(q.z), w(q.w) {}
};


union mat3 {
	float E[9];
	float M[3][3];
	float& operator()(int x, int y) {
		return E[x+3*y];
	}
};

union mat4 {
	struct {
	   vec3 right;
	   float _1;
	   vec3 up;
	   float _2;
	   vec3 forward;
	   float _3;
	   vec3 position;
	   float _4;
	};

	float E[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float M[4][4];

	struct {
		float m00;
		float m10;
		float m20;
		float m30;
		float m01;
		float m11;
		float m21;
		float m31;
		float m02;
		float m12;
		float m22;
		float m32;
		float m03;
		float m13;
		float m23;
		float m33;
	};

	mat4(){}
	mat4(float _m00, float _m10, float _m20, float _m30, float _m01, float _m11, float _m21, float _m31,
		 float _m02, float _m12, float _m22, float _m32, float _m03, float _m13, float _m23, float _m33) :
		 m00(_m00), m10(_m10), m20(_m20), m30(_m30), m01(_m01), m11(_m11), m21(_m21), m31(_m31),
		 m02(_m02), m12(_m12), m22(_m22), m32(_m32), m03(_m03), m13(_m13), m23(_m23), m33(_m33) {}

	mat4(const mat4& m) {
		m00 = m.m00;
		m10 = m.m10;
		m20 = m.m20;
		m30 = m.m30;
		m01 = m.m01;
		m11 = m.m11;
		m21 = m.m21;
		m31 = m.m31;
		m02 = m.m02;
		m12 = m.m12;
		m22 = m.m22;
		m32 = m.m32;
		m03 = m.m03;
		m13 = m.m13;
		m23 = m.m23;
		m33 = m.m33;
	}
	mat4(vec3 _right, vec3 _up, vec3 _forward) {
		right = _right;
		up = _up;
		forward = _forward;
		_1 = 0;
		_2 = 0;
		_3 = 0;
		_4 = 1;
	}

	float& operator()(int i, int j)
	{
		return E[4*j + i];
	}


	quaternion to_quaternion();
};




bool operator==(const vec3& a, const vec3& b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
bool operator!=(const vec3& a, const vec3& b) { return !(a == b); }
vec3 operator*(float a, const vec3& b) { return vec3(a*b.x, a*b.y, a*b.z); }
vec3 operator*(const vec3& a, const vec3& b) { return vec3(a.x*b.x, a.y*b.y, a.z*b.z); }

vec3 operator*(const vec3& b, float a) { return a*b; }
vec3& operator*=(vec3 &b, float a) {
	b = b * a;
	return b;
}

vec3 operator/(const vec3& b, float a) { return vec3(b.x/a, b.y/a, b.z/a); }
vec3& operator/=(vec3 &b, float a) {
	b = b / a;
	return b;
}

vec3 operator+(const vec3& a, const vec3& b) { return vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
vec3 operator+=(vec3 &a, const vec3& b) {
	a = a + b;
	return a;
}

vec3 operator-(const vec3& a, const vec3& b) { return vec3(a.x-b.x, a.y-b.y, a.z-b.z); }
vec3 operator-(const vec3& a) { return vec3(-a.x, -a.y, -a.z); }

vec3 & operator-=(vec3 &a, const vec3& b) {
	a = a - b;
	return a;
}

mat3 operator*(float A, const mat3& B) {
	mat3 result;
	for (int i = 0; i < 9; ++i)
		result.E[i] = A*B.E[i];
	return result;
}
mat3 operator*(const mat3& B, float A) { return A*B; }

mat4 operator*(float A, const mat4& B) {
	mat4 result;
	for (int i = 0; i < 16; ++i)
		result.E[i] = A*B.E[i];
	return result;
}
mat4 operator*(const mat4& B, float A) { return A*B; }
mat4 operator*(const mat4& A, const mat4& B) {
	mat4 result;

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			for (int k = 0; k < 4; ++k)
				result.M[i][j] += A.M[k][j] * B.M[i][k];

	return result;
}

vec3 operator*(const mat4& a, const vec3& b) {
	vec3 result;

	result.x = a.M[0][0]*b.x + a.M[0][1]*b.y + a.M[0][2]*b.z + a.M[0][3];
	result.y = a.M[1][0]*b.x + a.M[1][1]*b.y + a.M[1][2]*b.z + a.M[1][3];
	result.z = a.M[2][0]*b.x + a.M[2][1]*b.y + a.M[2][2]*b.z + a.M[2][3];

	return result;
}

float dot(const vec3& a, const vec3& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

float length(const vec3& v) {
	return sqrt(dot(v,v));
}

float length(const quaternion& q) {
	return sqrt(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);
}

vec3 normalize(const vec3& v) {
	float k = 1.0 / length(v);

	return vec3(v.x * k, v.y * k, v.z * k);
}

quaternion normalize(const quaternion& q) {
	float k = 1.0f / length(q);

	return quaternion(q.x*k, q.y*k, q.z*k, q.w*k);
}

vec3 reflect(vec3& v, vec3& n) {
    return v - 2*dot(v,n)*n;
}

vec3 refract(vec3& in, vec3& normal, double refraction_ratio) {
    auto cos_theta = fmin(dot(-in, normal), 1.0);
    vec3 r_out_perp =  refraction_ratio * (in + cos_theta*normal);
    vec3 r_out_parallel = -sqrt(fabs(1.0 - dot(r_out_perp, r_out_perp))) * normal;
    return r_out_perp + r_out_parallel;
}

thread_local static uint32_t xor_shift_state = 1;
static uint32_t xor_shift_32() {
    uint32_t x = xor_shift_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    xor_shift_state = x;
    return x;
}

float random_float_01() { return (xor_shift_32() & 0xFFFFFF) / 16777216.0f; }
float random_float(float min, float max) {
	return min + (max-min) * random_float_01();
}
vec3 vec3_random() {
	return vec3(random_float_01(), random_float_01(), random_float_01());
}

vec3 vec3_random(float min, float max) {
	return vec3(random_float(min, max), random_float(min, max), random_float(min, max));

}

vec3 random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_float(-1,1), random_float(-1,1), 0);
        if (dot(p,p) >= 1) continue;
        return p;
    }
}

vec3 random_on_unit_sphere() {
 	auto phi = random_float(0, 2*M_PI);
	auto theta = random_float(0, 2 * M_PI);

	auto sin_phi = sin(phi);
	auto sin_theta = sin(theta);
	auto cos_phi = cos(phi);
	auto cos_theta = cos(theta);

	return vec3(sin_phi * cos_theta, sin_phi * sin_theta, cos_phi);
}

vec3 random_in_unit_sphere() {
	auto radius = random_float(-1,1);

	auto phi = random_float(0, 2*M_PI);
	auto theta = random_float(0, 2 * M_PI);

	auto sin_phi = sin(phi);
	auto sin_theta = sin(theta);
	auto cos_phi = cos(phi);
	auto cos_theta = cos(theta);

	return vec3(radius * sin_phi * cos_theta, radius * sin_phi * sin_theta, radius * cos_phi);
}

vec3 random_in_hemisphere(vec3& normal) {
    vec3 in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0)
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

quaternion make_axis_angle_rotation(const vec3& v, float theta) {
	quaternion result;
	result.xyz = normalize(v) * sin(theta/2);
	result.w = cos(theta/2);
	return result;
}

vec3 make_euler(const quaternion& q) {
	double ysqr = q.y * q.y;
	double t0 = -2.0f * (ysqr + q.z * q.z) + 1.0f;
	double t1 = +2.0f * (q.x * q.y + q.w * q.z);
	double t2 = -2.0f * (q.x * q.z - q.w * q.y);
	double t3 = +2.0f * (q.y * q.z + q.w * q.x);
	double t4 = -2.0f * (q.x * q.x + ysqr) + 1.0f;

	t2 = t2 > 1.0f ? 1.0f : t2;
	t2 = t2 < -1.0f ? -1.0f : t2;

	vec3 result;
	result.x = asin(t2);
	result.y = atan2(t3, t4);
	result.z = atan2(t1, t0);
	return result;
}

quaternion conjugate(const quaternion& q) {
	return quaternion(-q.x, -q.y, -q.z, q.w);
}

void swap(float& l, float& r) {
	float tmp = l;
	l = r;
	r = tmp;
}

// i  j  k
// ux uy uz
// vx vy vz
vec3 cross(const vec3& u, const vec3& v) {
	return vec3(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x);
}

void transpose_matrix(mat4& m) {
	swap(m.m01, m.m10);
	swap(m.m02, m.m20);
	swap(m.m03, m.m30);

	swap(m.m12, m.m21);
	swap(m.m13, m.m31);

	swap(m.m23, m.m32);
}

mat4 make_identity_matrix() {
	mat4 result = {};
	result.M[0][0] = result.M[1][1] = result.M[2][2] = result.M[3][3] = 1;
	return result;
}

mat4 make_perspective_matrix(float fov, float aspect, float n, float f) {
	mat4 result;

	result.M[0][0] = 1.0f/tan(fov / 2) / aspect; /*(fov/2)*2*M_PI*/
	result.M[1][1] = aspect * result.M[0][0];
	result.M[2][2] = (f+n)/(f-n);
	result.M[3][2] = -2*f*n/(f-n);
	result.M[2][3] = 1;

	return result;
}

mat4 make_scale_matrix(float scaleX, float scaleY, float scaleZ) {
	mat4 result;

	result(0,0) = scaleX;
	result(1,1) = scaleY;
	result(2,2) = scaleZ;
	result(3,3) = 1;

	return result;
}

mat4 make_scale_matrix(const vec3& v) {
	return make_scale_matrix(v.x, v.y, v.z);
}

mat4 make_scale_matrix(float s) {
	return make_scale_matrix(s,s,s);
}

mat4 make_translation(float x, float y, float z) {
	mat4 result = make_identity_matrix();

	result.M[3][0] = x;
	result.M[3][1] = y;
	result.M[3][2] = z;

	return result;
}

mat4 make_translation(const vec3& v)
{
	return make_translation(v.x, v.y, v.z);
}

mat4 make_rotation_matrix(const quaternion& q) {
	mat4 result;

	result.E[0] = 1 - 2*(q.y*q.y + q.z*q.z);
	result.E[1] = 2*(q.x*q.y - q.w*q.z);
	result.E[2] = 2*(q.x*q.z + q.w*q.y);
	result.E[4] = 2*(q.x*q.y + q.w*q.z);
	result.E[5] = 1 - 2*(q.x*q.x + q.z*q.z);
	result.E[6] = 2*(q.y*q.z - q.w*q.x);
	result.E[8] = 2*(q.x*q.z - q.w*q.y);
	result.E[9] = 2*(q.y*q.z + q.w*q.x);
	result.E[10] = 1 - 2*(q.x*q.x + q.y*q.y);
	result.E[15] = 1;

	transpose_matrix(result);

	return result;
}

quaternion mat4::to_quaternion() {
	quaternion result;

	float sum = m00 + m11 + m22;
	if (sum > 0.0F) {
		result.w = sqrt(sum + 1.0F) * 0.5F;
		float f = 0.25F / result.w;

		result.x = (m21 - m12) * f;
		result.y = (m02 - m20) * f;
		result.z = (m10 - m01) * f;
	} else if ((m00 > m11) && (m00 > m22)) {
		result.x = sqrt(m00 - m11 - m22 + 1.0F) * 0.5F;
		float f = 0.25F / result.x;

		result.y = (m10 + m01) * f;
		result.z = (m02 + m20) * f;
		result.w = (m21 - m12) * f;
	} else if (m11 > m22) {
		result.y = sqrt(m11 - m00 - m22 + 1.0F) * 0.5F;
		float f = 0.25F / result.y;

		result.x = (m10 + m01) * f;
		result.z = (m21 + m12) * f;
		result.w = (m02 - m20) * f;
	} else {
		result.z = sqrt(m22 - m00 - m11 + 1.0F) * 0.5F;
		float f = 0.25F / result.z;

		result.x = (m02 + m20) * f;
		result.y = (m21 + m12) * f;
		result.w = (m10 - m01) * f;
	}

	return result;
}

void decompose_matrix(mat4* m, vec3* position, quaternion* rotation, vec3* scale) {
	mat4 tmp;
	memcpy(&tmp, m, 16 * sizeof(float));

	scale->x = length(tmp.right);
	scale->y = length(tmp.up);
	scale->z = length(tmp.forward);

	tmp.right /= scale->x;
	tmp.up /= scale->y;
	tmp.forward /= scale->z;

	*position = tmp.position;
	*rotation = tmp.to_quaternion();
}

#define quat_identity (quaternion(0,0,0,1))

quaternion operator*(const quaternion& l, const quaternion& r) {
	quaternion result;
	vec3 lv = {l.x, l.y, l.z};
	vec3 rv = {r.x, r.y, r.z};
	vec3 lvrvc = cross(lv, rv);
	result.x = lvrvc.x + l.w*rv.x + r.w*lv.x;
	result.y = lvrvc.y + l.w*rv.y + r.w*lv.y;
	result.z = lvrvc.z + l.w*rv.z + r.w*lv.z;
	result.w = l.w*r.w - dot(lv, rv);
	return result;
}

vec3 operator*(const quaternion& q, const vec3& v) {
	vec3 t = 2 * cross(q.xyz, v);
	return v + q.w * t + cross(q.xyz, t);
}

quaternion from_to_rotation(const vec3& from, const vec3& to) {
	quaternion q;
	vec3 a = cross(from, to);
	q.x = a.x;
	q.y = a.y;
	q.z = a.z;
	q.w = sqrt(dot(from,from) * dot(to,to)) + dot(from, to);
	q = normalize(q);

	return  q;
}

mat4 look_at(vec3 forward, vec3 up) {
	auto right = normalize(cross(forward, up));
	return mat4(right, up, forward);
}

quaternion look_rotation(vec3 forward, vec3 up) {
	if (forward == v3_zero) {
		return quat_identity;
	}

	if (up != forward) {
		up = normalize(up);
		vec3 v = forward - up * dot(up, forward);

		return from_to_rotation(v, forward) * from_to_rotation(v3_forward, v);
	}
	else {
		return from_to_rotation(v3_forward, forward);
	}
}

