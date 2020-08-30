#pragma once

namespace __details__ {
template <typename scalar_t>
struct vec2 {
	scalar_t x, y;

	vec2() {}
	explicit vec2(scalar_t s) : x(s), y(s) {}
	explicit vec2(scalar_t x, scalar_t y) : x(x), y(y) {}
	vec2(vec2 const& v) : x(v.x), y(v.y) {}

	vec2& operator+=(scalar_t s) {
		x += s;
		y += s;
		return *this;
	}

	vec2& operator-=(scalar_t s) {
		x -= s;
		y -= s;
		return *this;
	}

	vec2& operator*=(scalar_t s) {
		x *= s;
		y *= s;
		return *this;
	}

	vec2& operator/=(scalar_t s) {
		scalar_t const i = 1 / s;
		x *= i;
		y *= i;
		return *this;
	}

	vec2& operator+=(vec2 const& v) {
		x += v.x;
		y += v.y;
		return *this;
	}

	vec2& operator-=(vec2 const& v) {
		x -= v.x;
		y -= v.y;
		return *this;
	}

	vec2& operator*=(vec2 const& v) {
		x *= v.x;
		y *= v.y;
		return *this;
	}

	vec2& operator/=(vec2 const& v) {
		x /= v.x;
		y /= v.y;
		return *this;
	}
};

template <typename scalar_t>
vec2<scalar_t> operator+(vec2<scalar_t> const& v, scalar_t s) {
	return vec2<scalar_t>(v.x + s, v.y + s);
}

template <typename scalar_t>
vec2<scalar_t> operator-(vec2<scalar_t> const& v, scalar_t s) {
	return vec2<scalar_t>(v.x - s, v.y - s);
}

template <typename scalar_t>
vec2<scalar_t> operator*(vec2<scalar_t> const& v, scalar_t s) {
	return vec2<scalar_t>(v.x * s, v.y * s);
}

template <typename scalar_t>
vec2<scalar_t> operator/(vec2<scalar_t> const& v, scalar_t s) {
	scalar_t const i = 1 / s;
	return vec2<scalar_t>(v.x * i, v.y * i);
}

template <typename scalar_t>
vec2<scalar_t> operator+(vec2<scalar_t> const& a, vec2<scalar_t> const& b) {
	return vec2<scalar_t>(a.x + b.x, a.y + b.y);
}

template <typename scalar_t>
vec2<scalar_t> operator-(vec2<scalar_t> const& a, vec2<scalar_t> const& b) {
	return vec2<scalar_t>(a.x - b.x, a.y - b.y);
}

template <typename scalar_t>
vec2<scalar_t> operator*(vec2<scalar_t> const& a, vec2<scalar_t> const& b) {
	return vec2<scalar_t>(a.x * b.x, a.y * b.y);
}

template <typename scalar_t>
vec2<scalar_t> operator/(vec2<scalar_t> const& a, vec2<scalar_t> const& b) {
	return vec2<scalar_t>(a.x / b.x, a.y / b.y);
}

template <typename scalar_t>
vec2<scalar_t> inverse(vec2<scalar_t> const& v) {
	return vec2<scalar_t>(-v.x, -v.y);
}

template <typename scalar_t>
vec2<scalar_t> min(vec2<scalar_t> const& a, vec2<scalar_t> const& b) {
	return vec2<scalar_t>(min(a.x, b.x), min(a.y, b.y));
}

template <typename scalar_t>
vec2<scalar_t> max(vec2<scalar_t> const& a, vec2<scalar_t> const& b) {
	return vec2<scalar_t>(max(a.x, b.x), max(a.y, b.y));
}

template <typename scalar_t>
vec2<scalar_t> abs(vec2<scalar_t> const& v) {
	return vec2<scalar_t>(abs(v.x), abs(v.y));
}

template <typename scalar_t>
scalar_t dot(vec2<scalar_t> const& a, vec2<scalar_t> const& b) {
	return a.x * b.x + a.y * b.y;
}

template <typename scalar_t>
scalar_t length(vec2<scalar_t> const& v) {
	return sqrt(dot(v, v));
}

template <typename scalar_t>
scalar_t distance(vec2<scalar_t> const& a, vec2<scalar_t> const& b) {
	return length(a - b);
}

template <typename scalar_t>
vec2<scalar_t> normalize(vec2<scalar_t> const& v) {
	return v / length(v);
}

template <typename scalar_t>
vec2<scalar_t> mix(vec2<scalar_t> const& a, vec2<scalar_t> const& b,
									 scalar_t t) {
	return vec2<scalar_t>(mix(a.x, b.x, t), mix(a.y, b.y, t));
}

template <typename scalar_t>
vec2<scalar_t> reflect(vec2<scalar_t> const& i, vec2<scalar_t> const& n) {
	return i - n * dot(n, i) * 2;
}

template <typename scalar_t>
vec2<scalar_t> refract(vec2<scalar_t> const& i, vec2<scalar_t> const& n,
											 scalar_t eta) {
	scalar_t const dni = dot(n, i);
	scalar_t const k = 1 - eta * eta * (1 - dni * dni);
	return k < 0 ? vec2<scalar_t>(0) : (i * eta - n * (eta * dni + sqrt(k)));
}

template <typename scalar_t>
struct vec3 {
	scalar_t x, y, z;

	vec3() {}
	explicit vec3(scalar_t s) : x(s), y(s), z(s) {}
	explicit vec3(scalar_t x, scalar_t y, scalar_t z) : x(x), y(y), z(z) {}
	explicit vec3(vec2<scalar_t> const& v, scalar_t z) : x(v.x), y(v.y), z(z) {}
	vec3(vec3 const& v) : x(v.x), y(v.y), z(v.z) {}

	vec2<scalar_t> to_vec2() const { return vec2<scalar_t>(x, y); }

	vec3& operator+=(scalar_t s) {
		x += s;
		y += s;
		z += s;
		return *this;
	}

	vec3& operator-=(scalar_t s) {
		x -= s;
		y -= s;
		z -= s;
		return *this;
	}

	vec3& operator*=(scalar_t s) {
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}

	vec3& operator/=(scalar_t s) {
		scalar_t const i = 1 / s;
		x *= i;
		y *= i;
		z *= i;
		return *this;
	}

	vec3& operator+=(vec3 const& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	vec3& operator-=(vec3 const& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	vec3& operator*=(vec3 const& v) {
		x *= v.x;
		y *= v.y;
		z *= v.z;
		return *this;
	}

	vec3& operator/=(vec3 const& v) {
		x /= v.x;
		y /= v.y;
		z /= v.z;
		return *this;
	}
};

template <typename scalar_t>
vec3<scalar_t> operator+(vec3<scalar_t> const& v, scalar_t s) {
	return vec3<scalar_t>(v.x + s, v.y + s, v.z + s);
}

template <typename scalar_t>
vec3<scalar_t> operator-(vec3<scalar_t> const& v, scalar_t s) {
	return vec3<scalar_t>(v.x - s, v.y - s, v.z - s);
}

template <typename scalar_t>
vec3<scalar_t> operator*(vec3<scalar_t> const& v, scalar_t s) {
	return vec3<scalar_t>(v.x * s, v.y * s, v.z * s);
}

template <typename scalar_t>
vec3<scalar_t> operator/(vec3<scalar_t> const& v, scalar_t s) {
	scalar_t const i = 1 / s;
	return vec3<scalar_t>(v.x * i, v.y * i, v.z * i);
}

template <typename scalar_t>
vec3<scalar_t> operator+(vec3<scalar_t> const& a, vec3<scalar_t> const& b) {
	return vec3<scalar_t>(a.x + b.x, a.y + b.y, a.z + b.z);
}

template <typename scalar_t>
vec3<scalar_t> operator-(vec3<scalar_t> const& a, vec3<scalar_t> const& b) {
	return vec3<scalar_t>(a.x - b.x, a.y - b.y, a.z - b.z);
}

template <typename scalar_t>
vec3<scalar_t> operator*(vec3<scalar_t> const& a, vec3<scalar_t> const& b) {
	return vec3<scalar_t>(a.x * b.x, a.y * b.y, a.z * b.z);
}

template <typename scalar_t>
vec3<scalar_t> operator/(vec3<scalar_t> const& a, vec3<scalar_t> const& b) {
	return vec3<scalar_t>(a.x / b.x, a.y / b.y, a.z / b.z);
}

template <typename scalar_t>
vec3<scalar_t> inverse(vec3<scalar_t> const& v) {
	return vec3<scalar_t>(-v.x, -v.y, -v.z);
}

template <typename scalar_t>
vec3<scalar_t> min(vec3<scalar_t> const& a, vec3<scalar_t> const& b) {
	return vec3<scalar_t>(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
}

template <typename scalar_t>
vec3<scalar_t> max(vec3<scalar_t> const& a, vec3<scalar_t> const& b) {
	return vec3<scalar_t>(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z));
}

template <typename scalar_t>
vec3<scalar_t> abs(vec3<scalar_t> const& v) {
	return vec3<scalar_t>(abs(v.x), abs(v.y), abs(v.z));
}

template <typename scalar_t>
scalar_t dot(vec3<scalar_t> const& a, vec3<scalar_t> const& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

template <typename scalar_t>
scalar_t length(vec3<scalar_t> const& v) {
	return sqrt(dot(v, v));
}

template <typename scalar_t>
scalar_t distance(vec3<scalar_t> const& a, vec3<scalar_t> const& b) {
	return length(a - b);
}

template <typename scalar_t>
vec3<scalar_t> normalize(vec3<scalar_t> const& v) {
	return v / length(v);
}

template <typename scalar_t>
vec3<scalar_t> cross(vec3<scalar_t> const& a, vec3<scalar_t> const& b) {
	return vec3<scalar_t>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
												a.x * b.y - a.y * b.x);
}

template <typename scalar_t>
vec3<scalar_t> mix(vec3<scalar_t> const& a, vec3<scalar_t> const& b,
									 scalar_t t) {
	return vec3<scalar_t>(mix(a.x, b.x, t), mix(a.y, b.y, t), mix(a.z, b.z, t));
}

template <typename scalar_t>
vec3<scalar_t> reflect(vec3<scalar_t> const& i, vec3<scalar_t> const& n) {
	return i - n * dot(n, i) * 2;
}

template <typename scalar_t>
vec3<scalar_t> refract(vec3<scalar_t> const& i, vec3<scalar_t> const& n,
											 scalar_t eta) {
	scalar_t const dni = dot(n, i);
	scalar_t const k = 1 - eta * eta * (1 - dni * dni);
	return k < 0 ? vec3<scalar_t>(0) : (i * eta - n * (eta * dni + sqrt(k)));
}

}  // namespace __details__
