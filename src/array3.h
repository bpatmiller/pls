#pragma once

#include <cstring> // used for memset
#include <glm/glm.hpp>

template <class S, class T>
inline S lerp(const S &value0, const S &value1, T f) {
  return (1 - f) * value0 + f * value1;
}

template <class S, class T>
inline S bilerp(const S &v00, const S &v10, const S &v01, const S &v11, T fx,
                T fy) {
  return lerp(lerp(v00, v10, fx), lerp(v01, v11, fx), fy);
}

template <class T> struct Array3 {
  int sx, sy, sz;
  int size;
  T *data;

  Array3() {
    sx = 0;
    sy = 0;
    sz = 0;
    size = 0;
  }

  Array3(int sx_, int sy_, int sz_)
      : sx(sx_), sy(sy_), sz(sz_), size(0), data(0) {
    init();
  }

  ~Array3() {
    delete[] data;
    data = 0;
    sx = 0;
    sy = 0;
    sz = 0;
    size = 0;
  }

  void init() {
    size = sx * sy * sz;
    data = new T[size];
    clear();
  }

  void init(int sx_, int sy_, int sz_) {
    sx = sx_;
    sy = sy_;
    sz = sz_;
    size = sx * sy * sz;
    data = new T[size];
    clear();
  }

  void set(T val) {
    for (int i = 0; i < size; i++) {
      data[i] = val;
    }
  }

  // returns largest absolute value in the array
  T infnorm() const {
    T n = 0;
    for (int i = 0; i < size; i++) {
      if (std::fabs(data[i]) > n)
        n = std::fabs(data[i]);
    }
    return n;
  }

  double dot(Array3 &other) {
    double d = 0;
    for (int i = 0; i < size; i++)
      d += data[i] * other.data[i];
    return d;
  }

  void inc(double s, Array3 &other) {
    for (int i = 0; i < size; i++) {
      data[i] += s * other.data[i];
    }
  }

  void scale(float s) {
    for (int i = 0; i < size; i++) {
      data[i] *= s;
    }
  }

  void scale_inc(double s, Array3 &other) {
    for (int i = 0; i < size; i++) {
      data[i] = s * data[i] + other.data[i];
    }
  }

  void copy_to(Array3 &other) const {
    std::memcpy(other.data, data, other.size * sizeof(T));
  }

  // index (i,j,k cell index), coords (0-1f uv coords)
  T trilerp(glm::ivec3 index, glm::vec3 coords) {
    return (1 - coords.x) * (1 - coords.y) * (1 - coords.z) *
               (*this)(index.x, index.y, index.z) +
           (1 - coords.x) * (1 - coords.y) * (coords.z) *
               (*this)(index.x, index.y, index.z + 1) +
           (1 - coords.x) * (coords.y) * (1 - coords.z) *
               (*this)(index.x, index.y + 1, index.z) +
           (1 - coords.x) * (coords.y) * (coords.z) *
               (*this)(index.x, index.y + 1, index.z + 1) +
           (coords.x) * (1 - coords.y) * (1 - coords.z) *
               (*this)(index.x + 1, index.y, index.z) +
           (coords.x) * (1 - coords.y) * (coords.z) *
               (*this)(index.x + 1, index.y, index.z + 1) +
           (coords.x) * (coords.y) * (1 - coords.z) *
               (*this)(index.x + 1, index.y + 1, index.z) +
           (coords.x) * (coords.y) * (coords.z) *
               (*this)(index.x + 1, index.y + 1, index.z + 1);
  }

  glm::vec3 gradlerp(glm::ivec3 index, glm::vec3 coords) {
    int i = index.x;
    int j = index.y;
    int k = index.z;

    T v000 = (*this)(i, j, k);
    T v001 = (*this)(i, j, k + 1);
    T v010 = (*this)(i, j + 1, k);
    T v011 = (*this)(i, j + 1, k + 1);
    T v100 = (*this)(i + 1, j, k);
    T v101 = (*this)(i + 1, j, k + 1);
    T v110 = (*this)(i + 1, j + 1, k);
    T v111 = (*this)(i + 1, j + 1, k + 1);

    T ddx00 = (v100 - v000);
    T ddx10 = (v110 - v010);
    T ddx01 = (v101 - v001);
    T ddx11 = (v111 - v011);
    T dv_dx = bilerp(ddx00, ddx10, ddx01, ddx11, coords.y, coords.z);

    T ddy00 = (v010 - v000);
    T ddy10 = (v110 - v100);
    T ddy01 = (v011 - v001);
    T ddy11 = (v111 - v101);
    T dv_dy = bilerp(ddy00, ddy10, ddy01, ddy11, coords.x, coords.z);

    T ddz00 = (v001 - v000);
    T ddz10 = (v101 - v100);
    T ddz01 = (v011 - v010);
    T ddz11 = (v111 - v110);
    T dv_dz = bilerp(ddz00, ddz10, ddz01, ddz11, coords.x, coords.y);

    return glm::vec3(dv_dx, dv_dy, dv_dz);
  }

  void clear() { std::memset(data, 0, size * sizeof(T)); }

  T &operator()(int i, int j, int k) {
    return data[i + (sx * j) + (sx * sy * k)];
  }

  T &operator()(glm::ivec3 index) {
    return data[index.x + (sx * index.y) + (sx * sy * index.z)];
  }
};

typedef Array3<double> Array3d;
typedef Array3<float> Array3f;
typedef Array3<int> Array3i;
typedef Array3<glm::vec3> Array3v3;
typedef Array3<glm::dvec4> Array3v4;