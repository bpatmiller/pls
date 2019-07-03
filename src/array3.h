#pragma once

#include <cstring> // used for memset
#include <glm/glm.hpp>
#include <iostream>
#include <vector>

// lerp and bilerp inlines

template <class S, class T>
inline S lerp(const S &value0, const S &value1, T f) {
  return (1 - f) * value0 + f * value1;
}

template <class S, class T>
inline S bilerp(const S &v00, const S &v10, const S &v01, const S &v11, T fx,
                T fy) {
  return lerp(lerp(v00, v10, fx), lerp(v01, v11, fx), fy);
}

// class definition

template <class T> struct Array3 {
  int sx = 0;
  int sy = 0;
  int sz = 0;
  int size = 0;
  std::vector<T> data;

  Array3() {}

  Array3(int sx_, int sy_, int sz_) : sx(sx_), sy(sy_), sz(sz_), size(0) {
    init();
  }

  void init() {
    size = sx * sy * sz;
    data = std::vector<T>(size);
    clear();
  }

  void init(int sx_, int sy_, int sz_) {
    sx = sx_;
    sy = sy_;
    sz = sz_;
    init();
  }

  void clear() { set(static_cast<T>(0)); }

  void set(T val) { std::fill(data.begin(), data.end(), val); }

  // returns largest absolute value in the array
  T infnorm() const {
    T n = 0;
    for (int i = 0; i < size; i++) {
      if (std::fabs(data[i]) > n)
        n = std::fabs(data[i]);
    }
    return n;
  }

  T mean() const {
    T n = 0;
    for (int i = 0; i < size; i++) {
      n += data[i];
    }
    return n / (T)size;
  }

  void copy_to(Array3 &rhs) const {
    assert(size == rhs.size);
    rhs.data.assign(data.begin(), data.end());
  }

  void print() {
    for (int k = 0; k < sz; k++) {
      std::printf("z = %i\n", k);
      for (int j = 0; j < sy; j++) {
        for (int i = 0; i < sx; i++) {
          std::cout << (*this)(i, j, k) << "  ";
        }
        std::cout << "\n\n";
      }
    }
  }

  T &operator()(int i, int j, int k) {
    assert(i >= 0 and i <= sx);
    assert(j >= 0 and j <= sy);
    assert(k >= 0 and k <= sz);
    return data[i + (sx * j) + (sx * sy * k)];
  }

  T &operator()(glm::ivec3 index) { return (*this)(index.x, index.y, index.z); }

  // copy constructor
  Array3<T>(const Array3<T> &rhs) {
    // require they be of the same size if already initialized
    assert(size == rhs.size);

    for (int i = 0; i < size; i++) {
      data[i] = rhs.data[i];
    }
  }

  // ARRAY3<T> OPERATORS

  Array3<T> operator+(const Array3<T> &rhs) {
    assert(size == rhs.size);
    Array3<T> n;
    n.init(sx, sy, sz);
    for (int i = 0; i < size; i++) {
      n.data[i] = data[i] + rhs.data[i];
    }
    return n;
  }
  Array3<T> operator-(const Array3<T> &rhs) {
    assert(size == rhs.size);
    Array3<T> n;
    n.init(sx, sy, sz);
    for (int i = 0; i < size; i++) {
      n.data[i] = data[i] - rhs.data[i];
    }
    return n;
  }

  // component-wise multiplication
  Array3<T> operator*(const Array3<T> &rhs) {
    assert(size == rhs.size);
    Array3<T> n;
    n.init(sx, sy, sz);
    for (int i = 0; i < size; i++) {
      n.data[i] = data[i] * rhs.data[i];
    }
    return n;
  }

  Array3<T> operator/(const Array3<T> &rhs) {
    assert(size == rhs.size);
    Array3<T> n;
    n.init(sx, sy, sz);
    for (int i = 0; i < size; i++) {
      n.data[i] = data[i] / rhs.data[i];
    }
    return n;
  }
  // T OPERATORS
  Array3<T> operator+(const T rhs) {
    Array3<T> n;
    n.init(sx, sy, sz);
    for (int i = 0; i < size; i++) {
      n.data[i] = data[i] + rhs;
    }
    return n;
  }

  Array3<T> operator-(const T rhs) {
    Array3<T> n;
    n.init(sx, sy, sz);
    for (int i = 0; i < size; i++) {
      n.data[i] = data[i] - rhs;
    }
    return n;
  }
  Array3<T> operator*(const T rhs) {
    Array3<T> n;
    n.init(sx, sy, sz);
    for (int i = 0; i < size; i++) {
      n.data[i] = data[i] * rhs;
    }
    return n;
  }

  Array3<T> operator/(const T rhs) {
    Array3<T> n;
    n.init(sx, sy, sz);
    for (int i = 0; i < size; i++) {
      n.data[i] = data[i] / rhs;
    }
    return n;
  }

  // += and -=
  Array3<T> &operator+=(const T &rhs) {
    for (int i = 0; i < size; i++) {
      data[i] += rhs;
    }
    return *this;
  }

  Array3<T> &operator-=(const T &rhs) {
    for (int i = 0; i < size; i++) {
      data[i] -= rhs;
    }
    return *this;
  }
  // MATH FUNCTIONS

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

  // first-order upwind differencing [Osher and Fedkiw 2003]
  // for example, use forward differencing on dx when v.x < 0
  // since the flow is right-to-left
  glm::vec3 upwind_gradient(int i, int j, int k, float h, glm::vec3 v) {
    glm::vec3 grad(0.0f, 0.0f, 0.0f);
    float c = (*this)(i, j, k);
    // x
    if (i == 0) {
      grad.x = ((*this)(i + 1, j, k) - c);
    } else if (i == sx) {
      grad.x = (c - (*this)(i - 1, j, k));
    } else if (v.x < 0) {
      grad.x = ((*this)(i + 1, j, k) - c);
    } else {
      grad.x = (c - (*this)(i - 1, j, k));
    }
    // y
    if (j == 0) {
      grad.y = ((*this)(i, j + 1, k) - c);
    } else if (j == sy) {
      grad.y = (c - (*this)(i, j - 1, k));
    } else if (v.y < 0) {
      grad.y = ((*this)(i, j + 1, k) - c);
    } else {
      grad.y = (c - (*this)(i, j - 1, k));
    }
    // z
    if (k == 0) {
      grad.z = ((*this)(i, j, k + 1) - c);
    } else if (k == sz) {
      grad.z = (c - (*this)(i, j, k - 1));
    } else if (v.z < 0) {
      grad.z = ((*this)(i, j, k + 1) - c);
    } else {
      grad.z = (c - (*this)(i, j, k - 1));
    }

    return grad / h;
  }

  // central difference gradient at a given grid node
  glm::vec3 central_grad(int i, int j, int k, float h) {
    float d_dx = ((*this)(i + 1, j, k) - (*this)(i - 1, j, k));
    float d_dy = ((*this)(i, j + 1, k) - (*this)(i, j - 1, k));
    float d_dz = ((*this)(i, j, k + 1) - (*this)(i, j, k - 1));

    return glm::vec3(d_dx, d_dy, d_dz) / (2.0f * h);
  }

  // interpolated forward difference gradient at a position in space
  glm::vec3 gradlerp(glm::ivec3 index, glm::vec3 coords, float h) {
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

    return glm::vec3(dv_dx, dv_dy, dv_dz) / h;
  }
};

typedef Array3<double> Array3d;
typedef Array3<float> Array3f;
typedef Array3<int> Array3i;
typedef Array3<glm::vec3> Array3v3;
typedef Array3<glm::dvec4> Array3v4;