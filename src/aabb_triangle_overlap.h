/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Möller                              */
/* Function: int triBoxOverlap(float boxcenter[3],      */
/*          float boxhalfsize[3],float triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/

#ifndef _AABB_TRIANGLE_OVERLAP_H_
#define _AABB_TRIANGLE_OVERLAP_H_

#include <math.h>
#include <stdio.h>

#define XA 0
#define YA 1
#define ZA 2

#define CROSS(dest, v1, v2)                                                    \
  dest[0] = v1[1] * v2[2] - v1[2] * v2[1];                                     \
  dest[1] = v1[2] * v2[0] - v1[0] * v2[2];                                     \
  dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2)                                                      \
  dest[0] = v1[0] - v2[0];                                                     \
  dest[1] = v1[1] - v2[1];                                                     \
  dest[2] = v1[2] - v2[2];

#define FINDMINMAX(x0, x1, x2, min, max)                                       \
  min = max = x0;                                                              \
  if (x1 < min)                                                                \
    min = x1;                                                                  \
  if (x1 > max)                                                                \
    max = x1;                                                                  \
  if (x2 < min)                                                                \
    min = x2;                                                                  \
  if (x2 > max)                                                                \
    max = x2;

inline int planeBoxOverlap(float normal[3], float vert[3],
                           float maxbox[3]) // -NJMP-
{
  int q;
  float vmin[3], vmax[3], v;
  for (q = XA; q <= ZA; q++) {
    v = vert[q]; // -NJMP-
    if (normal[q] > 0.0f) {
      vmin[q] = -maxbox[q] - v; // -NJMP-
      vmax[q] = maxbox[q] - v;  // -NJMP-
    } else {
      vmin[q] = maxbox[q] - v;  // -NJMP-
      vmax[q] = -maxbox[q] - v; // -NJMP-
    }
  }
  if (DOT(normal, vmin) > 0.0f)
    return 0; // -NJMP-
  if (DOT(normal, vmax) >= 0.0f)
    return 1; // -NJMP-

  return 0;
}

/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)                                             \
  p0 = a * v0[YA] - b * v0[ZA];                                                \
  p2 = a * v2[YA] - b * v2[ZA];                                                \
  if (p0 < p2) {                                                               \
    min = p0;                                                                  \
    max = p2;                                                                  \
  } else {                                                                     \
    min = p2;                                                                  \
    max = p0;                                                                  \
  }                                                                            \
  rad = fa * boxhalfsize[YA] + fb * boxhalfsize[ZA];                           \
  if (min > rad || max < -rad)                                                 \
    return 0;

#define AXISTEST_X2(a, b, fa, fb)                                              \
  p0 = a * v0[YA] - b * v0[ZA];                                                \
  p1 = a * v1[YA] - b * v1[ZA];                                                \
  if (p0 < p1) {                                                               \
    min = p0;                                                                  \
    max = p1;                                                                  \
  } else {                                                                     \
    min = p1;                                                                  \
    max = p0;                                                                  \
  }                                                                            \
  rad = fa * boxhalfsize[YA] + fb * boxhalfsize[ZA];                           \
  if (min > rad || max < -rad)                                                 \
    return 0;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)                                             \
  p0 = -a * v0[XA] + b * v0[ZA];                                               \
  p2 = -a * v2[XA] + b * v2[ZA];                                               \
  if (p0 < p2) {                                                               \
    min = p0;                                                                  \
    max = p2;                                                                  \
  } else {                                                                     \
    min = p2;                                                                  \
    max = p0;                                                                  \
  }                                                                            \
  rad = fa * boxhalfsize[XA] + fb * boxhalfsize[ZA];                           \
  if (min > rad || max < -rad)                                                 \
    return 0;

#define AXISTEST_Y1(a, b, fa, fb)                                              \
  p0 = -a * v0[XA] + b * v0[ZA];                                               \
  p1 = -a * v1[XA] + b * v1[ZA];                                               \
  if (p0 < p1) {                                                               \
    min = p0;                                                                  \
    max = p1;                                                                  \
  } else {                                                                     \
    min = p1;                                                                  \
    max = p0;                                                                  \
  }                                                                            \
  rad = fa * boxhalfsize[XA] + fb * boxhalfsize[ZA];                           \
  if (min > rad || max < -rad)                                                 \
    return 0;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)                                             \
  p1 = a * v1[XA] - b * v1[YA];                                                \
  p2 = a * v2[XA] - b * v2[YA];                                                \
  if (p2 < p1) {                                                               \
    min = p2;                                                                  \
    max = p1;                                                                  \
  } else {                                                                     \
    min = p1;                                                                  \
    max = p2;                                                                  \
  }                                                                            \
  rad = fa * boxhalfsize[XA] + fb * boxhalfsize[YA];                           \
  if (min > rad || max < -rad)                                                 \
    return 0;

#define AXISTEST_Z0(a, b, fa, fb)                                              \
  p0 = a * v0[XA] - b * v0[YA];                                                \
  p1 = a * v1[XA] - b * v1[YA];                                                \
  if (p0 < p1) {                                                               \
    min = p0;                                                                  \
    max = p1;                                                                  \
  } else {                                                                     \
    min = p1;                                                                  \
    max = p0;                                                                  \
  }                                                                            \
  rad = fa * boxhalfsize[XA] + fb * boxhalfsize[YA];                           \
  if (min > rad || max < -rad)                                                 \
    return 0;

inline int triBoxOverlap(float boxcenter[3], float boxhalfsize[3],
                         float triverts[3][3]) {

  /*    use separating axis theorem to test overlap between triangle and box */
  /*    need to test for overlap in these directions: */
  /*    1) the {x,y,z}-directions (actually, since we use the AABB of the
   * triangle */
  /*       we do not even need to test these) */
  /*    2) normal of the triangle */
  /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
  /*       this gives 3x3=9 more tests */
  float v0[3], v1[3], v2[3];
  //   float axis[3];
  float min, max, p0, p1, p2, rad, fex, fey,
      fez; // -NJMP- "d" local variable removed
  float normal[3], e0[3], e1[3], e2[3];

  /* This is the fastest branch on Sun */
  /* move everything so that the boxcenter is in (0,0,0) */
  SUB(v0, triverts[0], boxcenter);
  SUB(v1, triverts[1], boxcenter);
  SUB(v2, triverts[2], boxcenter);

  /* compute triangle edges */
  SUB(e0, v1, v0); /* tri edge 0 */
  SUB(e1, v2, v1); /* tri edge 1 */
  SUB(e2, v0, v2); /* tri edge 2 */

  /* Bullet 3:  */
  /*  test the 9 tests first (this was faster) */
  fex = fabsf(e0[XA]);
  fey = fabsf(e0[YA]);
  fez = fabsf(e0[ZA]);
  AXISTEST_X01(e0[ZA], e0[YA], fez, fey);
  AXISTEST_Y02(e0[ZA], e0[XA], fez, fex);
  AXISTEST_Z12(e0[YA], e0[XA], fey, fex);

  fex = fabsf(e1[XA]);
  fey = fabsf(e1[YA]);
  fez = fabsf(e1[ZA]);
  AXISTEST_X01(e1[ZA], e1[YA], fez, fey);
  AXISTEST_Y02(e1[ZA], e1[XA], fez, fex);
  AXISTEST_Z0(e1[YA], e1[XA], fey, fex);

  fex = fabsf(e2[XA]);
  fey = fabsf(e2[YA]);
  fez = fabsf(e2[ZA]);
  AXISTEST_X2(e2[ZA], e2[YA], fez, fey);
  AXISTEST_Y1(e2[ZA], e2[XA], fez, fex);
  AXISTEST_Z12(e2[YA], e2[XA], fey, fex);

  /* Bullet 1: */
  /*  first test overlap in the {x,y,z}-directions */
  /*  find min, max of the triangle each direction, and test for overlap in */
  /*  that direction -- this is equivalent to testing a minimal AABB around */
  /*  the triangle against the AABB */

  /* test in X-direction */
  FINDMINMAX(v0[XA], v1[XA], v2[XA], min, max);
  if (min > boxhalfsize[XA] || max < -boxhalfsize[XA])
    return 0;

  /* test in Y-direction */
  FINDMINMAX(v0[YA], v1[YA], v2[YA], min, max);
  if (min > boxhalfsize[YA] || max < -boxhalfsize[YA])
    return 0;

  /* test in Z-direction */
  FINDMINMAX(v0[ZA], v1[ZA], v2[ZA], min, max);
  if (min > boxhalfsize[ZA] || max < -boxhalfsize[ZA])
    return 0;

  /* Bullet 2: */
  /*  test if the box intersects the plane of the triangle */
  /*  compute plane equation of triangle: normal*x+d=0 */
  CROSS(normal, e0, e1);
  // -NJMP- (line removed here)
  if (!planeBoxOverlap(normal, v0, boxhalfsize))
    return 0; // -NJMP-

  return 1; /* box and triangle overlaps */
}

#endif