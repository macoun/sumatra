#ifndef vec3_h
#define vec3_h

#include <math.h>

#define SMALL 1E-4

typedef struct vec3 vec3;
struct vec3 { float x, y, z; };

static inline vec3 vec3_mul(vec3 v, float s) {
  return (vec3){v.x*s, v.y*s, v.z*s};
}

static inline vec3 vec3_div(vec3 v, float s) {
  return (vec3){v.x/s, v.y/s, v.z/s};
}

static inline vec3 vec3_sub(vec3 a, vec3 b) {
  return (vec3){a.x-b.x, a.y-b.y, a.z-b.z};
}

static inline vec3 vec3_add(vec3 a, vec3 b) {
  return (vec3){a.x+b.x, a.y+b.y, a.z+b.z};
}

static inline float vec3_dot(vec3 a, vec3 b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline float vec3_mag2(vec3 v) {
  return v.x*v.x + v.y*v.y + v.z*v.z;
}

static inline float vec3_mag(vec3 v) {
  return sqrt(vec3_mag2(v));
}

static inline float vec3_dist(vec3 a, vec3 b) {
  return vec3_mag(vec3_sub(a, b));
}

static inline float vec3_dist2(vec3 a, vec3 b) {
  return vec3_mag2(vec3_sub(a, b));
}

static inline vec3 vec3_norm(vec3 v) {
  float mag = vec3_mag(v);
//  if (mag < 0.000001) mag = 0.00001;
  return (vec3){v.x/mag, v.y/mag, v.z/mag};
}

static inline vec3 vec3_cross(vec3 a, vec3 b) {
  return (vec3){
    a.y*b.z - a.z*b.y,
    a.z*b.x - a.x*b.z,
    a.x*b.y - a.y*b.x
  };
}

static inline vec3 vec3_parallel(vec3 v, vec3 axis) {
  float scale = vec3_dot(v, axis) / vec3_mag2(axis);
  return vec3_mul(axis, scale);
}

static inline vec3 vec3_perpendicular(vec3 v, vec3 axis) {
  return vec3_sub(v, vec3_parallel(v, axis));
}

//// Dihedral angle calculations

/// utility
static inline float match_sign(float a, float b) {
  return (b >= 0.) ? fabsf(a) : -fabsf(a);
}

// Angle between two unit vectors
static inline float vec3_ang_unit(vec3 u, vec3 v) {
  float dot = vec3_dot(u, v);
  float arg = fminf(fabsf(dot),1.0e0);
  if (dot < 0.)
    arg *= -1.;
  
  return acosf(arg);
}

// Angle between two vectors (p1-p2, p2-p3)
static inline float vec3_pos_ang(vec3 p1, vec3 p2, vec3 p3)
{
  vec3 v21 = vec3_norm(vec3_sub(p1, p2));
  vec3 v23 = vec3_norm(vec3_sub(p3, p2));
  
  return vec3_ang_unit(v21, v23);
}

// Angle between two vectors (a, b)
static inline float vec3_ang(vec3 a, vec3 b) {
  float amag = vec3_mag(a);
  float bmag = vec3_mag(b);
  
  if (amag * bmag < SMALL) return 0.0;
  
  float c = vec3_dot(a, b) / amag / bmag;
  
  if (c >=  1.0) return 0.0;
  if (c <= -1.0) return (float)M_PI;
  
  return acosf(c);
}

// dihedral angle between a and c through axis
static inline float vec3_dih(vec3 a, vec3 axis, vec3 c) {
  vec3 ap = vec3_perpendicular(a, axis);
  vec3 cp = vec3_perpendicular(c, axis);
  
  float ang = vec3_ang(ap, cp);
  
  if (vec3_dot(vec3_cross(ap, cp), axis) > 0)
    return -ang;
  else
    return ang;
}

static inline float vec3_pos_dih(vec3 p1, vec3 p2, vec3 p3, vec3 p4) {
  return vec3_dih(vec3_sub(p1, p2), vec3_sub(p2, p3), vec3_sub(p4, p3));
}

// Dihedral angle of three vectors, defined as an exterior spherical angle.
// Used in tripep_closure.c (Don't wanted to change the code)
static inline float vec3_dih_spherical(vec3 v1, vec3 v2, vec3 v3)
{
  float arg;
  vec3 p = vec3_cross(v1, v2);
  vec3 q = vec3_cross(v2, v3);
  vec3 s = vec3_cross(v3, v1);
  
  arg = vec3_dot(p,q)/sqrtf(vec3_mag2(p)*vec3_mag2(q));
  arg = match_sign(fminf(fabsf(arg),1.0e0), arg);
  return match_sign(acosf(arg), vec3_dot(s,v2));
}

// Dihedral angle defined by three bond vectors (p2-p1, p2-p3, p3-p4) connecting four atoms.
static inline float vec3_pos_dih_spherical(vec3 p1, vec3 p2, vec3 p3, vec3 p4)
{
  vec3 r12 = vec3_sub(p1, p2);
  vec3 r23 = vec3_sub(p3, p2);
  vec3 r34 = vec3_sub(p4, p3);
  
  return vec3_dih_spherical(r12, r23, r34);
}

// Benchmarks showed that this approach is faster than `vec3_dih`
static inline float vec3_dih_fast(vec3 v1, vec3 v2, vec3 v3)
{
  float arg;
  vec3 p = vec3_cross(v1, v2);
  vec3 q = vec3_cross(v3, v2);
  vec3 s = vec3_cross(v3, v1);
  
  arg = vec3_dot(p,q)/sqrtf(vec3_mag2(p)*vec3_mag2(q));
  arg = match_sign(fminf(fabsf(arg),1.0e0), arg);
  return match_sign(acosf(arg), vec3_dot(s,v2));
}


static inline float vec3_pos_dih_slow(vec3 p1, vec3 p2, vec3 p3, vec3 p4)
{
  vec3 r12 = vec3_norm(vec3_sub(p2, p1));
  vec3 r23 = vec3_norm(vec3_sub(p2, p3));
  vec3 r34 = vec3_norm(vec3_sub(p4, p3));
  
  vec3 n1 = vec3_cross(r12, r23);
  vec3 n2 = vec3_cross(r23, r34);
  vec3 m = vec3_cross(n1, r23);
  
  float x = vec3_dot(n1, n2);
  float y = vec3_dot(m, n2);
  
  return atan2f(y, x);
}

static inline int vec3_isnan(vec3 v)
{
  return isnan(v.x) || isnan(v.y) || isnan(v.z);
}

#undef SMALL

#endif
