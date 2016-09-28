//
//  smtr.h
//  Sumatra
//
//  Created by Ferhat Ayaz on 05/06/16.
//  Copyright © 2016 Ferhat Ayaz. All rights reserved.
//

// Model as described in
//
// JMB, 2003 H. Kaya, H. Chan: Solvation Effects and Driving Forces for
//    Protein Thermodynamic and Kinetic Cooperativity:
//      How Adequate is Native-centric Topological Modeling? (p. 911-931)
//

#ifndef smtr_h
#define smtr_h

#include "vec3.h"

// alpha: length scale (in angstron)
static const float kAlpha = 4.0; 

// epsilon: reference energy scale (in kcal/mol)
static const float kEpsilon = 1.0;

// mass: (in ?)
static const float kMass = 1.0;

// tau: time scale
extern float kTimeUnit;

// delta: integration time step
extern float kIntegrationStep;

// velocity based verlet algorithm coeeficients
// used in equation of motion
extern float kVelocityScaleFactor;

#define MAX_FORCE_COUNT 10
#define MAX_SUBSCRIBER_COUNT 10
#define MAX_NEIGHBOUR_COUNT 20
#define MAX_INTERACTING_DIST 30
#define DIST_CALC_INTERVAL 10

typedef struct SmtrContext SmtrContext;
typedef struct SmtrForceClass SmtrForceClass;
typedef struct SmtrSubscriber SmtrSubscriber;

typedef void (*SmtrUpdateForceFunc)(void *userData, vec3 *results);

struct SmtrForceClass
{
  vec3 *forces;
  void *userData;
  SmtrUpdateForceFunc update;
};

typedef int (*SmtrCallbackFunc)(void *userData);

struct SmtrSubscriber
{
  void *userData;
  SmtrCallbackFunc callback;
};

struct SmtrContext
{
  vec3 *particles;
  int particleCount;
  float *mass;
  float *distances;
  int *neighbours;
  long currentTimeStep;
  vec3 *velocities;
  float *vforceScalars;
  
  vec3 *forces;
  int forceCount;
  SmtrForceClass forceClasses[MAX_FORCE_COUNT];

  int subscriberCount;
  SmtrSubscriber subscribers[MAX_SUBSCRIBER_COUNT];
};

extern SmtrContext *smtr_ctx;

int smtr_init(vec3 *particles, float *mass,
              int particleCount, float temperature);
void smtr_run_loop(long steps);

void smtr_add_force(SmtrUpdateForceFunc func, void *userData);
void smtr_subscribe_event(SmtrCallbackFunc func, void *userData);

// Extensions to vec3
static inline float svec3_dist(int i, int j)
{
  return (i < j) ? smtr_ctx->distances[i*smtr_ctx->particleCount + j] :
  smtr_ctx->distances[j*smtr_ctx->particleCount + i];
}

static inline vec3 svec3_unit(int i, int j)
{
  vec3 *p = smtr_ctx->particles;
  return vec3_div(vec3_sub(p[i], p[j]), svec3_dist(i, j));
}


#endif /* smtr_ctx_h */

// time scale
// alpha: length scale
// epsilon: reference energy scale
// gamma: friction coefficient
// eta: random number generator
// delta: time step
