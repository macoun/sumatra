//
//  forces.c
//  Sumatra
//
//  Created by Ferhat Ayaz on 03/06/16.
//  Copyright © 2016 Ferhat Ayaz. All rights reserved.
//

#include "scg.h"
#include "smtr.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// -------------------------------------------
// Stretching Force
// -------------------------------------------
void update_stretching_forces(void *userData, vec3 *results)
{
  static const float kWeighted = 100;

  vec3 force, normv, sub;
  int i;
  
  float dEnergy;
  float dist;
  float *nativedists = userData;
  vec3 *p = smtr_ctx->particles;
  
  for (i = 0; i < smtr_ctx->particleCount - 1; i++)
  {
    dist = svec3_dist(i, i + 1);
    dEnergy = -2*kWeighted*(dist - nativedists[i]);
    sub = vec3_sub(p[i], p[i + 1]);
    normv = vec3_mul(sub , 1/dist);
    force = vec3_mul(normv, dEnergy);
    results[i] = vec3_add(results[i], force);
    results[i + 1] = vec3_sub(results[i + 1], force);
  }
}

void scg_add_streching_force()
{
  int i;
  float *nativeDistance;
  
  nativeDistance = calloc(smtr_ctx->particleCount - 1, sizeof(float));
  for (i = 0; i < smtr_ctx->particleCount - 1; i++)
  {
    nativeDistance[i] = vec3_dist(smtr_ctx->particles[i],
                                  smtr_ctx->particles[i+1]);
  }
  smtr_add_force(update_stretching_forces, nativeDistance);
}

// -------------------------------------------
// Bending Force
// -------------------------------------------
struct bending_force_context
{
  float *nativeTetha;
  float *currentTetha;
  float ck_theta;
};

float calculate_bending_energy(int i, int j, int k, void *userData)
{
  struct bending_force_context *bf = userData;
  float delta = bf->currentTetha[i] - bf->nativeTetha[i];
  return bf->ck_theta*delta*delta;
}

static
void update_bending_force(int i, int j, int k, vec3 *results, void *userData)
{
  static const float epsTheta = 1.0 - 1.0e-6;
  struct bending_force_context *bf = userData;
  float r_ij = svec3_dist(i, j);
  float r_kj = svec3_dist(j, k);
  vec3 u = svec3_unit(i, j);
  vec3 v = svec3_unit(k, j);
  vec3 dr_i, dr_k, tmp;
  float theta, sinTheta, cosTheta;
  float dEnergy;
  
  cosTheta = vec3_dot(u, v);
  if (fabs(cosTheta) > epsTheta)
  {
    cosTheta = cosTheta >= 0 ? fabs(epsTheta) : -fabs(epsTheta);
  }
  theta = acos(cosTheta);
  sinTheta = sin(theta);
  bf->currentTetha[i] = theta;
  dEnergy = -2*bf->ck_theta*(theta - bf->nativeTetha[i]);
  
  dr_i = vec3_div(vec3_sub(vec3_mul(u, cosTheta), v), r_ij*sinTheta);
  dr_k = vec3_div(vec3_sub(vec3_mul(v, cosTheta), u), r_kj*sinTheta);
  
  tmp = vec3_mul(vec3_add(dr_i, dr_k), -1);
  results[i] = vec3_add(results[i], vec3_mul(dr_i, dEnergy));
  results[j] = vec3_add(results[j], vec3_mul(tmp, dEnergy));
  results[k] = vec3_add(results[k], vec3_mul(dr_k, dEnergy));

  if (vec3_isnan(results[i]) || vec3_isnan(results[j]) || vec3_isnan(results[k]))
  {
    // Not good man. Not good.
  }
}

void update_bending_forces(void *userData, vec3 *results)
{
  int i;

  for(i = 0; i < smtr_ctx->particleCount - 2; i++)
    update_bending_force(i, i + 1, i + 2, results, userData);
}

void scg_add_bending_force(float ck_theta)
{
  int i;
  struct bending_force_context *bf;
  
  bf = calloc(1, sizeof(struct bending_force_context));
  bf->nativeTetha = calloc(smtr_ctx->particleCount - 2, sizeof(float));
  bf->currentTetha = calloc(smtr_ctx->particleCount - 2, sizeof(float));
  for (i = 0; i < smtr_ctx->particleCount - 2; i++)
  {
    bf->nativeTetha[i] = vec3_pos_ang(smtr_ctx->particles[i],
                                  smtr_ctx->particles[i + 1],
                                  smtr_ctx->particles[i + 2]);
  }
  
  bf->ck_theta = kEpsilon * ck_theta;
  smtr_add_force(update_bending_forces, bf);
}

// -------------------------------------------
// Torsion Force
// -------------------------------------------
struct torsion_force_context
{
  float *nativePhi;
  float ck_phi1;
  float ck_phi3;
};

void update_torsion_forces(void *userData, vec3 *results)
{
  int i, j, k, l;
  static const float eps = 1.0e-3;
  
  float s;
  float r_kj, r_nk,r_mj;
  float rijrkj,rklrkj;
  
  float phi, phi_0, diff_phi,
    cosphi, sinphi, cosdphi, sindphi,
    cos3dphi, sinp3p0, cos3p0, sin3p0;
  
  float derivative;
  float derivative1;
  float derivative2;
  
  vec3 *p = smtr_ctx->particles;
  
  vec3 v_ij, v_kj, v_kl, v_mj, v_nk, v_il;
  vec3 dr_i, dr_l, dr_j, dr_k;
  float comp;
  
  struct torsion_force_context *tf = userData;
  
  for(i = 0; i < smtr_ctx->particleCount - 3; i++)
  {
    j = i + 1;
    k = i + 2;
    l = i + 3;
    
    v_ij = vec3_sub(p[i], p[j]);
    v_kj = vec3_sub(p[k], p[j]);
    v_kl = vec3_sub(p[k], p[l]);

    v_mj = vec3_cross(v_ij, v_kj);
    v_nk = vec3_cross(v_kj, v_kl);
    v_il = vec3_cross(v_mj, v_nk);
    
    r_kj = svec3_dist(j, k);
    r_nk = vec3_mag(v_nk);
    r_mj = vec3_mag(v_mj);
    
    if(((r_nk * r_nk) < eps) || ((r_mj * r_mj) < eps))
      continue;
    
    cosphi = vec3_dot(v_nk, v_mj)/(r_nk*r_mj);
    
    if (cosphi < -1.0) cosphi = -1;
    if (cosphi >  1.0) cosphi = 1;
    
    comp = vec3_dot(v_kj, v_il);
    s = (fabs(comp) <= 0.0000001 || comp > 0) ? 1.0 : -1.0;
    phi = (fabs(cosphi) <= 0.000001) ? s * M_PI_2 : s * acos(cosphi);
    
    phi_0 = tf->nativePhi[i];
    
    diff_phi = phi - phi_0;
    sinphi = sin(phi);
    cosdphi = cos(diff_phi);
    sindphi = sin(diff_phi);
    cos3dphi = cos(3*diff_phi);
    sinp3p0 = sin(phi-3*phi_0);
    cos3p0 = cos(3*phi_0);
    sin3p0 = sin(3*phi_0);
    
    derivative1 = tf->ck_phi1*sindphi;
    derivative2 = 3*tf->ck_phi3*(4*cosphi*cosphi*sinp3p0+3*cosphi*sin3p0-sinphi*cos3p0);
    derivative = -(derivative1+derivative2);
    
    rijrkj = vec3_dot(v_ij, v_kj)/(r_kj*r_kj);
    rklrkj = vec3_dot(v_kl, v_kj)/(r_kj*r_kj);
    
    dr_i = vec3_mul(v_mj, r_kj/(r_mj*r_mj));
    dr_l = vec3_mul(v_nk, -r_kj/(r_nk*r_nk));

    dr_j = vec3_sub(vec3_mul(dr_i, rijrkj - 1), vec3_mul(dr_l, rklrkj));
    dr_k = vec3_sub(vec3_mul(dr_l, rklrkj - 1), vec3_mul(dr_i, rijrkj));
    
    results[i] = vec3_add(results[i], vec3_mul(dr_i, derivative));
    results[j] = vec3_add(results[j], vec3_mul(dr_j, derivative));
    results[k] = vec3_add(results[k], vec3_mul(dr_k, derivative));
    results[l] = vec3_add(results[l], vec3_mul(dr_l, derivative));
    
    if (vec3_isnan(results[i]) || vec3_isnan(results[j]) ||
        vec3_isnan(results[k]) || vec3_isnan(results[l])) {
      // Breakpoint before armageddon
    }
    
 }
}

void scg_add_torsion_force(float ck_phi1, float ck_phi3)
{
  int i;
  vec3 *p = smtr_ctx->particles;
  struct torsion_force_context *tf;
  
  tf = calloc(1, sizeof(struct torsion_force_context));
  tf->nativePhi = calloc(smtr_ctx->particleCount - 3, sizeof(float));
  
  for (i = 0; i < smtr_ctx->particleCount - 3; i++)
    tf->nativePhi[i] = vec3_pos_dih(p[i], p[i + 1], p[i + 2], p[i + 3]);
  
  tf->ck_phi1 = kEpsilon * ck_phi1;
  tf->ck_phi3 = kEpsilon * ck_phi3;

  smtr_add_force(update_torsion_forces, tf);
}

// -------------------------------------------
// Torsion Force
// -------------------------------------------
struct nonbond_force_context
{
  int *contacts;
  int contactCount;
  float *nativeDists;
  float *epsilon;
};

void update_nonbond_forces(void *userData, vec3 *results)
{
  struct nonbond_force_context *nf = userData;
  int i, j, k, n = nf->contactCount * 2;
  float r_nat, r_ij;
  float rij12, rij10, frac, dEnergy;
  vec3 u;
  
  for (k = 0; k < n; k += 2)
  {
    r_nat = nf->nativeDists[k/2];
    i = nf->contacts[k];
    j = nf->contacts[k + 1];
    r_ij = svec3_dist(i, j);
    
    if (r_ij > 1.2 * r_nat)
      continue;
    
    frac = r_nat/r_ij;
    rij10 = powf(frac, 10);
    rij12 = rij10 * frac*frac;
    dEnergy = 60 * nf->epsilon[k/2]*(rij12 - rij10)/r_ij;

    assert(!isnan(dEnergy));
  
    u = vec3_sub(smtr_ctx->particles[i], smtr_ctx->particles[j]);
    u = vec3_div(u, r_ij);

    results[i] = vec3_add(results[i], vec3_mul(u, dEnergy/r_ij));
    results[j] = vec3_sub(results[j], vec3_mul(u, dEnergy/r_ij));
    
    if (vec3_isnan(results[i]) || vec3_isnan(results[j])) {
      // Time to blow up
    }
  }
}

void scg_add_nonbond_force()
{
  static const float cutoff = 6.0f;
  int i, k = 0, *n;
  float dist;
  struct nonbond_force_context *nf;
  
  nf = calloc(1, sizeof(struct nonbond_force_context));
  nf->contacts = calloc(smtr_ctx->particleCount, 2 * sizeof(int) * MAX_NEIGHBOUR_COUNT);
  nf->nativeDists = calloc(smtr_ctx->particleCount, sizeof(float) * MAX_NEIGHBOUR_COUNT);
  nf->epsilon = calloc(smtr_ctx->particleCount, sizeof(float) * MAX_NEIGHBOUR_COUNT);
  
  for (i = 0; i < smtr_ctx->particleCount; i++)
  {
    n = smtr_ctx->neighbours + i * MAX_NEIGHBOUR_COUNT;
    for (; *n ; n++)
    {
      if (*n - i < 3)
        continue;
      dist = svec3_dist(i, *n);
      if (dist < cutoff)
      {
        nf->nativeDists[k/2] = dist;
        nf->epsilon[k/2] = kEpsilon;
        nf->contacts[k++] = i;
        nf->contacts[k++] = *n;
      }
    }
  }
  
  nf->contactCount = k/2;
  // printf("%d native contacts\n", nf->contactCount);
  smtr_add_force(update_nonbond_forces, nf);
}


// r, u, f
// r: Ca – Ca virtual bond length between successive residues
// u: Ca – Ca virtual bond angles
// f: Ca – Ca virtual torsion angle
