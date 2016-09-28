//
//  main.c
//  Sumatra
//
//  Created by Ferhat Ayaz on 02/06/16.
//  Copyright © 2016 Ferhat Ayaz. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "smtr.h"
#include "scg.h"

const char *resNames[] =
{
  "ILE", "GLN", "GLY", "GLU", "CYS", "ASP", "SER", "LYS", "PRO", "ASN",
  "VAL", "THR", "HIS", "TRP", "PHE", "ALA", "MET", "LEU", "ARG", "TYR"
};

float massByType[] =
{
  71.0779, 156.1857, 114.1026, 115.0874, 103.1429,
  129.114, 128.1292, 57.0513, 137.1393, 113.1576,
  113.1576, 128.1723, 131.1961, 147.1739, 97.1152,
  87.0773,  101.1039, 186.2099, 163.1733, 99.1311
};

static void normalize_mass(float mass[20]);
static int print_callback(void *data);
static int load_from_pdb(const char *pdbfile, vec3 **particles, float **mass);

int main(int argc, const char * argv[])
{
  int i, length;
  vec3 *particles;
  float *mass;
  float temperature = 0.5;

  if (argc < 2)
  {
    fprintf(stderr, "Usage %s <pdb-file>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  normalize_mass(massByType);
  
  length = load_from_pdb(argv[1], &particles, &mass);

  smtr_init(particles, mass, length, temperature);
  smtr_subscribe_event(print_callback, (void *)10);

  scg_add_streching_force();
  scg_add_bending_force(20.0);
  scg_add_torsion_force(1.0, 0.5);
  scg_add_nonbond_force();
  
  
//  print_particles(particles, chain->nres, stdout);
  for (i = 0; i < 10; i++)
  {
    printf("Running...\n");
    smtr_run_loop(100000);
  }
//  print_particles(particles, chain->nres, stdout);
}


static void normalize_mass(float mass[20])
{
  int i;
  float avg, total = 0;
  
  for (i = 0; i < 20; i++)
    total += mass[i];
  
  avg = total/20;
  
  for (i = 0; i < 20; i++)
    mass[i] /= avg;
}

void print_particles(vec3 *p, int n, FILE *o)
{
  while (n--)
  {
    printf("%3.3f %3.3f %3.3f\n", p->x, p->y, p->z);
    p++;
  }
}


static int print_callback(void *data)
{
  vec3 *p;
  if (smtr_ctx->currentTimeStep % 10000 == 0)
  {
    p = smtr_ctx->particles + (long)data;
    printf("%3.7f %3.3f %3.3f\n", p->x, p->y, p->z);
  }
  
  return 0;
}


static int load_from_pdb(const char *pdbfile, vec3 **particles, float **mass)
{
  char line[128];
  vec3 coord;
  FILE* f;
  int length, capacity, i;
  
  f = fopen(pdbfile, "r");
  if (!f) return 0;
  
  length = 0;
  capacity = 0;
  *particles = NULL;
  *mass = NULL;
  
  while (fgets(line, sizeof(line), f))
  {
    if (!strncmp(line, "ATOM  ", 6))
    {
      // This statement restricts loading to coarse grain only
      if (strncmp(line + 12, " CA", 3))
        continue;
      
      if (length == capacity)
      {
        capacity += 100;
        *particles = realloc(*particles, capacity * sizeof(vec3));
        *mass = realloc(*mass, capacity * sizeof(float));
      }
      
      coord.x = atof(line + 30);
      coord.y = atof(line + 38);
      coord.z = atof(line + 46);
      
      (*particles)[length] = coord;
      
      for (i = 0; i < 20; i++)
      {
        if (!strncmp(line + 17, resNames[i], 3))
        {
          (*mass)[length] = massByType[i];
          break;
        }
      }
      
      assert(i < 20);
      
      length++;
    }
    else if (!strncmp(line, "TER ", 4))
    {
      break;
    }
  }

  fclose(f);
  return length;
}
