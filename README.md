#SUMATRA
-
A simple molecular dynamics engine driven by Langevin dynamics.

Sumatra is inspired by the molecular dynamic model as described in the paper:

> JMB, 2003 H. Kaya, H. Chan (p. 911-931): 
> 
> Solvation Effects and Driving Forces for Protein Thermodynamic and Kinetic Cooperativity: How Adequate is Native-centric Topological Modeling?

The emphasis is put on simplicity and robustness of developing a molecular dynamics simulation written in C. A full functional simulation can be expressed by a few core functions of Sumatra.

## Building and Running

It's not rocket science.

	git clone https://github.com/macoun/sumatra.git
	cd sumatra
	make

I've also included a demo pdb file to start right away. Run the demo simulation as follows.

	./sumatra res/1cqu.pdb


## Architecture

Module | Description
------ | -----------
smtr.h/c | The Sumatra engine. The *verlet algorithm* and a small framework, which enables subscribing to events on trajectories and defining force fields, comprise the engine.
scg.h/c | The Sumatra Coarse Grain module implements the force fields as described in the reference paper.
vec3.h | You guessed it correctly! Another vector algebra library.

A schematic representation of a typical simulation application:

![Sumatra Diagram](https://github.com/macoun/sumatra/raw/master/res/diagram.png)

You can consider the SCG module as a helper or wrapper to `smtr_add_force()`. It also contains the *update functions* to the supported forces: *stretching, bending, torsion*, and *nonbonded*. It is seperated from the core Sumatra engine and can be replaced easily by your own force implementations.

## Example

The following example shows how to setup a simluation and run it for 100k time steps at 0.5 simulation temperature.
		
	void main()
	{
		vec3 *particles = //...
		float *mass = //...
		int length = // length of particles and mass arrays
	
		smtr_init(particles, mass, length, 0.5);
		smtr_subscribe_event(print_event, NULL);
		
		scg_add_streching_force();
		scg_add_bending_force(20.0);
		scg_add_torsion_force(1.0, 0.5);
		scg_add_nonbond_force();
		  
		smtr_run_loop(100000);
	}

For the sake of simplicity `print_event()` will print the coordinates of the 10th particle every 10000 steps.

	int print_event(void *userData)
	{
	  vec3 *p;
	  if (smtr_ctx->currentTimeStep % 10000 == 0)
	  {
	    p = smtr_ctx->particles + 10;
	    printf("%3.7f %3.3f %3.3f\n", p->x, p->y, p->z);
	  }
	  
	  // 0 means don't stop the simulation
	  return 0;
	}

## What else

- A random force as defined in Langevin dynamics will be added automatically if the simulation temperature is above zero. 

- Initial particle velocities are given randomly during initialization.
- This is an ongoing project and I discourage using this in productive researches.

## License

This version is licensed under MIT.
