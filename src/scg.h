//
//  forces.h
//  Sumatra
//
//  Created by Ferhat Ayaz on 03/06/16.
//  Copyright © 2016 Ferhat Ayaz. All rights reserved.
//

#ifndef forces_h
#define forces_h

void scg_add_streching_force();
void scg_add_bending_force(float ck_theta);
void scg_add_torsion_force(float ck_phi1, float ck_phi3);
void scg_add_nonbond_force();

#endif /* forces_h */
