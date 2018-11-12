/*
 * Consts.h
 *
 *  Created on: 13 cze 2018
 *      Author: oramus
 */

#ifndef CONSTS_H_
#define CONSTS_H_

const double DR_CORRECTION = 0.2;
const double DR_MAX = 0.3;
const double DR_MIN = 0.1;
const double DR_INITIAL_RATIO = 0.5 * ( DR_MAX + DR_MIN );

const double BOX_SIZE = 80.0;
const double DISTANCE = 0.9;

const int PARTICLES = 5000;
const int PARTICLES_PER_TEMP_MULTI = 3;
const int TEMPERATURES = 50;

// włączamy periodyczne warunki brzegowe
#define PBC

#endif /* CONSTS_H_ */
