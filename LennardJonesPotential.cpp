/*
 * LennardJonesPotential.cpp
 *
 *  Created on: 12 cze 2018
 *      Author: oramus
 */

#include <iostream>
#include <math.h>
#include "LennardJonesPotential.h"
#include "Particles.h"
#include "MonteCarlo.h"

LennardJonesPotential::LennardJonesPotential() :
		epsilon(1.0), rm(1.0) {
	epsilon12 = epsilon * 12;
	rmSQ = rm * rm;
}

LennardJonesPotential::~LennardJonesPotential() {
}

double LennardJonesPotential::getPotentialEnergy(double distance) {
	double h = rm / distance;

	h = h * h * h; // h^3
	h = h * h; // h^6

	return epsilon * h * (h - 2.0);
}

double LennardJonesPotential::getPotentialEnergyDistanceSQ(double distanceSQ) {
	double hSQ = rmSQ / distanceSQ;

	hSQ = hSQ * hSQ * hSQ; // h^6

	return epsilon * hSQ * (hSQ - 2.0);
}

double LennardJonesPotential::getForce(double distance) {
	double h = rm / distance;

	h = h * h * h; // h^3
	h = h * h; // h^6

	return epsilon12 * h * (1.0 - h) / distance;
}

double LennardJonesPotential::getForceDistanceSQ(double distanceSQ) {
	double distance = sqrt(distanceSQ);
	double h = rm / distance;
	double hSQ = rmSQ / distanceSQ;

	h = h * hSQ; // h^3
	h = h * h; // h^6

	return epsilon12 * h * (1.0 - h) / distance;
}
