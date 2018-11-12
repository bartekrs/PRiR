/*
 * LennardJonesPotential.h
 *
 *  Created on: 12 cze 2018
 *      Author: oramus
 */

#include"PotentialEnergy.h"

#ifndef LENNARDJONESPOTENTIAL_H_
#define LENNARDJONESPOTENTIAL_H_

class LennardJonesPotential : public PotentialEnergy {
private:
	const double epsilon;
	const double rm;
	 double rmSQ;
	 double epsilon12;
public:
	LennardJonesPotential();
	virtual ~LennardJonesPotential();

	double getPotentialEnergy( double distance );
	double getPotentialEnergyDistanceSQ( double distanceSQ );
	double getForce( double distace );
	double getForceDistanceSQ( double distaceSQ );
};

#endif /* LENNARDJONESPOTENTIAL_H_ */
