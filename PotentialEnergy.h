/*
 * PotentialEnergy.h
 *
 *  Created on: 12 cze 2018
 *      Author: oramus
 */

#ifndef POTENTIALENERGY_H_
#define POTENTIALENERGY_H_

class PotentialEnergy {
public:
	PotentialEnergy();
	virtual ~PotentialEnergy();

	virtual double getPotentialEnergy( double distance ) = 0;
	virtual double getPotentialEnergyDistanceSQ( double distanceSQ ) = 0;
	virtual double getForce( double distance ) = 0;
	virtual double getForceDistanceSQ( double distanceSQ ) = 0;
};

#endif /* POTENTIALENERGY_H_ */
