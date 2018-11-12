/*
 * MonteCarlo.h
 *
 *  Created on: 13 cze 2018
 *      Author: oramus
 */

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include"Particles.h"
#include "PotentialEnergy.h"
#include "MyMPI.h"

class MonteCarlo {
private:
	double dx, dy, dr;
	double kBTinv;
	double MAX_RANDOM;
	Particles *particles;
	Particles *processParticles;
	PotentialEnergy *energy;
	MyMPI *myMPI;
	double totalEp;
	void calcInitialDr();
	double calcContribution( int idx, double xx, double yy );
	double deltaEp( int idx, double oldX, double oldY, double newX, double newY );
	double rnd() {
		return random() * MAX_RANDOM;
	}
public:
	MonteCarlo();
	virtual ~MonteCarlo();
	void setParticles( Particles *particles );
	void setPotential( PotentialEnergy *energy );
	void calcMC( int draws );
	double calcAvrMinDistance();
	double calcTotalPotentialEnergy();

	double getTotalPotentialEnergy() {
		return totalEp;
	}
	void setMyMPI( MyMPI *myMPI ) {
		this->myMPI = myMPI;
	}
	void setKBTinv( double kBTinv ) {
		this->kBTinv = -kBTinv;
	}

	void shareParticles();
	void gatherParticles();
};

#endif /* MONTECARLO_H_ */
