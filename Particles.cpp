/*
 * Particles.cpp
 *
 *  Created on: 12 cze 2018
 *      Author: oramus
 */

#include "Particles.h"
#include <stdlib.h>
#include <iostream>

Particles::Particles( int size ) : size( size ), MAX_RANDOM( 1.0 / ( RAND_MAX + 1.0 ) ) {
	x = new double[ size ];
	y = new double[ size ];
	m = new double[ size ];
	Fx = new double[ size ];
	Fy = new double[ size ];
	Vx = new double[ size ];
	Vy = new double[ size ];

	for ( int i = 0; i < size; i++ )
		m[ i ] = 1.0;
}

Particles::Particles( Particles *pointer ) : size( pointer->size), MAX_RANDOM(1.0/(MAX_RANDOM+1.0)) {
	x = new double[ pointer->size ];
	y = new double[ pointer->size ];
	m = new double[ pointer->size ];
	Fx = new double[ pointer->size ];
	Fy = new double[ pointer->size ];
	Vx = new double[ pointer->size ];
	Vy = new double[ pointer->size ];

    for ( int i = 0; i < pointer->size; i++ ) {
		this->x[ i ]  = pointer->x[ i ];
		this->y[ i ]  = pointer->y[ i ];
		this->Fx[ i ] = pointer->Fx[ i ];
		this->Fy[ i ] = pointer->Fy[ i ];
		this->Vx[ i ] = pointer->Vx[ i ];
		this->Vy[ i ] = pointer->Vy[ i ];
		this->m[ i ]  = pointer->m[ i ];
	}

	boxSize = pointer->boxSize;
	revBoxSize = pointer->revBoxSize;
}

Particles::~Particles() {
	delete[] x;
	delete[] y;
	delete[] m;
	delete[] Fx;
	delete[] Fy;
	delete[] Vx;
	delete[] Vy;
}

double Particles::getDistanceSQToClosest( int idx ) {
	double min = 100000.0;
	double tmp;
	for ( int i = 0; i < size; i++ ) {
		if ( i != idx ) {
			tmp = getDistanceBetweenSQ( idx, i );
			if ( tmp < min ) {
				min = tmp;
			}
		}
	}
	return min;
}

void Particles::initializePositions( double length, double minDistance ) {

	double minDistanceSQ = minDistance * minDistance;

	x[ 0 ] = length * rnd();
	y[ 0 ] = length * rnd();

	int idx = 1;
	double xtmp, ytmp, dx, dy;
	bool isOK;
	do {
		xtmp = length * rnd();
		ytmp = length * rnd();

		isOK = true;
		for ( int i = 0; i < idx; i++ ) {
			dx = xtmp - x[ i ];
			dy = ytmp - y[ i ];

			if ( dx * dx + dy * dy < minDistanceSQ ) {
				isOK = false;
				break;
			}
		}

		if ( isOK ) {
			x[ idx ] = xtmp;
			y[ idx ] = ytmp;
			idx++;
			// std::cout << " idx = " << idx << " [ " << xtmp << ", " << ytmp << " ]" << std::endl;
		}

	} while ( idx < size );
}


double Particles::getKineticEnergy() {
	double sum = 0;

	for ( int i = 0; i < size; i++ ) {
		sum += m[ i ] * ( Vx[ i ] * Vx[ i ] + Vy[ i ] * Vy[ i ] );
	}

	return sum * 0.5;
}
