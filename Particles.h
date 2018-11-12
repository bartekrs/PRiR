/*
 * Particles.h
 *
 *  Created on: 12 cze 2018
 *      Author: oramus
 */

#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <stdlib.h>
#include <math.h>
#include "Consts.h"

class Particles {
private:
	const int size;
	const double MAX_RANDOM;
	double *x;
	double *y;
	double *m; // mass
	double *Vx;
	double *Vy; // velocity
	double *Fx;
	double *Fy; // force
	double boxSize;
	double revBoxSize;

	double rnd() {
		return random() * MAX_RANDOM;
	}
public:
	Particles( int size );
	Particles( Particles *pointer );

	double getX( int pos ) {
		return x[ pos ];
	}
	double getY( int pos ) {
		return y[ pos ];
	}
	double getM( int idx ) {
		return m[idx];
	}
	double getVx( int idx ) {
		return Vx[ idx ];
	}
	double getVy( int idx ) {
		return Vy[ idx ];
	}
	double getFx( int idx ) {
		return Fx[ idx ];
	}
	double getFy( int idx ) {
		return Fy[ idx ];
	}

	double pbc_pos( double x ) {
		return x - boxSize * floor(x/boxSize);
	}

	double pbc_dpos( double dx ) {
		return dx - boxSize * nearbyint(dx * revBoxSize);
	}

	double pbc_dpos_usign( double dx ) {
		dx = fabs( dx );
		return dx - static_cast<int>(dx * revBoxSize +  0.5) * boxSize;
	}

	void updateX( int idx, double value ) {
#ifdef PBC
		x[ idx ] = pbc_pos( x[ idx ] + value );
#else
		x[ idx ] += value;
#endif
	}
	void updateY( int idx, double value ) {
#ifdef PBC
		y[ idx ] = pbc_pos( y[ idx ] + value );
#else
		y[ idx ] += value;
#endif
	}
	void updateFx( int idx, double value ) {
		Fx[ idx ] += value;
	}
	void updateFy( int idx, double value ) {
		Fy[ idx ] += value;
	}
	void updateVx( int idx, double value ) {
		Vx[ idx ] += value;
	}
	void updateVy( int idx, double value ) {
		Vy[ idx ] += value;
	}
	void setVx( int idx, double value ) {
		Vx[ idx ] = value;
	}
	void setVy( int idx, double value ) {
		Vy[ idx ] = value;
	}
	void setXY( int pos, double newx, double newy ) {
#ifdef PBC
		x[pos] = pbc_pos( newx );
		y[pos] = pbc_pos( newy );
#else
		x[pos] = newx;
		y[pos] = newy;
#endif
	}

	void setBoxSize( double size ) {
		this->boxSize = size;
		this->revBoxSize = 1.0 / size;
	}

	int getNumberOfParticles() {
		return size;
	}

	inline double getDx( int idx1, int idx2 ) {
#ifdef PBC
		return pbc_dpos( x[ idx1 ] - x[ idx2 ] );
#else
		return x[ idx1 ] - x[ idx2 ];
#endif
	}

	inline double getDy( int idx1, int idx2 ) {
#ifdef PBC
		return pbc_dpos( y[ idx1 ] - y[ idx2 ] );
#else
		return y[ idx1 ] - y[ idx2 ];
#endif
	}

	inline double getDxUsign( int idx1, int idx2 ) {
#ifdef PBC
		return pbc_dpos_usign( x[ idx1 ] - x[ idx2 ] );
#else
		return x[ idx1 ] - x[ idx2 ];
#endif
	}

	inline double getDyUsign( int idx1, int idx2 ) {
#ifdef PBC
		return pbc_dpos_usign( y[ idx1 ] - y[ idx2 ] );
#else
		return y[ idx1 ] - y[ idx2 ];
#endif
	}

	double getDistanceSQ( double dx, double dy ) {
		return dx * dx + dy * dy;
	}

	double getDistance( double dx, double dy ) {
		return sqrt( dx * dx + dy * dy );
	}

	double getDistanceBetweenSQ( int idx1, int idx2 ) {
		double dx = getDxUsign( idx1, idx2 );
		double dy = getDyUsign( idx1, idx2 );
		return dx * dx + dy * dy;
	}

	double getDistanceBetweenSQ( int idx, double xx, double yy ) {
#ifdef PBC
		double dx = pbc_dpos_usign( x[ idx ] - xx );
		double dy = pbc_dpos_usign( y[ idx ] - yy );
#else
		double dx = x[ idx ] - xx;
		double dy = y[ idx ] - yy;
#endif
		return dx * dx + dy * dy;
	}

	void clearForce( int idx ) {
		Fx[ idx ] = 0;
		Fy[ idx ] = 0;
	}

	void stop( int idx ) {
		Vx[ idx ] = 0;
		Vy[ idx ] = 0;
	}

	double getKineticEnergy();

	double getDistanceSQToClosest( int idx );

	void initializePositions(double length, double minDistance  );

	virtual ~Particles();
};

#endif /* PARTICLES_H_ */
