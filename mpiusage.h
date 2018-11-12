/*
 * mpiusage.h
 *
 *  Created on: 22.05.2017
 *      Author: oramus
 */

#ifndef MPIUSAGE_H_
#define MPIUSAGE_H_

class mpiusage {
public:
	mpiusage();

	long counter;
	long bytes;

	virtual ~mpiusage();
};

#endif /* MPIUSAGE_H_ */
