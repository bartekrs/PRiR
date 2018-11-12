/*
 * MyMPI.cpp
 *
 *  Created on: 01.05.2017
 *      Author: oramus
 */

#include "MyMPI.h"
#include <iostream>

using namespace std;

MyMPI::MyMPI() {
	myRank = -1;
	initStruct(&send);
	initStruct(&recv);
	initStruct(&bcastS);
	initStruct(&bcastR);
	initStruct(&scatterS);
	initStruct(&scatterR);
	initStruct(&gatherS);
	initStruct(&gatherR);
	initStruct(&reduceS);
	initStruct(&reduceR);
}

int MyMPI::MPI_Sendrecv( void* sendbuf, int sendcount,
		MPI_Datatype sendtype, int dest, int sendtag, void* recvbuf,
		int recvcount, MPI_Datatype recvtype, int source, int recvtag,
		MPI_Comm comm, MPI_Status* status) {
	int result = ::MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,
			recvcount, recvtype, source, recvtag, comm, status );
	update(&send, sendcount, sendtype);
	update(&recv, status2count(status, recvtype), recvtype);

	return result;
}

MyMPI::~MyMPI() {
}

void MyMPI::sendTo0(mpiusage *mu) {
	long * tmp = new long[2];
	tmp[0] = mu->bytes;
	tmp[1] = mu->counter;
	::MPI_Send(tmp, 2, MPI_LONG, 0, 0, MPI_COMM_WORLD);
	delete tmp;
}

void MyMPI::sendTo0() {
	sendTo0(&send);
	sendTo0(&recv);
	sendTo0(&bcastS);
	sendTo0(&bcastR);
	sendTo0(&scatterS);
	sendTo0(&scatterR);
	sendTo0(&gatherS);
	sendTo0(&gatherR);
	sendTo0(&reduceS);
	sendTo0(&reduceR);
}

void MyMPI::recvFromN(int src) {
	recvFromN(&send, src);
	recvFromN(&recv, src);
	recvFromN(&bcastS, src);
	recvFromN(&bcastR, src);
	recvFromN(&scatterS, src);
	recvFromN(&scatterR, src);
	recvFromN(&gatherS, src);
	recvFromN(&gatherR, src);
	recvFromN(&reduceS, src);
	recvFromN(&reduceR, src);
}

void MyMPI::showAll(int proc) {
	cout << "------- process " << proc << endl;
	showReport(&send, "Simple Send....");
	showReport(&recv, "Simple Recv....");
	showReport(&bcastS, "Broadcast Send.");
	showReport(&bcastR, "Broadcast Recv.");
	showReport(&scatterS, "Scatter Send...");
	showReport(&scatterR, "Scatter Recv...");
	showReport(&gatherS, "Gather Send....");
	showReport(&gatherR, "Gather Recv....");
	showReport(&reduceS, "Reduce Send....");
	showReport(&reduceR, "Reduce Recv....");
}

void MyMPI::showAll() {
	showAll(0);
	for (int i = 1; i < size; i++) {
		recvFromN(i);
		showAll(i);
	}
}

void MyMPI::recvFromN(mpiusage *mu, int src) {
	long * tmp = new long[2];
	MPI_Status status;
	::MPI_Recv(tmp, 2, MPI_LONG, src, 0, MPI_COMM_WORLD, &status);
	mu->bytes = tmp[0];
	mu->counter = tmp[1];
	delete tmp;
}

int MyMPI::MPI_Init(int *argc, char ***argv) {
	int result = ::MPI_Init(argc, argv);
	::MPI_Comm_rank( MPI_COMM_WORLD, &myRank);
	::MPI_Comm_size( MPI_COMM_WORLD, &size);
	return result;
}

void MyMPI::showStats() {
	if (myRank == 0) {
		showAll();
	} else {
		sendTo0();
	}
}


int MyMPI::MPI_Finalize() {
        showStats();
	return ::MPI_Finalize();
}

int MyMPI::MPI_Comm_size(MPI_Comm comm, int *size) {
	return ::MPI_Comm_size(comm, size);
}
int MyMPI::MPI_Comm_rank(MPI_Comm comm, int *rank) {
	return ::MPI_Comm_rank(comm, rank);
}

int MyMPI::MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest,
		int tag, MPI_Comm comm) {
	update(&send, count, datatype);
	return ::MPI_Send(buf, count, datatype, dest, tag, comm);
}

int MyMPI::MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source,
		int tag, MPI_Comm comm, MPI_Status *status) {
	int result = ::MPI_Recv(buf, count, datatype, source, tag, comm, status);

	update(&recv, status2count(status, datatype), datatype);

	return result;
}

int MyMPI::MPI_Reduce(void *sendbuf, void *recvbuf, int count,
		MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {

	if (myRank == root) {
		update(&reduceR, count, datatype);
		update(&reduceS, count, datatype);
	} else {
		update(&reduceS, count, datatype);
	}

	return ::MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
}

int MyMPI::MPI_Get_count(MPI_Status *status, MPI_Datatype datatype,
		int *count) {
	return ::MPI_Get_count(status, datatype, count);
}

int MyMPI::MPI_Barrier(MPI_Comm comm) {
	return ::MPI_Barrier(comm);
}

int MyMPI::MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
		MPI_Comm comm) {

	if (myRank == root) {
		update(&bcastS, count, datatype);
	} else {
		update(&bcastR, count, datatype);
	}

	return ::MPI_Bcast(buffer, count, datatype, root, comm);
}

int MyMPI::MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
		MPI_Comm comm) {

	if (myRank == root) {
		update(&scatterS, sendcount, sendtype);
		update(&scatterR, recvcount, recvtype);
	} else {
		update(&scatterR, recvcount, recvtype);
	}

	return ::MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
			recvtype, root, comm);
}

int MyMPI::MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
		MPI_Comm comm) {

	if (myRank == root) {
		update(&gatherS, sendcount, sendtype);
		update(&gatherR, recvcount, recvtype);
	} else {
		update(&gatherS, sendcount, sendtype);
	}

	return ::MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
			recvtype, root, comm);
}

int MyMPI::MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
	update(&gatherR, recvcount, recvtype);
	update(&gatherS, sendcount, sendtype);
	return ::MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
			recvtype, comm);
}

int MyMPI::MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
		MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {

	update(&reduceR, count, datatype);
	update(&reduceS, count, datatype);

	return ::MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}


int MyMPI::MPI_Scatterv(void *sendbuf, int *sendcounts, int *displs,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount,
                 MPI_Datatype recvtype,
                 int root, MPI_Comm comm) {
        return ::MPI_Scatterv(sendbuf, sendcounts, displs, sendtype, 
                 recvbuf, recvcount, recvtype, root, comm);
}

int MyMPI::MPI_Gatherv( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf,  int *recvcounts,  int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm) {
	return ::MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, 
	        displs, recvtype, root, comm);
}

int MyMPI::MPI_Allgatherv( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf,  int *recvcounts,  int *displs,
                MPI_Datatype recvtype, MPI_Comm comm) {
                
        return ::MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, 
                displs, recvtype, comm);
}

//////////// private

int MyMPI::getDataTypeSize(MPI_Datatype datatype) {
	int size;
	MPI_Type_size(datatype, &size);
	return size;
}

void MyMPI::update(mpiusage* mu, long counts, MPI_Datatype datatype) {
	mu->bytes += counts * getDataTypeSize(datatype);
	mu->counter++;
}

long MyMPI::status2count(MPI_Status *status, MPI_Datatype datatype) {
	int count;

	::MPI_Get_count(status, datatype, &count);

	return count;
}

void MyMPI::setStruct(mpiusage* mu, long bytes, long count) {
	mu->bytes = bytes;
	mu->counter = count;
}

void MyMPI::initStruct(mpiusage* mu) {
	setStruct(mu, 0, 0);
}

void MyMPI::showReport(mpiusage *mu, const char *txt) {
	cout << txt << " " << mu->bytes << "B, executed " << mu->counter << " avg: "
			<< mu->bytes / (mu->counter ? mu->counter : 1) << endl;
	;
}

double MyMPI::MPI_Wtime() {
	return ::MPI_Wtime();
}
