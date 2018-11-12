/*
 * MyMPI.h
 *
 *      Author: oramus
 */

#ifndef MYMPI_H_
#define MYMPI_H_

#include "mpiusage.h"
#include <mpi.h>

/**
 * Klasa jest posrednikiem pomiedzy programem liczacym
 * a MPI. Jest podzbiorem funkcji dostepnych w MPI.
 * Zawarte w niej zostaly wylacznie najczesciej
 * uzywane funkcje.
 * Dodatkowo, klasa pozwala generowac komunikaty i zbierac
 * statystyki uzycia. W tej chwili zaimplementowana
 * jest podstawowa funkcjonalnosc, tj. wywolywanie
 * metod prowadzi do uruchamiania fukcji MPI o
 * tej samej nazwie.
 *
 * @brief Klasa posrednik pomiedzy kodem liczacym a MPI.
 *
 * @version 0.17
 * @author Piotr Marek Oramus
 */



using namespace std;

class MyMPI {
private:
	mpiusage send;
	mpiusage recv;
	mpiusage bcastS;
	mpiusage bcastR;
	mpiusage scatterS;
	mpiusage scatterR;
	mpiusage gatherS;
	mpiusage gatherR;
	mpiusage reduceS;
	mpiusage reduceR;

	int myRank;
	int size;

	int getDataTypeSize(MPI_Datatype datatype);
	void update(mpiusage *mu, long counts, MPI_Datatype datatype);
	long status2count(MPI_Status *status, MPI_Datatype datatype);
	void initStruct(mpiusage *mu);
	void setStruct(mpiusage *mu, long bytes, long count);
	void showReport(mpiusage *mu, const char *txt);
	void showAll(int proc);
	void showAll();
	void sendTo0(mpiusage *mu);
	void sendTo0();
	void recvFromN(mpiusage *mu, int src);
	void recvFromN(int src);
public:
	MyMPI();

	int MPI_Init(int *argc, char ***argv);
	int MPI_Finalize();
	int MPI_Comm_size(MPI_Comm comm, int *size);
	int MPI_Comm_rank(MPI_Comm comm, int *rank);
	int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
			MPI_Comm comm);
	int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source,
			int tag, MPI_Comm comm, MPI_Status *status);
	int MPI_Reduce(void *sendbuf, void *recvbuf, int count,
			MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

	int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);

	int MPI_Barrier(MPI_Comm comm);

	int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
			MPI_Comm comm);

	int MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype,
			void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
			MPI_Comm comm);

	int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
			void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
			MPI_Comm comm);

	int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
			void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);

	int MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
			MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

	int MPI_Sendrecv( void *sendbuf, int sendcount, MPI_Datatype sendtype,
	                int dest, int sendtag,
	                void *recvbuf, int recvcount, MPI_Datatype recvtype,
	                int source, int recvtag,
	                MPI_Comm comm, MPI_Status *status);

// implementacja ponizszych metod nie produkuje danych statystycznych
        int MPI_Scatterv( void *sendbuf,  int *sendcounts,  int *displs,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount,
                 MPI_Datatype recvtype,
                 int root, MPI_Comm comm);

        int MPI_Gatherv( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 void *recvbuf,  int *recvcounts,  int *displs,
                 MPI_Datatype recvtype, int root, MPI_Comm comm);

	int MPI_Allgatherv( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf,  int *recvcounts,  int *displs,
                   MPI_Datatype recvtype, MPI_Comm comm);
                   
	double MPI_Wtime();

	virtual ~MyMPI();
	
	void showStats();
};

#endif /* MYMPI_H_ */
