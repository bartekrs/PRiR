#include<iostream>
#include<time.h>
#include"LennardJonesPotential.h"
#include"MonteCarlo.h"
#include"Particles.h"
#include"MyMPI.h"

using namespace std;

void calc( Particles *particles, int myRank, MyMPI *mmpi ) {
	LennardJonesPotential *p = new LennardJonesPotential();
	MonteCarlo *mc = new MonteCarlo();
	mc->setParticles(particles);
	mc->setPotential(p);
	mc->setMyMPI( mmpi );

	if ( myRank == 0 ) {
  	    double Etot = mc->calcTotalPotentialEnergy();
	    cout << "Etot          = " << Etot << endl;
	}

     	
	double kBT = 0.9;

// start pomiaru czasu
    double tStart = mmpi->MPI_Wtime();

	mc->shareParticles();
	for (int i = 0; i < TEMPERATURES; i++) {
		mc->setKBTinv(kBT); // ustalenie parametrów temperatury
		mc->calcMC( PARTICLES * PARTICLES_PER_TEMP_MULTI );
		kBT += 0.1;
	}
	mc->gatherParticles();

    double tStop = mmpi->MPI_Wtime();
// koniec pomiaru czasu 

	if ( myRank == 0 ) {
		double Etot = mc->getTotalPotentialEnergy();
		cout << "Etot get  = " << Etot << endl;
		Etot = mc->calcTotalPotentialEnergy();
		cout << "Etot calc = " << Etot << endl;

		double avrMinDist = mc->calcAvrMinDistance();

		cout << "Srednia odleglosc do najblizszego sasiada: " << avrMinDist << endl;

		cout << "Obliczenie wykonano w czasie: " << ( tStop - tStart ) << " sekund";
	}
}

int main(int argc, char **argv) {
	srandom( time( NULL ) );

	MyMPI *mmpi = new MyMPI();
	mmpi->MPI_Init( &argc, &argv );

	int myRank;
	mmpi->MPI_Comm_rank( MPI_COMM_WORLD, &myRank );

	Particles *pa = new Particles(PARTICLES);
	Particles *paBackup;
	pa->setBoxSize( BOX_SIZE );

	if ( myRank == 0 ) {
		// inicjalizacja położeń tylko w jednym z procesów
		pa->initializePositions(BOX_SIZE, DISTANCE);
		paBackup = new Particles( pa ); // tu tworzona jest kopia położeń cząstek
		// przyda się w testach...
	}

	calc( pa, myRank, mmpi );

	mmpi->MPI_Finalize();
}
