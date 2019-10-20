#include<iostream>
#include<iomanip>
#include<bitset>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<vector>
#include<ctime>
#define random_number (double)rand()/(RAND_MAX +1.0)
	
	using namespace std;

const int nx = 20;
double Energy(int spin[nx][nx]);
void MonteCarloStep(int p,int q,double T,int (&spin)[nx][nx]);
int Magnet(int spin[nx][nx]);
main(){

	srand(time(NULL));
	//Initialize the spin system
	int nthermal = 20000;
	int nmeasure = 10000;
	//The lattice size
	int spin[nx][nx];
	for (int i = 0; i < nx; i++){
		for (int j = 0; j < nx; j++){
   			spin[i][j] = (random_number > 0.5)?1:-1;
			}
}
	double dT = 0.02;
	ofstream outfile("20x20.dat");
	for(double T = 0.2; T <= 5.0; T += dT){
		double hsum  = 0;
		double hsum2 = 0;
		double msum  = 0;
		double msum2 = 0;	
		double h;	
		double m;
		//Measure the thermal average of physical quantities
		for(int i = 0; i < nmeasure; i++){
			//Thermalize the system
			for(int l = 0; l < nthermal; l++){
				int i = rand()%nx;
				int j = rand()%nx;
				MonteCarloStep(i,j,T,spin);
			}	
			h = Energy(spin);
			hsum  += h;
			hsum2 += h*h;
			m = (double)Magnet(spin);
			msum  += m;
		
		}
		hsum  /= nmeasure;
		hsum2 /= nmeasure;
		msum  /= nmeasure;

		//Print out the measured quntities	
		outfile << T << setprecision(9) << " " << hsum/(nx*nx) <<" " <<(hsum2 - hsum*hsum)/(T*T*nx*nx) << " " << msum/(nx*nx) << endl;
	}
	outfile.close();
}

double Energy(int spin[nx][nx]){
	int sum = 0;
	for(int p = 0; p < nx; p++){
        	for(int q = 0; q < nx; q++){
                	p += nx;
                        q += nx;
                        sum += -(spin[p%nx][(q-1)%nx]+spin[p%nx][(q+1)%nx]+spin[(p-1)%nx][q%nx]+spin[(p+1)%nx][q%nx])*spin[p%nx][q%nx];	
                	p -= nx;
                        q -= nx;
			
		}
	}
	return (double)sum/2;
}
	
void MonteCarloStep(int p, int q, double T, int (&spin)[nx][nx]){
//	double dE,oldE,newE;
//	oldE = Energy(spin);
	double dE;
	spin[p][q] *= -1; 		
//	newE = Energy(spin);
//	dE = newE-oldE;
	p += nx;
        q += nx;
        dE = -2*(spin[p%nx][(q-1)%nx]+spin[p%nx][(q+1)%nx]+spin[(p-1)%nx][q%nx]+spin[(p+1)%nx][q%nx])*spin[p%nx][q%nx];
        p -= nx;
	q -= nx;
	if(dE <=0){ spin[p][q]*=1;}
	else if (random_number <= exp(-dE/T)){spin[p][q]*=1;}
	else {spin[p][q]*=-1;}				
}
int Magnet(int spin[nx][nx]){
	int m = 0;
	for(int i = 0; i < nx; i++){
		for(int j = 0; j < nx; j++){
			m += spin[i][j];
		}
	}
	return m;
}

