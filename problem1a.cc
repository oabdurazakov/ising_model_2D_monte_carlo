#include<iostream>
#include<bitset>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<vector>
#include<ctime>
#define random (double)rand()/(RAND_MAX +1.0)
	
	using namespace std;

main(){
	
	srand(time(NULL));
	//produce all the M possible microstates of spins
	const int M = 65536;
	vector< vector<int> > microstate(M);
	for(int i = 0; i < M; i++){
		bitset <16> bin(i);
		for(int j = 0; j < bin.size(); j++)
			microstate[i].push_back(pow(-1,int(bin[j])));
	}
	ofstream outfile("problem1a.dat");
	for(double T = 0.2; T < 5.0; T += 0.01){
		int m,n, k = 4;
		double h = 0;
		double h2 = 0;
		double part =0;
		for(int i = 0; i < M; i++){
			int site[k][k];
			int sum = 0;	
			for(int j = 0; j < k*k; j++){
				m = int(j/k);
				n = int(j%k);
				site[m][n] = microstate[i][j];
			}
			for(int p = 0; p < k; p++){
				for(int q = 0; q < k; q++){
					p += k;
					q += k; 
					sum += (site[p%k][(q-1)%k]+site[p%k][(q+1)%k]+site[(p-1)%k][q%k]+site[(p+1)%k][q%k])*site[p%k][q%k];	
					p -= k;
					q -= k;
				}
			}			
			h += -0.5*(double)sum*exp((double)sum/(2*T)); 
			h2+= 0.25*(double)(sum*sum)*exp((double)sum/(2*T));
			part += exp((double)sum/(2*T));
			//cout << h <<"\t" << h2 << endl;
		}
		h /= part;
		h2 /= part;
		outfile << T << "\t " << h/(k*k) << "\t" << (h2-h*h)/(T*T*k*k) << endl;
	}
	outfile.close();	

}
