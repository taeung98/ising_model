#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define N	7
#define Ns (1<<N)

int save_txt(char*,double*,int);

typedef struct Ising {
    double radius;
    double latvec[2][2];
    double supcoeff[2][2];
	double supvec[2][2];
	double subs[N][2];
	int nbs[N][6];
    int basis[Ns][N];
} Ising;

void find_subs(Ising* self){
	double lim[3] = {-1,0,1};
	int count = 0;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			double r[2] = {lim[i] * self->latvec[0][0] + lim[j] * self->latvec[1][0],
						lim[i] * self->latvec[0][1] + lim[j] * self->latvec[1][1]};
			if(sqrt(r[0]*r[0] + r[1]*r[1]) < self->radius + 1e-6){
				self ->subs[count][0] = r[0];
				self ->subs[count][1] = r[1];
				count++;
			}
		}
	}
}

void find_nbs(Ising* self){
	for(int i =0;i<N; i++){
		int count = 0;
		for(int j =0;j<N;j++){
			if(i==j) continue;
			for(int n1=-1;n1<=1;n1++){
				for(int n2=-1;n2<=1;n2++){
					double R[2] = {n1*self->supvec[0][0] + n2*self->supvec[1][0] + self->subs[j][0] - self->subs[i][0],
						n1*self->supvec[0][1] + n2*self->supvec[1][1] + self->subs[j][1] - self->subs[i][1]};
					if(fabs(sqrt(R[0]*R[0] + R[1]*R[1]) - self->radius) < 1e-6){
						self->nbs[i][count] = j;
						count++;
					}
				}
			}
		}
		if(count != 6){
		printf("6 neighbors are expected but %d obtained\n",count);
		}
	}
}

void mkbasis(Ising* self){
	for(int i=0;i<pow(2,N);i++){
		for(int j=0;j<N;j++){
			self->basis[i][j] = 2 * ((i>>j) & 0x01)-1;
		}
	}
}

double H(Ising* self, int* state){
	double energy = 0;
	for(int i=0;i<N;i++){
		int* nb = self->nbs[i];
		double sum = 0;
		for(int j=0;j<6;j++){
			sum += state[nb[j]];
		}
		energy -= state[i] * sum;
	}
	return energy/2;
}

double monte(Ising* self,int* state,double t){
	int s = rand() % N;
	double old_E,new_E,de;
	old_E = H(self,state);
	state[s] *= (-1);
	new_E = H(self,state);
	de = new_E - old_E;
	if (de < 0){
		return new_E;
	}
	else if( (rand()/(double)RAND_MAX) < exp(-de/t) ){
		return new_E;
	}
	else {
		state[s] *= (-1);
		return old_E;
	}
}
		
int main() {
	Ising ising;
	ising.radius = 1.0;
	double latvec[2][2] = { {1.0,0.0}, {0.5,0.5*sqrt(3.0)} };
	double supcoeff[2][2] = { {1,2} , {-3,1} };
	memcpy(ising.latvec, latvec, sizeof(latvec));
	memcpy(ising.supcoeff, supcoeff, sizeof(supcoeff));
	
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<2;k++){
				ising.supvec[i][j] += supcoeff[i][k] * latvec[k][j];
			}	
		}
	}
	find_subs(&ising);
	find_nbs(&ising);
	mkbasis(&ising);

#if 0
	for(int i=0;i<N;i++){
		for(int j=0;j<2;j++){
			printf("sub[%d][%d] = %f ",i,j,ising.subs[i][j]);
		}
		printf("\n");
	}

	for(int i=0;i<N;i++){
		for(int j=0;j<6;j++){
			printf("nbs[%d][%d] = %d ", i,j,ising.nbs[i][j]);
		}
		printf("\n");
	}
	
	for(int i=0;i<Ns;i++){
		for(int j=0;j<N;j++){
			printf("%d ",ising.basis[i][j]);
		}
		printf("\n");
	}
#endif
	
	double T[50];
	for(int i=49;i>=0;i--){
		T[i] = 0.1 + 4.9 / 49 * i;
//		printf("%f\n",T[i]);
	}
	int N_MC = 10000000;
	double Energy[50];
	srand(time(NULL));	
	for(int i=0;i<50;i++){
		double E = 0;
		int R = (int)( rand() % Ns );
//		printf("%d\n",R);
		int* state = ising.basis[ R ];
		for(int j=0;j<N_MC;j++){
			E += monte(&ising, state, T[i]);
		}
		E /= N_MC;
		Energy[i] = E;
//		printf("%f ",Energy[i]);
//		printf("\n");
	}
	
	save_txt("E_MN", Energy, N_MC);

	return 0;
}

	
int save_txt(char* head, double* arr, int N_MC){
	char fname[1024];
	sprintf(fname, "%s_%d.txt", head, N_MC);
	FILE* fp = fopen(fname,"w");
	for(int i=0;i<50;i++){
		fprintf(fp,"%lf\n",arr[i]);
	}
	fclose(fp);
	
	return 0;
}
#if 0

void print_array(double arr[], int n) {
    printf("[ ");
    for (int i = 0; i < n; i++) {
        printf("%.3f ", arr[i]);
    }
    printf("]\n");
}

void print_matrix(double mat[][2], int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.3f ", mat[i][j]);
        }
        printf("\n");
    }
}

void print_int_matrix(int mat[][6], int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", mat[i][j]);
        }
        printf("\n");
    }
}
#endif
