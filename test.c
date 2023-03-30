#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N		7

typedef struct {
    double x;
    double y;
} Vec;

typedef struct {
    double radius;
    Vec latvec[2];
    Vec supvec;
    Vec subs[N];
    int nbs[N][6];
    int basis[128][N];
//    int N;
} Ising;

double norm(Vec vec) {
    return sqrt(vec.x*vec.x + vec.y*vec.y);
}

Vec add(Vec vec1, Vec vec2) {
    Vec vec;
    vec.x = vec1.x + vec2.x;
    vec.y = vec1.y + vec2.y;
    return vec;
}

Vec sub(Vec vec1, Vec vec2) {
    Vec vec;
    vec.x = vec1.x - vec2.x;
    vec.y = vec1.y - vec2.y;
    return vec;
}

Vec scale(double scalar, Vec vec) {
    Vec result;
    result.x = scalar * vec.x;
    result.y = scalar * vec.y;
    return result;
}

double dot(Vec vec1, Vec vec2) {
    return vec1.x * vec2.x + vec1.y * vec2.y;
}

Ising* Ising_new(double radius, Vec* latvec, Vec svec) {
    Ising* self = malloc(sizeof(Ising));
    self->radius = radius;
    self->latvec[0] = latvec[0];
    self->latvec[1] = latvec[1];
    self->supvec = add(scale(svec.x, latvec[0]), scale(svec.y, latvec[1]));
    int n = 0;
    for (int n1=-1; n1<=1; n1++) {
        for (int n2=-1; n2<=1; n2++) {
            Vec r = add(scale(n1, latvec[0]), scale(n2, latvec[1]));
            if (norm(r) < radius + 1e-6) {
                self->subs[n] = r;
                n++;
            }
        }
    }
     
    for (int i=0; i<N; i++) {
        int k = 0;
        for (int N1=-1; N1<=1; N1++) {
            for (int N2=-1; N2<=1; N2++) {
                for (int j=0; j<N; j++) {
                    Vec R = add(add(scale(N1, self->supvec), self->subs[j]), scale(-1, self->subs[i]));
                    if (fabs(norm(R) - 1) < 1e-6) {
                        self->nbs[i][k] = j;
                        k++;
                    }
                }
            }
        }
        if (k != 6) {
            printf("6 neighbors are expected but %d obtained\n", k);
            return NULL;
        }
    }
    for (int i=0; i<128; i++) {
        int temp = i;
        for (int j=0; j<N; j++) {
            self->basis[i][j] = 2*(temp % 2) - 1;
            temp = temp >> 1;
        }
    }
    return self;
}

void Ising_free(Ising* self) {
    free(self);
}

double Ising_H(Ising* self, int* state) {
    double energy = 0;
	for(int i=0; i<N; i++) {
		int* nb = self->nbs[i];
		double spin_energy = 0;
		for(int j=0; j<6; j++) {
			spin_energy += state[nb[j]];
		}
		energy -= state[i]*spin_energy;
	}
	return energy/2.0;
}

	int Ising_M(Ising* self, int* state) {
	int M = 0;
	for(int i=0; i<N; i++) {
	M += state[i];
	}
	return M;
	}

int Ising_AM(Ising* self, int* state) {
	int M = 0;
	for(int i=0; i<N; i++) {
	M += abs(state[i]);
	}
	return M;
}

void Ising_compute_expect(Ising* self, double (func)(Ising*,int*)) {
	double T[50];
	double ens[N];
	double myexp[1<<N];
	double Eexp[50];
	for(int i=0; i<50; i++) {
	T[i] = 0.1 + 0.09*i;
	}
	for(int i=0; i<(1<<N); i++) {
		int state[N];
		for(int j=0; j<N; j++) {
			state[j] = 2*((i>>j)&0x01) - 1;
		}
		myexp[i] = func(self, state);
		ens[i] = Ising_H(self, state);
	}

	for(int i=0; i<50; i++) {
		double t = T[i];
		double Z = 0;
		for(int j=0; j<(1<<N); j++) {
			Z += exp(-ens[j]/t);
		}

		Eexp[i] = 0;
		for(int j=0; j<(1<<N); j++) {
			double ep = myexp[j];
			double en = ens[j];
			Eexp[i] += ep*exp(-en/t)/Z;
		}
	}

	for(int i=0; i<50; i++) {
		printf("%.4lf\t%.4lf\n", T[i], Eexp[i]);
	}
}

double Ising_monte(Ising* self, int* state, double t) {
	int s = rand() % N;
	double Ep = Ising_H(self, state);
	state[s] *= -1;
	double Ea = Ising_H(self, state);
	double de = Ep - Ea;
	if(de < 0 || exp(-de/t) > (double)rand()/RAND_MAX) {
		return Ep+de;
	} 
	else {
		state[s] *= -1;
		return Ep;
	}
}

int main() {
	Ising tri;	
	int nbs = (int)malloc(N * 6 * sizeof(int));
	for (int i=0; i<(1<<N); i++) {
		for (int N1=-1; N1<=1; N1++) {
			for (int N2=-1; N2<=1; N2++) {
				for (int j=0; j<N; j++) {
					double R[2];
					R[0] = N1 * tri.supvec[0] + N2 * tri.supvec[1] + tri.subs[j][0] - tri.subs[i][0];
					R[1] = N1 * tri.supvec[0] + N2 * tri.supvec[1] + tri.subs[j][1] - tri.subs[i][1];
					double dist = sqrt(R[0]*R[0] + R[1]*R[1]);
					if (fabs(dist - 1) < 1e-6) {
						tri.nbs[i6 + j] = j;
					}
				}
			}
		}
	}

	for (int i=0; i<(1<<N); i++) {
			for (int j=0; j<N; j++) {
				tri.basis[i][j] = 2 * ((i >> j) & 0x01) - 1;
			}
	}

	double ens = (double)malloc((1<<N) * sizeof(double));
	for (int i=0; i<(1<<N); i++) {
			double energy = 0;
			for (int j=0; j<N; j++) {
				energy -= basis[i][j] * np_sum(basis[i], nbs+j*6, 6);
			}
		ens[i] = energy / 2;
	}

	double myexp = (double)malloc((1<<N) * sizeof(double));
	for (int i=0; i<(1<<N); i++) {
			myexp[i] = tri.M(basis[i]);
	}

	double Eexp = (double)malloc(50 * sizeof(double));
	double T = (double)malloc(50 * sizeof(double));
	for (int i=0; i<50; i++) {
		T[i] = 0.1 + 4.9 / 49 * i;
		double Z = 0;
		for (int j=0; j<(1<<N); j++) {
			Z += exp(-ens[j] / T[i]);
		}	
		double E = 0;
		for (int j=0; j<(1<<N); j++) {
				E += myexp[j] * exp(-ens[j] / T[i]) / Z;
		}
		Eexp[i] = E;
	}

	return 0;
}
