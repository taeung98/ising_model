#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#define L (2)  // LxL격자

struct ising {	
	double E;       // <E>
	double E2;	//2<-square
	double M;	//<M>
        double absM;	// abs <- absolute
	double M2;	// 2<- square
	double C;	// capacity
	double K;	// Kai <- susceptability
	double Kp;	// Kai' <- susceptability'
	double M4;	// for Curie temperature	
};


double Energy(double arr[][L]) {

        double e = 0;
        for (int i =0; i<L; i++)
       		for (int j =0; j< L; j++)
                	e += -arr[i][j] * (arr[(i+1+L)%L][j]+arr[i][(j+1+L)%L]);
        return e;
}


double Magnetization(double arr[][L]) {
        double m = 0;
        for (int x =0; x<L; x++)
                for (int y =0; y<L; y++)
                	m += arr[x][y];
        return m;
}

int monte(struct ising *is, double arr[][L], double t, double n) //논문에서 주어진 알고리즘 구현
{
    double E=0,M=0,absM=0,dE=0,m2=0,m4=0,c=0,k=0,ct=0,e2=0,E_i=0,M_i=0;
	E = Energy(arr);
	M = Magnetization(arr);
 
    int x,y;
        
	srand(time(NULL));
	E_i = Energy(arr);
    M_i = Magnetization(arr);	
	
	for (int i =0; i<n*(L*L); i++){
        x= rand()%L;
        y= rand()%L;
		
        dE = 2*(arr[x][y])*(arr[(x-1+L)%L][y]+arr[x][(y-1+L)%L]+arr[(x+1+L)%L][y]+arr[x][(y+1+L)%L]);

		if (dE <= 0 || (rand()/(double)RAND_MAX) < exp(-dE/t) ){
            arr[x][y] *= (-1);
			E_i += dE;
			M_i += 2*arr[x][y];
		}

		M += M_i;
		E += E_i;
		m2 += pow(M_i,2);
		e2 += pow(Energy(arr),2);
        absM += fabs(M_i);
		m4 += pow(M_i,4);
     }

    E /= (n*(L*L));
	M /= (n*(L*L));
	absM /= (n*(L*L));
	m2 /= (n*(L*L)); 
	e2 /= (n*(L*L));
	m4 /= (n*(L*L));

	is->E = E;
	is->M = M;
	is->absM = absM;	
	is->M2 = m2;
	is->E2 = e2;
	is->M4 = m4;

    return 0;
}

int save_txt(char *head, double *arr){
	char fname[1024];
	sprintf(fname, "%s%d.txt", head, L);
	FILE *fp = fopen(fname, "w");
	for(int i=0; i<46; i++){
		fprintf(fp, "%lf\n", arr[i]);
	}

	fclose(fp);

	return 0;
}

int main()
{

	clock_t start = clock();

	struct ising i1;

    double n = 25000;                 // mcs/(L^2)

    double arr[L][L] = {0,};          // preparation of initial spin_arr

    double T;
    double Energy_per_spin[46];
	double Magnetization_per_spin[46];
	double absMagnetization_per_spin[46];
	double squrEnergy_per_spin[46];
	double squrMagnetization_per_spin[46];
	double capacity[46];
	double susceptability[46];
	double susceptabilityp[46];
	double Curie_temperature[46];

    for (int t=5; t<51; t++) {
		for (int i =0; i<L;i++)                 // initialize spin_arr
			for (int j=0; j<L; j++)
				arr[i][j]=1;
        
		T = t*0.1;
        monte(&i1, arr,T,n);
		
		Energy_per_spin[t-5] = i1.E/(L*L);
       	Magnetization_per_spin[t-5] = i1.M/(L*L);		
		absMagnetization_per_spin[t-5] = i1.absM/(L*L); 		
		squrEnergy_per_spin[t-5] = i1.E2/(L*L);
		squrMagnetization_per_spin[t-5] = i1.M2/pow(L,4);	// <M>^2 하고 차원 맞춰야되나까 L^4
		capacity[t-5] = (i1.E2-pow(i1.E,2))/((T*T)*(L*L));		
		susceptability[t-5] = (i1.M2-pow(i1.M,2))/(T*(L*L));
		susceptabilityp[t-5] = (i1.M2-pow(i1.absM,2))/(T*(L*L));		
       	Curie_temperature[t-5] = 1 - i1.M4/(3*pow(i1.M2,2));
       	}
	
	save_txt("energy",Energy_per_spin);
/*	save_txt("absm",absMagnetization_per_spin);
	save_txt("C",capacity);
	save_txt("M2_",squrMagnetization_per_spin);
	save_txt("M",Magnetization_per_spin);
	save_txt("K",susceptability);
	save_txt("Kp",susceptabilityp);
	save_txt("UL",Curie_temperature);
*/
	clock_t end = clock();

	printf("Time: %lf\n", (double)(end - start)/CLOCKS_PER_SEC);

	return 0;
}
