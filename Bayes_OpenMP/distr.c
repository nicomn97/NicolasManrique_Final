#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    double dom;
    dom=0.3989422804014327;
    //Inicia el espacio en paralelo
    #pragma omp parallel
    {
        //A cada hilo le manda el comando
        int thread_id = omp_get_thread_num();
	int thread_count = omp_get_num_threads();
        int i;
        double cadena[1000];
        i=0;
        char fname[50];
        sprintf(fname, "data%d", thread_id);
        srand((unsigned int)time(NULL));
        double a;
        double b;
        double posV,posN;
        posV=0;
        FILE *fp;
        fp=fopen(fname,"w");
        fprintf(fp, "A character: %d\n", thread_id);
        fprintf(fp, "Rand: %lf\n", a);
        while(i<1000){
        a=((double)rand()/(double)(RAND_MAX));
        b=((double)rand()/(double)(RAND_MAX));
        posN=posV+a;
        double norm;
        norm=dom*exp(-0.5*posN*posN);
        if(a<norm){
        cadena[i]=posN;
        posV=posN;
        i+=1;
        }
        else if(a>b){
        cadena[i]=posN;
        posV=posN;
        i+=1;
        }
        }
        fclose(fp);
    }
    return 0;

}
