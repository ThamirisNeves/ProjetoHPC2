#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<mpi.h>
#define	nsp 16 
#define mcstp 10000
#define	J 1
#define temp 0.5
//variaveis globais
int NP,myid;

//void init_it(int *argc, char **argv)
//{	
//	MPI_Init(argc,&argv); //inicia o ambiente paralelo
//     	MPI_Comm_size(MPI_COMM_WORLD, &NP); //retorna o numero de processos
//    	MPI_Comm_rank(MPI_COMM_WORLD, &myid); //retorna o id dos processos
//	
//	return;
//}

void inicializa(double **spin)
{	int i, j;
	double aleat;
	
	for(i=0;i<nsp;i++)
	{	for(j=0;j<nsp+2;j++)
		{	//aleat=rand()%101;
			//aleat=aleat/100;
			//printf("%lf\n" ,aleat);
			//if(aleat>0.5)
				spin[i][j]=1.0;
			//else
				//spin[i][j]=-1.0;
		}
		
	}
		
	return;
}
void imprimematriz(double **spin)
{	int i, j;

	for(i=0;i<nsp;i++)
        {       for(j=0;j<nsp+2;j++)
                {       printf("%lf \t ",spin[i][j]);
                }
                printf("\n");
        }
        printf("\n\n");	
	return;
}
void condicaoborda(double **spin, int l)
{       int i;

        for(i=1;i<l;i++)
        {       spin[i][0]=spin[i][nsp];
                spin[i][nsp+1]=spin[i][1];
        }

        return;

}
void condicaocore(double **matriz, double *vect, double *vecr, int n, int next, int preview, MPI_Status status)
{	int i, j, tag;
	tag=1234;
	
	MPI_Send(vect,nsp+2,MPI_DOUBLE,next,tag,MPI_COMM_WORLD);
        MPI_Recv(vecr,nsp+2,MPI_DOUBLE,preview,tag,MPI_COMM_WORLD,&status);
	
	j=0; //a linha recebida entra na primeira linha
	for(i=0;i<nsp+2;i++)
	{	matriz[j][i]=vecr[i];
	}

	return;
}
double energia(double **spin, int n)
{	int i, j;
	double ene;
	
	ene = 0.;
	for(i=1;i<n+1;i++)
	{	for(j=1;j<nsp;j++)
			ene = ene - J*spin[i][j]*(spin[i+1][j] + spin[i][j+1]);
	}
	//printf("%lf \n", ene);
	
	return(ene);
}
void distribui(double **matriz, double **spin, int n)
{	int i, j;
	
	for(i=1;i<n+1;i++)
	{	for(j=0;j<nsp+2;j++)
			matriz[i][j]=spin[myid*n+i-1][j];
	}
	return;
}
void main(int argc, char *argv[])
{	int i, j, k, id, next, preview, n=nsp/NP;
	double ene, ene1, ene2, enet, de, mag, expmet, **spin, *enmc, *mgmc, aleat, **matriz, *vect, *vecr, redmag=0, redene=0;
	FILE *fp, *ft, *fr;
	clock_t inicio, fim;
	MPI_Status status;

	inicio=clock();

	//init_it(&argc,argv);
	MPI_Init(&argc,&argv); //inicia o ambiente paralelo
        MPI_Comm_size(MPI_COMM_WORLD, &NP); //retorna o numero de processos
        MPI_Comm_rank(MPI_COMM_WORLD, &myid); //retorna o id dos processos

	next=(myid+1)%NP;
        preview=(myid+NP-1)%NP;

	if(myid==0)
	{	fp=fopen("energia.dat", "w");
		ft=fopen("tempo.dat", "w");
		fr=fopen("media.dat","w");
	}	

	//alocando memoria para spin, enmc, mgmc, matriz
	spin=(double**)malloc(nsp*sizeof(double*));
        for(i=0;i<nsp;i++)
                spin[i]=(double*)malloc((nsp+2)*sizeof(double));
	vect=(double*)malloc((nsp+2)*sizeof(double));
	vecr=(double*)malloc((nsp+2)*sizeof(double));
	enmc=(double*)malloc(mcstp*sizeof(double));
	mgmc=(double*)malloc(mcstp*sizeof(double));
	matriz=(double**)malloc((n+1)*sizeof(double*));
        for(i=0;i<n+1;i++)
                matriz[i]=(double*)malloc((nsp+2)*sizeof(double));

	inicializa(spin);
	distribui(matriz,spin,n);
	condicaocore(matriz,vect,vecr,n,next,preview,status);
	//imprimematriz(spin, nsp+2);
	ene=energia(matriz, n);	
	ene1=ene;
	
	//Loop sob os passos de MC
	for(i=0;i<mcstp;i++)
	{	for(j=0;j<nsp;j++)
		{	for(k=0;k<nsp;k++)
			{	//Flipando o spin da particula
				//ene1=ene;
				matriz[j][k]=-matriz[j][k];
				condicaoborda(matriz,n+1);
				condicaocore(matriz,vect,vecr,n,next,preview,status);
				ene=energia(matriz, n); //so vou calcular a energia nos processos
				//printf("%d \t %lf\n", j, ene);
				ene2=ene;
				de=ene2-ene1;
				//printf("ok1");
				if(de<=0.)
				{	ene1=ene2;	 
				}
				else
				{	//Algoritmo de Metropolis
					aleat=rand()%101;
					aleat=aleat/100;
					expmet=exp(-de/temp);
					//printf("%lf \t %lf\n", aleat, expmet);
					if(aleat>expmet)
					{	matriz[j][k]=-matriz[j][k];
	                	                condicaoborda(matriz,n+1);
        		                        condicaocore(matriz,vect,vecr,n,next,preview,status);
					}
					else
					{	ene1=ene2; //Troca aceita
					}
				
				}
				//printf("ok2");
			
			}
		}
		//printf("ok3");
		//Calculo da Magnetizacao
		redmag=0.; 
		redene=0.;
		if(myid==0)
		{	MPI_Reduce(&mag, &redmag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&ene1, &redene, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		//for(j=1;j<n+1;j++)
		//{	for(k=1;k<nsp+1;k++)
		//		mag= mag + spin[j][k];
		//}

			mgmc[i] = mag/(n*nsp);
			enmc[i]= ene1/(n*nsp); 
			fprintf(fp, "%d \t %lf \t %lf \n", i, mgmc[i], enmc[i]);
		}
		
	}
	//printf("ok4");
	// Calculo da magnetizacao e energia media:
	if(myid==0)
	{	mag = 0.;
		ene = 0.;

		for(i=0;i<mcstp;i++)
		{ 	mag = mag + mgmc[i];
			ene = ene + enmc[i];
		}
	
		ene = ene/(double)mcstp;
		mag = mag/(double)mcstp;
	
	//imprimematriz(spin);
		printf("Energia final média = %lf \n", ene);
		printf("Magnetização final média = %lf \n", mag); 
		fprintf(fr, "%lf \t %lf \t %lf \n", temp, ene, mag);
	
		fim=clock();	
		printf("Tempo gasto: %lf\n\n",((fim-inicio)/(double)CLOCKS_PER_SEC));
		fprintf(ft,"%lf \n", ((fim-inicio)/(double)CLOCKS_PER_SEC));
		fclose(ft);
		fclose(fp);
	}

	MPI_Finalize();
}

