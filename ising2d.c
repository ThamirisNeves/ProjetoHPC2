#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#define	nsp 102 //para a condição periódica de contorno o número de spin deve ser n+2
#define mcstp 100000
#define	J 1
#define temp 0.1

void inicializa(double **spin)
{	int i, j;
	double aleat;
	
	for(i=1;i<nsp-1;i++)
	{	for(j=1;j<nsp-1;j++)
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
void condicaoperiodica(double **spin)
{	int i;

	for(i=1;i<nsp-1;i++)
        {       spin[0][i]=spin[nsp-2][i];
                spin[i][0]=spin[i][nsp-2];
                spin[nsp-1][i]=spin[1][i];
                spin[i][nsp-1]=spin[i][1];
        }
	
	return;	
	
}
void imprimematriz(double **spin)
{	int i, j;

	for(i=0;i<nsp;i++)
        {       for(j=0;j<nsp;j++)
                {       printf("%lf \t ",spin[i][j]);
                }
                printf("\n");
        }
        printf("\n\n");	
	return;
}
double energia(double **spin)
{	int i, j;
	double ene;
	
	ene = 0.;
	for(i=1;i<nsp-1;i++)
	{	for(j=1;j<nsp-1;j++)
			ene = ene - J*spin[i][j]*(spin[i+1][j] + spin[i][j+1]);
		//printf("%lf \n", ene);
	}
	//printf("%lf \n", ene);
	
	return(ene);
}

void main()
{	int i, j, k;
	double ene, ene1, ene2, enet, de, mag, expmet, **spin, *enmc, *mgmc, aleat;
	FILE *fp, *ft;
	clock_t inicio, fim;
		
	inicio=clock();

	fp=fopen("energia.dat", "w");
	ft=fopen("time.dat", "w");
	
	//alocando memoria para spin, enmc, mgmc
	spin=(double**)malloc(nsp*sizeof(double*));
        for(i=0;i<nsp;i++)
                spin[i]=(double*)malloc(nsp*sizeof(double));
	enmc=(double*)malloc(mcstp*sizeof(double));
	mgmc=(double*)malloc(mcstp*sizeof(double));
	
	for(j=0;j<nsp;j++)
	{	for(k=0;k<nsp;k++)
		{	spin[j][k]=0;
		}
	}
		
	inicializa(spin);
	condicaoperiodica(spin);
	//imprimematriz(spin);
	ene=energia(spin);	
	ene1=ene;
	
	//Loop sob os passos de MC
	for(i=0;i<mcstp;i++)
	{	for(j=1;j<nsp-1;j++)
		{	for(k=1;k<nsp-1;k++)
			{	//Flipando o spin da particula
				//ene1=ene;
				spin[j][k]=-spin[j][k];
				condicaoperiodica(spin);
				ene=energia(spin);
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
					{	spin[j][k]=-spin[j][k]; //Troca rejeitada
						condicaoperiodica(spin);
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
		mag=0.; 
		for(j=1;j<nsp-1;j++)
		{	for(k=1;k<nsp-1;k++)
				mag= mag + spin[j][k];
		}
		
		mgmc[i] = mag/((nsp-2)*(nsp-2));
		enmc[i]= ene1/((nsp-2)*(nsp-2)); //ver depois
			
		fprintf(fp, "%d \t %lf \t %lf \n", i, mgmc[i], enmc[i]);
		
	}
	//printf("ok4");
	// Calculo da magnetizacao e energia media:
	mag = 0.;
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
	
	
	fclose(fp);	
	fim=clock();

	printf("Tempo gasto: %lf\n\n",((fim-inicio)/(double)CLOCKS_PER_SEC));
	fprintf(ft,"%lf \n", ((fim-inicio)/(double)CLOCKS_PER_SEC));
	fclose(ft);
}

