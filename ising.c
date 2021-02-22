#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define	nsp 25
#define mcstp 1000000
#define	J 1
#define temp 0.1

double *inicializa()
{	int i, j;
	double aleat, *spin;
	
	spin=(double*)malloc(nsp*sizeof(double));
	
	for(i=0;i<nsp;i++)
	{	aleat=rand()%101;
		aleat=aleat/100;
		//printf("%lf\n" ,aleat);
		if(aleat>0.5)
			spin[i]=1.0;
		else
			spin[i]=-1.0;
	}
	/*for(i=0;i<nsp;i++)
		printf("%lf \t %d \n",spin[i],i);
	*/
		
	return(spin);
}
double energia(double *spin, double ene)
{	int i;
	
	ene = 0.;
	for(i=0;i<nsp-1;i++)
	{	ene = ene - J*spin[i]*spin[i+1]*spin[i-1];
	}
	ene = ene - spin[0]*spin[i];  //Condicao de contorno Periodica
	//printf("%lf \n", ene);
	
	return(ene);
}


void main()
{	int i, j, k;
	double ene, ene1, ene2, de, mag, expmet, *spin, *enmc, *mgmc, aleat;
	FILE *fp;
	
	fp=fopen("energia.dat", "w");
	
	//alocando memoria para enmc, mgmc
	enmc=(double*)malloc(mcstp*sizeof(double));
	mgmc=(double*)malloc(mcstp*sizeof(double));
	
	spin=inicializa();
	ene=energia(spin,ene);
	ene1=ene;
	
	
	//Loop sob os passos de MC
	for(i=0;i<mcstp;i++)
	{	for(j=0;j<nsp;j++)
		{	//Flipando o spin da particula
			spin[j]=-spin[j];
			ene=energia(spin,ene);
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
				{	spin[j]=-spin[j]; //Troca rejeitada
				}
				else
				{	ene1=ene2; //Troca aceita
				}
				
			}
			//printf("ok2");
			
		}
		//printf("ok3");
		//Calculo da Magnetizacao
		mag=0.; 
		for(j=0;j<nsp;j++)
		{	mag= mag + spin[j];
		}
		mgmc[i] = mag/(double)nsp;
		enmc[i]= ene1/(double)nsp; //ver depois
			
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

	printf("Energia final média = %lf \n", ene);
	printf("Magnetização final média = %lf \n", mag); 
	
	
	fclose(fp);	
}

