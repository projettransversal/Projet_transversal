#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <unistd.h>

#define LOG2 10

#define SWAP(a,b) ctmp=(a); (a)=(b); (b)=ctmp
#define ARRAY_LEN(x) ((int) (sizeof (x) / sizeof (x [0])))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

// Taille de lecture, nombre de points
#define BUFFER_LEN 1024
// Nombre de bits qui sont traités par unité de temps qui est de 128kbits en moyenne
#define BITRATE 512

/*
Commande de compilation : gcc -std=c99 fftread.c -lm -lsndfile -o fftread.o
*/


/* 
Pour utiliser les modules permettant l'importation des fichiers son :
*/
SNDFILE *infile, *outfile;
SF_INFO sfinfo;


// fonction sfx_mix_mono_read_double
sf_count_t sfx_mix_mono_read_double (SNDFILE * file, double * data, sf_count_t datalen){
  SF_INFO info;
  static double multi_data [2048] ;
  int k, ch, frames_read ;
  sf_count_t dataout = 0 ;
  sf_command (file, SFC_GET_CURRENT_SF_INFO, &info, sizeof (info)) ;
  if (info.channels == 1)
    	return sf_read_double (file, data, datalen) ;
  while (dataout < datalen){   
   		int this_read ;
		this_read = MIN (ARRAY_LEN (multi_data) / info.channels, datalen) ;
		frames_read = sf_readf_double (file, multi_data, this_read) ;
		if (frames_read == 0)
				break ;
		for (k = 0 ; k < frames_read ; k++){       
				double mix = 0.0 ;
				for (ch = 0 ; ch < info.channels ; ch++)
						mix += multi_data [k * info.channels + ch] ;
		        data [dataout + k] = mix / info.channels ;
        } ;
    dataout += this_read ;
    } ;
  	return dataout ;
}

// fonction bitrev 
int bitrev(int inp, int numbits){
  int rev=0;
  for (int i=0; i<numbits;i++){
      rev=(rev << 1) | (inp & 1);
      inp >>=1;
    }
  return rev;
}


/************* Programmation des FFT et DFT  *************/

/* Rôle : Calcul des Twiddles Factors */
complex TW[BUFFER_LEN];
void twiddle(complex *TW, unsigned int size){
	complex phi = cexp(-2*I*3.14150/size);

	TW[0]=1;
	for(int i=1; i<size; i++)
		TW[i]=TW[i-1]*phi;
}

/* Rôle : Calcul de la FFT itérative */
void fftiterTW(complex *data, unsigned int size, int log2n){
	int j, N2, Bpair, Bimpair, Bp=1, N=size;
	complex impair, pair, ctmp;

	for(int k=0; k<log2n;k++){
		N2=N/2;
		Bpair=0;
		for(int b=0; b<Bp;b++){
			Bimpair=Bpair+N2;
			for(int n=0;n<N2;n++){
				impair = data[Bpair+n] + data[Bimpair+n];
				pair = (data[Bpair+n] - data[Bimpair+n])*TW[n*size/N];
				data[Bpair+n] = pair;
				data[Bimpair+n] = impair;
			}
			Bpair = Bpair+N;
		}
		Bp=Bp*2;
		N=N2;
	}
	for(int i=0;i<size;i++){
		j=bitrev(i,log2n);
		if(j>i){
			SWAP(data[j],data[i]);
		}
	}
	for(int i=size-1;i>0;i--){
		data[i]=data[i-1];
	}
	data[0]=ctmp;
	return;
}

/* Rôle : Calcul de la FFT récursive */
void fftrec(complex *data, complex *result, unsigned int size, int log2n){
	complex ypair[size], yimpair[size], Fimpair[size], Fpair[size];
	int n,N2;
	if(size>1){
		N2=size/2;
		for(n=0;n<N2;n++){
			ypair[n] = data[n]+data[n+N2];
			yimpair[n] = (data[n]-data[n+N2])*cexp(-2*I*3.14150*n/size);
		}
		fftrec(ypair,Fpair,N2,log2n);
		fftrec(yimpair,Fimpair,N2,log2n);
		for(n=0;n<N2;n++){
			result[2*n]=Fpair[n];
			result[2*n+1]=Fimpair[n];
		}
	}
	else{
		result[0]=data[0];
		return;
	}

}

/* Rôle : Calcul de la DFT du signal */
void DFT(complex *data, complex *result, unsigned int buffer){
	for(int k=0; k<BUFFER_LEN;++k){
		for(int n=0; n<BUFFER_LEN;++n){
			result[k]=result[k]+data[n]*cexp(-2*I*3.14150*k*n/buffer);
		}
	}
}

/*************** fonction utiles ***************/

/* Rôle : Création d'un nombre complexe */
complex conversioncomplexe(double valeur){
	return valeur + I*0.0;
}

/* Rôle : Calcul du module d'un nombre complexe */
double module_complexe(complex c){
	return sqrt(creal(c)*creal(c)+cimag(c)*cimag(c));
}

/* Rôle : Effectue un delai avec le temps t passé en parametre */
void delais(int t){
	clock_t time =clock();
	while( (clock()-time) <= (t*CLOCKS_PER_SEC/1000));
}


/*************** FIN fonction utiles ***************/

/*************** Fonction time de la fft ***************/

int time_fft(){
	sf_count_t readcount;
	clock_t start, stop;
	double data[BUFFER_LEN];
	double tempsDFT, tempsFFTiterTW, tempsFFTrec;
	double complex resultDFT[BUFFER_LEN];
	double complex datac[BUFFER_LEN];
	double complex spectre[BUFFER_LEN];

	infile = sf_open("lapur.wav",SFM_READ, &sfinfo);
	if(infile == NULL){
		printf("Not available to open input file \n");
		sf_perror(NULL);
		return 1; 
	}

	while((readcount=sfx_mix_mono_read_double(infile,data,BUFFER_LEN))>0){
		for(int i=0; i<BUFFER_LEN; i++){
			datac[i] = conversioncomplexe(data[i]);
		}
		start=clock();
		twiddle(TW,BUFFER_LEN);
		fftiterTW(datac,BUFFER_LEN,LOG2);
		stop = clock();

		tempsFFTiterTW = (double)(stop-start)/(double)(CLOCKS_PER_SEC);

		start=clock();
		fftrec(datac,spectre,BUFFER_LEN,LOG2);
		stop = clock();

		tempsFFTrec = (double)(stop-start)/(double)(CLOCKS_PER_SEC);

		start=clock();
		DFT(datac,resultDFT,BUFFER_LEN);
		stop = clock();

		tempsDFT = (double)(stop-start)/(double)(CLOCKS_PER_SEC);
	}
	sf_close(infile);


	printf("Temps pour une FFT itérative avec Twiddle Factor : %lf ms\n",tempsFFTiterTW*1000);
	printf("Temps pour une FFT récursive : %lf ms\n",tempsFFTrec*1000);
	printf("Temps pour une DFT : %lf ms\n",tempsDFT*1000);
	return 0;


}
/*************** Fin Fonction time de la fft ***********/

/*************** Fonction affichage du spectre ***************/

void Spectre_affichage(double spectre[], double max){
	int hauteur = 30; // Hauteur de la fenetre voulu pour l'affichage
	int m = max; 
	while(hauteur >= 0){
		while(m<hauteur){
			printf("\n");
			hauteur--;
		}
		for(int i=0; i<BITRATE/2; i++){
			if(spectre[i] >= m){
				printf("*");
				spectre[i]=m-1;
			}
			else{
				printf(" ");
			} 
		}
		printf("\n");
		hauteur--;
		m--;
	}
}

int affichage_fft(){
	sf_count_t readcount;
	double data[BUFFER_LEN];
	double complex datac[BUFFER_LEN];
	double spectre[BITRATE];
	int maximun, min;

	

	twiddle(TW, BUFFER_LEN);
	infile = sf_open("Ca va mieux en le disant 24.11.2015.wav",SFM_READ, &sfinfo);
	if(infile==NULL){
		fprintf(stderr, "Not available to open input file \n");
		exit(1);
	}

	double samplerate = sfinfo.samplerate;
	double channels = sfinfo.channels;

	while((readcount=sfx_mix_mono_read_double(infile,data,BUFFER_LEN)) >0){
		for(int i=0; i<BUFFER_LEN; i++){
			datac[i] = conversioncomplexe(data[i]);
		}
	fftiterTW(datac,BUFFER_LEN,LOG2);

		int echantillon = 0; 
		double moyenne = 0;
		int place_spectre = 0;
		for(int i=0; i<BUFFER_LEN;i++){
			moyenne = moyenne + module_complexe(datac[i]);
			echantillon++;

			if(echantillon==BUFFER_LEN/BITRATE){
				spectre[place_spectre] = 20*log10(moyenne/echantillon);

				echantillon = 0;
				place_spectre++;	
				moyenne=0;
			}
			if(place_spectre==BITRATE){
				min=spectre[0];
				for(int j=0; j<BITRATE;j++){
					if(spectre[j]<min){

						min = spectre[j];
					}
				}
			}
		}
		double max =spectre[0];

		for(int j=0; j<BITRATE;j++){
			spectre[j]=spectre[j]-min;
			if(max<spectre[j]){
				max=spectre[j];
				maximun = j;
			}
		}	
		system("clear");
		Spectre_affichage(spectre,max);
		printf("%f\n",(maximun*samplerate)/128);
		delais(80);
	}
	sf_close(infile);
	return EXIT_SUCCESS;
}
int main(int argc, char*argv[]){
	//printf("%d",time_fft());
	printf("%d",affichage_fft());
	
}
