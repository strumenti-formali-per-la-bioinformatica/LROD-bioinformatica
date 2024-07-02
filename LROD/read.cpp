#ifndef Read_CPP_INCLUDED 
#define Read_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "read.h"

using namespace std;

long int max (long int x, long int y)
{
	long int z;
	if (x > y) z = x;
	else z = y;
	return (z);
}

long int min (long int x, long int y)
{
	long int z;
	if (x < y) z = x;
	else z = y;
	return (z);
}

//@param StrLine -> array di size maxsize
//@param maxSize -> taglia dell'array (1000000)
//@param fp -> puntatore al file che contiene la read (long_read.fa)
ReadSetHead * GetReadSetHead(char *filename,char *StrLine, long int maxSize){
	
	printf("\nfunzione GetReadSetHead\n");
	ReadSetHead * readSetHead = (ReadSetHead * )malloc(sizeof(ReadSetHead));
	
	FILE *fp; 
	long int m=1;
    if((fp = fopen(filename,"r")) == NULL){ 
        printf("error!"); 
        return NULL; 
    } 
	
    readSetHead->readCount = 0;

	// conta il nuemero delle reads
	while((fgets(StrLine, maxSize, fp)) != NULL){	//legge maxSize caratteri in ogni ciclo
		if(StrLine[0]=='>'){ //se il primo carattere è > -> è una read -> incrementa il contantore
			readSetHead->readCount++;
		}
	}
	fclose(fp);
	
	printf("Number of long reads: %ld;\n",readSetHead->readCount);
	//alloca un readSet per ogni read (contante al passo precedente)
	readSetHead->readSet = (ReadSet *)malloc(sizeof(ReadSet)*readSetHead->readCount);
	for(long int i = 0; i <readSetHead->readCount; i++){
		readSetHead->readSet[i].readLength = 0;
	}//inizializza a 0 la lunghezza di ogni read
	
	FILE *fp1; //apre il file della long read
	if((fp1 = fopen(filename,"r")) == NULL){ 
        printf("error!"); 
        return NULL; 
    } 
	long int readIndex = -1;
	//long int j = 0;
	long int allReadLen = 0;

	//printf("stampa della lista di reads\n");
	while((fgets(StrLine, maxSize, fp1)) != NULL){

		if(StrLine[0]=='>'){
			readIndex++;
			continue;
		}
		//elimina il carattere di terminazione /n o /r 
		readSetHead->readSet[readIndex].readLength = strlen(StrLine); //salva lunghezza della read
		if(StrLine[readSetHead->readSet[readIndex].readLength - 1]=='\n' || StrLine[readSetHead->readSet[readIndex].readLength - 1]=='\r'){
			readSetHead->readSet[readIndex].readLength--;
		}
		
		readSetHead->readSet[readIndex].read = (char *)malloc(sizeof(char)*(readSetHead->readSet[readIndex].readLength + 1));
		strncpy(readSetHead->readSet[readIndex].read, StrLine, readSetHead->readSet[readIndex].readLength);
		readSetHead->readSet[readIndex].read[readSetHead->readSet[readIndex].readLength] = '\0';
	}
	return readSetHead;
}


#endif