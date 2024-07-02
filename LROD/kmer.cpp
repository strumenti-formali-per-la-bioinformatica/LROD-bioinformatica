#ifndef Kmer_CPP_INCLUDED 
#define Kmer_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
//#include <stdint.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include<assert.h> 
#include <time.h>

#include "kmer.h"
#include "bitarray.h"

char * strup(char * str)     
{  
    assert(str);                  
    char *ret = str;              
    while(*str != '\0')           
    {     
           if((*str >= 'a')&&(*str <= 'z'))
        {  
            *str = *str -32;          
            str++;  
        }  
        else  
            str++;  
    }  
    return ret;             
}  
/*
	@param unsigned int key -> il kmer
	restituisce il valore hash del kmer sfruttando operazioni bitwise
*/
unsigned int hash32shift(unsigned int key) 
{ 
  key = ~key + (key << 15); 
  key = key ^ (key >> 12); 
  key = key + (key << 2); 
  key = key ^ (key >> 4); 
  key = key * 2057; 
  key = key ^ (key >> 16); 
  return key; 
}

/*
	@param unsigned int kmer -> il kmer visto come valori binari
	@param unsigned int max -> numero di kmer (che rientrano nell'intevllo min max) moltiplicato per 2. Allocation count
		rappresenta il valore massimo che il kmer puo assumere come indice dell'array
	
	@return il valore hash del kmer che servira da indice per memorizzare e trovare il kmer stesso 
	nota: %max ci assicura che l'indice sia compreso tra 0 e max-1 (max-1 è il numero di kmer allocati)
*/
unsigned int Hash(unsigned int kmer, unsigned int max)  
{  
	//per il primo kmer HASH INDEX: 102
	//per il secondo kmer HASH INDEX: 100
	//per il terzo kmer HASH INDEX: 0
	//per il quarto kmer HASH INDEX: 176
	//per il quinto kmer HASH INDEX: 112 e cosi via
	//printf("HASH INDEX: %d\n", ((hash32shift(kmer)*kmer) % max));
    return (hash32shift(kmer)*kmer) % max;  
} 



//Questa funzione prende in input un k-mer rappresentato come una stringa di caratteri (char * kmer) e la sua lunghezza (long int kmerLength).
//Essenzialmente la funzione esegue queste due operazioni:
//1. Reverse: Inverte il k-mer. Ciò significa che il primo carattere diventa l'ultimo, il secondo diventa il penultimo e così via fino alla metà della lunghezza del k-mer.
//2. Complement: Sostituisce ciascun carattere del k-mer con il suo complemento nucleotidico. Ad esempio, A diventa T, T diventa A, G diventa C e viceversa. Il carattere N rimane invariato.
//kmerLength = 15
void ReverseComplementKmer(char * kmer, long int kmerLength){
	for(int k = 0; k < kmerLength/2; k++){ //Il ciclo for scorre metà della lunghezza del k-mer
		char temp = kmer[k];
		kmer[k] = kmer[kmerLength -1 - k]; //Viene scambiato il carattere alla posizione k con il carattere corrispondente alla posizione simmetrica all'interno del k-mer.
		kmer[kmerLength -1 - k] = temp;
	}
	//printf("KMER1 DOPO IL ReverseKmer = %s\n", kmer);

	//Successivamente, un altro ciclo for scorre attraverso ogni carattere del k-mer invertito.
	//Durante ogni iterazione, il carattere corrente viene sostituito con il suo complemento nucleotidico:
	for(int i = 0; i < kmerLength; i++){
		if(kmer[i] == 'A' || kmer[i] == 'a'){
			kmer[i] = 'T';
		}else if(kmer[i] == 'T' || kmer[i] == 't'){
			kmer[i] = 'A';
		}else if(kmer[i] == 'G' || kmer[i] == 'g'){
			kmer[i] = 'C';
		}else if(kmer[i] == 'C' || kmer[i] == 'c'){
			kmer[i] = 'G';
		}else if(kmer[i] == 'N' || kmer[i] == 'n'){ 
			kmer[i] = 'N'; //Il carattere 'N' o 'n' viene lasciato invariato, poiché rappresenta una posizione non definita.
		}
	}
}

//Determine whether a kmer exists in the kmer hash table, if it exists, return the hash value, otherwise return -1
//cerca un determinato k-mer all'interno della tabella hash dei k-mer.
//Prende due parametri: un puntatore alla testa della tabella hash dei k-mer e il valore del k-mer da cercare.
long int SearchKmerHashTable(KmerHashTableHead * kmerHashTableHead, unsigned int kmer){
	//Utilizza una funzione di hash per calcolare l'indice hash (hashIndex) del k-mer all'interno della tabella hash. Questo indice viene calcolato utilizzando la funzione di hash e il numero di elementi nell'array (dimensione della tabella hash).
	long int hashIndex = Hash(kmer, kmerHashTableHead->allocationCount);
	//printf("VALORE IN KMER.CPP DI kmer = %u e di kmerHashTableHead->kmerHashNode[%ld].kmer = %u\n", kmer, hashIndex, kmerHashTableHead->kmerHashNode[hashIndex].kmer);
	while(true){ //Esegue un ciclo infinito che cerca il k-mer all'interno della tabella hash.
		if(kmerHashTableHead->kmerHashNode[hashIndex].kmer == 0){ //controlla se il k-mer nella posizione indicata dall'indice hash è uguale a 0. 
			return -1;;  //indica che il k-mer non è presente nella tabella hash.
		}

		///printf ("kmerHashTableHead");
		//Se il k-mer nella posizione indicata dall'indice hash è uguale al k-mer cercato più uno (kmer + 1), restituisce l'indice hash, indicando che il k-mer è stato trovato.
		if(kmerHashTableHead->kmerHashNode[hashIndex].kmer == kmer + 1){
			//printf ("SONO NELL'IF E RITORNO PROPRIO L'HASHINDEX = %ld\n", hashIndex);
			//printf("VALORE DI kmer = %u e di kmerHashTableHead->kmerHashNode[%ld].kmer = %u\n", kmer, hashIndex, kmerHashTableHead->kmerHashNode[hashIndex].kmer);
			return hashIndex;
		}else{
			//Se il k-mer nella posizione indicata dall'indice hash non è né 0 né kmer + 1, il ciclo incrementa l'indice hash di 5
			hashIndex = (hashIndex + 5)%kmerHashTableHead->allocationCount;
			//Questo è un metodo di risoluzione delle collisioni chiamato "linear probing", che cerca la prossima posizione disponibile nell'array per inserire il k-mer.
		}
	}
	return -1; //Il ciclo continua finché non viene trovato il k-mer o fino a quando non viene determinato che non esiste nella tabella hash. In quest'ultimo caso, viene restituito -1.
}

//Questa funzione implementa l'algoritmo di ordinamento noto come QuickSort, il quale suddivide ricorsivamente l'array in due sotto-array, ordinandoli separatamente.
//*a sarebbe kmerReadNodeHead->kmerReadNode, left = 0, right = 3
void sort(KmerReadNode * a, long int left, long int right)
{
    if(left >= right){ //caso base
        return ; //significa che l'array ha dimensione zero o uno e quindi è già ordinato, quindi si esce dalla funzione.
    }
    long int i = left;
    long int j = right;
	//si memorizza il valore del k-mer e le sue informazioni corrispondenti
    unsigned int key = a[left].kmer; //key è il valore del k-mer dell'elemento nella posizione iniziale (left) dell'array.
	unsigned int position = a[left].position;
	unsigned int readIndex = a[left].readIndex;
	bool orientation = a[left].orientation;

	//printf("PRIMA DEL SORT: i=%ld, j=%ld, key=%u, position=%u, readIndex=%u, orientation=%b\n", i, j, key, position, readIndex, orientation);
	
	//Si esegue il partizionamento dell'array: si sposta l'indice j verso sinistra fino a trovare un elemento minore della chiave, quindi si sposta l'indice i verso destra fino a trovare un elemento maggiore della chiave. Se i e j non si sono ancora incrociati, si scambiano gli elementi in posizioni i e j.
	//Due cicli while scorrono l'array in direzioni opposte, cercando di trovare elementi più piccoli e più grandi del valore di riferimento (key).
	//Quando trovano tali elementi, li scambiano tra loro.
	//Questo processo continua finché i e j non si incrociano.
    while(i < j)                       
    {
        while(i < j && key <= a[j].kmer){
            j--;
        }
		
		if(i < j){
			a[i].kmer = a[j].kmer;
			a[i].readIndex = a[j].readIndex;
			a[i].position = a[j].position;
			a[i].orientation = a[j].orientation;
			i++;
		}
         
        while(i < j && key >= a[i].kmer){
            i++;
        }
		
		if(i < j){
			a[j].kmer = a[i].kmer;
			a[j].readIndex = a[i].readIndex;
			a[j].position = a[i].position;
			a[j].orientation = a[i].orientation;
			j--;
		} 
    }
    
	//Restoring Key (Ripristino della chiave):
	//Alla fine del partizionamento, i indica la posizione corretta della chiave nell'array. Quindi, il valore della chiave e le sue informazioni corrispondenti 
	//vengono posizionati nell'array in posizione i, garantendo che tutti gli elementi a sinistra di i siano minori o uguali a key, mentre tutti quelli a destra di i siano maggiori o uguali a key.
    a[i].kmer = key;
	a[i].readIndex = readIndex;
	a[i].position = position;
	a[i].orientation = orientation;

	//Chiamate ricorsive:
	//Dopo il partizionamento e il ripristino della chiave, la funzione sort viene richiamata ricorsivamente sui sotto-array a sinistra e a destra dell'elemento i.
    sort(a, left, i - 1);
    sort(a, i + 1, right);  
	//Questo processo continua finché l'array non è completamente ordinato. Alla fine, l'intero array sarà stato ordinato in base al valore del k-mer.                
	//printf("ALLA FINE DEL SORT: i=%ld, j=%ld, key=%u, position=%u, readIndex=%u, orientation=%d\n", i, j, key, position, readIndex, orientation);
}

//restituisce true se almeno un numero della frequenza è diverso
//false se tutti i numeri della frequenza sono uguali
//prende come parametri la frequenza di un singolo kmer e la lunghezza del kmer stesso
bool DetectSameKmer(char * kmerf, long int kmerLength){ 
	long int i = 1;
	//printf ("Stampa di kmerf in DetectSameKmer %s\n", kmerf);
	//printf ("Stampa di kmer length in DetectSameKmer %ld\n", kmerLength);
	for(; i < kmerLength; i++){
		if(kmerf[0] != kmerf[i]){
			//printf("stampa di kmerf[0] %c\n", kmerf[0]);
			//printf("stampa di kmerf[i] %c\n", kmerf[i]);
			break;
		}
	}
	if(i < kmerLength){
		//printf("HA RITORNATO TRUE e il valore della i è questo: %d\n");
		return true; //almeno un numero all'interno della frequnza del kmer è diverso
	}
	//printf("Tutti i numeri della frequenza sono uguali, quindi ritorno false. Questo è il valore finale della i: %d\n", i);
	return false;//tutti i numeri all'interno della stessa frequenza del kmer sono uguali

}

KmerReadNodeHead * InitKmerReadNodeHead(char * address, ReadSetHead * readSetHead, long int kmerLength, long int step , KmerHashTableHead * kmerHashTableHead, int frequencyCutOff){

	FILE * fp; 
    if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	long int hashIndex = 0;
	long int maxSize = 100000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
	char * kmer = (char *)malloc(sizeof(char)*(kmerLength + 1));
	char * kmerC = (char *)malloc(sizeof(char)*(10 + 1));
	long int kmerCount = 0;
	long int length = 0;
	int frequency = 0;
	
	long int arrayCount = 1000000;
	int * freArray = (int *)malloc(sizeof(int)*arrayCount);
	for(long int i = 0; i < arrayCount; i++){
		freArray[i] = 0;
	}
	cout<<"bb"<<endl;
	long int allKmerFrequency = 0;
	while((fgets(line, maxSize, fp)) != NULL){
		
		length = strlen(line);
		
		strncpy(kmer, line + kmerLength + 1, length - kmerLength - 1);
		kmer[length-kmerLength-1] = '\0';
		if(DetectSameKmer(kmer, kmerLength) != true){
			continue;
		}
		
		frequency = atoi(kmer);	
		if(frequency > arrayCount - 10){
			continue;
		}
		freArray[frequency]++;
		kmerCount++;
		allKmerFrequency = allKmerFrequency + frequency;
	}
	
	printf("The number of kmer types: %ld;\n",kmerCount);
	printf("Sum of frequencies for different kmer: %ld;\n",allKmerFrequency);
	
	float acc = 0;
	long int max = 0;
	long int min = 0;
	
	for(long int i = 0; i < arrayCount; i++){
		if(freArray[i] != 0){
			acc = acc + (float)(freArray[i]*i)/allKmerFrequency;
			if(acc > 0.9 && max == 0){
				max = i;
				break;
			}
		}
	}
	
	if(min < 2){
		min = 2;
	}
	min = 2;
	if(min > max){
		cout<<"min is larger than max!"<<endl;
		exit(0);
	}
	
	printf("The range of kmer frequency is:[%ld,%ld];\n",min,max);
	
	fclose(fp);

	kmerCount = 0;
	allKmerFrequency = 0;
	

    if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	while((fgets(line, maxSize, fp)) != NULL){
		strncpy(kmer, line, kmerLength);
		kmer[kmerLength]='\0';
		if(DetectSameKmer(kmer, kmerLength) != true){
			continue;
		}
		
		length = strlen(line);
		strncpy(kmerC, line + kmerLength + 1, length - kmerLength - 1);
		kmerC[length-kmerLength-1] = '\0';
		frequency = atoi(kmerC);
		
		if(frequency > arrayCount - 10){
			continue;
		}
		if(frequency <= max && frequency >= min){
			kmerCount++;
			allKmerFrequency = allKmerFrequency + frequency;
		}
	}
	
	printf("There are %ld kmer for overlap detection;\n",kmerCount);
	
	fclose(fp);
	
	kmerHashTableHead->allocationCount = kmerCount*1.2;
	kmerHashTableHead->kmerHashNode = (KmerHashNode *)malloc(sizeof(KmerHashNode)*kmerHashTableHead->allocationCount);
	
	for(unsigned long int i = 0; i < kmerHashTableHead->allocationCount; i++){
		kmerHashTableHead->kmerHashNode[i].kmer = 0;
		kmerHashTableHead->kmerHashNode[i].startPositionInArray = -1;
	}
	
	if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	
	unsigned long int kmerInteger = 0;

	while((fgets(line, maxSize, fp)) != NULL){

		strncpy(kmer, line, kmerLength);
		kmer[kmerLength]='\0';
		if(DetectSameKmer(kmer, kmerLength) != true){
			continue;
		}
		length = strlen(line);
		strncpy(kmerC, line + kmerLength + 1, length - kmerLength - 1);
		kmerC[length-kmerLength-1] = '\0';
		frequency = atoi(kmerC);	
		if(frequency > arrayCount - 10){
			continue;
		}
		if(!(frequency <= max && frequency >= min)){
			continue;
		}
		
		
		SetBitKmer(&kmerInteger, kmerLength, kmer);
		
		hashIndex = Hash(kmerInteger, kmerHashTableHead->allocationCount);
		while(true){
			if(kmerHashTableHead->kmerHashNode[hashIndex].kmer == 0){
				kmerHashTableHead->kmerHashNode[hashIndex].kmer = kmerInteger + 1;
				break;
			}else{
				hashIndex = (hashIndex + 1) % kmerHashTableHead->allocationCount;
			}
		}
	}
	fclose(fp);
	
	KmerReadNodeHead * kmerReadNodeHead = (KmerReadNodeHead *)malloc(sizeof(KmerReadNodeHead));
	kmerReadNodeHead->realCount = 0;
	kmerReadNodeHead->allocationCount = allKmerFrequency*1.1;
	kmerReadNodeHead->kmerReadNode = (KmerReadNode *)malloc(sizeof(KmerReadNode)*kmerReadNodeHead->allocationCount);
	
	for(long int i = 0; i < kmerReadNodeHead->allocationCount; i++){
		kmerReadNodeHead->kmerReadNode[i].kmer = 0;
	}
	cout<<"bb3333:"<<kmerReadNodeHead->allocationCount<<endl;

	long int readLength = 0;

	char * kmer1 = (char *)malloc(sizeof(char)*(kmerLength + 1));
	

	
	for(long int i = 0; i < readSetHead->readCount; i++){
		readLength = readSetHead->readSet[i].readLength;
		for(int j = 0; j < readLength - kmerLength + 1-step ; j = j+step){
			
			strncpy(kmer1,readSetHead->readSet[i].read + j, kmerLength);
			kmer1[kmerLength] = '\0';
			if(DetectSameKmer(kmer1, kmerLength) != true){
				continue;
			}
			SetBitKmer(&kmerInteger, kmerLength, kmer1);
			hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
			if(hashIndex!=-1){
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].kmer = kmerInteger;
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].readIndex = i+1; 
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].position= j; 
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].orientation = true;
				kmerReadNodeHead->realCount++;
			}else{
				ReverseComplementKmer(kmer1, kmerLength);
				SetBitKmer(&kmerInteger, kmerLength, kmer1);
				hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
				if(hashIndex!=-1){
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].kmer = kmerInteger;
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].readIndex = i+1; 
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].position= j; 
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].orientation = false;
					kmerReadNodeHead->realCount++;
				}
			}
			
		}	
	}
	
	sort(kmerReadNodeHead->kmerReadNode, 0, kmerReadNodeHead->realCount - 1);

	kmerInteger = kmerReadNodeHead->kmerReadNode[0].kmer + 1;
	for(long int i = 0; i < kmerReadNodeHead->realCount; i++){
		if(kmerReadNodeHead->kmerReadNode[i].kmer != kmerInteger){
			kmerInteger = kmerReadNodeHead->kmerReadNode[i].kmer;
			hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
			kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray = i;
		}
	}
	free(line);
	free(kmer);
	free(kmerC);
	free(kmer1);
	kmerReadNodeHead->kmerLength = kmerLength;
	
	return kmerReadNodeHead;
}
// According to the frequency of kmer, extract useful kmer and store it in the hash table
// crea la tabella hash dei kmer
// @param  char* address -> kmer file
// @param ReadSetHead * readSetHead -> lista delle reads
// @param long int kmerLenght -> lunghezza dei kmer (15)
// @param int step -> lunghezza dello step per estrazione kkmer (1)
// @param int min -> frequenza minima del kmer (2)
// @param float maxRatio -> tasso di frequenza massima del kmer (0,9)
KmerHashTableHead * GetKmerHashTableHead(char * address, ReadSetHead * readSetHead, long int kmerLength, long int step, long int min, float maxRatio){
	printf("\n funzione GetKmerHashTableHead\n");
	KmerHashTableHead * kmerHashTableHead = (KmerHashTableHead *)malloc(sizeof(KmerHashTableHead));
	FILE * fp; 
    if((fp = fopen(address, "r")) == NULL){ //apro il kmer file
        printf("%s, does not exist!", address);
        exit(0);
    }
	long int hashIndex = 0;
	long int maxSize = 100000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
	char * kmer = (char *)malloc(sizeof(char)*(kmerLength + 1));
	char * kmerC = (char *)malloc(sizeof(char)*(10 + 1));
	long int kmerCount = 0;
	long int length = 0;
	int frequency = 0;
	long int arrayCount = 1000000;
	int * freArray = (int *)malloc(sizeof(int)*arrayCount);
	for(long int i = 0; i < arrayCount; i++){
		freArray[i] = 0;
	}
	
	long int allKmerFrequency = 0;

	while((fgets(line, maxSize, fp)) != NULL){
		length = strlen(line);
		//copia dal file dei kmer, le frequnze dei kmer e le inserisce nell'array kmer
		strncpy(kmer, line + kmerLength + 1, length - kmerLength - 1); //si copia solo la frequenza del kmer (%ld) e non la stringa che c'era prima
		kmer[length-kmerLength-1] = '\0';
		//printf("%s\n",kmer);

		//analizza solo i kmer che differiscono almeno in un carattere

		//errore riferito nel paper

		if(DetectSameKmer(kmer, kmerLength) == false){ //se tutte le cifre all'interno della frequnza del kmer sono uguali
			//printf("frequenza del kmer con tutte cifre uguali: %s\n",kmer);
			continue; //ritorna al while -> passa al kmer successivo
		}
		

		//se esiste almeno un carattere diverso
		//printf ("Valore di arrayCount: %ld\n", arrayCount);
		frequency = atoi(kmer);	//frequenza del kmer
		//arrayCount = 1 milione
		if(frequency > arrayCount - 10){ //capire perche
			//printf("ho skippato il kmer con la frquenza superiore a 1 milione meno 10 %d\n",k);
			continue;
			
		}

		//freArray è il contatore delle frequenze. Cioè freArray[1] conta quanti kmer ci sono con frequenza 1
		//freArray[1945] conta quanti kmer con frequenza 1945 ci sono
		//ad esempio se la frequenza del kmer che legge è 1945, va in freArray[1945] e incrementa il contatore a 1
		//somma cumulativa (sommatoria nel paper)
		freArray[frequency]++;
		
		printf("frequency: %ld, freeArray: %ld\n",frequency, freArray[frequency]); //UTILE
		
		//kmerCount inizialmente è zero
		kmerCount++;  // The number of kmer types
		//allKmerFrequency inizialmente è 0. In allKmerFrequency somma le frequenze totali di tutti i kmer
		//ad esempio legge tre kmer con frequenza 1945, 25, 30. allKmerFrequency avrà il valore 2000
		allKmerFrequency = allKmerFrequency + frequency;  //The sum of the frequency of different types of kmer
	}
	printf("all kmer frequency: %ld\n",allKmerFrequency);
	printf("The number of kmer types: %ld;\n",kmerCount);

	
	float acc = 0;
	long int max = 0;
	
	if(maxRatio > 1){ //maxRatio = 0.9
		maxRatio = 1;
	}
	
	//arrayCount = 1 milione
	for(long int i = 0; i < arrayCount; i++){
		//freArray[0] è il contatore dei kmer con frequenza 0
		//freArray[1945] è il contatore dei kmer con frequenza è 1945
		//salta i valori delle frequenze che non sono collegati ad alcun kmer
		//ad esempio nessun kmer ha frequenza 60, quindi salta freArray[60]
		if(freArray[i] != 0){ 
			//il valore di allKmerFrequency è 4587
			//acc inizialmente è zero ed è di tipo float
			//prima di arrivare a freArray[1945], acc=0.151951
			//consideriamo freArray[1945] e consideriamo che ci siano due kmer con frequenza 1945
			//quindi freArray[1945] = 2
			//andiamo ad aggiungere ad acc il valore freArray[1945] moltiplicato per i che è 1945
			//quindi freArray[i]*1 = 2*1945 = 3890
			//acc = 0.151951 + (3890 / 4587) = 1.0
			acc = acc + (float)(freArray[i]*i)/allKmerFrequency; //in free array ci sono tutti 0 e 1 -> moltiplica 1 
													//per il valore corrente di i e lo divide per la frequenza totale dei kmers
			
			//printf ("valori delle i: %ld\n", i);

			//non abbiamo kmer che hanno frequenze da 51 a 1944.
			//per questo motivo si passa da freArray[50] a freArray[1945]
			if (i==50) {
				printf ("Valore di acc arrivati alla casella 50 che è quella che precede 1945: %lf\n", acc);
			}
		
			//maxRatio = 0.9
			//nel nostro caso acc = 1.0, quindi acc è maggiore di maxRatio e max = 0
			//allora max = i
			if(acc > maxRatio && max == 0){
				//max = 1945
				max = i;  //max assume la frequenza del kmer maggiore. Nostro caso il kmer maggiore ha frequenza 1945 e ce ne sono due
				break;
			}
		}
	}

	printf("il max di acc è: %ld\n",max);
	
	if(min < 2){ //min = 2 di default
		//printf("minimumKmerFrequency < 2, LROD sets minimumKmerFrequency = 2!");
		min = 2;
	}

	if(min > max){	
		cout<<"The paramter minimumKmerFrequency is larger than maxKmerFrequencyRatio. Please increas the value of maxKmerFrequencyRatio or decrease the value of minimumKmerFrequency!"<<endl;
		exit(0);
	}

	//ricapitolando min = 2 e max = 1945
		
	fclose(fp);

	kmerCount = 0;
	allKmerFrequency = 0;
	

    if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	//Get the number of available kmer types
	while((fgets(line, maxSize, fp)) != NULL){
		//mentre prima si copiava solo la frequenza, ora si copia solo il kmer
		strncpy(kmer, line, kmerLength);
		kmer[kmerLength]='\0';
		//printf("KMER COPIATO: %s\n", kmer);
		
		//qui effettua il confronto sui caratteri del kmer
		//se c'è un kmer con tutti caratteri uguali, lo salta
		if(DetectSameKmer(kmer, kmerLength) != true){
			//printf("kmer con tutti caratteri uguali: %s\n",kmer);
			continue;
		}
		//considera quelli che hanno almeno un carattere diverso
		length = strlen(line); //length = 15
		//si copia le frequenze in kmerC
		strncpy(kmerC, line + kmerLength + 1, length - kmerLength - 1);
		//quindi ha copiato il kmer vero e proprio in kmer e ha copiato la frequenza in kmerC
		kmerC[length-kmerLength-1] = '\0';
		frequency = atoi(kmerC);
		
		if(frequency > arrayCount - 10){ //arrayCount = 1 milione
			continue;
		}
		if(frequency <= max && frequency >= min){ //verifica che la frequenza sia compresa tra min e max
			kmerCount++; //si conta i kmer che hanno la frequenza compresa tra min e max. 
			//inizialmente allKmerFrequency=0
			//calcola la somma cumulativa
			allKmerFrequency = allKmerFrequency + frequency; //fa la sommatoria delle frequenze dei kmer
		}
	}
	//alla fine in kmerCount avremo il numero di kmer con frequenza compresa tra min e max
	//ci sono x candidati per calcolare la regione di overlap -> gia ha costruito la chain di kmer consistenti(?)
	//nel nostro caso abbiamo 99 kmer perchè solo il primo è stato escluso, avendo 15 caratteri tutte A
	printf("There are %ld kmer for overlap detection;\n",kmerCount);
	if(kmerCount <= 0){
		exit(0); //non c'è overlap
	}
	fclose(fp);
	
	//allocationCount = 99*2 = 198
	kmerHashTableHead->allocationCount = kmerCount*2;
	//printf ("ALLOCATION COUNT: %ld\n", kmerHashTableHead->allocationCount);
	kmerHashTableHead->kmerHashNode = (KmerHashNode *)malloc(sizeof(KmerHashNode)*kmerHashTableHead->allocationCount);
	
	//inizializza i nodi della kmer hashtable per i che va da 0 a 198
	for(unsigned long int i = 0; i < kmerHashTableHead->allocationCount; i++){
		kmerHashTableHead->kmerHashNode[i].kmer = 0; //kmer (visto come intero)  viene posto a 0
		kmerHashTableHead->kmerHashNode[i].startPositionInArray = -1; //posizionne del kmer nell'array = -1
	}
	
	if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	
	unsigned long int kmerInteger = 0; //variabile che contiene il kmer visto come intero (operazioni bit a bit)

//Convert kmer with a frequency between max and min into bytes and store it in a hash table
	printf("\n");
	printf("costruiamo la tabella: key= hash(kmerInteger) value = kmerInteger\n");

	while((fgets(line, maxSize, fp)) != NULL){

		strncpy(kmer, line, kmerLength); //copia il singolo kmer (senza la frequenza) in kmer
		kmer[kmerLength]='\0';
		//considera solo quelli che hanno almeno un carattere diverso nel kmer (non considera quelli tutti uguali)
		if(DetectSameKmer(kmer, kmerLength) != true){
			continue;
		}
		length = strlen(line);
		strncpy(kmerC, line + kmerLength + 1, length - kmerLength - 1); //copia solo la frequenza in kmerC
		kmerC[length-kmerLength-1] = '\0';
		frequency = atoi(kmerC);	
		if(frequency > arrayCount - 10){
			continue;
		}
		
		if(!(frequency <= max && frequency >= min)){
			continue;
		}

		//H(x) -> valore hash |||| x=3 -> ABCDERFGKEJKLKFS   -> H(kmerInteger) = Indice_tabella_hash
		SetBitKmer(&kmerInteger, kmerLength, kmer); //converte i kmer in valori binari e li salva in kmerInteger
		
		hashIndex = Hash(kmerInteger, kmerHashTableHead->allocationCount); //calcola gli indici sul kmerHash
		
		while(true){
			//dimensione della tabella hash
			//hashIndex = 3 per il kmer = AAA ||||| kmer = TTT  hashIndex=2 
			if(kmerHashTableHead->kmerHashNode[hashIndex].kmer == 0){ //se l'hashIndex non è stato gia usato in una precedente iterazione ->mettici il kmer integer+1
				//printf ("hashIndex = %ld ,kmerInteger = %ld, kmerHashTableHead-> kmerHashNode[hashindex].kmer = %u\n", hashIndex,kmerInteger, kmerHashTableHead->kmerHashNode[hashIndex].kmer); UTILE
				kmerHashTableHead->kmerHashNode[hashIndex].kmer = kmerInteger + 1;//pone il valore del kmer a kmer integer e poi si muove tramite l'hash di kmer integer
				break;
			}else{
					
				//semplicemente incrementa l'hashIndex di 5
				hashIndex = (hashIndex + 5) % kmerHashTableHead->allocationCount; //vai avanti
				//printf("HashIndex incrementato di 5 = %ld\n",hashIndex); UTILE
			}
		} 
		//ottine la tabella hash con key= hash(kmerInteger) value = kmerInteger
	
	}
	fclose(fp);
	printf("\n");

	printf("stampa della kmerHashtable: \n");
	for (int i = 0; i <  kmerHashTableHead->allocationCount; i++){
		printf("elemento %d = %ld\n",i, kmerHashTableHead->kmerHashNode[i].kmer);
	}
	
	free(line);
	free(kmer);
	free(kmerC);
	kmerHashTableHead->min = min;
	kmerHashTableHead->max = max;

	printf("kmer = %u e startPositionInArray = %d\n", kmerHashTableHead->kmerHashNode->kmer, kmerHashTableHead->kmerHashNode->startPositionInArray);

	return kmerHashTableHead;
}
//check sul numero di kmer
int GetKmerHashTableHead_UnitTest(KmerHashTableHead * kmerHashTableHead){
	long int kmerCount = 0;
	for(unsigned long int i = 0; i < kmerHashTableHead->allocationCount; i++){
		//se la tabella hash dei k-mer ha una dimensione non valida
		if(kmerHashTableHead->kmerHashNode[i].kmer < 0){
			return 0;
		//Se trova almeno un k-mer valido, incrementa il contatore 
		}else if(kmerHashTableHead->kmerHashNode[i].kmer > 0){
			kmerCount++;
		}  
	}

	printf("KMERCOUNT = %lu\n", kmerCount);

	if(kmerCount <= 0){
		return 0;
	}
	return 1;
}

/*
@param readSetHead ->lista di reads
@param kmerLenght -> lunghezza del kmer 15
@param step -> 1
@param intervalCount -> intervallo di 50000 reads
@return KmerReadNodeHead ->
*/
KmerReadNodeHead * GetKmerReadNodeHeadSub(ReadSetHead * readSetHead, long int kmerLength, long int step, long int intervalCount){
	
	printf("\n Funzione GetKmerReadNodeHeadSub\n");

	long int max = 0;
	long int allKmerFrequency = 0;
	long int startReadIndex = 0;
	//intervalCount=50.000
	long int endReadIndex = startReadIndex + intervalCount;

	printf("readSetHead->readCount: %ld\n", readSetHead->readCount);

	while(true){
		if(startReadIndex >= readSetHead->readCount){ //se lo start read index supera il numero di letture, interrompi
			break;
		}
		if(endReadIndex >= readSetHead->readCount){
			endReadIndex = readSetHead->readCount - 1; //decrementa l'indice di lettura
		}
		allKmerFrequency = 0; //ad ogni intervallo, si pone la frequenza a zero e la si compara a quella precedente per calcolare il max
		for(long int i = startReadIndex; i <= endReadIndex; i++){
			allKmerFrequency = allKmerFrequency + readSetHead->readSet[i].readLength - kmerLength;
		}
		if(allKmerFrequency > max){
			max = allKmerFrequency;
		}
		startReadIndex = endReadIndex + 1; //startReadIndex = 4 (comincia dalla fine)
		endReadIndex = endReadIndex + intervalCount; //endReadIndex = 50.003
	}
	//step = 1
	printf("End allkmerfrequency: %ld\n", allKmerFrequency);
	max = max/step; //max = 23189
	
	
	KmerReadNodeHead * kmerReadNodeHead = (KmerReadNodeHead *)malloc(sizeof(KmerReadNodeHead)); //alloca la testa della lista di kmer
	kmerReadNodeHead->realCount = 0;
	kmerReadNodeHead->allocationCount = max*1.1; //si vuole aumentare la dimensione di circa il 10%
	kmerReadNodeHead->kmerReadNode = (KmerReadNode *)malloc(sizeof(KmerReadNode)*kmerReadNodeHead->allocationCount); //alloca gli altri nodi kmer
	
	//inizializzazione
	for(long int i = 0; i < kmerReadNodeHead->allocationCount; i++){
		kmerReadNodeHead->kmerReadNode[i].kmer = 0;
	}
	
	return kmerReadNodeHead; //ritorna il puntatore alla testa
}

// Initialize the kmer hash table node data
//@param int startIndex -> indice della prima read del batch di reads
//@param int endIndex -> indice dell'ultima read del batch di reads
void InitKmerReadNodeHeadSub(ReadSetHead * readSetHead, KmerReadNodeHead * kmerReadNodeHead, KmerHashTableHead * kmerHashTableHead, long int kmerLength, long int step, long int startReadIndex, long int endReadIndex){
	printf("funzione: InitKmerReadNodeHeadSub\n");
	printf("\n startIndex = %ld endReadIndex = %ld\n ",startReadIndex, endReadIndex);
	
	kmerReadNodeHead->realCount = 0; //inizializzazione contatore dei k-mer
	
	for(unsigned long int i = 0; i < kmerHashTableHead->allocationCount; i++){
		kmerHashTableHead->kmerHashNode[i].startPositionInArray = -1; //inizializzazione delle posizioni dei k-mer nella tabella hash dei k-mer a -1 per indicare che non sono ancora stati assegnati:
	}
	
	for(long int i = 0; i < kmerReadNodeHead->allocationCount; i++){
		kmerReadNodeHead->kmerReadNode[i].kmer = 0; //Inizializzazione dei k-mer
	}
	
	long int readLength = 0; //conterrà la lunghezza della read corrente

	char * kmer1 = (char *)malloc(sizeof(char)*(kmerLength + 1)); //conterrà il k-mer estratto da ogni posizione nel read.
	
	long int hashIndex = 0; //utilizzato per memorizzare l'indice del k-mer nella tabella hash dei k-mer.
	unsigned long int kmerInteger = 0; //memorizzerà la rappresentazione intera del k-mer estratto.
	//Extract available kmer in reads
	//Il ciclo esterno scorre attraverso tutti i read nell'intervallo specificato da startReadIndex a endReadIndex.
	//startReadIndex = 0 e endReadIndex = 3

	printf("\n stampa della struttura kmerReadNode\n");

	for(long int i = startReadIndex; i <= endReadIndex; i++){
		//printf("readSetHead->readSet[i].readLength = %ld\n", readSetHead->readSet[i].readLength);
		readLength = readSetHead->readSet[i].readLength;

		//All'interno del ciclo esterno, un altro ciclo scorre attraverso ogni posizione nel read corrente
		//step = 1, quindi j viene incrementato di 1
		//j viene incrementato di 1

		//quando i=0, sta leggendo la prima read
		for(int j = 0; j < readLength - kmerLength + 1-step ; j = j+step){

		// Estrae un k-mer dalla read di lunghezza k-mer length
        // Verifica se il k-mer esiste nella tabella hash dei k-mer
        // Aggiunge il k-mer alla struttura KmerReadNodeHead con l'indice del read, la posizione nel read e l'orientamento (true indica l'orientamento originale del k-mer).
        // Gestisce anche il complemento inverso del k-mer

			//copia al massimo kmerLength caratteri dalla stringa di origine readSetHead->readSet[i].read + j nella stringa di destinazione kmer1
			//in pratica quando j=0, si copia i primi 15 caratteri da 0 a 14
			//quando j=1, si copia i 15 caratteri da 1 a 15
			//poi quelli da 2 a 16 e cosi via
			//in questo modo si copia tutti i possibili kmer di lunghezza 15
			strncpy(kmer1,readSetHead->readSet[i].read + j, kmerLength);
			kmer1[kmerLength] = '\0';
			//printf("kmer1 = %s\n", kmer1);
			
			//chiama setBitKmer sul kmer1 di lunghezza 15 che sta analizzando in questo momento
			SetBitKmer(&kmerInteger, kmerLength, kmer1);

			//hashIndex è quasi sempre -1 tranne con hashIndex = 94, 1, 46
			hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
			//printf("hashIndex in kmer.cpp = %ld\n", hashIndex);
			
			if(hashIndex!=-1){
				printf("\n kmer: %s kmerInteger: %ld readIndex: %ld\n",kmer1, kmerInteger,i);
				//printf("\n realcount = %ld\n",kmerReadNodeHead->realCount);

				//printf("\n SONO NEL PRIMO IF -> il kmer è presente nella tabella hash\n");
				//sta considerando i tre indici 94, 1 e 46 (controllare questa cosa)
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].kmer = kmerInteger;
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].readIndex = i+1; 
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].position= j;  //punto di partenza del kmer nella read
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].orientation = true;
				kmerReadNodeHead->realCount++; //incremento il realCount
				
 				//printf("HASHINDEX=%ld, KMER=%lu, READINDEX=%lu, POSITION=%d\n", hashIndex,kmerInteger,i+1,j);
			}else{
				

				//printf("\n SONO NEL ELSE -> il kmer non è presente nella tabella hash,verifico la presenza del kmer reverse\n");

				//Se il k-mer non esiste nella tabella hash, viene generato il suo complemento inverso, e se esiste nella tabella hash, viene aggiunto alla struttura KmerReadNodeHead con le stesse informazioni e l'orientamento impostato su false per indicare che è stato estratto dal complemento inverso del read.
				ReverseComplementKmer(kmer1, kmerLength);  //Reverse complementary kmer, and then check whether it exists
				SetBitKmer(&kmerInteger, kmerLength, kmer1); //Il complemento inverso viene convertito in un intero kmerInteger usando la funzione SetBitKmer.
				hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger); //Viene quindi verificato se il suo complemento inverso esiste nella tabella hash dei k-mer
				if(hashIndex!=-1){ //Se il suo complemento inverso esiste nella tabella hash
					
					printf("\n kmer: %s kmerInteger: %ld readIndex: %ld\n",kmer1, kmerInteger,i);
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].kmer = kmerInteger;
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].readIndex = i+1; 
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].position= j; 
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].orientation = false; //impostando l'orientamento su false per indicare che il k-mer è stato estratto dal complemento inverso del read.
					kmerReadNodeHead->realCount++; //viene incrementato per tener traccia del numero totale di k-mer estratti.
					//printf("HASHINDEX=%ld, KMER=%lu, READINDEX=%lu, POSITION=%d\n", hashIndex,kmerInteger,i+1,j);
				}
			}
		}	
	}
	sort(kmerReadNodeHead->kmerReadNode, 0, kmerReadNodeHead->realCount - 1);//Ordina la struttura KmerReadNodeHead in base al valore del k-mer per facilitare l'accesso in seguito.

	kmerInteger = kmerReadNodeHead->kmerReadNode[0].kmer + 1;
	//Indica la posizione iniziale di ogni k-mer nella tabella hash dei k-mer, in modo che possa essere facilmente accessibile in seguito poichè saranno organizzati in ordine crescente
	for(long int i = 0; i < kmerReadNodeHead->realCount; i++){
		// Individua il primo k-mer con un valore diverso dagli altri e registra la sua posizione
		//Questo viene fatto per garantire che la posizione iniziale di ciascun k-mer nella tabella hash dei k-mer sia facilmente accessibile in seguito. 
		if(kmerReadNodeHead->kmerReadNode[i].kmer != kmerInteger){
			kmerInteger = kmerReadNodeHead->kmerReadNode[i].kmer;
			hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
			kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray = i; //La posizione viene memorizzata nella struttura KmerHashTableHead nella posizione corrispondente al k-mer nella tabella hash.
		}
	}
	//Viene aggiornata la lunghezza del k-mer, l'indice di inizio e l'indice di fine dei read associati alla struttura KmerReadNodeHead.
	kmerReadNodeHead->kmerLength = kmerLength;
	kmerReadNodeHead->startReadIndex = startReadIndex;
	kmerReadNodeHead->endReadIndex = endReadIndex;
}

int GetKmerReadNodeHeadSub_UnitTest(KmerReadNodeHead * kmerReadNodeHead){
	if(kmerReadNodeHead->realCount <= 0){
		cout<<kmerReadNodeHead->realCount<<endl;
		return 0;
	}
	for(long int i = 0; i < kmerReadNodeHead->realCount - 1; i++){
		cout<<kmerReadNodeHead->kmerReadNode[i].kmer<<endl;
		cout<<kmerReadNodeHead->kmerReadNode[i + 1].kmer<<endl;
		if(kmerReadNodeHead->kmerReadNode[i].kmer < kmerReadNodeHead->kmerReadNode[i + 1].kmer){
			return 0;
		}
	}
	return 1;
}


#endif