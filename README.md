# LROD-bioinformatica
## Problema 
Le tecnologie di sequenziamento di terza generazione frammentano il genoma in un gran numero di long reads e il processo di ricombinazione di queste letture in una sequenza di DNA completa è chiamato assemblaggio del genoma. L’overlap detection tra due long reads è utile per il processo di assemblaggio del genoma, allineando e unendo i frammenti di DNA in un’unica sequenza continua. A differenza delle short reads, le long reads permettono un assemblaggio più accurato, evitando regioni ripetute e regioni più complesse. Tuttavia l’alto numero di errori che derivano dal sequenziamento di terza generazione (TGS) implica che ottenere l’overlap detection è ancora un compito impegnativo. 
## Obiettivo
Questo lavoro presenta l'analisi di un algoritmo di overlap detection a lettura lunga (LROD) che può migliorare l’accuratezza dei risultati delle sovrapposizioni tra long reads. L’obiettivo è quello di trovare sovrapposizioni in base alla distribuzione dei k-mers, ovvero delle sottostringhe del genoma di lunghezza k. Diversamente da altri algoritmi che utilizzano i k-mers, LROD conserva innanzitutto solo i k-mer solidi comuni tra le long reads che, rispetto ai k-mers tradizionali, hanno una frequenza più elevata. Dato che le tecniche TGS hanno un alto tasso di errore e le regioni ripetitive complicano il processo di overlap detection, LROD tenta di risolvere il problema sfruttando i solid k-mers. LROD utilizza una strategia che si adopera in varie fasi:
(1) innanzitutto trova l’insieme di k-mers comuni solidi;
(2) in secondo luogo trova una catena che include i k-mers comuni consistenti. La consistenza della catena è determinata da alcune condizioni che vedremo durante lo studio degli algoritmi;
(3) infine, tramite la catena, valuta la regione di overlap candidata e la restituisce.

# Passi per utilizzare il progetto
1) Introduzione
```
    I dati in input ad LROD sono le long reads (in formato fasta).
```
2) Prima dell'installazione e dell'esecuzione
```
    E' necessario installare DSK, un programma di conteggio dei k-mers, per calcolare la frequenza di ciascun k-mer nel dataset.
```
3) Installazione.
```
    LROD dovrebbe essere eseguito sul sistema operativo Linux con gcc. Testiamo LROD utilizzando gcc 4.6.3 su Ubuntu.
    Creare una main directory (esempio:LROD). Copiare tutto il codice sorgente in questa directory.
	cd LROD
	make all
```
4) Esecuzione.
```
    Step 1: Usa DSK per creare il file con le frequenze dei k-mers.
        dsk -file <long-read-file.fa>  -kmer-size 15
	dsk2ascii -file <long-read-file.h5> -out kmer-frequency-file.txt
    Step 2: LROD -r <long-read-file> -c <kmer-frequency-file> -o result-file [options]
    	-r long-read-file: file in input in formato fasta;
	-c kmer-frequency-file: ogni linea nel kmer-frequency-file dovrebbe essere formata dalla coppia "kmer kmer-frequency";
	-o result-file: result file;
	-t count: conteggio dei thread (default 1);
	-k kmerLength: lunghezza del kmer (default 15);
	-q smallKmerLength: lunghezza kmer più piccoli (default 9);
	-f minimumKmerFrequency: frequenza minima del kmer (default 2);
	-m maxKmerFrequencyRatio: massimo rapporto di frequenza del kmer (default 0.9);
	-s kmerStep: kmer step (default 1);
	-d distance: il valore piccolo della distanza usato per determinare se due kmer sono consistenti (default 400);
	-e distance: il valore grande della distanza usato per determinare se due kmer sono consistenti (default 1500);
	-a min-overlap-length: la lunghezza minima di sovrapposizione tra due long reads (default 500);
	-b length-ratio: il massimo rapporto di lunghezza tra due regioni allineate (default 0.3); 
	-h -help: mostra le regole di utilizzo per LROD.
	
    Nota:
    	ogni linea nel kmer-frequency-file dovrebbe essere formata dalla coppia "kmer kmer-frequency".
	Per esempio:
	AGTCCAGGCCGGGAA 3
	GAAATCCAGCCGCCG 6
	AACCGGCGAATCGGA 3
	TATTTTAACATTCTC 2
	TATGGCCGATGAATT 4
	AAAGCCGAAGCCTAG 3
	CATCTTCACATCAGA 2
	ATAAGTGATAGCTTC 4
	TCGGCCATATTACCA 4
	ATTATTGCAATACTT 6
     Di seguito è mostrato un esempio di riga di comando. Il file di output è result.csv.
	LROD -r sra.fasta -c kmer-cout -o result
```
5) Output.
```
    Il file di output "output-file-name.csv" è il risultato della sovrapposizione.
    La prima colonna è il primo numero letto.
    La seconda colonna è il secondo numero letto.
    La terza colonna  è l'orientamento dell'allineamento. 0 rappresenta l'allineamento in avanti (forward). 1 rappresenta l'allineamento inverso (reverse).
    La quarta colonna è la posizione di partenza nella prima read.
    La quinta colonna è la posizione finale nella prima read.
    La sesta colonna è la lunghezza della prima read.
    La settima colonna è la posizione di partenza nella seconda read.
    L'ottava colonna è la posizione finale nella seconda read.
    La nona colonna è la lunghezza della seconda read.
    
    Di seguito è mostrato un esempio nel file dei risultati.
    36423,1,0,5326,9923,9923,1,4364,10479
    Significa che la regione [5326,9923] nella prima read è sovrapposta alla regione [1,4364] nella seconda read.
```
6) Test.
```
    Nella directory Test, sono presenti il ​​file di long read e il relativo file di frequenza del k-mer. Gli utenti possono usare questo dataset per eseguire LROD.
```
