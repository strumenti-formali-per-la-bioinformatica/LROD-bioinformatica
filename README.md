# LROD-bioinformatica
Questo lavoro presenta l'analisi di un algoritmo di overlap detection a lettura lunga (LROD) che può migliorare l’accuratezza dei risultati delle sovrapposizioni tra long reads. 
Le tecnologie di sequenziamento di terza generazione frammentano il genoma in un gran numero di long reads e il processo di ricombinazione di queste letture in una sequenza di DNA completa è chiamato assemblaggio del genoma. L’overlap detection tra due long reads è utile per il processo di assemblaggio del genoma, allineando e unendo i frammenti di DNA in un’unica sequenza
continua. A differenza delle short reads, le long reads permettono un assemblaggio più accurato, evitando regioni ripetute e regioni più complesse. Tuttavia l’alto numero di errori che derivano dal sequenziamento di terza generazione (TGS) implica che ottenere l’overlap detection è ancora un compito impegnativo. L’obiettivo è quello di trovare sovrapposizioni in base alla distribuzione dei k-mers, ovvero delle sottostringhe del genoma di lunghezza k. Diversamente da altri algoritmi che utilizzano i k-mers, LROD conserva innanzitutto solo i k-mer solidi comuni tra le long reads che, rispetto ai k-mers tradizionali, hanno una frequenza più elevata. Dato che le tecniche TGS hanno un alto tasso di errore e le regioni ripetitive complicano il processo di overlap detection, LROD tenta di risolvere il problema sfruttando i solid k-mers. LROD utilizza una strategia che si adopera in varie fasi:
(1) innanzitutto trova l’insieme di k-mers comuni solidi;
(2) in secondo luogo trova una catena che include i k-mers comuni consistenti. La consistenza della catena è determinata da alcune condizioni che vedremo durante lo studio degli algoritmi;
(3) infine, tramite la catena, valuta la regione di overlap candidata e la restituisce.
# Passi per utilizzare il progetto
## Prerequisiti
E' necessario installare DSK, un programma di conteggio dei k-mers, per calcolare la frequenza di ciascun k-mer nel dataset.
## Table of Contents
<a name="general-info"></a>
### General Info
## General Info
***
Write down the general informations of your project. It is worth to always put a project status in the Readme file. This is where you can add it. 
### Screenshot
![Image text](/path/to/the/screenshot.png)
## Technologies
***
A list of technologies used within the project:
* [Technology name](https://example.com): Version 12.3 
* [Technology name](https://example.com): Version 2.34
* [Library name](https://example.com): Version 1234
