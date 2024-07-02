# LROD
## Contenuto della cartella
La cartella contiene:
* la presentazione per l'esame e un riassunto del paper
* il paper originale, sottolineato
* il paper scritto da Venturino Silvio e Staiano Catello.

Il link del progetto originale è: https://github.com/luojunwei/LROD
Il link del tool dsk è: https://github.com/GATB/dsk
Sui rispettivi siti è descritto come installare i tools
Sul sito di DSK è descritto come ottenere il file dei k-mer con le frequenze associate, nel formato che LROD richiede. Anche nel nostro paper.
Sul sito di LROD è descrtto come avviare LROD e il comando da eseguire. Anche nel nostro paper.
## Note:
* il codice contiene i commenti per indicare un **riferimento** agli algoritmi spiegati nel paper. Per trovarli nel codice basta cercare "Riferimento algoritmo";
* i commenti nel codice che sono in inglese sono quelli originali degli autori, presenti nel codice originale.
* i punti critici individuati da noi, sono segnalati nel paper. Il riferimento alle righe di codice (riga x)  che abbiamo inserito per indicare dove si trovano gli algoritmi oppure i punti critici, fanno riferimento al progetto originale e NON a quello commentato da noi. 
* il paper non fornisce il modo in cui ha valutato F1-score, precision e recall di LROD. Generalmente questi indicatori vengono utilizzati per valutare un classificatore, nel contesto dell'intelligenza artificiale.
LROD invece è un algoritmo (probabilistico) che sfrutta la catena di kmer consistenti per trovare regioni di overlap. Gli indicatori F1-score, precision e recall nel contesto di LROD, valutano con quanta "precisione" l'algoritmo rilevi le regioni di overlap, sfruttando la catena di k-mer consistenti. Si valutano i risultati dell'overlap, contenuti nel file result.csv.

* non è chiaro dai link su quali **dataset** (reali) i ricercatori abbiano effettuato i test. 
* non sono disponibili i dataset simulati, ovvero quelli creati da loro, ma c'è scritto il tool che hanno usato.
* LROD su fedora 38 dava errore (illegal instruction) e non poteva essere eseguito. Per ovviare al problema si è usato virtual box con un'altra distribuzione linux (mint)
* Il codice commentato da noi, non è compilato. Si consiglia di scaricare il progetto originale da github e compilarlo nuovamente col proprio compilatore.
* potrebbero esserci altri punti critici che non siamo riusciti ad individuare..

SV

