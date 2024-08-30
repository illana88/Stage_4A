# Pour exécuter le pipeline original
Pour commencer, il faut installer conda (par exemple via ce lien : https://docs.anaconda.com/free/miniconda/miniconda-other-installer-links/), puis créer l'environnement comme indiqué sur le pipeline original. Enfin, exéctuer les commandes suivantes :  
conda config --add channels bioconda  
conda config --add channels conda-forge  

## Partie A
Les fichiers nécessaires à l'exécution de cette partie sont les suivant (à avoir dans le dossier 'pgp-a') :  
- abundant_tx.R  
- all_bams.tsv  
- Auto_CoverV4_layered_intronV3.R  
- gencode.v38.annotation.gtf  
- ggsashimi_txV3.py  
- GRCh38_appris_data.principal.txt  
- Homo_sapiens.GRCh38.103.chr.sorted_new.gtf  
- palette.txt  
- pgp-a.sh  
- run_sashimiV1.sh  
- selected_events.csv  
- TxEnsDB103_layeredV6.R  

### Option 0  
```bash
bash pgp-a.sh 0 > pgp-a-0.txt 2> pgp-a-0.error.txt
```

### Option 1  
Dans cette option, il y avoir plusieurs points à revoir par rapport au pipeline original :  
- Installation des librairies suivantes (lancer R dans l'environnement conda) : BiocManager, RCurl, GenomeInfoDB,
et GenomicRanges à l'aide des commandes suivantes :
```R  
install.packages("BiocManager")
```
Puis :  
```bash
conda install -c r r-devtools r-bitops r-rcurl
```
Et enfin :  
```R
BiocManager::install("GenomicRanges")
BiocManager::install("GenomeInfoDB")
```

- Suppression des gènes TNNI3 et DNAAF3 de LGR 432 dans le fichier ”principal txs.csv” (en s'aidant de la commande suivante :)
```bash
grep -w -v TNNI3 principal_txs.csv
```
  
- Dans Tx_EnsDB103_layeredV6.R, remplacer les '\n\t+' et les '\n +' par des '\n' (utiliser BBEdit par exemple)
- Exécuter Tx_EnsDB103_layeredV6.R directement dans le terminal (20 lignes par 20 lignes)

  
Pour finir, on peut exécuter la commande suivante :  
```bash
bash pgp-a.sh 1 selected_events.csv > pgp-a-1.txt 2> pgp-a-1.error.txt
```

### Option 2
On peut exécuter la commande suivante :  
```bash
bash pgp-a.sh 2 selected_events.csv > pgp-a-2.txt 2> pgp-a-2.error.txt
```

### Option 3
On peut exécuter la commande suivante (environ 45min de temps d'exécution) :  
```bash
bash pgp-a.sh 3 > pgp-a-3.txt 2> pgp-a-3.error.txt
```

### Option 4
On peut exécuter la commande suivante (entre 10 et 15 jours de temps d'exécution) :  
```bash
bash pgp-a.sh 4 selected_events.csv > pgp-a-4.txt 2> pgp-a-4.error.txt
```


## Partie B
Les fichiers nécessaires à l'exécution de cette partie sont les suivant (à avoir dans le dossier 'pgp-b') :  
- all_bams.tsv  
- Auto_CoverV4_layered_intronV3.R
- CE_inclusion.csv (à récupérer dans la partie A après exécution)  
- check_aaV4_allFrames.R  
- clean_selected_events.csv (à récupérer dans la partie A après exécution)  
- esV5_layered_CDSV3.sh  
- gencode.v38.annotation.gtf  
- get_orf_cds.R  
- ggsashimi_txV3.py
- GRCh38_appris_data.principal.txt
- GRCh38.p13.genome.fa
- Homo_sapiens.GRCh38.103.gtf
- merge_sashimis.py
- palette.txt
- pgp_b_ce_ir.sh
- pgp_b.sh
- principal_txs.csv (à récupérer dans la partie A après exécution)
- run_sashimiV1.sh
- selected_events.csv
- TxEnsDB103_layeredV6.R
  
### Première exécution
On peut exécuter la commande suivante :  
```bash
bash pgp_b.sh selected_events.csv > pgp_b_es.txt 2>pgp_b_es.error.txt
```
Ou :
```bash
bash pgp_b.sh sorted_selected_events.csv > pgp_b_es.txt 2>pgp_b_es.error.txt
```

<br><br><br>

# Pour exécuter mon pipeline
On commence de la même manière que précédemment en installant l'environnement conda et en s'assurant que la version de R sur cet environnement est bien '4.3.3'.  

## Partie A
