# Pour exécuter le pipeline original
Lien du pipeline original : https://github.com/NIH-CARD/proteogenomic-pipeline  
<br>
Pour commencer, il faut installer conda (par exemple via ce lien : https://docs.anaconda.com/free/miniconda/miniconda-other-installer-links/), puis créer l'environnement comme indiqué sur le pipeline original. Enfin, exéctuer les commandes suivantes :  
```bash
conda config --add channels bioconda  
conda config --add channels conda-forge
```

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
On peut exécuter la commande suivante :  
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

<br><br>

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

<br><br><br><br>

# Pour exécuter mon pipeline
On commence de la même manière que précédemment en installant l'environnement conda et en s'assurant que la version de R sur cet environnement soit bien '4.3.3'.  

## Partie A
Les fichiers nécessaires à l'exécution de cette partie sont les suivant (à avoir dans le dossier 'pgp-a') :  
- abundant_tx.py  
- all_bams.tsv  
- Auto_CoverV4_layered_intronV3.py  
- gencode.v38.annotation.gtf  
- ggsashimi_txV3.py  
- GRCh38_appris_data.principal.txt  
- Homo_sapiens.GRCh38.103.chr.sorted_new.gtf  
- palette.txt  
- pgp-a.py  
- run_sashimiV1.py  
- selected_events.csv  
- TxEnsDB103_layeredV6.py

### Options
On peut exécuter la commande suivante :  
```python
python3 pgp-a.py
```

Puis la ligne suivante va apparaître :  
'Option : '  
 <br>
Il faudra alors indiquer l'entier correspondant à l'option que l'on souhaite exécuter (de 0 à 5).  
<br>
De plus, si on exécute les options 0, 2, 3, 4, ou 5, il apparaîtra la ligne suivante :  
'protein : '  
<br>
Il faudra ici indiquer le nom de la protéine que l'on souhaite étudier (par exemple : 'TDP43').  

### Vérifications
La partie A a entièrement été retranscrite en Python (tous les fichiers de code), et a été testée en grande partie (toutes les options sauf la fin de l'option 4 ainsi que la fin du fichier TxEnsDB103_layeredV6.py - chercher en ctrl+F "##### Teste jusqu'ici et tout fonctionne #####").  
<br>
Il y a dans la fin de l'option 4 des appels à des variables nommées "usd_abs" ainsi que "dsd_abs" qui ne sont pas définies (j'ai simplement suivi le code du pipeline original, donc c'est peut-être une erreur ou pas le bon nom des variables).  

<br><br>

## Partie B
Les fichiers nécessaires à l'exécution de cette partie sont les suivant (à avoir dans le dossier 'pgp-b') :  
- all_bams.tsv  
- Auto_CoverV4_layered_intronV3.py
- CE_inclusion.csv (à récupérer dans la partie A après exécution)  
- check_aaV4_allFrames.R  
- clean_selected_events.csv (à récupérer dans la partie A après exécution)  
- esV5_layered_CDSV3.py  
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
- pgp_b.py (seulement pour le début)
- principal_txs.csv (à récupérer dans la partie A après exécution)
- run_sashimiV1.py
- selected_events.csv
- TxEnsDB103_layeredV6.py  

### Première exécution
On peut exécuter la commande suivante :  
```python
python3 pgp_b.py
```

### Vérifications  
Pour cette partie, je n'ai eu le temps de retranscrire que le début de pgp_b.sh ainsi que tout le fichier esV5_layered_CDSV3.sh en Python, et je n'ai eu le temps que de tester le début de pgp_b.py.  
