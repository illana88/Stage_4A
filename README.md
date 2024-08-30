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
- Suppression des gènes TNNI3 et DNAAF3 de LGR 432 dans le fichier ”principal txs.csv”
- Dans Tx_EnsDB103_layeredV6.R, remplacer les '\n\t+' et les '\n +' par des '\n'
- Exécuter Tx_EnsDB103_layeredV6.R directement dans le terminal (20 lignes par 20 lignes)




```bash
bash pgp-a.sh 1 selected_events.csv > pgp-a-1.txt 2> pgp-a-1.error.txt
```

