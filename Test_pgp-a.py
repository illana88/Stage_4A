######## Proteogenomic Pipeline for Biomarker Discovery in iPSC neurons
######### Inline Readme
# 1. Get List of most abundant Transcripts from KD iPSC samples
# This step assumes that all bam files are in current folder and follow have
# change directory to part-a and issue run part-a.sh with appropriate option and input file

import os
import glob
import csv
import subprocess
import re
import pandas as pd
import sys

arg1 = int(input("Argument 1 :"))
if os.access("Summary_stats.txt",os.F_OK):
    os.remove(("Summary_stats.txt"))

if arg1==0 or arg1==2 or arg1==5 :
    with open("Summary_stats.txt","w") as fichier :
        fichier.write("STARTING pgp-a.sh WITH FLAG 0\n")
    phrase = "STARTED SringTie Calculations - Will take a while depending on number of samples\n"
    print(phrase)
    with open("Summary_stats.txt","a") as fichier :
        fichier.write(phrase)
    os.makedirs("iPSC_gtfs",exist_ok=True)
    if len(os.listdir("iPSC_gtfs"))!=0 :
        for f in glob.glob("iPSC_gtfs/*.*") :
            os.remove(f)
    # First get bam files
    file = open("all_bams.tsv")
    reader = csv.reader(file,delimiter='\t') 
    samples = []
    for row in reader :
        if "TDP43" in row[2] :
            samples.append(row[1])
    # Check if we got any samples
    if "samples" in globals() and len(samples)!=0 :
        for sample in samples :
            print(f"STARTED PROCESSING SAMPLE {sample}")
            with open("Summary_stats.txt","a") as fichier :
                fichier.write(f"STARTED PROCESSING SAMPLE {sample}\n")
            samp = os.path.basename(sample)
            samp = samp.split('.')[0]
            commande_shell = f"stringtie {sample} -p 8 -G gencode.v38.annotation.gtf -o iPSC_gtfs/{samp}.gtf"
            subprocess.run(commande_shell, shell=True, check=True)
        if len(os.listdir("iPSC_gtfs"))!=0 :
            print("DONE WITH SringTie Calculations - NOW GENERATING LIST of ABUNDANT Txs")
            with open("Summary_stats.txt","a") as fichier :
                fichier.write("DONE WITH SringTie Calculations - NOW GENERATING LIST of ABUNDANT Txs\n")
            # Now get unique transcripts in csv files for each sample
            for i in glob.glob("iPSC_gtfs/*.gtf") :
                if (os.path.isfile(i)) :
                    samp = i.split('/')[1].split('.')[0]
                    print(f"PROCESSING SAMPLE {i}")
                    with open("Summary_stats.txt","a") as fichier :
                        fichier.write(f"PROCESSING SAMPLE {i}\n")
                    # Step 1. Get all Txs having reference_id (ENST), ref_gene_id (ENSG) and ref_gene_name (HUGO SYMBOL), these are 24 column lines in stringtie's gtf file
                    with open(i, 'r') as file:
                        for line in file:
                            tab = line.replace(' ', '\t')
                            tab = tab.strip().split('\t')
                            if len(tab)==24 :
                                col = [tab[13], tab[15], tab[17], tab[19], tab[21], tab[23]]
                                new_col = re.sub(r';', '\t', "\t".join(col))
                                with open(f"iPSC_gtfs/{samp}.csv","a") as fichier :
                                    fichier.write(new_col + "\n")
                else :
                    break
            # Now select most abundant transcripts from all samples
            print("Now calling abundant_tx.R")
            with open("Summary_stats.txt","a") as fichier :
                fichier.write("Now calling abundant_tx.R\n")
                
                
            ########## abundant_tx.R codé en Pyhton ##########
            if os.access("principal_tx.csv",os.F_OK):
                # Delete file if it exists
                os.remove(("principal_tx.csv"))
            # Create a list of the files from your target directory
            file_list = [f for f in os.listdir("iPSC_gtfs/") if f.endswith('.csv')]
            ddf = pd.DataFrame({
                'TxID': pd.Series(dtype='str'),
                'GeneID': pd.Series(dtype='str'),
                'Gene_Name': pd.Series(dtype='str'),
                'cov': pd.Series(dtype='float'),
                'FPKM': pd.Series(dtype='float'),
                'TPM': pd.Series(dtype='float')
            })
            all_ddf = pd.DataFrame({
                'TxID': pd.Series(dtype='str'),
                'GeneID': pd.Series(dtype='str'),
                'Gene_Name': pd.Series(dtype='str'),
                'cov': pd.Series(dtype='float'),
                'FPKM': pd.Series(dtype='float'),
                'TPM': pd.Series(dtype='float')
            })
            for i in range(len(file_list)):
                print(f"reading file : {file_list[i]}")
                rec = pd.read_csv(f"iPSC_gtfs/{file_list[i]}", sep="\t", header=None, dtype=str, encoding='utf-8')
                print("rec : ", rec)
                rec = rec.drop(columns=[7])
                rec.columns = ["TxID","GeneID","Gene_Name", "cov","FPKM","TPM"] + list(rec.columns[6:])
                all_ddf = pd.concat([all_ddf, rec], ignore_index=True) ##### Testé jusqu'ici et tout fonctionne #####
#             # Remove spaces
#             print("all_ddf before : ", all_ddf)
#             def remove_spaces(s):
#                 if isinstance(s, str):
#                     return re.sub(r' ', '', s)
#             all_ddf = all_ddf.apply(lambda col: col.map(remove_spaces)) # all_ddf.applymap(remove_spaces)
#             print("all_ddf after revoming space : ", all_ddf)
#             # And sort for fast retrieval
#             all_ddf = all_ddf.sort_values(by='Gene_Name')
#             print("all_ddf after sorted : ", all_ddf)
#             # Get unique gene names
#             unique_genes = pd.unique(all_ddf['Gene_Name'])
#             print("unique_genes : ", unique_genes)
#             Tx_ddf = pd.DataFrame({
#                 'TxID': pd.Series(dtype='str'),
#                 'GeneID': pd.Series(dtype='str'),
#                 'Gene_Name': pd.Series(dtype='str') ##### Testé jusqu'ici et à vérifier si ça fonctionne #####
#             })
#             # Select row with max cov for each gene
#             print('Now Generating Txs Table, Will take a while !!!!!!!')
#             for i in range (len(unique_genes)):
#                 TxSubset = all_ddf[all_ddf['Gene_Name'] == unique_genes[i]]
#                 tx_gene = TxSubset[TxSubset['cov'] == TxSubset['cov'].max()][0:3]
#                 Tx_ddf = pd.concat([Tx_ddf, tx_gene], ignore_index=True)
#             # Write in file
#             # Also remove trailing version numbers
#             Tx_ddf['GeneID'] = Tx_ddf['GeneID'].str.split('.', n=1, expand=True)[0]
#             Tx_ddf['TxID'] = Tx_ddf['TxID'].str.split('.', n=1, expand=True)[0]
#             Tx_ddf.to_csv("principal_txs.csv", index=False, quoting=pd.QUOTE_NONE, sep=',', header=False)
#             print('Done With Txs Table: principal_txs.csv')
#             ########## FIN abundant_tx.R codé en Pyhton ##########
            
            
#             print("DONE WITH ABUNDANT Txs file Generation (principal_txs.csv) for input BAM SAMPLES")
#             with open("Summary_stats.txt","a") as fichier :
#                 fichier.write("DONE WITH ABUNDANT Txs file Generation (principal_txs.csv) for input BAM SAMPLES\n")
#                 fichier.write("FINISHED pgp-a.sh WITH FLAG 0\n")
#             if arg1 == 0 :
#                 print("++++++NOW Exiting+++++, PLEASE CHECK Summary_stats.txt file for details")
#                 with open("Summary_stats.txt","a") as fichier :
#                     fichier.write("++++++NOW Exiting+++++, PLEASE CHECK Summary_stats.txt file for details\n")
#                 sys.exit(1)
#             else : # ABORT if we got no samples
#                 print("No .gtf FILE GENERATED BY StringTie2 (found empty iPSC_gtfs folder) ABORTING!!!! Please check BAM samples path and Rerun again")
#                 with open("Summary_stats.txt","a") as fichier :
#                     fichier.write("No .gtf FILE GENERATED BY StringTie2 (found empty iPSC_gtfs folder) ABORTING!!!! Please check BAM samples path and Rerun again\n")
#         else :
#             print("!!!!!!!!!!!Something went wrong with SringTie Calculations (Please check all_bams.tsv for BAM FILE PATHS) - Please correct the error and re-run again !!!!!!!!!!!")
#             with open("Summary_stats.txt","a") as fichier :
#                 fichier.write("!!!!!!!!!!!Something went wrong with SringTie Calculations (Please check all_bams.tsv for BAM FILE PATHS) - Please correct the error and re-run again !!!!!!!!!!!\n")
#             print("++++++NOW Exiting+++++")
#             with open("Summary_stats.txt","a") as fichier :
#                 fichier.write("++++++NOW Exiting+++++\n")
#             sys.exit(1)
            
            
            
            
# # NEXT SORT FILE ON COLUMN 5 - THIS IS IMPORTANT AS PIPELINE DEPENDS ON IT
# if arg1==1 or arg1==2 or arg1==5 :
#     with open("Summary_stats.txt","a") as fichier :
#         fichier.write("STARTING pgp-a.sh WITH FLAG 1 for all events sashimi plots\n")
#     arg2 = input("Argument 2 :") # arg2 = "sorted_selected_events.csv
#     if os.access(arg2 ,os.F_OK):
#         splicing_events_file = arg2.split('.')[0]
        