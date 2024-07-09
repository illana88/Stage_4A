######## Proteogenomic Pipeline for Biomarker Discovery in iPSC neurons
######### Inline Readme
# 1. Get List of most abundant Transcripts from KD iPSC samples
# This step assumes that all bam files are in current folder and follow have
# change directory to part-a and issue run part-a.sh with appropriate option and input file

arg1 = int(input("Argument 1 :"))

import os
import glob
import csv
import subprocess
import re
import pandas as pd
import sys
from pyensembl import EnsemblRelease
from bioservices import BioMart
from pybedtools import BedTool

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
            
            print('type de cov dans all_ddf à la création : ', all_ddf['cov'].dtype)
            
            for i in range(len(file_list)):
                print(f"reading file : {file_list[i]}")
                rec = pd.read_csv(f"iPSC_gtfs/{file_list[i]}", sep="\s+", header=None, dtype=str)
                rec.columns = ["TxID","GeneID","Gene_Name","cov","FPKM","TPM"]
                
                rec = rec.astype({
                    'TxID': 'str',
                    'GeneID': 'str',
                    'Gene_Name': 'str',
                    'cov': 'float',
                    'FPKM': 'float',
                    'TPM': 'float'
                })
                
                all_ddf = pd.concat([all_ddf, rec], ignore_index=True)
                
            print('type de cov dans all_ddf après concatenate : ', all_ddf['cov'].dtype)
            
            # # Remove spaces
            # def remove_spaces(s):
            #     if isinstance(s, str):
            #         return re.sub(r' ', '', s)
                
            # all_ddf = all_ddf.apply(lambda col: col.map(remove_spaces))

            # print('type de cov dans all_ddf après remove spaces : ', all_ddf['cov'].dtype)
            # pd.set_option('display.max_columns', None)
            # print('all_ddf après remove spaces : ', all_ddf.head())
            
            # And sort for fast retrieval
            all_ddf = all_ddf.sort_values(by='Gene_Name')
            print('type de cov dans all_ddf après sorted : ', all_ddf['cov'].dtype)
            print('all_ddf après sorted : ', all_ddf.head())
            
            # Get unique gene names
            unique_genes = pd.unique(all_ddf['Gene_Name']) ##### Testé jusqu'ici et tout fonctionne #####
            
            Tx_ddf = pd.DataFrame({
                'TxID': pd.Series(dtype='str'),
                'GeneID': pd.Series(dtype='str'),
                'Gene_Name': pd.Series(dtype='str')
            })
            
            # Select row with max cov for each gene
            print('Now Generating Txs Table, Will take a while !!!!!!!')
            
            print('type de cov dans all_ddf final : ', all_ddf['cov'].dtype)
            print('all_ddf final : ', all_ddf.head())
            
            for i in range (len(unique_genes)):
                TxSubset = all_ddf[all_ddf['Gene_Name'] == unique_genes[i]]
                #print('toutes les cov : ', TxSubset['cov'])
                #print('cov max : ', TxSubset['cov'].max())
                tx_gene = TxSubset[TxSubset['cov'] == TxSubset['cov'].max()][0:3]
                Tx_ddf = pd.concat([Tx_ddf, tx_gene], ignore_index=True)
                
            print('Tx_ddf : ', Tx_ddf)
            
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
#     arg2 = "selected_events.csv" # arg2 = input("Argument 2 :")
#     if os.access(arg2 ,os.F_OK):
#         splicing_events_file = arg2.split('.')[0]
#         sorted_selected_events = pd.read_csv(arg2, delimiter=',')
#         print('verification of column 5 : ', sorted_selected_events.columns[4])
#         sorted_sorted_selected_events = sorted_selected_events.sort_values(by=sorted_selected_events.columns[4])
#         sorted_sorted_selected_events.to_csv(sorted_sorted_selected_events.csv, index=False)
#     else :
#         print('PLEASE PROVIDE SPLICING_EVENTS csv file  bash pgp_0.sh flag splicing_events.csv and RERUN')
#         with open("Summary_stats.txt","a") as fichier :
#             fichier.write('PLEASE PROVIDE SPLICING_EVENTS csv file  bash pgp_0.sh flag splicing_events.csv and RERUN\n')
#             fichier.write('NOW EXITING\n')
#         sys.exit(1)
#     with open('sorted_sorted_selected_events.csv', 'r') as file:
#         lines = file.readlines()
#         total_events = len(lines)
#     print('CAME IN FOR SASHIMI PLOTS FOR ALL $total_events EVENTS')
#     with open("Summary_stats.txt","a") as fichier :
#         fichier.write('CAME IN FOR SASHIMI PLOTS FOR ALL $total_events EVENTS\n')
#     # ALSO get bed file for sashimi plots for all events.
#     # FORMAT IS: chr, start of up_ex,end of ds_ex, 1, 0, strand, gene_name, TxID
#     generate_all_events_sashimi_beds = 1
#     if generate_all_events_sashimi_beds == 1 :
#         os.makedirs("temp_all_events_sashimi",exist_ok=True)
#         print('################################ GENERATING BED FILE FOR SASHIMI PLOTS FOR ALL EVENTS')
#         with open("Summary_stats.txt","a") as fichier :
#             fichier.write('################################ GENERATING BED FILE FOR SASHIMI PLOTS FOR ALL EVENTS\n')
#         # SHOULD ADD CHECKS ON EXON_NUM READING FROM FILE
#         if len(os.listdir("temp_all_events_sashimi"))!=0 :
#             for f in glob.glob("temp_all_events_sashimi/*.*") :
#                 os.remove(f)
#         generate_event_beds = 1
#         if generate_event_beds == 1 :
#             os.makedirs("event_bedfiles",exist_ok=True)
#             if len(os.listdir("event_bedfiles"))!=0 :
#                 for f in glob.glob("event_bedfiles/*.*") :
#                     os.remove(f)
#             with open("Summary_stats.txt","a") as fichier :
#                 fichier.write('************************ From pgp-a.sh, CALLING TxEnsDB103_layeredV6.R \n')
            
            
            
            
#             ########## TxEnsDB103_layeredV6.R codé en Pyhton ##########
#             # Using ENSG (commented out this - junction start and end range) to get transcripts
#             # Starting from previous version (TxEnsDB103_layeredV1.R), it has 41 events for which Tx selected does not encapsulate the event
#             # An example is event: chr12	22515151	22517985	AC053513.1
#             # So coordinate based search is back on to see if it makes difference
#             # THIS SCRIPT HAS BEEN MODIFIED - SO REPLACE IT IN ALL VERSIONS
#             # 04/20/2022 - removing 5'utr
#             ensembl = EnsemblRelease(103)
#             transcripts = ensembl.transcripts()
#             transcript_data = []
            
#             for tx in transcripts:
#                 transcript_data.append({
#                     'transcript_id': tx.transcript_id,
#                     'gene_id': tx.gene_id,
#                     'exon_lengths': tx.exon_lengths,
#                     '5_utr_length': tx.five_prime_utr_length,
#                     '3_utr_length': tx.three_prime_utr_length
#                 })
            
#             tx_lens = pd.DataFrame(transcript_data)
#             GeneIDField = 6
#             # Read Peaks File
#             SpliceData = pd.read_csv('sorted_sorted_selected_events.csv', header=None)
#             # Also read Tx list
#             Tx_list = pd.read_csv('principal_txs.csv', header=None)
#             # Also read appris annoation data to get principal 1 isoform
#             appr_anno1 = pd.read_csv("GRCh38_appris_data.principal.txt", sep="\t", header=None)
#             # V1 - hugo_symbol, V2 - ENSG_ID, V3 - TX_ID, V4 - , V2 - PRINCIPAL/ALTERNATE
#             appr_anno1.columns = ['V1', 'V2', 'V3', 'V4', 'V5']
#             # Select only PRINCIPAL.1 ISOFORMS
#             appr_anno = appr_anno1[appr_anno1['V5'].str.contains("PRINCIPAL:1")]
            
#             if os.access("all_tx_events.csv",os.F_OK):
#                 os.remove(("all_tx_events.csv"))
                
#             if os.access("all_events_bed_sashimi.tab",os.F_OK):
#                 os.remove(("all_events_bed_sashimi.tab"))
                
#             if os.access("events_to_tx_mapping_valid.csv",os.F_OK):
#                 os.remove(("events_to_tx_mapping_valid.csv"))
                
#             if os.access("events_to_tx_mapping_invalid.csv",os.F_OK):
#                 os.remove(("events_to_tx_mapping_invalid.csv"))
            
#             with open("temp_all_events_sashimi/FINAL_STATS_ALL_SASHIMIS.txt","a") as fichier :
#                 fichier.write('                                 \n')
#                 fichier.write(f"Starting From TxEnsDB103_layeredV6.R: --------------- Processing file: 'sorted_sorted_selected_events.csv' with: {SpliceData.shape[0]} events to generate each event .bed files in event_bedfiles/ folder: \n")
            
#             print("Started Generating BED files for Splicing Events in folder event_bedfiles/ from File: 'sorted_sorted_selected_events.csv'")
            
#             trackj = 1
#             temp_gene = ""
#             current_gene = ""
#             tx_lengths = []
            
#             df_notfound = pd.DataFrame({
#                 'seqnames': pd.Series(dtype='str'),
#                 'start': pd.Series(dtype='numeric'),
#                 'end': pd.Series(dtype='numeric'),
#                 'strand': pd.Series(dtype='str'),
#                 'genename': pd.Series(dtype='str'),
#                 'junc_type': pd.Series(dtype='str')
#             })
            
#             df_zeroutr = pd.DataFrame({
#                 'seqnames': pd.Series(dtype='str'),
#                 'start': pd.Series(dtype='numeric'),
#                 'end': pd.Series(dtype='numeric'),
#                 'strand': pd.Series(dtype='str'),
#                 'genename': pd.Series(dtype='str'),
#                 'junc_type': pd.Series(dtype='str')
#             })
            
#             repeated_entries = 0
#             iPSC_events = 0
#             appris_events = 0
#             principalTx_events = 0
#             events_xTx = 0
#             Tx_str = 0 # 0 for iPSC, 1 for APPRIS and 2 for EnsDB
#             Tx_valid = 0
#             Total_Events = SpliceData.shape[0]
#             probable_noise_events = 0
#             probable_noncoding_events = 0
#             utr5_events = 0
            
#             for i in range(SpliceData.shape[0]) :
#                 # Step 0 - get transcripts for each gene (MAJIQ only reports for some)
#                 # Deal with multiple events for the same gene
                
#                 if "Tx_name" in globals() :
#                     del Tx_name
                
#                 if i == 1 :
#                     temp_gene = SpliceData[i,5]
#                     trackj = 1
                
#                 elif temp_gene == SpliceData[i,5] :
#                     trackj = trackj+1
                    
#                 else :
#                     temp_gene = SpliceData[i,5]
#                     trackj = 1
                
#                 # Get gene name and gene_id (from granges filter) using event coordinates
                                
#             ########## FIN TxEnsDB103_layeredV6.R codé en Pyhton ##########