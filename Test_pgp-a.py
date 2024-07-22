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
import pybedtools
from io import StringIO

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
                rec = pd.read_csv(f"iPSC_gtfs/{file_list[i]}", sep=r"\s+", header=None, dtype=str)
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
            
            # # Remove spaces -> vraiment nécessaire ???
            # def remove_spaces(s):
            #     if isinstance(s, str):
            #         return re.sub(r' ', '', s)
                
            # all_ddf = all_ddf.apply(lambda col: col.map(remove_spaces))

            # print('type de cov dans all_ddf après remove spaces : ', all_ddf['cov'].dtype)
            # pd.set_option('display.max_columns', None)
            # print('all_ddf après remove spaces : ', all_ddf.head())
            
            # And sort for fast retrieval
            all_ddf = all_ddf.sort_values(by='Gene_Name')
            
            # Get unique gene names
            unique_genes = pd.unique(all_ddf['Gene_Name'])
            
            Tx_ddf = pd.DataFrame({
                'TxID': pd.Series(dtype='str'),
                'GeneID': pd.Series(dtype='str'),
                'Gene_Name': pd.Series(dtype='str')
            })
            
            # Select row with max cov for each gene
            print('Now Generating Txs Table, Will take a while !!!!!!!')
            
            for i in range (len(unique_genes)):
                TxSubset = all_ddf[all_ddf['Gene_Name'] == unique_genes[i]]
                tx_gene = TxSubset[TxSubset['cov'] == TxSubset['cov'].max()]
                tx_gene = tx_gene.iloc[:, :3]
                Tx_ddf = pd.concat([Tx_ddf, tx_gene], ignore_index=True)
            
            # Write in file
            # Also remove trailing version numbers
            Tx_ddf['GeneID'] = Tx_ddf['GeneID'].str.split('.', n=1, expand=True)[0]
            Tx_ddf['TxID'] = Tx_ddf['TxID'].str.split('.', n=1, expand=True)[0]
            Tx_ddf.to_csv("principal_txs.csv", index=False, sep=',', header=False)
            print('Done With Txs Table: principal_txs.csv')

            ########## FIN abundant_tx.R codé en Pyhton ##########
            
            
            
            
            print("DONE WITH ABUNDANT Txs file Generation (principal_txs.csv) for input BAM SAMPLES")
            with open("Summary_stats.txt","a") as fichier :
                fichier.write("DONE WITH ABUNDANT Txs file Generation (principal_txs.csv) for input BAM SAMPLES\n")
                fichier.write("FINISHED pgp-a.sh WITH FLAG 0\n")
            if arg1 == 0 :
                print("++++++NOW Exiting+++++, PLEASE CHECK Summary_stats.txt file for details")
                with open("Summary_stats.txt","a") as fichier :
                    fichier.write("++++++NOW Exiting+++++, PLEASE CHECK Summary_stats.txt file for details\n")
                sys.exit(1)
            else : # ABORT if we got no samples
                print("No .gtf FILE GENERATED BY StringTie2 (found empty iPSC_gtfs folder) ABORTING!!!! Please check BAM samples path and Rerun again")
                with open("Summary_stats.txt","a") as fichier :
                    fichier.write("No .gtf FILE GENERATED BY StringTie2 (found empty iPSC_gtfs folder) ABORTING!!!! Please check BAM samples path and Rerun again\n")
        else :
            print("!!!!!!!!!!!Something went wrong with SringTie Calculations (Please check all_bams.tsv for BAM FILE PATHS) - Please correct the error and re-run again !!!!!!!!!!!")
            with open("Summary_stats.txt","a") as fichier :
                fichier.write("!!!!!!!!!!!Something went wrong with SringTie Calculations (Please check all_bams.tsv for BAM FILE PATHS) - Please correct the error and re-run again !!!!!!!!!!!\n")
            print("++++++NOW Exiting+++++")
            with open("Summary_stats.txt","a") as fichier :
                fichier.write("++++++NOW Exiting+++++\n")
            sys.exit(1) ##### Testé jusqu'ici et tout fonctionne #####
            
            
            
            
            
            
# NEXT SORT FILE ON COLUMN 5 - THIS IS IMPORTANT AS PIPELINE DEPENDS ON IT
if arg1==1 or arg1==2 or arg1==5 :
    with open("Summary_stats.txt","a") as fichier :
        fichier.write("STARTING pgp-a.sh WITH FLAG 1 for all events sashimi plots\n")
        
    arg2 = "selected_events.csv" # arg2 = input("Argument 2 :")
    
    if os.access(arg2 ,os.F_OK):
        splicing_events_file = arg2.split('.')[0]
        selected_events = pd.read_csv(arg2, delimiter=',')
        sorted_selected_events = selected_events.sort_values(by=selected_events.columns[4])
        sorted_selected_events.to_csv("sorted_selected_events.csv", index=False)
        
    else :
        print('PLEASE PROVIDE SPLICING_EVENTS csv file  bash pgp_0.sh flag splicing_events.csv and RERUN')
        
        with open("Summary_stats.txt","a") as fichier :
            fichier.write('PLEASE PROVIDE SPLICING_EVENTS csv file  bash pgp_0.sh flag splicing_events.csv and RERUN\n')
            fichier.write('NOW EXITING\n')
            
        sys.exit(1)
        
    with open('sorted_selected_events.csv', 'r') as file:
        lines = file.readlines()
        total_events = len(lines)
        
    print(f"CAME IN FOR SASHIMI PLOTS FOR ALL {total_events} EVENTS")
    
    with open("Summary_stats.txt","a") as fichier :
        fichier.write(f"CAME IN FOR SASHIMI PLOTS FOR ALL {total_events} EVENTS\n")
        
    # ALSO get bed file for sashimi plots for all events.
    # FORMAT IS: chr, start of up_ex,end of ds_ex, 1, 0, strand, gene_name, TxID
    generate_all_events_sashimi_beds = 1
    
    if generate_all_events_sashimi_beds == 1 :
        os.makedirs("temp_all_events_sashimi",exist_ok=True)
        print('################################ GENERATING BED FILE FOR SASHIMI PLOTS FOR ALL EVENTS')
        
        with open("Summary_stats.txt","a") as fichier :
            fichier.write('################################ GENERATING BED FILE FOR SASHIMI PLOTS FOR ALL EVENTS\n')
            
        # SHOULD ADD CHECKS ON EXON_NUM READING FROM FILE
        
        if len(os.listdir("temp_all_events_sashimi"))!=0 :
            for f in glob.glob("temp_all_events_sashimi/*.*") :
                os.remove(f)
        generate_event_beds = 1
        if generate_event_beds == 1 :
            os.makedirs("event_bedfiles",exist_ok=True)
            if len(os.listdir("event_bedfiles"))!=0 :
                for f in glob.glob("event_bedfiles/*.*") :
                    os.remove(f)
            with open("Summary_stats.txt","a") as fichier :
                fichier.write('************************ From pgp-a.sh, CALLING TxEnsDB103_layeredV6.R \n')
            
            
            
            
            ########## TxEnsDB103_layeredV6.R codé en Pyhton ##########

            command = [
                "Rscript",
                "txens.R",
                "sorted_selected_events.csv",
                "principal_txs.csv",
                "temp_all_events_sashimi/FINAL_STATS_ALL_SASHIMIS.txt"
            ]

            subprocess.run(command, capture_output=True, text=True) ##### Testé jusqu'ici et tout fonctionne #####

            # Using ENSG (commented out this - junction start and end range) to get transcripts
            # Starting from previous version (TxEnsDB103_layeredV1.R), it has 41 events for which Tx selected does not encapsulate the event
            # An example is event: chr12	22515151	22517985	AC053513.1
            # So coordinate based search is back on to see if it makes difference
            # THIS SCRIPT HAS BEEN MODIFIED - SO REPLACE IT IN ALL VERSIONS
            # 04/20/2022 - removing 5'utr

            # def download_ftp_file(ftp_url, local_file):
            #     ftp = ftplib.FTP("ftp.ensembl.org")
            #     ftp.login()
            #     ftp.cwd(ftp_url)
            #     with open(local_file, "wb") as f:
            #         ftp.retrbinary("RETR " + local_file, f.write)
            #     ftp.quit()

            # files = {
            #     "gtf": "ftp://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz",
            #     "transcript_fasta": "ftp://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
            #     "protein_fasta": "ftp://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
            # }

            # paths = {
            #     "gtf": "Homo_sapiens.GRCh38.103.gtf.gz",
            #     "transcript_fasta": "Homo_sapiens.GRCh38.cdna.all.fa.gz",
            #     "protein_fasta": "Homo_sapiens.GRCh38.pep.all.fa.gz"
            # }

            # for file_type, ftp_path in files.items():
            #     local_path = paths[file_type]
            #     if not os.path.exists(local_path):
            #         print(f"Téléchargement de {ftp_path}...")
            #         download_ftp_file(os.path.dirname(ftp_path), os.path.basename(local_path))
            #         print(f"{local_path} téléchargé avec succès.")
            #     else:
            #         print(f"{local_path} existe déjà. Téléchargement ignoré.")

            # GeneIDField = 6
            
            # # Read Peaks File
            # SpliceData = pd.read_csv('sorted_sorted_selected_events.csv', header=None)
            
            # # # Also read Tx list
            # Tx_list = pd.read_csv('principal_txs.csv', header=None)
            
            # # Also read appris annoation data to get principal 1 isoform
            # appr_anno1 = pd.read_csv("GRCh38_appris_data.principal.txt", sep="\t", header=None)
            # # V1 - hugo_symbol, V2 - ENSG_ID, V3 - TX_ID, V4 - , V2 - PRINCIPAL/ALTERNATE
            # appr_anno1.columns = ['V1', 'V2', 'V3', 'V4', 'V5']
            # # Select only PRINCIPAL.1 ISOFORMS
            # appr_anno = appr_anno1[appr_anno1['V5'].str.contains("PRINCIPAL:1")]
            
            # if os.access("all_tx_events.csv",os.F_OK):
            #     os.remove(("all_tx_events.csv"))
                
            # if os.access("all_events_bed_sashimi.tab",os.F_OK):
            #     os.remove(("all_events_bed_sashimi.tab"))
                
            # if os.access("events_to_tx_mapping_valid.csv",os.F_OK):
            #     os.remove(("events_to_tx_mapping_valid.csv"))
                
            # if os.access("events_to_tx_mapping_invalid.csv",os.F_OK):
            #     os.remove(("events_to_tx_mapping_invalid.csv"))
            
            # with open("temp_all_events_sashimi/FINAL_STATS_ALL_SASHIMIS.txt","a") as fichier :
            #     fichier.write('                                 \n')
            #     fichier.write(f"Starting From TxEnsDB103_layeredV6.R: --------------- Processing file: 'sorted_sorted_selected_events.csv' with: {SpliceData.shape[0]} events to generate each event .bed files in event_bedfiles/ folder: \n")
            
            # print("Started Generating BED files for Splicing Events in folder event_bedfiles/ from File: 'sorted_sorted_selected_events.csv'")
            
            # trackj = 1
            # temp_gene = ""
            # current_gene = ""
            # tx_lengths = []
            
            # df_notfound = pd.DataFrame({
            #     'seqnames': pd.Series(dtype='str'),
            #     'start': pd.Series(dtype='numeric'),
            #     'end': pd.Series(dtype='numeric'),
            #     'strand': pd.Series(dtype='str'),
            #     'genename': pd.Series(dtype='str'),
            #     'junc_type': pd.Series(dtype='str')
            # })
            
            # df_zeroutr = pd.DataFrame({
            #     'seqnames': pd.Series(dtype='str'),
            #     'start': pd.Series(dtype='numeric'),
            #     'end': pd.Series(dtype='numeric'),
            #     'strand': pd.Series(dtype='str'),
            #     'genename': pd.Series(dtype='str'),
            #     'junc_type': pd.Series(dtype='str')
            # })
            
            # repeated_entries = 0
            # iPSC_events = 0
            # appris_events = 0
            # principalTx_events = 0
            # events_xTx = 0
            # Tx_str = 0 # 0 for iPSC, 1 for APPRIS and 2 for EnsDB
            # Tx_valid = 0
            # Total_Events = SpliceData.shape[0]
            # probable_noise_events = 0
            # probable_noncoding_events = 0
            # utr5_events = 0
            
            # for i in range(SpliceData.shape[0]) :
            #     # Step 0 - get transcripts for each gene (MAJIQ only reports for some)
            #     # Deal with multiple events for the same gene
                
            #     if "Tx_name" in globals() :
            #         del Tx_name
                
            #     if i == 1 :
            #         temp_gene = SpliceData[i,5]
            #         trackj = 1
                
            #     elif temp_gene == SpliceData[i,5] :
            #         trackj = trackj+1
                    
            #     else :
            #         temp_gene = SpliceData[i,5]
            #         trackj = 1
                
            #     # Get gene name and gene_id (from granges filter) using event coordinates
            #     intervals = []
                
            #     for index, row in SpliceData.iterrows():
            #         chrom = row['col1'][3:]
            #         start = row['col2']
            #         end = row['col3']
            #         strand = row['col4']
                    
            #         interval = Interval(chrom, start, end, strand=strand)
            #         intervals.append(interval)
                    
            #     event_bed = BedTool(intervals)
            #     gene_intervals = []
                
            #     for index, row in gene_df.iterrows():
            #         interval = Interval(row['chrom'], row['start'], row['end'], strand=row['strand'], name=row['gene_id'])
            #         gene_intervals.append(interval)
                
            #     gene_bed = BedTool(gene_intervals)
            #     gn = event_bed.intersect(gene_bed, wao=True)
                
            #     # First check if gene_name exactly matches upto gene_version and then get all its Txs
            #     genes_data = gn.to_dataframe()
            #     gn_id = SpliceData[i,6]
                
            #     # Reset Tx flag
            #     Tx_flg = 0
                
            #     if genes_data.shape[1]>0 :
            #         # First try iPSC Tx
            #         if genes_data.shape[1]>1 :
            #             # Some events (coordinates) are mapped to multiple gene_id's, so select one that has same gene_id as majiq and spans
            #             # The genomic range or the one which spans the genomic range
            #             flg_ex = 0
                        
            #             for ti in range(genes_data.shape[1]) :
            #                 if (SpliceData[i,2]>=genes_data.iloc[ti]['start'] and SpliceData[i,3]<=genes_data.iloc[ti]['end'] and len(Tx_list[Tx_list['V2'] == genes_data.iloc[ti]['gene_id']][1])>0) :
            #                     Tx_name = Tx_list[Tx_list['V2'] == genes_data.iloc[ti]['gene_id']][1]
            #                     flg_ex = 1
            #                     break
            #                 if flg_ex==1 :
            #                     break
                            
            #         else :
            #             if (SpliceData[i,2]>=genes_data.iloc['start'] and SpliceData[i,3]<=genes_data.iloc['end'] and len(Tx_list[Tx_list['V2'] == genes_data.iloc['gene_id']][1])>0) :
            #                 Tx_name = Tx_list[Tx_list['V2'] == genes_data.iloc['gene_id']][1]
                            
            #         if ("Tx_name" in globals() and len(Tx_name)>0):                   
            #             filtered_transcripts = [edb.transcript_by_id(tx_id) for tx_id in Tx_name if edb.transcript_by_id(tx_id) is not None]
            #             tl1 = {transcript.transcript_id: len(transcript) for transcript in filtered_transcripts}
            #             # Make sure that we have transcript in EnsDB
                        
            #             if len(tl1)>0 :
            #                 # And event coordinates lies within the Tx
            #                 exon_data = []
                            
            #                 for tx_id in Tx_name :
            #                     transcript = edb.transcript_by_id(tx_id)
                                
            #                     for exon in transcript.exons :
            #                         exon_data.append({
            #                             'transcript_id': transcript.transcript_id,
            #                             'gene_id': transcript.gene_id,
            #                             'chromosome': exon.contig,
            #                             'start': exon.start,
            #                             'end': exon.end,
            #                             'strand': exon.strand
            #                         })
                                    
            #                 dat = pd.DataFrame(exon_data)
                            
            #                 if (SpliceData[i,2]>=min(dat.iloc[:, 1]) and SpliceData[i,3]<=max(dat.iloc[:, 2])) : # Make sure that event lies within the transcript
            #                     iPSC_events = iPSC_events+1
            #                     Tx_str = 'iPSC'
            #                     Tx_flg = 1
                    
            #         if (Tx_flg==0 and len(pd.merge(appr_anno, genes_data, left_on='V2', right_on='gene_id')['V3']))
            
            # ########## FIN TxEnsDB103_layeredV6.R codé en Pyhton ##########




            with open("Summary_stats.txt","a") as fichier :
                fichier.write("************************ BACK FROM TxEnsDB103_layeredV6.R, CONTINUING pgp-a.sh \n")
            
            csv_data = pd.read_csv('all_tx_events.csv')
            csvi = 0
            samples = glob.glob('event_bedfiles/temp_*.bed')
    
            print("csv_data : ", csv_data.head())
            print("samples : ", samples[5])
    
            with open("Summary_stats.txt","a") as fichier :
                fichier.write("GENERATING BED FILES FOR EACH EVENT\n")
    
            for sample in samples :
                # Read csv entry
                csv_ln = csv_data[csvi]
                csvi = csvi + 1
    
                print("csv_ln : ", csv_ln)
                print("csvi : ", csvi)
            
                allexons = sample.split('/')[1].split('_')[1]
                gene_name1 = allexons.split('.')[0]
                gene_name = gene_name1.split('.')[0]
    
                print("gene_name : ", gene_name)

        #     # First sort the bed
        #     sorted_bed = f'event_bedfiles/{allexons}'.sort()
        #     sorted_bed.saveas(f'event_bedfiles/t{allexons}')

        #     # Also read Tx Files to retrieve selected Tx - should find better ways
        #     df = pd.read_csv(f'event_bedfiles/TxID{allexons}', delim_whitespace=True, header=None)
        #     TxID = df.iloc[0, 6]

        #     df = pd.read_csv(sample, delim_whitespace=True, header=None)
        #     strnd = df.iloc[0, 5]

        #     # Get distance to downstream exon (for ties, report first) from current reference and pick start, end and d
        #     a = pybedtools.BedTool(sample)
        #     b = pybedtools.BedTool(allexons)
        #     closest = a.closest(b, s=True, D="a", iu=True, d=True, t="first")
        #     ds = closest.to_dataframe()

        #     # Also get distance to upstream exon from current reference and pick start, end and d
        #     closest = a.closest(b, s=True, D="a", id=True, d=True, t="last")
        #     us = closest.to_dataframe()

        #     # Get up and down stream exon numbers
        #     upexon = us.iloc[:, 10]
        #     dnexon = ds.iloc[:, 10]

        #     # Events star and end
        #     event_st = us.iloc[:, 1]
        #     event_end = us.iloc[:, 2]

        #     diff_exon = upexon - dnexon

        #     # Take absolute value
        #     diff_exon_abs = abs(diff_exon)

        #     if diff_exon_abs>=1 : #ALL EVENTS THAT SPANS 2 OR MORE EXONS
        #         if strnd == '+' :
        #             start = us.iloc[:, 7]
        #             end = ds.iloc[:, 8]
        #         else :
        #             start = ds.iloc[:, 7]
        #             end = us.iloc[:, 8]
                
        #         # Also save
        #         # FORMAT IS: chr, start of up_ex,end of ds_ex, 1, 0, strand, gene_name, TxID
		# 		# First check if event lies between selected exons
        #         if (start<=event_st and end>=event_end) :
        #             input_data = us + [start, end, strnd] + gene_name + TxID

        #             output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
        #             with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.bed", "a") as f:
        #                 f.write("\t".join(output_data) + "\n")

        #             with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.bed", "a") as f:
        #                 f.write("\t".join(output_data) + "\n")
                    
        #             with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.csv", "a") as f:
        #                 f.write(csv_ln)

        #             with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.csv", "a") as f:
        #                 f.write(csv_ln)
                    
        #         else :
        #             if strnd == '+' :
        #                 start = us.iloc[:, 7]
        #                 end = ds.iloc[:, 8] # Both are same
        #                 # First check if star > event_start, then select upstream exon

        #                 if start>=event_st :
        #                     # Get exon (line number in bed file) to read
        #                     exon = ds.iloc[:, 10]
        #                     exon = exon - 2

        #                     with open(f"event_bedfiles/{allexons}", 'r') as file:
        #                         bed_data = [line.strip() for line in file]
                            
        #                     bed_ln = bed_data[exon]

        #                     # Update start
        #                     start = bed_ln.iloc[:, 1]
                        
        #                 # Now check if end <event_end
        #                 if end<=event_end :
        #                     # Get exon (line number in bed file) to read
        #                     exon = us.iloc[:, 10]

        #                     with open(f"event_bedfiles/{allexons}", 'r') as file:
        #                         bed_data = [line.strip() for line in file]
                            
        #                     print('bed_data ',bed_data)
        #                     bed_ln = bed_data[exon]
        #                     print('exon ', exon, ' bed_ln ', bed_ln)

        #                     # Update end
        #                     end = bed_ln.iloc[:, 2]
                        
        #                 # Now one more time check if event lies between selected exons
        #                 if (start<=event_st and end>=event_end):
        #                     input_data = us + [start, end, strnd] + gene_name + TxID

        #                     output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.bed", "a") as f:
        #                         f.write("\t".join(output_data) + "\n")

        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.bed", "a") as f:
        #                         f.write("\t".join(output_data) + "\n")

        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.csv", "a") as f:
        #                         f.write(csv_ln)

        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.csv", "a") as f:
        #                         f.write(csv_ln)
                            
        #                 else :
        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_progress2.txt", "a") as f:
        #                         f.write(f"ds 1 {ds}\n")

        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_progress2.txt", "a") as f:
        #                         f.write(f"us 1 {us}\n")
                            
        #                     print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                        
        #             else : # THIS IS FOR NEGATIVE STRAND
        #                 start = ds.iloc[:, 7]
        #                 end = us.iloc[:, 8]

        #                 # First check if star > event_start, then select upstream exon
        #                 if start>=event_st :
        #                     # Get exon (line number in bed file) to read
        #                     exon = ds.iloc[:, 10]
        #                     exon = int(exon) # Tx on -ve strand has exons listed from bottom to top in increasing order

        #                     with open(f"event_bedfiles/{allexons}", 'r') as file:
        #                         bed_data = [line.strip() for line in file]
                            
        #                     bed_ln = bed_data[exon]

        #                     # Update start
        #                     start = bed_ln.iloc[:, 1]
                        
        #                 # Now check if end <event_end
        #                 if end<=event_end :
        #                     # Get exon (line number in bed file) to read
        #                     exon = us.iloc[:, 10]
        #                     exon = exon - 2 # Tx on -ve strand has exons listed from bottom to top in increasing order

        #                     with open(f"event_bedfiles/{allexons}", 'r') as file:
        #                         bed_data = [line.strip() for line in file]

        #                     bed_ln = bed_data[exon]

        #                     # Update end
        #                     end = bed_ln.iloc[:, 2]

        #                 # Now one more time check if event lies between selected exons
        #                 if (start<=event_st and end>=event_end):
        #                     input_data = us + [start, end, strnd] + gene_name + TxID

        #                     output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.bed", "a") as f:
        #                         f.write("\t".join(output_data) + "\n")

        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.bed", "a") as f:
        #                         f.write("\t".join(output_data) + "\n")

        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.csv", "a") as f:
        #                         f.write(csv_ln)

        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.csv", "a") as f:
        #                         f.write(csv_ln)
                            
        #                 else :
        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_progress2.txt", "a") as f:
        #                         f.write(f"ds 2 {ds}\n")

        #                     with open(f"temp_all_events_sashimi/{splicing_events_file}_progress2.txt", "a") as f:
        #                         f.write(f"us 2 {us}\n")
                            
        #                     print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                        
        #     elif diff_exon_abs==0 :
        #         if strnd == '+' :
        #             start = us.iloc[:, 7]
        #             end = ds.iloc[:, 8] # Both are same

		# 			# First check if star > event_start, then select upstream exon
        #             if start>=event_st :
        #                 # Get exon (line number in bed file) to read
        #                 exon = us.iloc[:, 10]
        #                 exon = exon - 2

        #                 with open(f"event_bedfiles/{allexons}", 'r') as file:
        #                     bed_data = [line.strip() for line in file]

        #                 bed_ln = bed_data[exon]

        #                 # Update start
        #                 start = bed_ln.iloc[:, 1]

        #                 # Now go on the other side
        #                 if start>=event_st : # Get the other exon
        #                     # Get exon (line number in bed file) to read
        #                     exon = us.iloc[:, 10]

        #                     with open(f"event_bedfiles/{allexons}", 'r') as file:
        #                         bed_data = [line.strip() for line in file]

        #                     bed_ln = bed_data[exon]

        #                     # Update start
        #                     start = bed_ln.iloc[:, 1]

        #             # Now check if end <event_end
        #             if end<=event_end :
        #                 # Get exon (line number in bed file) to read
        #                 exon = us.iloc[:, 10]

        #                 with open(f"event_bedfiles/{allexons}", 'r') as file:
        #                     bed_data = [line.strip() for line in file]

        #                 print(f"bed_data {bed_data}")
        #                 bed_ln = bed_data[exon]
        #                 print(f"exon {exon} bed_ln {bed_ln}")

        #                 # Update end
        #                 end = bed_ln.iloc[:, 2]

        #                 if end<=event_end :
        #                     # Get exon (line number in bed file) to read
        #                     exon = us.iloc[:, 10]
        #                     exon = exon - 2 # Reading line for readarray starts from 0

        #                     with open(f"event_bedfiles/{allexons}", 'r') as file:
        #                         bed_data = [line.strip() for line in file]

        #                     bed_ln = bed_data[exon]

        #                     # Update end
        #                     end = bed_ln.iloc[:, 2]

        #             # Now one more time check if event lies between selected exons
        #             if (start<=event_st and end>=event_end) :
        #                 input_data = us + [start, end, strnd] + gene_name + TxID

        #                 output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi01.bed", "a") as f:
        #                     f.write("\t".join(output_data) + "\n")

        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.bed", "a") as f:
        #                     f.write("\t".join(output_data) + "\n")

        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi01.csv", "a") as f:
        #                     f.write(csv_ln)

        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.csv", "a") as f:
        #                     f.write(csv_ln)
                        
        #             else :
        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_progress01.txt", "a") as f:
        #                     f.write(f"ds {ds}\n")

        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_progress01.txt", "a") as f:
        #                     f.write(f"us {us}\n")

        #                 print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")

        #         else : # THIS IS FOR NEGATIVE STRAND
        #             start = ds.iloc[:, 7]
        #             end = us.iloc[:, 8]

        #             # First check if star > event_start, then select upstream exon
        #             if start>=event_st :
        #                 # Get exon (line number in bed file) to read
        #                 exon = us.iloc[:, 10]

        #                 with open(f"event_bedfiles/{allexons}", 'r') as file:
        #                     bed_data = [line.strip() for line in file]

        #                 bed_ln = bed_data[exon]

        #                 # Update start
        #                 start = bed_ln.iloc[:, 1]
                    
        #             # Now check if end <event_end
        #             if end<=event_end :
        #                 # Get exon (line number in bed file) to read
        #                 exon =us.iloc[:, 10]
        #                 exon = exon - 2 #Tx on -ve strand has exons listed from bottom to top in increasing order
                        
        #                 with open(f"event_bedfiles/{allexons}", 'r') as file:
        #                     bed_data = [line.strip() for line in file]

        #                 bed_ln = bed_data[exon]

        #                 # Update end
        #                 end = bed_ln.iloc[:, 2]

        #             # Now one more time check if event lies between selected exons
        #             if (start<=event_st and end>=event_end) :
        #                 input_data = us + [start, end, strnd] + gene_name + TxID

        #                 output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi02.bed", "a") as f:
        #                     f.write("\t".join(output_data) + "\n")

        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.bed", "a") as f:
        #                     f.write("\t".join(output_data) + "\n")

        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi02.csv", "a") as f:
        #                     f.write(csv_ln)

        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.csv", "a") as f:
        #                     f.write(csv_ln)
                        
        #             else :
        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_progress02.txt", "a") as f:
        #                     f.write(f"ds {ds}\n")

        #                 with open(f"temp_all_events_sashimi/{splicing_events_file}_progress02.txt", "a") as f:
        #                     f.write(f"us {us}\n")
                        
        #                 print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                
        #         # First check if star > event_start, then select upstream exon

        #     else :
        #         with open(f"temp_all_events_sashimi/{splicing_events_file}_progress_all.txt", "a") as f:
        #             f.write(f"ds {ds}\n")

        #         with open(f"temp_all_events_sashimi/{splicing_events_file}_progress_all.txt", "a") as f:
        #             f.write(f"us {us}\n")

        # print("################################ DONE GENERATING BED AND OTHER RELATED FILES FOR ALL EVENTS TO CONTINUE FOR SASHIMI PLOTD")

        # with open("Summary_stats.txt", "a") as f:
        #     f.write("################################ DONE GENERATING BED AND OTHER RELATED FILES FOR ALL EVENTS TO CONTINUE FOR SASHIMI PLOTD\n")
        
        # if arg1==1 or arg1==2 or arg1==5 :
        #     print("################################ STARTED CREATING SASHIMI PLOTS FOR ALL EVENTS - MAY TAKE MANY HOURS ")

        #     with open("Summary_stats.txt", "a") as f:
        #         f.write("################################ STARTED CREATING SASHIMI PLOTS FOR ALL EVENTS - MAY TAKE MANY HOURS \n")

        #     # CALL run_sashimiV1.sh for SASHIMI PLOTS




            ########## run_sashimiV1.sh codé en Pyhton ##########

            # THIS SI FINAL SCRIPT FOR SASHIMI PLOTS FOR ALL PARTS OF PGP

            # THIS SCRIPT CONTAINS SASHIMI PLOT CODE FOR
            # 0. PLEASE NOTE THAT FOLLOWING 2 FILES (OR SOFT LINKS) SHOULD BE IN CURRENT FOLDER
                #01: ggsashimi_txV4.py
                #02: Homo_sapiens.GRCh38.103.chr.sorted_new.gtf
            # 1. Skiptic Events
            # 2. ALL MAJIQ EVENTS
            # 3. CE (INCLUDING INCLUSION, EXTENSION AND IR) events

            args = [f"temp_all_events_sashimi/{splicing_events_file}", f"temp_all_events_sashimi/{splicing_events_file}", 2]

            # First CHECK IF Called from pgp-a/b or pgp-c
            # CHECK IF 3 ARGUMENTS ARE PROVIDED

            if len(args)==3 :

                ######### NEW - Now read input csv and bed files and flag
                inp_csv = args[0]
                inp_bed = args[1]

                # Get folder
                inp_prefix = args[0].split('/')[0]

                if args[2]==1 :
                    os.makedirs(f"{inp_prefix}/sashimi_plots",exist_ok=True)
                    if len(os.listdir(f"{inp_prefix}/sashimi_plots/"))!=0 :
                        for f in glob.glob(f"{inp_prefix}/sashimi_plots/*.*") :
                            os.remove(f)
                    
                    bed = f"{inp_bed}"

                    with open(bed, 'r') as file:
                        all_bed_data = [line.strip() for line in file]
                        
                    with open(bed, 'r') as file:
                        nrecrds = sum(1 for line in file)

                    nrecrdst = nrecrds/2
                    print(f"read {nrecrdst} records")
                    csv = f"{inp_csv}"

                    with open(csv, 'r') as file:
                        all_csv_data = [line.strip() for line in file]
                    
                    i = 0
                    eventn = 0
                    
                    while i<nrecrds :
                        # Construct string for ggsashimi
                        line1 = all_bed_data[i]
                        line2 = all_bed_data[i+1]
                        i = i + 2

                        # Read strand
                        strnd = line1.iloc[:, 5]
                        # Also read TxID
                        TxID = line1.iloc[:, 7]

                        # us exon length
                        exon1 = 0
                        # ds exon length
                        exon2 = 0

                        combined_line = f"{line1} {line2}".split()
                        field1 = combined_line[0]
                        field2 = int(combined_line[1]) - 50
                        field11 = int(combined_line[10]) + 50
                        print(f"{field1}:{field2}-{field11}")

                        event = all_csv_data[eventn]
                        gene_name = event.split(',')[7]

                        fields = event.split(',')
                        field8 = fields[7]
                        field2 = fields[1]
                        field3 = fields[2]
                        field4 = fields[3]
                        fn = f"{field8}-{field2}_{field3}-{field4}"

                        # String for majiq event
                        chr_name = event.split(',')[1]
                        start = event.split(',')[2]
                        end = event.split(',')[3]

                        # For now using strand from ggsashimi
                        comb_line = f"{chr_name} {start} {end}".split()
                        majiq_event = f"{comb_line[0]}-{comb_line[1]}-{comb_line[2]}"

                        # Also get actual event identified
                        event_identified = f"PGPEvent-{combined_line[2]}-{combined_line[9]}"
                        eventn = eventn + 1

                        print("processing event num {eventn} and event {fn}") # And majiq event is $majiq_event and event identified is $event_identified

                        # Here removed -PGPTx flag
                        command = [
                            "./ggsashimi_txV3.py",
                            "-A", "median_j",
                            "-b", "all_bams.tsv",
                            "-c", line,
                            "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
                            "-GeneName", gene_name,
                            "-MajiqStrnd", strnd,
                            "-ORIG", "1",
                            "-UEX", exon1,
                            "-DEX", exon2,
                            "-MajiqTx", majiq_event,
                            "-Majiq", fn,
                            "-Tx", TxID,
                            "-M", "1",
                            "-C", "3",
                            "-o", f"{inp_prefix}/sashimi_plots/{fn}",
                            "-O", "3",
                            "--alpha", "0.25",
                            "--base-size=20",
                            "--ann-height=2.5",
                            "--height=2.5",
                            "--width=18",
                            "-P", "palette.txt"
                        ]

                        subprocess.run(command)

                    # Now merge all pdf's
                    command = [
                        "python", "merge_sashimis.py", f"{inp_prefix}/sashimi_plots/"
                    ]

                    subprocess.run(command)

                ###### THIS IS FOR ALL MAJIQ EVENTS
                if args[2]==2 :
                    os.makedirs("all_events_sashimi",exist_ok=True)
                    if len(os.listdir("all_events_sashimi/"))!=0 :
                        for f in glob.glob("all_events_sashimi/*.*") :
                            os.remove(f)
                
                    bed = f"{args[1]}_all_sashimi.bed"

                    with open(bed, 'r') as file:
                        all_bed_data = [line.strip() for line in file]
                    
                    with open(bed, 'r') as file:
                        nrecrds = sum(1 for line in file)

                    print(f"read {nrecrds} records")
                    csv = f"{args[1]}_all_sashimi.csv"

                    with open(csv, 'r') as file:
                        all_csv_data = [line.strip() for line in file]
                    
                    i = 0
                    eventn = 0

                    while i<nrecrds :
                        # Construct string for ggsashimi
                        line1 = all_bed_data[i]
                        i = i + 1

                        # Read strand
                        strnd = '+'

                        # Also read TxID
                        TxID = line1.iloc[7]
                        line = f"{line1.iloc[0]}:{line1.iloc[1]-50}-{line1.iloc[2]+50}" # This is the actual event

                        # Get majiq event
                        event = all_csv_data[eventn]
                        eventn = eventn + 1

                        # THIS IS TO KEEP DIFFERENT FILE NAMES FOR EVENTS WITH IDENTICAL COORDINATES AND GENE_NAMES
                        gene_name = event.split(',')[4]

                        if i==1 :
                            temp_gene = event.split(',')[4]
                            trackj=1

                        elif temp_gene==gene_name :
                            trackj = trackj + 1

                        else :
                            temp_gene = event.split(',')[4]
                            trackj=1
                        
                        # SHOULD NOT HAVE THIS BUT NEED TO KEEP AS PER REQUEST
                        fields = event.split(',')
                        fn = f"{fields[4]}-{fields[0]}-{fields[1]}-{fields[2]}-{trackj}"

                        # String for majiq event
                        chr_name = event.split(',')[0]
                        start = event.split(',')[1]
                        end = event.split(',')[2]
                        
                        comb_line = f"{chr_name} {start} {end}".split()
                        majiq_event = f"{comb_line[0]}-{comb_line[1]}-{comb_line[2]}"

                        exon1 = 0
                        exon2 = 0

                        # Also get actual event identified
                        event_identified = f"None -{line1.iloc[:,1]}-{line1.iloc[:,2]}"

                        print(f"processing event num {eventn} and event {fn}")

                        command = [
                            "./ggsashimi_txV3.py",
                            "-A", "median_j",
                            "-b", "all_bams.tsv",
                            "-c", line,
                            "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
                            "-GeneName", gene_name,
                            "-MajiqStrnd", strnd,
                            "-ORIG", "1",
                            "-UEX", exon1,
                            "-DEX", exon2,
                            "-MajiqTx", majiq_event,
                            "-Majiq", fn,
                            "-Tx", TxID,
                            "-M", "1",
                            "-C", "3",
                            "-o", f"all_events_sashimi/{fn}",
                            "-O", "3",
                            "--alpha", "0.25",
                            "--base-size=20",
                            "--ann-height=2.5",
                            "--height=2.5",
                            "--width=18",
                            "-P", "palette.txt"
                        ]

                        subprocess.run(command)

                    # Now merge all pdf's
                    command = ["python", "merge_sashimis.py", "all_events_sashimi/"]
                    subprocess.run(command)
                
                # THIS SECTION IS FOR CE_INCLUSION EVENTS
                # WILL MERGE INCLUSION AND EXTENSION EVENTS
                if args[2]==3 :
                    os.makedirs(f"{inp_prefix}/ce_incl_sashimi_plots",exist_ok=True)
                    if len(os.listdir(f"{inp_prefix}/ce_incl_sashimi_plots/"))!=0 :
                        for f in glob.glob(f"{inp_prefix}/ce_incl_sashimi_plots/*.*") :
                            os.remove(f)
                
                    bed = inp_bed

                    with open(bed, 'r') as file:
                        all_bed_data = [line.strip() for line in file]
                    
                    with open(bed, 'r') as file:
                        nrecrds = sum(1 for line in file)
                    
                    nrecrdst = nrecrds/3

                    print(f"read {nrecrdst} records")
                    csv = inp_csv

                    with open(csv, 'r') as file:
                        all_csv_data = [line.strip() for line in file]
                    
                    i = 0
                    eventn = 0

                    while i<nrecrds :
                        # Construct string for ggsashimi
                        line1 = all_bed_data[i]
                        line2 = all_bed_data[i+1]
                        line3 = all_bed_data[i+2]

                        i = i + 3

                        # Read strand
                        strnd = line1.iloc[:, 5]

                        # Also read TxID
                        TxID = line1.iloc[7]

                        if strnd=='+' :
                            strndflg = "plus"

                            # us exon length
                            exon1 = line1.iloc[:, 3]

                            # ds exon length
                            exon2 = line3.iloc[:, 3]

                            # Now modify to reflect whole up/dn exons
                            up = line1.iloc[:, 1]
                            upn = up - exon1
				             
                            # ds exon
                            dn = line3.iloc[:, 2]
                            dnn = dn + exon2

                            l123 = f"{line1} {line2} {line3}".split()
                            fulltitle = f"{l123[0]}-{l123[1]}:{l123[2]}-{l123[9]}:{l123[10]}-{l123[17]}:{l123[18]}"

                            comb = f"{line1} {upn} {dnn}".split()
                            line = f"{comb[0]}:{comb[8]}-{comb[9]}"

                        else :
                            strndflg = "minus"

                            # us exon length
                            exon1 = line3.iloc[:, 3]

                            # ds exon length
                            exon2 = line1.iloc[:, 3]

                            # Now modify to reflect whole up/dn exons
                            dn = line3.iloc[:, 1]
                            dnn = dn - exon1
				             
                            # ds exon
                            up = line1.iloc[:, 2]
                            upn = up + exon2

                            l123 = f"{line1} {line2} {line3}".split()
                            fulltitle = f"{l123[0]}-{l123[17]}:{l123[18]}-{l123[9]}:{l123[10]}-{l123[1]}:{l123[2]}"

                            comb = f"{line1} {dnn} {upn}".split()
                            line = f"{comb[0]}:{comb[8]}-{comb[9]}"

                        event = all_csv_data[eventn]

                        fields = event.split(',')
                        fn = f"{fields[7]}-{fields[1]}_{fields[2]}-{fields[3]}"

                        # String for majiq event
                        chr_name = event.split(',')[1]
                        start = event.split(',')[2]
                        end = event.split(',')[3]
                        gene_name = event.split(',')[7]
                        intron = end - start

                        #gene_name-start-intron-end #for now using strand from ggsashimi
                        comb_line = f"{chr_name} {start} {end}".split()
                        majiq_event = f"{comb_line[0]}-{comb_line[1]}-{comb_line[2]}"

                        # Also get actual event identified
                        event_identified = f"{line2.iloc[:,0]}-{line2.iloc[:,1]}-{line2.iloc[:,2]}"
                        eventn = eventn + 1

                        print(f"processing event num {eventn} and event {fn}")

                        command = [
                            "./ggsashimi_txV3.py",
                            "-A", "median_j",
                            "-b", "all_bams.tsv",
                            "-c", line,
                            "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
                            "-GeneName", gene_name,
                            "-MajiqStrnd", strnd,
                            "-ORIG", "1",
                            "-UEX", exon1,
                            "-DEX", exon2,
                            "-FullTitle", fulltitle,
                            "-MajiqTx", majiq_event,
                            "-Majiq", fn,
                            "-Tx", TxID,
                            "-M", "1",
                            "-C", "3",
                            "-o", f"{inp_prefix}/ce_incl_sashimi_plots/{fn}",
                            "-O", "3",
                            "--alpha", "0.25",
                            "--base-size=20",
                            "--ann-height=2.5",
                            "--height=2.5",
                            "--width=18",
                            "-P", "palette.txt"
                        ]

                        subprocess.run(command)

                    # Now merge all pdf's
                    command = ["python", "merge_sashimis.py", f"{inp_prefix}/ce_incl_sashimi_plots/"]
                    subprocess.run(command)

                # CE_EXTENSION
                if args[2]==4 :
                    os.makedirs(f"{inp_prefix}/ce_ext_sashimi_plots",exist_ok=True)
                    if len(os.listdir(f"{inp_prefix}/ce_ext_sashimi_plots/"))!=0 :
                        for f in glob.glob(f"{inp_prefix}/ce_ext_sashimi_plots/*.*") :
                            os.remove(f)
                
                    bed = inp_bed

                    with open(bed, 'r') as file:
                        all_bed_data = [line.strip() for line in file]
                    
                    with open(bed, 'r') as file:
                        nrecrds = sum(1 for line in file)
                    
                    nrecrdst = nrecrds/3

                    print(f"read {nrecrdst} records")
                    csv = inp_csv

                    with open(csv, 'r') as file:
                        all_csv_data = [line.strip() for line in file]
                    
                    i = 0
                    eventn = 0

                    while i<nrecrds :
                        # Construct string for ggsashimi
                        line1 = all_bed_data[i]
                        line2 = all_bed_data[i+1]
                        line3 = all_bed_data[i+2]

                        i = i + 3

                        # Read strand
                        strnd = line1.iloc[:, 5]

                        # Also read TxID
                        TxID = line1.iloc[7]

                        if strnd=='+' :
                            strndflg = "plus"

                            # us exon length
                            exon1 = line1.iloc[:, 3]

                            # ds exon length
                            exon2 = line3.iloc[:, 3]

                            # Now modify to reflect whole up/dn exons
                            up = line1.iloc[:, 1]
                            upn = up - exon1
				             
                            # ds exon
                            dn = line3.iloc[:, 2]
                            dnn = dn + exon2

                            l123 = f"{line1} {line2} {line3}".split()
                            fulltitle = f"{l123[0]}-{l123[1]}:{l123[2]}-{l123[9]}:{l123[10]}-{l123[17]}:{l123[18]}"

                            comb = f"{line1} {upn} {dnn}".split()
                            line = f"{comb[0]}:{comb[8]}-{comb[9]}"

                        else :
                            # us exon length
                            exon1 = line3.iloc[:, 3]

                            # ds exon length
                            exon2 = line1.iloc[:, 3]

                            # Now modify to reflect whole up/dn exons
                            dn = line3.iloc[:, 1]
                            dnn = dn - exon1
				             
                            # ds exon
                            up = line1.iloc[:, 2]
                            upn = up + exon2

                            l123 = f"{line1} {line2} {line3}".split()
                            fulltitle = f"{l123[0]}-{l123[17]}:{l123[18]}-{l123[9]}:{l123[10]}-{l123[1]}:{l123[2]}"

                            comb = f"{line1} {dnn} {upn}".split()
                            line = f"{comb[0]}:{comb[8]}-{comb[9]}"

                        event = all_csv_data[eventn]

                        fields = event.split(',')
                        fn = f"{fields[7]}-{fields[1]}_{fields[2]}-{fields[3]}"

                        # String for majiq event
                        chr_name = event.split(',')[1]
                        start = event.split(',')[2]
                        end = event.split(',')[3]
                        gene_name = event.split(',')[7]
                        intron = end - start

                        #gene_name-start-intron-end #for now using strand from ggsashimi
                        comb_line = f"{chr_name} {start} {end}".split()
                        majiq_event = f"{comb_line[0]}-{comb_line[1]}-{comb_line[2]}"

                        # Also get actual event identified
                        event_identified = f"{line2.iloc[:,0]}-{line2.iloc[:,1]}-{line2.iloc[:,2]}"
                        eventn = eventn + 1

                        print(f"processing event num {eventn} and event {fn}")

                        command = [
                            "./ggsashimi_txV3.py",
                            "-A", "median_j",
                            "-b", "all_bams.tsv",
                            "-c", line,
                            "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
                            "-GeneName", gene_name,
                            "-MajiqStrnd", strnd,
                            "-ORIG", "1",
                            "-UEX", exon1,
                            "-DEX", exon2,
                            "-FullTitle", fulltitle,
                            "-MajiqTx", majiq_event,
                            "-Majiq", fn,
                            "-Tx", TxID,
                            "-M", "1",
                            "-C", "3",
                            "-o", f"{inp_prefix}/ce_ext_sashimi_plots/{fn}",
                            "-O", "3",
                            "--alpha", "0.25",
                            "--base-size=20",
                            "--ann-height=2.5",
                            "--height=2.5",
                            "--width=18",
                            "-P", "palette.txt"
                        ]

                        subprocess.run(command)

                    # Now merge all pdf's
                    command = ["python", "merge_sashimis.py", f"{inp_prefix}/ce_ext_sashimi_plots/"]
                    subprocess.run(command)
                
                ######THIS IS FOR ALL IR EVENTS
                if args[2]==5 :
                    os.makedirs(f"{inp_prefix}/ir_sashimi_plots",exist_ok=True)
                    if len(os.listdir(f"{inp_prefix}/ir_sashimi_plots/"))!=0 :
                        for f in glob.glob(f"{inp_prefix}/ir_sashimi_plots/*.*") :
                            os.remove(f)
                
                    bed = inp_bed

                    with open(bed, 'r') as file:
                        all_bed_data = [line.strip() for line in file]
                    
                    with open(bed, 'r') as file:
                        nrecrds = sum(1 for line in file)

                    print(f"read {nrecrds} records")
                    csv = inp_csv

                    with open(csv, 'r') as file:
                        all_csv_data = [line.strip() for line in file]
                    
                    i = 0
                    eventn = 0

                    while i<nrecrds :
                        # Construct string for ggsashimi
                        line1 = all_bed_data[i]

                        i = i + 1

                        # Read strand
                        strnd = line1.iloc[:, 5]

                        # Also read TxID
                        TxID = line1.iloc[7]

                        line = f"{line1[0]}:{line1[1]-50}-{line1[2]+50}"

                        # Get majiq event
                        event = all_csv_data[eventn]
                        eventn = eventn + 1

                        # THIS IS TO KEEP DIFFERENT FILE NAMES FOR EVENTS WITH IDENTICAL COORDINATES AND GENE_NAMES
                        gene_name = event.split(',')[4]

                        if i==1 :
                            temp_gene = event.split(',')[4]
                            trackj = 1
                        
                        elif temp_gene == gene_name :
                            trackj = trackj + 1
                        
                        else :
                            temp_gene = event.split(',')[4]
                            trackj = 1
                        
                        # SHOULD NOT HAVE THIS BUT NEED TO KEEP AS PER REQUEST
                        fields = event.split(',')
                        fn = f"{fields[4]}-{fields[0]}_{fields[1]}-{fields[2]}-{trackj}" #also gene_id to avoid same file names for repeated events

                        # String for majiq event
                        chr_name = event.split(',')[0]
                        start = event.split(',')[1]
                        end = event.split(',')[2]

                        comb_line = f"{chr_name} {start} {end}".split()
                        majiq_event = f"{comb_line[0]}-{comb_line[1]}-{comb_line[2]}"
                        exon1 = 0
                        exon2 = 0

                        # Also get actual event identified
                        event_identified = f"{line1.iloc[:,0]}-{line1.iloc[:,1]}-{line1.iloc[:,2]}"

                        print(f"processing event num {eventn} and event {fn}")

                        command = [
                            "./ggsashimi_txV3.py",
                            "-A", "median_j",
                            "-b", "all_bams.tsv",
                            "-c", line,
                            "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
                            "-GeneName", gene_name,
                            "-MajiqStrnd", strnd,
                            "-ORIG", "1",
                            "-UEX", exon1,
                            "-DEX", exon2,
                            "-MajiqTx", majiq_event,
                            "-Majiq", fn,
                            "-Tx", TxID,
                            "-M", "1",
                            "-C", "3",
                            "-o", f"{inp_prefix}/ir_sashimi_plots/{fn}",
                            "-O", "3",
                            "--alpha", "0.25",
                            "--base-size=20",
                            "--ann-height=2.5",
                            "--height=2.5",
                            "--width=18",
                            "-P", "palette.txt"
                        ]

                        subprocess.run(command)

                    # Now merge all pdf's
                    command = ["python", "merge_sashimis.py", f"{inp_prefix}/ir_sashimi_plots/"]
                    subprocess.run(command)
                    
            else :
                print("came in for peaksbackmapping")
                
                ###### THIS PART IS FOR BACKMAPPING
                flg_sashimi_files_only = 1
                
                if flg_sashimi_files_only==1 :
                    step01_flg = 1
                    
                    if step01_flg==1 :
                        ################################ GENERATING BED FILE FOR SASHIMI PLOTS FOR ALL EVENTS
                        # IMPORTANT - THIS CODE RELIES ON THE OUTPUT OF BEDTOOLS CLOSEST FUNCTION
                        # GIVEN A MAJIQ EVENT (as a bed file) and BED FILE FOR THE TRANSCRIPT IT BELIEVES TO BE PART OF
                        # BEDTOOLS CLOSEST FUNCTION USES INPUT RANGE (AMJIQ JUNTION) and finds closest exons (and their ranges) from the TRANSCRIPT BED FILE
                        # BEDTOOLS CLOSEST RETURNS A BED FILE WITH FOLLOWING OUTPUT
                        # INPUT: chr,start,end,1,0,+ (majiq event with strand)
                        # OUTPUT: chr,start,end,1,0,+,chr,start,end,size,exon_num,+,distance
                        # Here first 6 entries are the original majiq input and next 7 entries are the resulting closest exon coordinates, its size (in bp),exon_num,strand and distance from reference
                        # SOULD ADD CHECKS ON EXON_NUM READING FROM FILE
                        
                        inp = args[0].split('.')[0]
                        print(f"run_sashimiV1.sh, got inp as {inp}")
                        
                        # Also get folder
                        folder = inp.split('/')[0]
                        print(f"run_sashimiV1.sh, got folder as {folder}")
                        
                        # Also get event TYPE
                        eventyp = args[1].split('/')[1].split('.')[0]
                        print(f"run_sashimiV1.sh, got eventype as {eventyp}")
                        
                        if os.path.exists(f"{inp}_majiq.bed"):
                            os.remove(f"{inp}_majiq.bed")
                        
                        if os.path.exists(f"{inp}_majiq.csv"):
                            os.remove(f"{inp}_majiq.csv")
                            
                        if os.path.exists("all_tx_events.csv"):
                            os.remove("all_tx_events.csv")
                            
                        if os.path.exists("majiq_events_sashimi2.bed"):
                            os.remove("majiq_events_sashimi2.bed")
                            
                        if os.path.exists("majiq_events_sashimi2.csv"):
                            os.remove("majiq_events_sashimi2.csv")
                            
                        if os.path.exists("majiq_events_sashimi01.bed"):
                            os.remove("majiq_events_sashimi01.bed")
                            
                        if os.path.exists("majiq_events_sashimi01.csv"):
                            os.remove("majiq_events_sashimi01.csv")
                            
                        if os.path.exists("majiq_events_sashimi02.bed"):
                            os.remove("majiq_events_sashimi02.bed")
                            
                        if os.path.exists("majiq_events_sashimi02.csv"):
                            os.remove("majiq_events_sashimi02.csv")
                            
                        if os.path.exists("majiq_events_progress2.txt"):
                            os.remove("majiq_events_progress2.txt")
                            
                        if os.path.exists("majiq_events_progress01.txt"):
                            os.remove("majiq_events_progress01.txt")
                            
                        if os.path.exists("majiq_events_progress02.txt"):
                            os.remove("majiq_events_progress02.txt")
                            
                        if os.path.exists("majiq_events_progress_all.txt"):
                            os.remove("majiq_events_progress_all.txt")
                            
                        if os.path.exists("majiq_events_progress1.txt"):
                            os.remove("majiq_events_progress1.txt")
                            
                        flg = 1
                        
                        if flg==1 :
                            os.makedirs("event_bedfiles",exist_ok=True)
                            if len(os.listdir("event_bedfiles/"))!=0 :
                                for f in glob.glob("event_bedfiles/*.*") :
                                    os.remove(f)
                            
                            command = [
                                "Rscript",
                                "txens.R",
                                args[0],
                                "principal_txs.csv"
                            ]

                            subprocess.run(command, capture_output=True, text=True)
                            
                        with open("all_tx_events.csv", 'r') as file:
                            csv_data = [line.strip() for line in file]
                        
                        csvi = 0
                        samples = glob.glob('event_bedfiles/temp_*.bed')
                        
                        for sample in samples :
                            # Read csv entry
                            csv_ln = csv_data[csvi]
                            csvi = csvi + 1
                            
                            print(f"processing {sample}")
                            allexons = sample.split('/')[1].split('_')[1]
                            gene_name1 = allexons.split('.')[0]
                            gene_name = gene_name1.split('-')[0]
                            
                            # First sort the bed
                            bed = pybedtools.BedTool(f"event_bedfiles/{allexons}")
                            sorted_bed = bed.sort()
                            sorted_bed.saveas(f"event_bedfiles/t{allexons}")
                            
                            # Also read Tx Files to retrieve selected Tx - should find better ways
                            with open(f"event_bedfiles/TxID{allexons}", 'r') as file:
                                first_line = file.readline().strip()
                                TxID = first_line.split()[6]
                                
                            samp = pd.read_csv(sample, sep='\t', header=None)
                            strnd = samp.iloc[:, 5]
                            
                            # Get distance to downstream exon (for ties, report first) from current reference and pick start, end and d
                            a = pybedtools.BedTool(sample)
                            b = pybedtools.BedTool(f'event_bedfiles/t{allexons}')
                            closest = a.closest(b, d=True, s=True, t='first')
                            ds = [line.strip() for line in closest]
                            ds = closest.to_dataframe()
                            
                            # Also get distance to upstream exon from current reference and pick start, end and d
                            closest = a.closest(b, d=True, s=True, t='last')
                            us = [line.strip() for line in closest]
                            us = closest.to_dataframe()
                            
                            # Get up and down stream exon numbers
                            upexon = us.iloc[:,10]
                            dnexon = ds.iloc[:,10]
                            
                            # Events star and end
                            event_st = us.iloc[:,1]
                            event_end = us.iloc[:,2]
                            diff_exon = upexon - dnexon
                            
                            # Take absolute value
                            diff_exon_abs = abs(diff_exon)
                            
                            if diff_exon_abs>=1 : # ALL EVENTS THAT SPANS 2 OR MORE EXONS
                                if strnd=="+" :
                                    start = us.iloc[:,7]
                                    end = ds.iloc[:,8]
                                    
                                else :
                                    start = ds.iloc[:,7]
                                    end = us.iloc[:,8]
                                
                                # Also save
                                # FORMAT IS: chr, start of up_ex,end of ds_ex, 1, 0, strand, gene_name, TxID
                                # First check if event lies between selected exons
                                if (start<=event_st and end>=event_end) :
                                    input_data = us + [start, end, strnd] + gene_name + TxID
                                    output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
                                    with open(f"{inp}_majiq.bed", "a") as f:
                                        f.write("\t".join(output_data) + "\n")
                                        
                                    with open(f"{inp}_majiq.csv", "a") as f:
                                        f.write(csv_ln)
                                    
                                else :
                                    if strnd=="+" :
                                        start = us.iloc[:,7]
                                        end = ds.iloc[:,8] # Both are same
                                        
                                        # First check if star > event_start, then select upstream exon
                                        if start>=event_st :
                                            # Get exon (line number in bed file) to read
                                            exon = ds.iloc[:,10]
                                            exon = exon - 2
                                            
                                            with open(f"event_bedfiles/{allexons}", 'r') as file:
                                                bed_data = [line.strip() for line in file]
                                            
                                            bed_ln = bed_data[exon]
                                            
                                            # Update start
                                            start = bed_ln.iloc[:,1]
                                            
                                        # Now check if end <event_end
                                        if end<=event_end :
                                            # Get exon (line number in bed file) to read
                                            exon = us.iloc[:,10]
                                            
                                            with open(f"event_bedfiles/{allexons}", 'r') as file:
                                                bed_data = [line.strip() for line in file]
                                            
                                            print(f"bed_data {bad_data}")
                                            bed_ln = bed_data[exon]
                                            print(f"exon {exon} bed_ln {bed_ln}")
                                            
                                            # Update end
                                            end = bed_ln.iloc[:,2]
                                            
                                        # Now one more time check if event lies between selected exons
                                        if (start<=event_st and end>=event_end):
                                            input_data = us + [start, end, strnd] + gene_name + TxID
                                            output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
                                            with open(f"{inp}_majiq.bed", "a") as f:
                                                f.write("\t".join(output_data) + "\n")
                                                
                                            with open(f"{inp}_majiq.csv", "a") as f:
                                                f.write(csv_ln)
                                            
                                        else :
                                            with open("majiq_events_progress2.txt", "a") as f:
                                                f.write(f"ds 1 {ds}\n")
                                                f.write(f"us 1 {us}\n")
                                                
                                            print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                                                
                                        else : # THIS IS FOR NEGATIVE STRAND
                                            start = ds.iloc[:,7]
                                            end = us.iloc[:,8]
                                            
                                            # First check if star > event_start, then select upstream exon
                                            if start>=event_st :
                                                # Get exon (line number in bed file) to read
                                                exon = ds.iloc[:,10]
                                                