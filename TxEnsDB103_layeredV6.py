########## TxEnsDB103_layeredV6.R codé en Pyhton ##########
import subprocess

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