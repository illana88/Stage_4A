########## TxEnsDB103_layeredV6.R code en Pyhton ##########
# -*- coding: utf-8 -*-

# Using ENSG (commented out this - junction start and end range) to get transcripts
# Starting from previous version (TxEnsDB103_layeredV1.R), it has 41 events for which Tx selected does not encapsulate the event
# An example is event: chr12	22515151	22517985	AC053513.1
# So coordinate based search is back on to see if it makes difference
# THIS SCRIPT HAS BEEN MODIFIED - SO REPLACE IT IN ALL VERSIONS
# 04/20/2022 - removing 5'utr

for name in list(globals().keys()):
    if not name.startswith('_'):
        del globals()[name]

# USE EnsDb V86
# USE EnsDB V99 from Annotation hub
        
import sys
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import gffutils
import os
import subprocess
import sqlite3

# Also load gtf file fron V86
# Get object of EnsDBV99

conn = sqlite3.connect("transcript_lengths.db")
cursor = conn.cursor()
cursor.execute("SELECT * FROM transcripts")
edb = cursor.fetchall()
conn.close()

# Read bed file from MAJIQ
args = sys.argv[1:]

GeneIDField = 6

# Read Peaks File
SpliceData = pd.read_csv(args[0], header=None)

# Also read Tx list
Tx_list = pd.read_csv(args[1], header=None)

# Also read appris annoation data to get principal 1 isoform
appr_anno1 = pd.read_csv("GRCh38_appris_data.principal.txt", sep="\t", header=None)

# V1 - hugo_symbol, V2 - ENSG_ID, V3 - TX_ID, V4 - , V2 - PRINCIPAL/ALTERNATE
appr_anno1.columns = ['V1', 'V2', 'V3', 'V4', 'V5']

# Select only PRINCIPAL.1 ISOFORMS
appr_anno = appr_anno1[appr_anno1['V5'].str.contains("PRINCIPAL:1")] #APPRIS has multiple P1 Tx for a gene- e.g RALYL
print("appr_anno : ", appr_anno)

if os.access("all_tx_events.csv",os.F_OK):
    # Delete file if it exists
    os.remove(("all_tx_events.csv"))
    
if os.access("all_events_bed_sashimi.tab",os.F_OK):
    # Delete file if it exists
    os.remove(("all_events_bed_sashimi.tab"))
    
if os.access("events_to_tx_mapping_valid.csv",os.F_OK):
    # Delete file if it exists
    os.remove(("events_to_tx_mapping_valid.csv"))
    
if os.access("events_to_tx_mapping_invalid.csv",os.F_OK):
    # Delete file if it exists
    os.remove(("events_to_tx_mapping_invalid.csv"))

with open(args[2],"a") as fichier :
    fichier.write('                                 \n')
    fichier.write(f"Starting From TxEnsDB103_layeredV6.R: --------------- Processing file: {args[0]} with: {SpliceData.shape[0]} events to generate each event .bed files in event_bedfiles/ folder: \n")

print("Started Generating BED files for Splicing Events in folder event_bedfiles/ from File: {args[0]}")

trackj = 1
temp_gene = ""
current_gene = ""
tx_lengths = []

df_notfound = pd.DataFrame({
    'seqnames': pd.Series(dtype='str'),
    'start': pd.Series(dtype='numeric'),
    'end': pd.Series(dtype='numeric'),
    'strand': pd.Series(dtype='str'),
    'genename': pd.Series(dtype='str'),
    'junc_type': pd.Series(dtype='str')
})

print("df_notfound : ", df_notfound)

df_zeroutr = pd.DataFrame({
    'seqnames': pd.Series(dtype='str'),
    'start': pd.Series(dtype='numeric'),
    'end': pd.Series(dtype='numeric'),
    'strand': pd.Series(dtype='str'),
    'genename': pd.Series(dtype='str'),
    'junc_type': pd.Series(dtype='str')
})

print("df_zeroutr : ", df_zeroutr)

repeated_entries = 0
iPSC_events = 0
appris_events = 0
principalTx_events = 0
events_xTx = 0
Tx_str = 0 # 0 for iPSC, 1 for APPRIS and 2 for EnsDB
Tx_valid = 0
Total_Events = SpliceData.shape[0]
probable_noise_events = 0
probable_noncoding_events = 0
utr5_events = 0

# # Get gene name and gene_id (from granges filter) using event coordinates
# def filter_genes_for_row(i):
#     chromosome = SpliceData.iloc[i, 0]
#     start = SpliceData.iloc[i, 1]
#     end = SpliceData.iloc[i, 2]
#     strand = SpliceData.iloc[i, 3]

#     gr_range = gr.GRanges(
#         seqnames=[chromosome],
#         ranges=[(start, end)],
#         strand=[strand]
#     )

#     genes = list(edb.features_of_type('gene'))
#     genes_data = {
#         'seqnames': [gene.chrom for gene in genes],
#         'ranges': [(gene.start, gene.end) for gene in genes],
#         'strand': [gene.strand for gene in genes]
#     }

#     gr_genes = gr.GRanges(
#         seqnames=genes_data['seqnames'],
#         ranges=genes_data['ranges'],
#         strand=genes_data['strand']
#     )

#     filtered_genes = gr_range.intersect(gr_genes, which="any")

#     return filtered_genes


# def parallel_filter_transcripts():
#     with ThreadPoolExecutor() as executor:
#         futures = [executor.submit(filter_genes_for_row, i) for i in range(len(SpliceData))]
#         results = [future.result() for future in futures]
    
#     return results

# gn = parallel_filter_transcripts()

# for i in range(SpliceData.shape[0]) :
#     # Step 0 - get transcripts for each gene (MAJIQ only reports for some)
#     # Deal with multiple events for the same gene
    
#     if "Tx_name" in globals() :
#         del Tx_name
    
#     if i == 1 :
#         temp_gene = SpliceData.iloc[i,5]
#         trackj = 1
    
#     elif temp_gene == SpliceData.iloc[i,5] :
#         trackj = trackj+1
        
#     else :
#         temp_gene = SpliceData.iloc[i,5]
#         trackj = 1

#     # First check if gene_name exactly matches upto gene_version and then get all its Txs
#     genes_data = gn.to_dataframe()
#     gn_id = SpliceData.iloc[i,6].split('.')[0]
    
#     # Reset Tx flag
#     Tx_flg = 0
    
#     if genes_data.shape[0]>0 :
#         # First try iPSC Tx
#         if genes_data.shape[0]>1 :
#             # Some events (coordinates) are mapped to multiple gene_id's, so select one that has same gene_id as majiq and spans
#             # The genomic range or the one which spans the genomic range
#             flg_ex = 0
            
#             for ti in range(genes_data.shape[0]) :
#                 if (SpliceData.iloc[i,2]>=genes_data.iloc[ti]['start'] and SpliceData.iloc[i,3]<=genes_data.iloc[ti]['end'] and len(Tx_list.loc[Tx_list['V2'].str.contains(genes_data.iloc[ti]['gene_id']), 'V1'].values)>0) :
#                     Tx_name = Tx_list.loc[Tx_list['V2'].str.contains(genes_data.iloc[ti]['gene_id']), 'V1'].values
#                     flg_ex = 1
#                     break
#                 if flg_ex==1 :
#                     break
                
#         else :
#             if (SpliceData.iloc[i,2]>=genes_data.iloc['start'] and SpliceData.iloc[i,3]<=genes_data.iloc['end'] and len(Tx_list.loc[Tx_list['V2'].str.contains(genes_data.iloc['gene_id']), 'V1'].values)>0) :
#                 Tx_name = Tx_list.loc[Tx_list['V2'].str.contains(genes_data.iloc['gene_id']), 'V1'].values
                
#         if ("Tx_name" in globals() and len(Tx_name)>0):                   
            # filtered_transcripts = [edb.transcript_by_id(tx_id) for tx_id in Tx_name if edb.transcript_by_id(tx_id) is not None]
            # tl1 = {transcript.id: len(transcript.sequence) for transcript in filtered_transcripts if transcript.sequence}
#             # Make sure that we have transcript in EnsDB
            
#             if len(tl1)>0 :
#                 # And event coordinates lies within the Tx
                # exon_data = []
                
                # for tx_id in Tx_name :
                #     transcript = edb.transcript_by_id(tx_id)
                    
                #     for exon in transcript.exons :
                #         exon_data.append({
                #             'transcript_id': transcript.transcript_id,
                #             'gene_id': transcript.gene_id,
                #             'chromosome': exon.contig,
                #             'start': exon.start,
                #             'end': exon.end,
                #             'strand': exon.strand
                #         })
                        
                # dat = pd.DataFrame(exon_data)
                
#                 if (SpliceData[i,2]>=min(dat.iloc[:, 1]) and SpliceData[i,3]<=max(dat.iloc[:, 2])) : # Make sure that event lies within the transcript
#                     iPSC_events = iPSC_events+1
#                     Tx_str = 'iPSC'
#                     Tx_flg = 1
        
        # if (Tx_flg==0 and len(pd.merge(appr_anno, genes_data, left_on='V2', right_on='gene_id')['V3'])>0) :
        #     Tx_name = appr_anno.loc[genes_data["gene_id"].isin(appr_anno["V2"]), "V3"].tolist()

        #     filtered_transcripts = [edb.transcript_by_id(tx_id) for tx_id in Tx_name if edb.transcript_by_id(tx_id) is not None]
        #     tltt = {transcript.id: len(transcript.sequence) for transcript in filtered_transcripts if transcript.sequence}

        #     # As APPRIS has multiple P1 Txs for some genes (like RALYL), make sure that the one with maximum exons and highest Tx length (in bp) is selected


########## FIN TxEnsDB103_layeredV6.R code en Pyhton ##########