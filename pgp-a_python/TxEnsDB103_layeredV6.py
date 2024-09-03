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
import os
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import rpy2.rinterface_lib.callbacks
from rpy2.robjects.packages import importr

pandas2ri.activate()
rpy2.rinterface_lib.callbacks.logger.setLevel('INFO')

# Also load gtf file fron V86
# Get object of EnsDBV99

base = importr('base')
GenomicFeatures = importr('GenomicFeatures')
Repitools = importr('Repitools')
GenomicRanges = importr('GenomicRanges')
IRanges = importr('IRanges')

ro.r('''
library(AnnotationHub)
library(ensembldb)

ah <- AnnotationHub()
edb <- query(ah, c("EnsDb", "Hsapiens", "103"))[[1]]
''')

pd.set_option('display.max_columns', None)

#Also get Tx lengths to remove 5'utr
tx_lens = ro.r['transcriptLengths'](ro.r['edb'], with_utr5_len=True, with_utr3_len=True)

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

print(f"Started Generating BED files for Splicing Events in folder event_bedfiles/ from File: {args[0]}")

trackj = 1
temp_gene = ""
current_gene = ""
tx_lengths = []

df_notfound = pd.DataFrame({
    'seqnames': pd.Series(dtype='str'),
    'start': pd.Series(dtype='int'),
    'end': pd.Series(dtype='int'),
    'strand': pd.Series(dtype='str'),
    'genename': pd.Series(dtype='str'),
    'junc_type': pd.Series(dtype='str')
})

df_zeroutr = pd.DataFrame({
    'seqnames': pd.Series(dtype='str'),
    'start': pd.Series(dtype='int'),
    'end': pd.Series(dtype='int'),
    'strand': pd.Series(dtype='str'),
    'genename': pd.Series(dtype='str'),
    'junc_type': pd.Series(dtype='str')
})

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

def create_granges_object(seqnames, ranges_start, ranges_end, strand):
    granges_func = ro.r['GRanges']
    iranges_func = ro.r['IRanges']
    
    iranges = iranges_func(ranges_start, ranges_end)
    granges = granges_func(seqnames, ranges=iranges, strand=strand)
    return granges

def create_granges_filter(granges):
    granges_filter_func = ro.r['GRangesFilter']
    return granges_filter_func(granges, type='any')

def get_genes(edb, granges_filter):
    genes_func = ro.r['genes']
    return genes_func(edb, filter=granges_filter)

def create_tx_id_filter(tx_name):
    TxIdFilter = ro.r['TxIdFilter']
    return TxIdFilter(tx_name)

def create_gene_name_filter(gene_name):
    GeneNameFilter = ro.r['GeneNameFilter']
    return GeneNameFilter(gene_name)

def get_exons_by_tx(edb, tx_filter):
    exonsBy = ro.r['exonsBy']
    return exonsBy(edb, by="tx", filter=tx_filter)

def get_unlisted_data(exons):
    return exons.do_slot("unlistData")

def convert_to_dataframe(r_object):
    data_frame = ro.r['data.frame']
    return data_frame(r_object)

for i in range(SpliceData.shape[0]) :
    # Step 0 - get transcripts for each gene (MAJIQ only reports for some)
    # Deal with multiple events for the same gene
    
    if "Tx_name" in globals() :
        del Tx_name
    
    if i == 1 :
        temp_gene = SpliceData.iloc[i,4]
        trackj = 1
    
    elif temp_gene == SpliceData.iloc[i,4] :
        trackj = trackj+1
        
    else :
        temp_gene = SpliceData.iloc[i,4]
        trackj = 1
    
    seqnames = ro.r('substr')(SpliceData.iloc[i, 0], 4, ro.r('nchar')(ro.r('as.character')(SpliceData.iloc[i, 0])))
    ranges_start = SpliceData.iloc[i, 1]
    ranges_end = SpliceData.iloc[i, 2]
    strand = SpliceData.iloc[i, 3]
    
    granges = create_granges_object(seqnames, ranges_start, ranges_end, strand)
    granges_filter = create_granges_filter(granges)
    gn = get_genes(ro.r['edb'], granges_filter)
    
    # First check if gene_name exactly matches upto gene_version and then get all its Txs
    genes_data = ro.r('annoGR2DF')(gn)
    gn_id = SpliceData.iloc[i,5].split('.')[0]
    
    # Reset Tx flag
    Tx_flg = 0
    
    if len(genes_data)>0 :
        # First try iPSC Tx
        if len(genes_data)>1 :
            # Some events (coordinates) are mapped to multiple gene_id's, so select one that has same gene_id as majiq and spans
            # The genomic range or the one which spans the genomic range
            flg_ex = 0
            
            for ti in range(len(genes_data)) :
                if (SpliceData.iloc[i,1]>=genes_data.rx(ti+1,'start')[0] and SpliceData.iloc[i,2]<=genes_data.rx(ti+1,'end')[0] and len(Tx_list.loc[Tx_list[1].str.contains(genes_data.rx(ti+1,'gene_id')[0]), 0].values)>0) :
                    Tx_name = Tx_list.loc[Tx_list[1].str.contains(genes_data.rx(ti+1,'gene_id')[0]), 0].values
                    flg_ex = 1
                    break
                if flg_ex==1 :
                    break
                
        else :
            if (SpliceData.iloc[i,1]>=genes_data.rx['start'] and SpliceData.iloc[i,2]<=genes_data.rx['end'] and len(Tx_list.loc[Tx_list[1].str.contains(genes_data.rx['gene_id']), 0].values)>0) :
                Tx_name = Tx_list.loc[Tx_list['V2'].str.contains(genes_data.iloc['gene_id']), 'V1'].values
                
        if ("Tx_name" in globals() and len(Tx_name)>0):
            tx_filter = create_tx_id_filter(Tx_name)
            tl1 = ro.r['transcriptLengths'](ro.r['edb'], filter = tx_filter)
            
            # Make sure that we have transcript in EnsDB
            if len(tl1)>0 :
                # And event coordinates lies within the Tx
                exons = get_exons_by_tx(ro.r['edb'], tx_filter)
                unlisted_data = get_unlisted_data(exons)
                dat = convert_to_dataframe(unlisted_data)
                
                if (SpliceData.iloc[i,1]>=min(dat[1]) and SpliceData.iloc[i,2]<=max(dat[2])) : # Make sure that event lies within the transcript
                    iPSC_events = iPSC_events + 1
                    Tx_str = 'iPSC'
                    Tx_flg = 1 ##### Teste jusqu'ici et tout fonctionne #####
        
        genes_data = pandas2ri.rpy2py(genes_data)
        
        if (Tx_flg==0 and len(pd.merge(appr_anno, genes_data, left_on='V2', right_on='gene_id')['V3'])>0) :
            Tx_name = appr_anno.loc[genes_data.rx["gene_id"].isin(appr_anno["V2"]), "V3"].tolist()

            tx_filter = create_tx_id_filter(Tx_name)
            tltt = ro.r['transcriptLengths'](ro.r['edb'], filter = tx_filter)

            # As APPRIS has multiple P1 Txs for some genes (like RALYL), make sure that the one with maximum exons and highest Tx length (in bp) is selected
            if len(tltt)>0 :
                # If multiple P1 Txs, then select the one with highet num of exons and largest size (in bp) and encapsulates event
                tflg = 0

                for ti in range(len(tltt)) :
                    tx_filter = create_tx_id_filter(tltt[ti,'tx_id'])
                    exons = get_exons_by_tx(ro.r['edb'], tx_filter)
                    unlisted_data = get_unlisted_data(exons)
                    dat = convert_to_dataframe(unlisted_data)
                    
                    if (SpliceData.iloc[i,1]>=min(dat[1]) and SpliceData.iloc[i,2]<=max(dat[2])) : # Make sure that event lies within the transcript
                        if tflg==0 :
                            tlt = tltt[ti,:]
                        
                        else :
                            tlt = pd.concat(tlt, tltt[ti].to_frame().T, ignore_index=True)
                        
                        tflg = 1
                
                ########
                # Now make sure that longest Tx is selected
                if tflg==1 :
                    if len(tlt)>0 :
                        jj1 = tlt.index[tlt['nexon'] == tlt['nexon'].max()].tolist() # For nexon
                        tl1t = tlt[jj1,:]
                        jj = tl1t.index[tl1t['tx_len'] == tl1t['tx_len'].max()].tolist() # Which(tl$nexon==max(tl$nexon)&&tl$tx_len==max(tl$tx_len)) # For nexon
                        tl1 = tl1t[jj,:]
                        Tx_name = tl1[1,2]
                        
                        if len(tl1)>0 :
                            ########
                            # And event coordinates lies within the Tx
                            tx_filter = create_tx_id_filter(Tx_name)
                            exons = get_exons_by_tx(ro.r['edb'], tx_filter)
                            unlisted_data = get_unlisted_data(exons)
                            dat = convert_to_dataframe(unlisted_data)
                            
                            if (SpliceData.iloc[i,1]>=min(dat[1]) and SpliceData.iloc[i,2]<=max(dat[2])) : # Make sure that event lies within the transcript
                                appris_events = appris_events + 1
                                Tx_str = 'APPRIS'
                                Tx_flg = 1
        
        gene_name_filter = create_gene_name_filter(gn['gene_name'])
        transcript_lengths = ro.r['transcriptLengths'](ro.r['edb'], filter=gene_name_filter)
        
        if (Tx_flg==0) and (len(transcript_lengths)>0) : #and Finally longest Tx
            if len(transcript_lengths)>1 :
                events_xTx = events_xTx + 1
            
            tltt = ro.r['transcriptLengths'](ro.r['edb'], filter = gene_name_filter)
            
            # First get list of Tx's which encapsulates the event
            tflg = 0
            
            for ti in len(tltt) :
                # Make sure that transcipt starting with ENST is used - was case with gene FH, where one of the Tx_id is LRG_504t1
                if 'ENST' in tltt.loc[ti, 'tx_id'] :
                    tx_filter = create_tx_id_filter(tltt[ti,'tx_id'])
                    exons = get_exons_by_tx(ro.r['edb'], tx_filter)
                    unlisted_data = get_unlisted_data(exons)
                    dat = convert_to_dataframe(unlisted_data)
                    
                    if (SpliceData.iloc[i,1]>=min(dat[1]) and SpliceData.iloc[i,2]<=max(dat[2])) : # Make sure that event lies within the transcript
                        if tflg==0 :
                            tlt = tltt.iloc[ti,:]
                        
                        else :
                            tlt = pd.concat(tlt, tltt.loc[ti].to_frame().T, ignore_index=True)
                        
                        tflg = 1
            
            # Now make sure that longest Tx is selected
            if len(tlt)>0 :
                jj1 = tlt.index[tlt['nexon'] == tlt['nexon'].max()].tolist() # For nexon
                tl1t = tlt.iloc[jj1,:]
                jj = tl1t.index[tl1t['tx_len'] == tl1t['tx_len'].max()].tolist() # Which(tl$nexon==max(tl$nexon)&&tl$tx_len==max(tl$tx_len)) # For nexon
                tl1 = tl1t.iloc[jj,:]
                Tx_name = tl1.iloc[1,2]
                
                if len(tl1)>0 :
                    # And event coordinates lies within the Tx
                    tx_filter = create_tx_id_filter(Tx_name)
                    exons = get_exons_by_tx(ro.r['edb'], tx_filter)
                    unlisted_data = get_unlisted_data(exons)
                    dat = convert_to_dataframe(unlisted_data)
                    
                    if (SpliceData.iloc[i,1]>=min(dat[1]) and SpliceData.iloc[i,2]<=max(dat[2])) : # Make sure that event lies within the transcript
                        principalTx_events = principalTx_events + 1
                        Tx_str = 'EnsDB'
                        Tx_flg = 1
                        
        ################3
        # Step 1 - for each transcript, get all exons for it
        ddf = pd.DataFrame({
            'seqnames': pd.Series(dtype='str'),
            'start': pd.Series(dtype='float'),
            'end': pd.Series(dtype='float'),
            'width': pd.Series(dtype='float'),
            'exon_rank': pd.Series(dtype='float'),
            'strand': pd.Series(dtype='str')
        })
        
        # Includes Tx id to use with sashimi plots
        ddf1 = pd.DataFrame({
            'seqnames': pd.Series(dtype='str'),
            'start': pd.Series(dtype='float'),
            'end': pd.Series(dtype='float'),
            'width': pd.Series(dtype='float'),
            'exon_rank': pd.Series(dtype='float'),
            'strand': pd.Series(dtype='str'),
            'TxID': pd.Series(dtype='str')
        })
        
        # Check if Tx is available from either of the three resources
        tx_filter = create_tx_id_filter(Tx_name)
        transcript_lengths = ro.r['transcriptLengths'](ro.r['edb'], filter=tx_filter)
        
        if (Tx_flg==1 and len(tl1)>0 and len(Tx_name)>0 and len(transcript_lengths)>0) :
            # Now use Tx to
            EnsGenes1 = get_exons_by_tx(ro.r['edb'], tx_filter)
            EnsGenes_cds = get_exons_by_tx(ro.r['edb'], tx_filter, columns=["tx_seq_start", "tx_cds_seq_start", "tx_seq_end", "tx_cds_seq_end", "tx_biotype"])
            
            # Get 5utr length
            utr5l1 = tx_lens.iloc['Tx_name', 'utr5_len']
            
            if len(utr5l1)>1 :
                utr5l = utr5l1[0]
            
            else :
                urt5l = utr5l1
            
            utr3l1 = tx_lens.iloc['Tx_name', 'utr3_len']
            
            if len(utr3l1)>1 :
                utr3l = utr3l1[0]
            
            # Get exons
            # Now get data frame and subtract 5utr and 3utr
            unlisted_data1 = get_unlisted_data(EnsGenes1)
            g_dat = convert_to_dataframe(unlisted_data1)
            
            unlisted_data_cds = get_unlisted_data(EnsGenes_cds)
            g_dat1 = convert_to_dataframe(unlisted_data_cds)
            
            if g_dat1[0,'tx_biotype']=='protein_coding' :
                if g_dat[0, 'strand']=='+' :
                    # First check if event lies within the cds
                    if SpliceData.iloc[i,'V2']>=g_dat1[0,'tx_cds_seq_start'] :
                        if g_dat[0,'width']>utr5l : #AARS1 has trl5l > exon 1 width
                            g_dat[0,'start'] = g_dat[0,'start'] + utr5l - 1
                            g_dat[0,'width'] = g_dat[0,'width'] - utr5l
                            
                    else :
                        probable_noncoding_events = probable_noncoding_events + 1
                        print(f"*** gene: {SpliceData.iloc[i,'V5']} @ line: {i} lies outside of 5 UTR")
                        
                        with open(args[2], "a") as f :
                            f.write(f"Probable non-cds events#: {probable_noncoding_events} ** gene: {SpliceData.iloc[i, 'V5']} @ line: {i} lies outside of 5 UTR\n")
                            
            if g_dat1[0, 'tx_biotype']=='protein_coding' :
                probable_noise_events = probable_noise_events + 1
                print(f"Noise Event#: {probable_noise_events}--- Transcript: {Tx_name} gene: {SpliceData.iloc[i,'V5']} @ line: {i} is: {g_dat1.iloc[0,'tx_biotype']}")
                
                with open(args[2], "a") as f :
                    f.write(f"Noise Event#: {probable_noise_events}--- Transcript: {Tx_name} gene: {SpliceData.iloc[i,'V5']} @ line: {i} is: {g_dat1.iloc[0,'tx_biotype']}\n")
            
            df = g_dat
            df['seqnames'] == f"chr{df.iloc['seqnames']}"
            
            # Include junction width as column 4 and exon rank as column 5
            cols = list(range(0, 4)) + [7, 4]
            ddf = pd.concat([ddf, df[:, cols]], ignore_index=True)
            cols1 = list(range(0, 4)) + [7, 4, 6]
            ddf1 = pd.concat([ddf1, df[:, cols1]], ignore_index=True)
            
            if len(ddf)>0 :
                # Save bedfile for current junction
                SJdat = pd.DataFrame({
                    'seqnames': [SpliceData.iloc[i,0]],
                    'start': [SpliceData.iloc[i,1]],
                    'end': [SpliceData.iloc[i,2]],
                    'name': [1],
                    'score': [0],
                    'strand': [SpliceData.iloc[i,3]]
                })
                
                SJdat1 = pd.DataFrame({
                    'seqnames': [SpliceData.iloc[i,0]],
                    'start': [SpliceData.iloc[i,1]],
                    'end': [SpliceData.iloc[i,2]],
                    'name': [1],
                    'score': [0],
                    'strand': [SpliceData.iloc[i,3]],
                    'strand': [SpliceData.iloc[i,4]],
                    'strand': [Tx_name]
                })
                
                if (SpliceData.iloc[i,1]>=min(ddf.iloc[:,1])) and (SpliceData.iloc[i,2]<=max(ddf.iloc[:,2])) : # Make sure that event lies within the transcript
                    Tx_valid = Tx_valid + 1
                    SJdat2 = pd.DataFrame({
                        'seqnames': [SpliceData.iloc[i,0]],
                        'start': [SpliceData.iloc[i,1]],
                        'end': [SpliceData.iloc[i,2]],
                        'gene_name': [SpliceData.iloc[i,4]],
                        'gene_id': [SpliceData.iloc[i,5]],
                        'Tx': [Tx_name],
                        'Tx_start': [min(ddf.iloc[:,1])],
                        'Tx_end': [max(ddf.iloc[:,2])],
                        'Tx_type': [Tx_str],
                        'valid': [1]
                    })
                
                else :
                    SJdat2 = pd.DataFrame({
                        'seqnames': [SpliceData.iloc[i,0]],
                        'start': [SpliceData.iloc[i,1]],
                        'end': [SpliceData.iloc[i,2]],
                        'gene_name': [SpliceData.iloc[i,4]],
                        'gene_id': [SpliceData.iloc[i,5]],
                        'Tx': [Tx_name],
                        'Tx_start': [min(ddf.iloc[:,1])],
                        'Tx_end': [max(ddf.iloc[:,2])],
                        'Tx_type': [Tx_str],
                        'valid': [0]
                    })
                
                SJdat.to_csv(f"event_bedfiles/temp_{SpliceData.iloc[i,4]}-{trackj}.bed", sep='\t', index=False, header=False, quoting=3)
                           
                # Also save all_events bed file with Tx selected for sashimi plots
                SJdat1.to_csv("all_events_bed_sashimi.tab", sep='\t', index=False, header=False, mode='a', quoting=3)
                ddf.to_csv(f"event_bedfiles/{SpliceData.iloc[i,4]}-{trackj}.bed", sep='\t', index=False, header=False, quoting=3)
                
                # TxID is used with sashimi plots
                ddf1.to_csv(f"event_bedfiles/TxID{SpliceData.iloc[i,4]}-{trackj}.bed", sep='\t', index=False, header=False, quoting=3)
                
                # Save for ce script
                SpliceData.iloc[i,:].to_csv("all_tx_events.csv", sep=',', index=False, header=False, mode='a', quoting=3)
                
                # Save for Tx mapping
                SJdat2.to_csv("events_to_tx_mapping_valid.csv", sep=',', index=False, header=False, mode='a', quoting=3)
                
            else :
                print(f"got null at record: {i} for gene: {SpliceData.iloc[i,4]}")
                df_notfound = pd.concat([df_notfound, SpliceData.iloc[i,:]], ignore_index=True) # Concatenate all
                
        else :
            print(f"got null at record: {i} for gene: {SpliceData.iloc[i,4]}")
            df_notfound = pd.concat([df_notfound, SpliceData.iloc[i,:]], ignore_index=True) # Concatenate all
        
    else :
        print(f"got null at record: {i} for gene: {SpliceData.iloc[i,4]}")
        df_notfound = pd.concat([df_notfound, SpliceData.iloc[i,:]], ignore_index=True) # Concatenate all
        
# Finally write it
if df_notfound.shape[0]>0 :
    df_notfound.to_csv("EnsDB_tx_not_found.csv", sep=',', index=False, header=False, quoting=3)
    
print(f"Out of a total of: {Total_Events} Events: Transcripts for: {df_notfound.shape[0]} Events are not_found")
print(f"Total : {iPSC_events} iPSC Princiapl Tx are selected")
print(f"Total : {appris_events} APPRIS Principal Tx are selected")
print(f"Total : {principalTx_events} EnsDB Princiapl Tx (maximum number of Exons and largest size (in bp)) are selected")
print(f"Out of Total : {Total_Events} A total of: {Tx_valid} Events have Valid Transcripts (event lies between Transcript ends) and {Total_Events-Tx_valid} Events have Transcripts that does not Encapsulate events")

with open(args[2], "a") as f :
    f.write(f"Finally total Noise Events detected are: {probable_noise_events}\n")
    f.write(f"Finally total Probable non-cds events#: {probable_noncoding_events}\n")
