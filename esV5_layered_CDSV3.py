import os
import glob
import sys
import subprocess
import shutil
import pybedtools
import io
import pandas as pd
import csv
from pybedtools import BedTool
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# THIS VERSION NOW ADDS EVENT_ID FLAG TO FINAL AA FILE FOR FUSED EVENTS TO MAP BACK SASHIMI PLOTS FOR THE PEAKS EVENTS

# STARTED USING TxEnsDB103_layeredV2.R which takes care of the Txs seleceted not encapsulating the event - 01-13-2022
# Adding up/dn stream exon length for sashimi plots
# Above change INTERFERE with removal of repeated entries, so CREATING NEW bed file (skiptics_uniq_sashimi.bed) that conatin exon lengths
# STARTING NEW CODE TO TIGHTLY ONLY EXON_SKIP EVENTS - Dec2, 2021
# ONLY EXON_SKIP EVENTS - USE

args = sys.argv[1:]

bed_flg = 1

os.makedirs("event_bedfiles",exist_ok=True)

if bed_flg==1 :
    if len(os.listdir("event_bedfiles/"))!=0 :
        for f in glob.glob("event_bedfiles/*.*") :
            os.remove(f)
    
if os.path.exists("events_to_tx_mapping_valid.csv"):
    os.remove("events_to_tx_mapping_valid.csv")

if os.path.exists("events_tx_mapping_invalid.csv"):
    os.remove("events_tx_mapping_invalid.csv")

print("CLEARING UP res_skiptics folder from previous run if any. SASHIMI FOLDERS WILL ONLY BE DELETED BEFORE CREATING NEW SASHIMI PLOTS.")

with open("res_skiptics/Summary_stats.txt", "a") as f :
    f.write("CLEARING UP res_skiptics folder from previous run if any. SASHIMI FOLDERS WILL ONLY BE DELETED BEFORE CREATING NEW SASHIMI PLOTS.\n")

if len(os.listdir("res_skiptics/"))!=0 :
    for f in glob.glob("res_skiptics/*.*") :
        os.remove(f)

inpfile = args[0] # RBP_all.csv
txfile = args[1]

with open("res_skiptics/Summary_stats.txt", "a") as f :
    f.write("######################################\n")
    f.write("CALLING TxEnsDB103_layeredV6.R FROM SKIPTICS SCRITP TO GENERATE BED FILES FOR EACH EVENT\n")

if bed_flg==1 :
    command = [
        "python",
        "TxEnsDB103_layeredV6.py",
        f"{inpfile}",
        f"{txfile}",
        "res_skiptics/Summary.Skiptics.txt"
    ]
    
    result = subprocess.run(command, capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)
    
    if os.path.exists("EnsDB_tx_not_found.csv") :
        shutil.copy("EnsDB_tx_not_found.csv", "res_skiptics/")
    
with open("res_skiptics/Summary_stats.txt", "a") as f :
    f.write("BACK FROM TxEnsDB103_layeredV6.R script\n")
    f.write("#######################################\n")

# Also read csv file to filter out non_es events for ce script
with open("all_tx_events.csv", "r") as file:
    csv_data = file.readlines()
    csv_data = [line.strip() for line in csv_data]

csvi = 0
flg_bed = 1

if bed_flg==1 :
    samples = glob.glob('event_bedfiles/temp_*.bed')
    
    for sample in samples :
        # Read csv entry
        csv_ln = csv_data[csvi]
        gene_id = csv_ln.split(',')[5]
        csvi = csvi + 1
        
        # Get  all exons bed
        allexons = sample.split('/')[1]
        allexons = allexons.split('_')[1]
        
        # Also get 5UTR and 3UTR bed files
        gene_name1 = allexons.split('.')[0]
        gene_name = gene_name1.split('-')[0]
        
        # First sort the bed
        bed = pybedtools.BedTool(f"event_bedfiles/{allexons}")
        sorted_bed = bed.sort()
        sorted_bed.saveas(f"event_bedfiles/t{allexons}")
        
        # Also read Tx Files to retrieve selected Tx - should find better ways
        with open(f"event_bedfiles/TxID{allexons}", 'r') as file:
            first_line = file.readline().strip()
        
        columns = first_line.split()
        TxID = columns[6]
        
        with open(sample, 'r') as file:
            strnd = [line.split()[5] for line in file if line.strip()]
        
        # Get distance to downstream exon (for ties, report first) from current reference and pick start, end and d
        a = pybedtools.BedTool(sample)
        b = pybedtools.BedTool(f'event_bedfiles/t{allexons}')
        
        closest = a.closest(b, s=True, D="a", iu=True, d=True, t="first")
        ds = closest.to_dataframe(names = [
            "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
            "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
        ])
    
        # Also get distance to upstream exon from current reference and pick start, end and d
        closest = a.closest(b, s=True, D="a", id=True, d=True, t="last")
        us = closest.to_dataframe(names = [
            "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
            "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
        ])
        
        # Get up/dn exon lengths
        upexonl = us.iloc[:,9].tolist()[0]
        dnexonl = ds.iloc[:,9].tolist()[0]
    
        # Get up and down stream exon numbers
        upexon = us.iloc[:,10].tolist()[0]
        dnexon = ds.iloc[:,10].tolist()[0]
        
        # Get overlap with up/dn exon - 0 means complete overlap which is assumed as exon skip event
        dsovlp = ds.iloc[:,12].tolist()[0]
        usovlp = us.iloc[:,12].tolist()[0]
        diff_exon = int(upexon) - int(dnexon)
    
        # Take absolute value
        diff_exon_abs = abs(diff_exon)
        
        if (upexon==dnexon) and (upexonl==dnexonl) :
            # THIS SECTION IS FOR EXON_SKIP events
            sample_bed = pybedtools.BedTool(sample)
            allexons_bed = pybedtools.BedTool(f"event_bedfiles/t{allexons}")
            dsi = sample_bed.intersect(allexons_bed, s=True, wo=True)
            dsi = dsi.to_dataframe()
            
            usi = sample_bed.intersect(allexons_bed, s=True, wo=True)
            usi = usi.to_dataframe()
            
            dsovlp = dsi.iloc[:, 12].tolist()[0]
            usovlp = usi.iloc[:, 12].tolist()[0]
            
            # Get up/dn exon lengths
            upexonl = dsi.iloc[:, 9].tolist()[0]
            dnexonl = dsi.iloc[:, 9].tolist()[0]
            
            diff1 = int(dnexonl) - int(dsovlp)
            diff2 = int(upexonl) - int(usovlp) # If single exon is intersected
            
            if (dsovlp==usovlp) and (upexonl==dnexonl) and (diff1==1) and (diff2==1) :
                # As algorithm returns skipped exon for exon_skip events, so determine coordinates of actual exons involved in the event
                us_df = pd.read_csv(io.StringIO(us), delim_whitespace=True, header=None)
                ds_df = pd.read_csv(io.StringIO(ds), delim_whitespace=True, header=None)

                formatted_udf = pd.DataFrame({
                    0: us_df[0],
                    1: us_df[1] - 10,
                    2: us_df[1],
                    3: us_df[3],
                    4: us_df[4],
                    5: us_df[5]
                })
                
                formatted_ddf = pd.DataFrame({
                    0: ds_df[0],
                    1: ds_df[1] - 10,
                    2: ds_df[1],
                    3: ds_df[3],
                    4: ds_df[4],
                    5: ds_df[5]
                })
                
                formatted_udf.to_csv("dump_us.bed", sep="\t", header=False, index=False)
                formatted_ddf.to_csv("dump_ds.bed", sep="\t", header=False, index=False)
                
                a = pybedtools.BedTool("dump_ds.bed")
                b = pybedtools.BedTool(f'event_bedfiles/t{allexons}')
                
                closest = a.closest(b, s=True, D="a", iu=True, d=True, t="last")
                dsn = closest.to_dataframe(names = [
                    "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
                    "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
                ])
            
                a = pybedtools.BedTool("dump_us.bed")
                
                closest = a.closest(b, s=True, D="a", id=True, d=True, t="first")
                usn = closest.to_dataframe(names = [
                    "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
                    "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
                ])
                
                # Get up and down stream exon numbers
                upexon = usn.iloc[:,10].tolist()[0]
                dnexon = dsn.iloc[:,10].tolist()[0]
                
                dsovlp = dsn.iloc[:,12].tolist()[0]
                usovlp = usn.iloc[:,12].tolist()[0]
                
                if (int(upexon)==1) or (int(dnexon)==1) or (int(dsovlp)!=0) or (int(usovlp!=0)) :
                    df = pd.read_csv(sample, sep='\t', header=None)
                    formatted_df = pd.DataFrame({
                        0: [df[0].astype(str) + ':' + df[1].astype(str) + '-' + df[2].astype(str)],
                        1: df[0],
                        2: df[1],
                        3: df[2],
                        4: df[3],
                        5: df[4],
                        6: df[5],
                        7: [gene_name]
                    })
                    formatted_df.to_csv("res_skiptics/IGV_all_others.csv", sep=',', header=False, index=False)
                    
                    # Also save csv file
                    with open("all_non_skiptics.csv", "a", newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow(csv_ln)
                    
                else :
                    u = usn.split('\t')
                    d = dsn.split('\t')
                    
                    formatted_df = pd.DataFrame({
                        0: u[0],
                        1: [u[7]-1], #-1 is to comply with bedtools getfasta
                        2: u[8],
                        3: u[3],
                        4: u[4],
                        5: u[5],
                        6: [gene_name]
                    })
                    upex = formatted_df.to_csv(sep='\t', header=False, index=False)
                    
                    formatted_df = pd.DataFrame({
                        0: d[0],
                        1: [d[7]-1], #-1 is to comply with bedtools getfasta
                        2: d[8],
                        3: d[3],
                        4: d[4],
                        5: d[5],
                        6: [gene_name]
                    })
                    dsex = formatted_df.to_csv(sep='\t', header=False, index=False)
                    
                    # Saving exon size instead for sashimi plots
                    formatted_df = pd.DataFrame({
                        0: u[0],
                        1: [u[7]-1], #-1 is to comply with bedtools getfasta
                        2: u[8],
                        3: u[9],
                        4: u[4],
                        5: u[5],
                        6: [gene_name]
                    })
                    upex1 = formatted_df.to_csv(sep='\t', header=False, index=False)
                    
                    formatted_df = pd.DataFrame({
                        0: d[0],
                        1: [d[7]-1], #-1 is to comply with bedtools getfasta
                        2: d[8],
                        3: d[9],
                        4: d[4],
                        5: d[5],
                        6: [gene_name]
                    })
                    dsex1 = formatted_df.to_csv(sep='\t', header=False, index=False)
                    
                    # Finally paste three segments to a bed file
                    # First remove if any such file exists
                    if os.path.isfile(f"{gene_name}_nt.bed"):
                        os.remove(f"{gene_name}_nt.bed")
                    
                    with open(f"{gene_name}_nt.bed", 'a') as file:
                        file.write(upex + '\n') # All_junctions.bed
                        file.write(dsex + '\n') # All_junctions.bed
                    
                    with open(f"{gene_name}_nt.bed", 'r') as src_file:
                        content = src_file.read()
                
                    with open("res_skiptics/skiptics_all_coordinates.bed", 'a') as dest_file:
                        dest_file.write(content)
                    
                    # And also save for sashimi plote
                    with open(f"{gene_name}_nt1.bed", 'a') as file:
                        file.write(upex1 + '\n') # All_junctions.bed
                        file.write(dsex1 + '\n') # All_junctions.bed
                    
                    with open(f"{gene_name}_nt1.bed", 'r') as src_file:
                        content1 = src_file.read()
                
                    with open("res_skiptics/skiptics_sashimi.bed", 'a') as dest_file:
                        dest_file.write(content1)
                    
                    if os.path.isfile(f"{gene_name}_nt.bed"):
                        os.remove(f"{gene_name}_nt.bed")
                    
                    if os.path.isfile(f"{gene_name}_nt1.bed"):
                        os.remove(f"{gene_name}_nt1.bed")
                    
                    with open(sample, 'r') as infile, open("res_skiptics/IGV_skiptics.csv", 'a') as outfile:
                        for line in infile:
                            fields = line.strip().split(',')
                            formatted_line = f"{fields[0]}:{fields[1]}-{fields[2]},{fields[0]},{fields[1]},{fields[2]},{fields[3]},{fields[4]},{fields[5]},{gene_name},{TxID},{gene_id}\n"
                            outfile.write(formatted_line)
                    
            else :
                with open(sample, 'r') as infile, open("res_skiptics/IGV_all_others.csv", 'a') as outfile:
                    for line in infile:
                        fields = line.strip().split(',')
                        formatted_line = f"{fields[0]}:{fields[1]}-{fields[2]},{fields[0]},{fields[1]},{fields[2]},{fields[3]},{fields[4]},{fields[5]},{gene_name}\n"
                        outfile.write(formatted_line)
                
                with open("all_non_skiptics.csv", "a", newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(csv_ln)
                
        elif (diff_exon_abs>=1) and (dsovlp==0) and (usovlp==0) : # THIS IS FOR MULTI-SKIP
            # As algorithm detects middle exon for exon_skip events, so determine coordinates of actual exons involved in the event
            us_df = pd.read_csv(io.StringIO(us), delim_whitespace=True, header=None)
            ds_df = pd.read_csv(io.StringIO(ds), delim_whitespace=True, header=None)

            formatted_udf = pd.DataFrame({
                0: us_df[0],
                1: us_df[1] - 10,
                2: us_df[2],
                3: us_df[3],
                4: us_df[4],
                5: us_df[5]
            })
            
            formatted_ddf = pd.DataFrame({
                0: ds_df[0],
                1: ds_df[1],
                2: ds_df[2] + 10,
                3: ds_df[3],
                4: ds_df[4],
                5: ds_df[5]
            })
            
            formatted_udf.to_csv("dump_us.bed", sep="\t", header=False, index=False)
            formatted_ddf.to_csv("dump_ds.bed", sep="\t", header=False, index=False)
            
            a = pybedtools.BedTool("dump_ds.bed")
            b = pybedtools.BedTool(f'event_bedfiles/t{allexons}')
            
            closest = a.closest(b, s=True, D="a", iu=True, d=True, t="last")
            dsn = closest.to_dataframe(names = [
                "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
                "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
            ])
        
            a = pybedtools.BedTool("dump_us.bed")
            
            closest = a.closest(b, s=True, D="a", id=True, d=True, t="first")
            usn = closest.to_dataframe(names = [
                "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
                "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
            ])
            
            u = usn.split('\t')
            d = dsn.split('\t')
            
            formatted_df = pd.DataFrame({
                0: u[0],
                1: [u[7]-1],
                2: u[8],
                3: u[3],
                4: u[4],
                5: u[5],
                6: [gene_name]
            })
            upex = formatted_df.to_csv(sep='\t', header=False, index=False)
            
            formatted_df = pd.DataFrame({
                0: d[0],
                1: [d[7]-1],
                2: d[8],
                3: d[3],
                4: d[4],
                5: d[5],
                6: [gene_name]
            })
            dsex = formatted_df.to_csv(sep='\t', header=False, index=False)
            
            # Saving exon size instead for sashimi plots
            formatted_df = pd.DataFrame({
                0: u[0],
                1: [u[7]-1], #-1 is to comply with bedtools getfasta
                2: u[8],
                3: u[9],
                4: u[4],
                5: u[5],
                6: [gene_name]
            })
            upex1 = formatted_df.to_csv(sep='\t', header=False, index=False)
            
            formatted_df = pd.DataFrame({
                0: d[0],
                1: [d[7]-1], #-1 is to comply with bedtools getfasta
                2: d[8],
                3: d[9],
                4: d[4],
                5: d[5],
                6: [gene_name]
            })
            dsex1 = formatted_df.to_csv(sep='\t', header=False, index=False)
            
            # Finally paste three segments to a bed file
            # First remove if any such file exists
            if os.path.isfile(f"{gene_name}_nt.bed"):
                os.remove(f"{gene_name}_nt.bed")
            
            with open(f"{gene_name}_nt.bed", 'a') as file:
                file.write(upex + '\n') # All_junctions.bed
                file.write(dsex + '\n') # All_junctions.bed
            
            with open(f"{gene_name}_nt.bed", 'r') as src_file:
                content = src_file.read()
        
            with open("res_skiptics/skiptics_all_coordinates.bed", 'a') as dest_file:
                dest_file.write(content)
            
            # And also save for sashimi plote
            with open(f"{gene_name}_nt1.bed", 'a') as file:
                file.write(upex1 + '\n') # All_junctions.bed
                file.write(dsex1 + '\n') # All_junctions.bed
            
            with open(f"{gene_name}_nt1.bed", 'r') as src_file:
                content1 = src_file.read()
        
            with open("res_skiptics/skiptics_sashimi.bed", 'a') as dest_file:
                dest_file.write(content1)
            
            if os.path.isfile(f"{gene_name}_nt.bed"):
                os.remove(f"{gene_name}_nt.bed")
            
            if os.path.isfile(f"{gene_name}_nt1.bed"):
                os.remove(f"{gene_name}_nt1.bed")
            
            with open(sample, 'r') as infile, open("res_skiptics/IGV_skiptics.csv", 'a') as outfile:
                for line in infile:
                    fields = line.strip().split(',')
                    formatted_line = f"{fields[0]}:{fields[1]}-{fields[2]},{fields[0]},{fields[1]},{fields[2]},{fields[3]},{fields[4]},{fields[5]},{gene_name},{TxID},{gene_id}\n"
                    outfile.write(formatted_line)
            
            with open(sample, 'r') as infile, open("res_skiptics/IGV_multi_skiptics.csv", 'a') as outfile:
                for line in infile:
                    fields = line.strip().split(',')
                    formatted_line = f"{fields[0]}:{fields[1]}-{fields[2]},{fields[0]},{fields[1]},{fields[2]},{fields[3]},{fields[4]},{fields[5]},{gene_name}\n"
                    outfile.write(formatted_line)
            
        else :
            with open(sample, 'r') as infile, open("res_skiptics/IGV_all_others.csv", 'a') as outfile:
                for line in infile:
                    fields = line.strip().split(',')
                    formatted_line = f"{fields[0]}:{fields[1]}-{fields[2]},{fields[0]},{fields[1]},{fields[2]},{fields[3]},{fields[4]},{fields[5]},{gene_name}\n"
                    outfile.write(formatted_line)
            
            with open("all_non_skiptics.csv", "a", newline='') as f:
                writer = csv.writer(f)
                writer.writerow(csv_ln)
            
####### NEW CODE
with open('res_skiptics/skiptics_all_coordinates.bed', 'r') as file:
    lines = file.readlines()
    nrecrds = len(lines)

nrecrdst = nrecrds/2

# PRINT STATS SO FAR
with open(inpfile, 'r') as file:
    lines = file.readlines()
    total_events = len(lines)

with open("res_skiptics/Summary_stats.txt", "a") as f :
    f.write(f"ALL EVENTS READ ARE: {total_events}\n")

ensdb_notfound = 0

if os.path.exists("res_skiptics/EnsDB_tx_not_found.csv"):
    with open("res_skiptics/EnsDB_tx_not_found.csv", 'r') as file:
        lines = file.readlines()
        ensdb_notfound = len(lines)

with open("res_skiptics/Summary_stats.txt", "a") as f :
    f.write(f"EVENTS NOT FOUND IN EnsDB are: {ensdb_notfound} , PLEASE SEE res_skiptics/EnsDB_tx_not_found.csv file\n")

all_others = 0

if os.path.exists("res_skiptics/IGV_all_others.csv"):
    with open("res_skiptics/IGV_all_others.csv", 'r') as file:
        lines = file.readlines()
        all_others = len(lines)

remaining_events_processed = total_events - ensdb_notfound

with open("res_skiptics/Summary_stats.txt", "a") as f :
    f.write(f"REMAINING EVENTS PROCESSED ARE: {remaining_events_processed}\n")
    f.write(f"OUT OF THESE {remaining_events_processed} events, non_skiptics are: {all_others}, PLEASE SEE res_skiptics/IGV_all_others.csv and all_non_skiptics.csv files\n")
    f.write(f"TOTAL SKIPTIC EVENTS ARE: {nrecrdst}\n")

if nrecrdst>0 :
    multi_es = 0
    
    if os.path.exists("res_skiptics/IGV_multi_skiptics.csv"):
        with open("res_skiptics/IGV_multi_skiptics.csv", 'r') as file:
            lines = file.readlines()
            multi_es = len(lines)
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write(f"OF THESE {nrecrdst} SKIPTIC EVENTS, MULTI_SKIPTIC EVENTS ARE: {multi_es}, PLEASE SEE res_skiptics/IGV_multi_skiptics.csv\n")
    
    print("NOW CHECKING FOR REPEATED EVENTS IF ANY")
    
    with open('res_skiptics/skiptics_all_coordinates.bed', 'r') as file:
        all_data = [line.strip() for line in file]

    # Also read csv for TxID
    with open('res_skiptics/IGV_skiptics.csv', 'r') as file:
        all_csv_data = [line.strip() for line in file]
    
    # Also read sashimi bed file to get unique entries
    with open('res_skiptics/skiptics_sashimi.bed', 'r') as file:
        all_sashimi_data = [line.strip() for line in file]
    
    i = 0
    eventn = 0
    
    while i<nrecrds :
        line11 = all_data[i]
        line12 = all_data[i+1]
        
        # Sashimi records
        lines1 = all_sashimi_data[i]
        lines2 = all_sashimi_data[i+1]
        
        txid_recrd = all_csv_data[eventn]
        txid = txid_recrd.split(',')[8]
        
        eventn = eventn + 1
        i = i + 2
        j = i
        
        # Now go through rest of the data
        flg = 0
        
        while j<nrecrds :
            line21 = all_data[j]
            line22 = all_data[j+1]
            
            if (line11==line21) and (line12==line22) :
                flg = 1
            
            j = j + 2
            
        if flg==0 :
            columns11 = line11.split("\t")[:7]
            columns12 = line12.split("\t")[:7]
            
            with open('res_skiptics/skiptics_unique.bed', 'a') as output_file:
                output_file.write("\t".join(columns11) + "\n")
                output_file.write("\t".join(columns12) + "\n")
            
            # This is used of sashimi plots, should do better
            columns1 = lines1.split("\t")[:7].append(txid)
            columns2 = lines2.split("\t")[:7].append(txid)
            
            with open('res_skiptics/skiptics_uniq_sashimi.bed', 'a') as output_file:
                output_file.write("\t".join(columns1) + "\n")
                output_file.write("\t".join(columns2) + "\n")
            
        else :
            columns11 = line11.split("\t")[:7]
            columns12 = line12.split("\t")[:7]
            
            with open('res_skiptics/skiptics_repeated.bed', 'a') as output_file:
                output_file.write("\t".join(columns11) + "\n")
                output_file.write("\t".join(columns12) + "\n")
            
            columns11 = line11.split(",")
            columns12 = line12.split(",")
            
            col11 = f"{columns11[0]}:{columns11[1]}-{columns11[2]}" + columns11[:7]
            col12 = f"{columns12[0]}:{columns12[1]}-{columns12[2]}" + columns12[:7]
            
            with open('res_skiptics/IGV_skiptics_repeated.bed', 'a') as output_file:
                output_file.write("\t".join(col11) + "\n")
                output_file.write("\t".join(col12) + "\n")
            
    uniq_eventst = 0
    
    if os.path.exists("res_skiptics/skiptics_unique.bed"):
        with open("res_skiptics/skiptics_unique.bed", 'r') as file:
            lines = file.readlines()
            uniq_eventst = len(lines)
    
    uniq_events = uniq_eventst/2
    rep_eventst = 0
    
    if os.path.exists("res_skiptics/skiptics_repeated.bed"):
        with open("res_skiptics/skiptics_repeated.bed", 'r') as file:
            lines = file.readlines()
            rep_eventst = len(lines)
    
    rep_events = rep_eventst/2
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write(f"OUT OF TOTAL {nrecrdst} SKIPTIC EVENTS, UNIQUE SKIPTIC EVENTS ARE: {uniq_events}, PLEASE SEE res_skiptics/skiptics_unique.bed\n")
        f.write(f"AND REPEATED SKIPTIC EVENTS ARE: {rep_events}, PLEASE SEE res_skiptics/skiptics_repeated.bed and res_skiptics/IGV_skiptics_repeated.csv\n")
        
    print("NOW STARTING CDS PHASE LIST")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("###########\n")
        f.write("NOW STARTING CDS PHASE LIST\n")
    
    # Also read sashimi bed file to get TxID
    # GET US_EXON FOR CDS PHASE CALCULATIONS
    if os.path.exists("res_skiptics/skiptics_us_exons.txt"):
        os.remove("res_skiptics/skiptics_us_exons.txt")
    
    with open('res_skiptics/skiptics_uniq_sashimi.bed', 'r') as file:
        all_sashimi_data = [line.strip() for line in file]
    
    # Also get gene_id as well
    with open('res_skiptics/IGV_skiptics.csv', 'r') as file:
        csv_data = [line.strip() for line in file]
    
    i = 0
    ev_id = 0
    
    with open("res_skiptics/skiptics_unique.bed", 'r') as file:
        all_sashimi_data = [line.strip() for line in file]
    
    while i < len(all_sashimi_data) - 1 : # Write for IGV
        lines1 = all_sashimi_data[i]
        lines2 = all_sashimi_data[i+1]
        txid = lines1.split('\t')[7]
        strnd = lines1.split('\t')[5]
        
        # Get gene_id
        csvln = csv_data[ev_id]
        ev_id = ev_id + 1
        gid = csvln.split(',')[9]
        i = i + 2
        
        r1 = lines1.split('\t')
        r2 = lines2.split('\t')
        comb = r1 + r2
        comb = comb.split(',')
        col = f"{comb[0]}:{comb[2]}-{comb[8]},{comb[0]},{comb[2]},{comb[8]}" + columns12[4:8] + txid
        
        with open('res_skiptics/IGV_unique_skiptics_translated.csv', 'a') as output_file:
            output_file.write("\t".join(col) + "\n")
        
        col = f"{comb[0]},{comb[2]},{comb[8]},{comb[5]},{comb[6]}" + gid
        
        with open('res_skiptics/IGV_skiptics_uniq.csv', 'a') as output_file:
            output_file.write("\t".join(col) + "\n")
        
        # Also save bed file containing up_stream exons to retrieve cds
        if strnd=='+' :
            with open("res_skiptics/skiptics_us_exons.txt", "a") as f :
                f.write(lines1 + '\n')
            
        else :
            with open("res_skiptics/skiptics_us_exons.txt", "a") as f :
                f.write(lines2 + '\n')
            
    print("EVENTS LIST FOR US EXONS FOR CDS PHASE LIST IS DONE, PLEASE SEE  res_skiptics/skiptics_us_exons.txt")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("EVENTS LIST FOR US EXONS FOR CDS PHASE LIST IS DONE, PLEASE SEE  res_skiptics/skiptics_us_exons.txt\n")
    
    ########## NEW CODE ENDS
    fasta_flg = 0
    
    if fasta_flg==1 : # FOR US AND DS FASTA
        print("STARTING NT and AA TRANSLATION FOR US and DS EXONS")
        
        # Now call bedtools getfasta function to get nt sequence from reference_genome.fa file
        bed = BedTool("res_skiptics/skiptics_unique.bed")
        sequences = bed.sequence(fi="GRCh38.p13.genome.fa", s=True)
        sequences.save_seqs("res_skiptics/skiptics_unique_nt.fasta")
        
        # Now remove >chr lines from the res_skipticsulting file
        with open("res_skiptics/skiptics_unique_nt.fasta", 'r') as infile, open("skiptics_unique_nt1.fasta", 'w') as outfile:
            for line in infile:
                if not line.startswith('>chr'):
                    outfile.write(line)
        
        # Combine the two files to get desired csv file - chrXX,start,end,genename
        with open("res_skiptics/skiptics_unique.bed", 'r') as file1, open("skiptics_unique_nt1.fasta", 'r') as file2, open("res_skiptics/skiptics_unique_nt.bed", 'w') as outfile :
            for line1, line2 in zip(file1, file2) :
                line1 = line1.rstrip()
                line2 = line2.rstrip()
                outfile.write(f"{line1}\t{line2}\n")
        
        # Remove temp file
        if os.path.exists("skiptics_unique_nt1.fasta") :
            os.remove("skiptics_unique_nt1.fasta")
        
        # And finally transeq compatible
        with open("res_skiptics/skiptics_unique_nt.bed", 'r') as infile, open("res_skiptics/skiptics_unique_transeq_in.fasta", 'w') as outfile :
            for line in infile:
                line = line.rstrip()
                fields = line.split('\t')
                
                strand = fields[5]
                col7 = fields[6]
                col1 = fields[0]
                col2 = fields[1]
                col3 = fields[2]
                col8 = fields[7]
                
                if strand == "+":
                    header = f">sp|{col7}_{col1}_{col2}_{col3}_plus"
                else:
                    header = f">sp|{col7}_{col1}_{col2}_{col3}_minus"
                
                outfile.write(f"{header}\n{col8}\n")
        
        if os.path.exists("res_skiptics/skiptics_unique_nt.bed") :
            os.remove("res_skiptics/skiptics_unique_nt.bed")
        
        # Now do the 3-frame translation
        TEMPFILE = re.sub(r'\.fasta$', '.temp', "res_skiptics/skiptics_unique_transeq_in.fasta")
        OUTPUTFILE = re.sub(r'\.fasta$', '.trans', "res_skiptics/skiptics_unique_transeq_in.fasta")
        
        # Translate sequence
        with open(TEMPFILE, "w") as out_handle :
            for record in SeqIO.parse("res_skiptics/skiptics_unique_transeq_in.fasta", "fasta"):
                seq = record.seq
                
                for frame in range(3):
                    translated_seq = seq[frame:].translate(to_stop=True)
                    new_record = SeqIO.SeqRecord(
                        translated_seq,
                        id=f"{record.id}_frame{frame+1}",
                        description=f"Translated frame {frame+1}"
                    )
                    SeqIO.write(new_record, out_handle, "fasta")
        
        # Rename sequence
        with open(TEMPFILE, "r") as infile, open(OUTPUTFILE, "w") as outfile :
            for line in infile :
                if line.startswith(">"):
                    line = line.replace(">", ">sp|", 1)
                outfile.write(line)
        
        # Also remove all newlines from the
        with open(OUTPUTFILE, "r") as infile, open("res_skiptics/SKIPTICS_AA.fasta", "w") as outfile :
            sequence = ""
            
            for line in infile:
                line = line.strip()
                
                if line.startswith(">"):
                    if sequence:
                        outfile.write(sequence + "\n")
                    
                    outfile.write(line + "\n")
                    sequence = ""
                else:
                    sequence += line
            
            if sequence:
                outfile.write(sequence + "\n")
        
        # Also remove "$OUTPUTFILE" file
        if os.path.exists(OUTPUTFILE) :
            os.remove(OUTPUTFILE)
        
        # Remove temp file
        if os.path.exists(TEMPFILE) :
            os.remove(TEMPFILE)
        
        my_sym = '>'
        
        with open("res_skiptics/SKIPTICS_AA.fasta", "r") as infile, open("res_skiptics/FINAL_SKIPTICS_AA.fasta", "w") as outfile :
            for line in infile:
                line = line.rstrip()
                
                a = line[0:1]
                
                if my_sym==a :
                    line1 = line
                    i = 0
                    
                else :
                    b = line.replace("*", "\n")
                    
                    for li in b.splitlines() :
                        size = len(li)
                        
                        if size>8 :
                            i = i + 1
                            print(f"{line1}_{i}")
                            print(li)
                
                outfile.write(line + "\n")
        
    # Now find cds for us_exon of each event
    orf_flg = 1
    if orf_flg==1 :
        print("NOW CALLING get_orf_cds.R TO IDENTIFY PROTEIN CODING GENES FROM THE LIST OF IDENTIFIED SKIPTICS AND CDS PHASE")
        
        with open("res_skiptics/Summary_stats.txt", "a") as f :
            f.write("NOW CALLING get_orf_cds.R TO IDENTIFY PROTEIN CODING GENES FROM THE LIST OF IDENTIFIED SKIPTICS AND CDS PHASE\n")
        
        command = [
            "Rscript",
            "get_orf_cds.R",
            "res_skiptics/skiptics_us_exons.txt",
            "res_skiptics/Summary_stats.txt",
            "skiptics"
        ]
        
        result = subprocess.run(command, capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)
        
        if os.path.exists("protein_coding.bed") :
            shutil.move("protein_coding.bed", "res_skiptics/protein_coding.bed")
        
        if os.path.exists("cds_unsuccessful_frames_list.csv") :
            shutil.move("cds_unsuccessful_frames_list.csv", "res_skiptics/cds_unsuccessful_frames_list.csv")
        
        if os.path.exists("cds_successful_frames_list.csv") :
            shutil.move("cds_successful_frames_list.csv", "res_skiptics/cds_successful_frames_list.csv")
        
        if os.path.exists("protein_coding_a.bed") :
            shutil.move("protein_coding_a.bed", "res_skiptics/protein_coding_a.bed") # This is to compare with skiptics_unique.bed
        
        print("DONE WITH CDS LIST and PHASE FOR SKIPTICS, Please see res_skiptics/protein_coding.bed")
        
        with open("res_skiptics/Summary_stats.txt", "a") as f :
            f.write("DONE WITH CDS LIST and PHASE FOR SKIPTICS, Please see res_skiptics/protein_coding.bed\n")
        
    else :
        print("SKIPPING ORF CALCULATIONS, PLEASE SET orf_flg flag to 1 to repeat these calculations")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("SKIPPING ORF CALCULATIONS, PLEASE SET orf_flg flag to 1 to repeat these calculations\n")
    
    # Get difference between two files
    with open("res_skiptics/protein_coding.bed", "r") as pc_file:
        protein_coding_lines = set(line.strip() for line in pc_file)
    
    with open("res_skiptics/skiptics_us_exons.txt", "r") as ue_file:
        us_exons_lines = ue_file.readlines()
    
    with open("res_skiptics/non_protein_coding.bed", "w") as out_file:
        for line in us_exons_lines:
            if line.strip() not in protein_coding_lines:
                out_file.write(line)
    
    with open("res_skiptics/non_protein_coding.bed", "r") as file:
        num_lines = sum(1 for line in file)
    
    print(f"!!!!!!!! GOT {num_lines}  ZERO_CDS OR NON-PROTEIN_CODING EVENTS")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write(f"!!!!!!!! GOT {num_lines}  ZERO_CDS OR NON-PROTEIN_CODING EVENTS\n")
    
    print("NOW UPDATING SKIPTICS LIST WITH CDS, Please see non-protein coding/zero_cds list for non cds events res_skiptics/non_protein_coding.bed")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("NOW UPDATING SKIPTICS LIST WITH CDS, Please see non-protein coding list for non cds events res_skiptics/non_protein_coding.bed\n")
    
    # Also modify skiptics_unique.bed list by removing non-coding genes
    ############
    if os.path.exists("res_skiptics/cds_skiptics_unique.bed"):
        os.remove("res_skiptics/cds_skiptics_unique.bed")
    
    with open("res_skiptics/protein_coding.bed", "r") as file:
        nrecrds = sum(1 for line in file)
    
    with open("res_skiptics/skiptics_unique.bed", "r") as file:
        nrecrds_u = sum(1 for line in file)
    
    with open("res_skiptics/protein_coding_a.bed", "r") as file:
        cds_data = [line.strip() for line in file]
    
    with open("res_skiptics/skiptics_unique.bed", "r") as file:
        uniq_data = [line.strip() for line in file]
    
    # Also MODIFY SASHIMI PLOTS
    if os.path.exists("res_skiptics/cds_skiptics_uniq_sashimi.bed"):
        os.remove("res_skiptics/cds_skiptics_uniq_sashimi.bed")
    
    if os.path.exists("res_skiptics/cds_IGV_unique_skiptics_translated.csv"):
        os.remove("res_skiptics/cds_IGV_unique_skiptics_translated.csv")
    
    with open("res_skiptics/skiptics_uniq_sashimi.bed", "r") as file:
        all_sashimi_bed = [line.strip() for line in file]
    
    with open("res_skiptics/IGV_unique_skiptics_translated.csv", "r") as file:
        all_sashimi_csv = [line.strip() for line in file]
    
    i = 0
    eventn = 0
    
    while i<nrecrds_u :
        line21 = uniq_data[i]
        line22 = uniq_data[i+1]
        j = 0
        
        while j<nrecrds :
            cds = cds_data[j]
            
            if (cds==line21) or (cds==line22) :
                col21 = line21.split("\t")[:7]
                col22 = line22.split("\t")[:7]
                
                with open('res_skiptics/cds_skiptics_unique.bed', 'a') as output_file:
                    output_file.write("\t".join(col21) + "\n")
                    output_file.write("\t".join(col22) + "\n")
                
                # Also copy sashimi bed and csv files
                bed21 = all_sashimi_bed[i]
                bed22 = all_sashimi_bed[i+1]
                
                colbed21 = bed21.split("\t")[:8]
                colbed22 = bed22.split("\t")[:8]
                
                with open('res_skiptics/cds_skiptics_uniq_sashimi.bed', 'a') as output_file:
                    output_file.write("\t".join(colbed21) + "\n")
                    output_file.write("\t".join(colbed22) + "\n")
                
                csv21 = all_sashimi_csv[eventn]
                
                with open("res_skiptics/cds_IGV_unique_skiptics_translated.csv", "a", newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(csv21)
                
                break
            
            j = j + 1
        
        i = i + 2
        eventn = eventn + 1
    
    print("DONE WITH GENERATION VALID CDS bed and csv files for FINAL NT AND AA TRANSLATION, Please see res_skiptics/cds_skiptics_unique.bed, res_skiptics/cds_skiptics_uniq_sashimi.bed and res_skiptics/cds_IGV_unique_skiptics_translated.csv")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("DONE WITH GENERATION VALID CDS bed and csv files for FINAL NT AND AA TRANSLATION, Please see res_skiptics/cds_skiptics_unique.bed, res_skiptics/cds_skiptics_uniq_sashimi.bed and res_skiptics/cds_IGV_unique_skiptics_translated.csv\n")
    
    ############ Done modifying skiptics_unique.bed list by removing non-coding genes
    print("AND NOW STARTING NT and AA TRANSLATION FOR CONCATENATED EXON RANGES")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("AND NOW STARTING NT and AA TRANSLATION FOR CONCATENATED EXON RANGES\n")
    
    # NOW concatenate nt sequence for exon_skip
    # Now form concatenated coordinates
    if os.path.exists("res_skiptics/cds_SKIPTICS_FUSED_AA.fasta"):
        os.remove("res_skiptics/cds_SKIPTICS_FUSED_AA.fasta")
    
    if os.path.exists("res_skiptics/cds_skiptics_fused_transeq_in.fasta"):
        os.remove("res_skiptics/cds_skiptics_fused_transeq_in.fasta")
    
    event_i = 1
    i = 0
    
    # Also read unsuccessful_frames_list.csv file to decide the frame of translation
    with open("res_skiptics/cds_successful_frames_list.csv", "r") as file:
        cds_data = [line.strip() for line in file]
    
    with open("res_skiptics/cds_skiptics_unique.bed", "r") as file:
        lines = file.readlines()
    
    for i in range(0, len(lines), 2):
        r1 = lines[i].strip()
        if i + 1 < len(lines):
            r2 = lines[i + 1].strip()
        
        strnd1 = r1[5]
        if strnd1=='+' :
            concatenated = f"{r1} {r2}"
            formatted = concatenated.replace(' ', '\n')
            
            with open("dump.bed", "w") as file:
                file.write(formatted)
            
        else :
            concatenated = f"{r2} {r1}"
            formatted = concatenated.replace(' ', '\n')
            
            with open("dump.bed", "w") as file:
                file.write(formatted)
           
            bed_df = pd.read_csv("dump.bed", sep='\t', header=None, names=['chrom', 'start', 'end'])
            genome = SeqIO.to_dict(SeqIO.parse("GRCh38.p13.genome.fa", "fasta"))
            
            with open("dump.fasta", "w") as out_file:
                for _, row in bed_df.iterrows():
                    chrom = row['chrom']
                    start = int(row['start'])
                    end = int(row['end'])
                    
                    if chrom in genome:
                        sequence = genome[chrom].seq[start:end]
                        out_file.write(f">{chrom}:{start}-{end}\n")
                        out_file.write(f"{sequence}\n")
    
        
        # Now remove >chr lines from the resulting file
        with open("dump.fasta", "r") as infile, open("dump1.fasta", "w") as outfile:
            for line in infile:
                if not line.startswith(">chr"):
                    outfile.write(line)
            
        # Combine the two files to get desired bed file - chrXX,start,end,genename
        if strnd1=='+' :
            with open("dump.bed", "r") as file :
                content = file.read().replace('\n', ' ')
            
            fields = content.split()
            title = f">sp|{fields[13]}_{fields[0]}_{fields[1]}_{fields[2]}_{fields[8]}_{fields[9]}_{event_i}_plus"
            
            with open("dump.fasta", "r") as infile :
                filtered_lines = [line.strip() for line in infile if not line.startswith(">chr")]
            
            concatenated = ''.join(filtered_lines)
            formatted = '\n'.join([f"{title} {char}" for char in concatenated])
            
            with open("res_skiptics/skiptics_fused_transeq_in1.fasta", "w") as outfile:
                outfile.write(formatted)
            
        else :
            # NOW CHANGIN
            with open("dump.bed", "r") as file :
                content = file.read().replace('\n', ' ')
            
            fields = content.split()
            title = f">sp|{fields[13]}_{fields[0]}_{fields[1]}_{fields[2]}_{fields[8]}_{fields[9]}_{event_i}_minus"
            
            with open("dump.fasta", "r") as infile :
                filtered_lines = [line.strip() for line in infile if not line.startswith(">chr")]
            
            concatenated = ''.join(filtered_lines)
            formatted = '\n'.join([f"{title} {char}" for char in concatenated])
            
            with open("res_skiptics/skiptics_fused_transeq_in1.fasta", "w") as outfile:
                outfile.write(formatted)
            
        with open("res_skiptics/skiptics_fused_transeq_in1.fasta", "r") as src, open("res_skiptics/cds_skiptics_fused_transeq_in.fasta", "a") as dest:
            dest.write(src.read())
    
        # Also get translation fram
        frame_line = cds_data[i]
        framen = frame_line.split(',')[9]
        framen = framen + 1
        
        # Now do the 3-frame translation
        TEMPFILE = re.sub(r'\.fasta$', '.temp', "res_skiptics/skiptics_fused_transeq_in1.fasta")
        OUTPUTFILE = re.sub(r'\.fasta$', '.trans', "res_skiptics/skiptics_fused_transeq_in1.fasta")
        
        # Translate sequence
        seq_records = SeqIO.parse("res_skiptics/skiptics_fused_transeq_in1.fasta", "fasta")
    
        translated_records = []
        for record in seq_records:
            protein_seq = record.seq.translate(to_stop=True)
            translated_record = SeqRecord(protein_seq, id=record.id, description=record.description)
            translated_records.append(translated_record)
        
        with open(TEMPFILE, "w") as outfile:
            SeqIO.write(translated_records, outfile, "fasta")
        
        # Rename sequence
        with open(TEMPFILE, "r") as infile, open(OUTPUTFILE, "w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    line = line.replace(">", ">sp|", 1)
                outfile.write(line)
        
        # Also remove all newlines from the
        with open(OUTPUTFILE, "r") as infile, open("res_skiptics/cds_SKIPTICS_FUSED_AA.fasta", "a") as outfile:
            sequence = ""
            for line in infile:
                if line.startswith(">"):
                    if sequence:
                        outfile.write(sequence + "\n")
                        sequence = ""
                    outfile.write(line)
                else:
                    sequence += line.strip()
    
            if sequence:
                outfile.write(sequence + "\n")
        
        # Also remove "$OUTPUTFILE" file
        if os.path.exists(OUTPUTFILE):
            os.remove(OUTPUTFILE)
        
        # Remove temp file
        if os.path.exists(TEMPFILE):
            os.remove(TEMPFILE)
        
        event_i = event_i + 1
        i = i + 1
        
    if os.path.exists("res_skiptics/skiptics_fused_transeq_in1.fasta"):
        os.remove("res_skiptics/skiptics_fused_transeq_in1.fasta")
    
    print("GOT FUSED NT and AA FASTA FILE FOR SKIPTICS EVENTS, PLEASE SEE res_skiptics/cds_skiptics_fused_transeq_in.fasta and res_skiptics/cds_SKIPTICS_FUSED_AA.fasta file")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("GOT FUSED NT and AA FASTA FILE FOR SKIPTICS EVENTS, PLEASE SEE res_skiptics/cds_skiptics_fused_transeq_in.fasta and res_skiptics/cds_SKIPTICS_FUSED_AA.fasta file\n")
    
    # Cleanup
    if os.path.exists("dump1.fasta"):
        os.remove("dump1.fasta")
    
    print("NOW CALL R SCRIPT TO MATCH cannonical frame to get final_aa.fasta")
    
    command = [
        "Rscipt",
        "check_aaV4_allFrames.R",
        "res_skiptics/aa.fasta",
        "res_skiptics/cds_SKIPTICS_FUSED_AA.fasta",
        "res_skiptics/Summary_stats.txt"
    ]
    
    result = subprocess.run(command, capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)
    
    stop_codon_flg = 1
    
    if stop_codon_flg==1 :
        # Now truncate lines at *
        with open("res_skiptics/final_aa.fasta", 'r') as infile, open("res_skiptics/cds_SKIPTICS_FUSED_AA1.fasta", 'w') as outfile :
            for line in infile:
                cleaned_line = re.sub(r'\*.*', '', line)
                outfile.write(cleaned_line)
        
        print("DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS, PLEASe SEE res_skiptics/cds_SKIPTICS_FUSED_AA1.fasta file")
        
        with open("res_skiptics/Summary_stats.txt", "a") as f :
            f.write("DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS, PLEASe SEE res_skiptics/cds_SKIPTICS_FUSED_AA1.fasta file\n")
        
        # Also truncate lines at * for allframes files
        with open("res_skiptics/all_frames_final_aa.fasta", 'r') as infile, open("res_skiptics/AllFrames_SKIPTICS_FUSED_AA1.fasta", 'w') as outfile:
            for line in infile:
                cleaned_line = re.sub(r'\*.*', '', line)
                outfile.write(cleaned_line)
    
        
        print("DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS FOR ALL FRAMES FASTA, PLEASe SEE res_skiptics/AllFrames_SKIPTICS_FUSED_AA1.fasta file")
        
        with open("res_skiptics/Summary_stats.txt", "a") as f :
            f.write("DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS FOR ALL FRAMES FASTA, PLEASe SEE res_skiptics/AllFrames_SKIPTICS_FUSED_AA1.fasta file\n")
        
    else :
        with open("res_skiptics/final_aa.fasta", 'r') as infile, open("res_skiptics/cds_SKIPTICS_FUSED_AA1.fasta", 'w') as outfile:
            outfile.write(infile.read())
        
        with open("res_skiptics/all_frames_final_aa.fasta", 'r') as infile, open("res_skiptics/AllFrames_SKIPTICS_FUSED_AA1.fasta", 'w') as outfile:
            outfile.write(infile.read())
        
    # Delete empty line if any
    print("NOW REMOVING EMPTY LINES LEFT BY EVENTS WITH STOP CODONS IN THE BEGINNING ""ARHGAP22_chr10"" is one EXAMPLE")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("NOW REMOVING EMPTY LINES LEFT BY EVENTS WITH STOP CODONS IN THE BEGINNING ""ARHGAP22_chr10"" is one EXAMPLE\n")
    
    if os.path.exists("res_skiptics/cds_SKIPTICS_FUSED_AA2.fasta"):
        os.remove("res_skiptics/cds_SKIPTICS_FUSED_AA2.fasta")
    
    if os.path.exists("res_skiptics/cds_SKIPTICS_FUSED_AA_Empty.fasta"):
        os.remove("res_skiptics/cds_SKIPTICS_FUSED_AA_Empty.fasta")
    
    i = 0
    
    with open("res_skiptics/cds_SKIPTICS_FUSED_AA1.fasta", 'r') as infile :
        lines = iter(infile)
        
        for r1 in lines:
            r2 = next(lines, None)
        
        if r2!="" :
            with open("res_skiptics/cds_SKIPTICS_FUSED_AA2.fasta", "a") as f :
                f.write(r1 + '\n')
                f.write(r2 + '\n')
            
        else :
            with open("res_skiptics/cds_SKIPTICS_FUSED_AA_Empty.fasta", "a") as f :
                f.write(r1 + '\n')
                
            i = i + 1
            
    with open('res_skiptics/cds_SKIPTICS_FUSED_AA_Empty.fasta', 'r') as file:
        lines = file.readlines()
        nb_lines = len(lines)
    
    if nb_lines>0 :
        print(f"DONE WITH REMOVING EMPTY LINES GOT {nb_lines} SKIPTICS HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_skiptics/cds_SKIPTICS_FUSED_AA_Empty.fasta file")
        
        with open("res_skiptics/Summary_stats.txt", "a") as f :
            f.write(f"DONE WITH REMOVING EMPTY LINES GOT {nb_lines} SKIPTICS HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_skiptics/cds_SKIPTICS_FUSED_AA_Empty.fasta file\n")
        
    else :
        with open("res_skiptics/cds_SKIPTICS_FUSED_AA1.fasta", 'r') as infile, open("res_skiptics/cds_SKIPTICS_FUSED_AA2.fasta", 'w') as outfile :
            outfile.write(infile.read())
        
    ########## ALLFRAMES FASTA FILE
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("NOW REMOVING EMPTY LINES LEFT BY EVENTS WITH STOP CODONS IN THE BEGINNING ""ARHGAP22_chr10"" is one EXAMPLE\n")
    
    if os.path.exists("res_skiptics/AllFrames_SKIPTICS_FUSED_AA2.fasta"):
        os.remove("res_skiptics/AllFrames_SKIPTICS_FUSED_AA2.fasta")
    
    if os.path.exists("res_skiptics/AllFrames_SKIPTICS_FUSED_AA_Empty.fasta"):
        os.remove("res_skiptics/AllFrames_SKIPTICS_FUSED_AA_Empty.fasta")
    
    i = 0
    
    with open("res_skiptics/AllFrames_SKIPTICS_FUSED_AA1.fasta", 'r') as infile :
        lines = iter(infile)
        
        for r1 in lines:
            r2 = next(lines, None)
        
        if r2!="" :
            with open("res_skiptics/AllFrames_SKIPTICS_FUSED_AA2.fasta", "a") as f :
                f.write(r1 + '\n')
                f.write(r2 + '\n')
            
        else :
            with open("res_skiptics/AllFrames_SKIPTICS_FUSED_AA_Empty.fasta", "a") as f :
                f.write(r1 + '\n')
            
            i = i + 1
        
    with open('res_skiptics/AllFrames_SKIPTICS_FUSED_AA_Empty.fasta', 'r') as file:
        lines = file.readlines()
        nb_lines = len(lines)
    
    if nb_lines>0 :
        print("DONE WITH REMOVING EMPTY LINES GOT {nb_lines} SKIPTICS HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_skiptics/AllFrames_SKIPTICS_FUSED_AA_Empty.fasta file")
        
        with open("res_skiptics/Summary_stats.txt", "a") as f :
            f.write("DONE WITH REMOVING EMPTY LINES GOT {nb_lines} SKIPTICS HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_skiptics/AllFrames_SKIPTICS_FUSED_AA_Empty.fasta file\n")
        
    else :
        with open("res_skiptics/AllFrames_SKIPTICS_FUSED_AA1.fasta", 'r') as infile, open("res_skiptics/AllFrames_SKIPTICS_FUSED_AA2.fasta", 'w') as outfile :
            outfile.write(infile.read())
        
    ##########
    print("STARTED REMOVING AA < 8")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("STARTED REMOVING AA < 8")
    
    my_sym = '>'
    
    with open("res_skiptics/cds_SKIPTICS_FUSED_AA2.fasta", "r") as infile, open("res_skiptics/cds_PEAKS_SKIPTICS_FUSED_AA1.fasta", "w") as outfile :
        for line in infile:
            line = line.rstrip()
            
            a = line[0:1]
            
            if my_sym==a :
                line1 = line
                i = 0
                
            else :
                b = line.replace("*", "\n")
                
                for li in b.splitlines() :
                    size = len(li)
                    
                    if size>8 :
                        i = i + 1
                        print(f"{line1}_{i}")
                        print(li)
            
            outfile.write(line + "\n")
    
    # Also remove trailing X in fasta files
    with open("res_skiptics/cds_PEAKS_SKIPTICS_FUSED_AA1.fasta", 'r') as infile, open("res_skiptics/cds_PEAKS_SKIPTICS_FUSED_AA.fasta", 'w') as outfile:
            for line_number, line in enumerate(infile, start=1):
                if line_number % 2 == 0:
                    line = line.rstrip('\n')
                    if line.endswith('X'):
                        line = line[:-1]
                    outfile.write(line + '\n')
                else:
                    outfile.write(line)
    
    if os.path.exists("res_skiptics/cds_PEAKS_SKIPTICS_FUSED_AA1.fasta"):
        os.remove("res_skiptics/cds_PEAKS_SKIPTICS_FUSED_AA1.fasta")
    
    ######### ALLFRAMES
    my_sym = '>'
    
    with open("res_skiptics/AllFrames_SKIPTICS_FUSED_AA2.fasta", "r") as infile, open("res_skiptics/AllFrames_PEAKS_SKIPTICS_FUSED_AA1.fasta", "w") as outfile :
        for line in infile:
            line = line.rstrip()
            
            a = line[0:1]
            
            if my_sym==a :
                line1 = line
                i = 0
                
            else :
                b = line.replace("*", "\n")
                
                for li in b.splitlines() :
                    size = len(li)
                    
                    if size>8 :
                        i = i + 1
                        print(f"{line1}_{i}")
                        print(li)
            
            outfile.write(line + "\n")
    
    # Also remove trailing X in fasta files
    
    
    if os.path.exists("res_skiptics/AllFrames_PEAKS_SKIPTICS_FUSED_AA1.fasta"):
        os.remove("res_skiptics/AllFrames_PEAKS_SKIPTICS_FUSED_AA1.fasta")
    
    ###################
    
    print("DONE WITH REMOVING AA "<" 8, FINAL FILE IS res_skiptics/cds_PEAKS_SKIPTICS_FUSED_AA.fasta")
    
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write("DONE WITH REMOVING AA "<" 8, FINAL FILE IS res_skiptics/cds_PEAKS_SKIPTICS_FUSED_AA.fasta\n")
    
else :
    # Final clean up
    with open("res_skiptics/Summary_stats.txt", "a") as f :
        f.write(f"AS SKIPTIC EVENTS ARE: {nrecrdst} , SO EXITING\n")
    
if os.path.exists("dump_ds.bed"):
    os.remove("dump_ds.bed")

if os.path.exists("dump_us.bed"):
    os.remove("dump_us.bed")

if os.path.exists("dump.bed"):
    os.remove("dump.bed")

if os.path.exists("dump.fasta"):
    os.remove("dump.fasta")

with open("res_skiptics/Summary_stats.txt", "a") as f :
    f.write("THOSE ARE ALL PERTINENT STATISTICS, PLEASE LET US KNOW IF ANY OTHER INFORMATION MIGHT BE USEFUL!!!\n")
    f.write("######## DONE STATS FROM esV5_layered_CDSV1.py")
