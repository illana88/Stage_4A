######## Proteogenomic Pipeline for Biomarker Discovery in iPSC neurons
######### Inline Readme
# 1. Get List of most abundant Transcripts from KD iPSC samples
# This step assumes that all bam files are in current folder and follow have
# change directory to part-a and issue run part-a.sh with appropriate option and input file

arg1 = int(input("Argument 1 : "))
arg2 = "selected_events.csv"

import os
import glob
import csv
import subprocess
import re
import pandas as pd
import sys
import pybedtools
import pysam
import shutil
import gffutils
import rpy2.robjects as ro

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
            result = subprocess.run(commande_shell, shell=True, check=True)
            print(result.stdout)
            
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
            
            command = [
                "python",
                "abundant_tx.py"
            ]

            result = subprocess.run(command, capture_output=True, text=True)
            print(result.stdout)

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
            
            
            
            
            
            
# # NEXT SORT FILE ON COLUMN 5 - THIS IS IMPORTANT AS PIPELINE DEPENDS ON IT
if arg1==1 or arg1==2 or arg1==5 :
    with open("Summary_stats.txt","a") as fichier :
        fichier.write("STARTING pgp-a.sh WITH FLAG 1 for all events sashimi plots\n")
    
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
            
            # gtf_file = "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf"
            
            # os.makedirs("database", exist_ok=True)
            # edb_file = os.path.join("database", gtf_file + ".db")
            
            # if not os.path.exists(edb_file):
            #     edb = gffutils.create_db(gtf_file, dbfn=edb_file, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)
            # else:
            #     edb = gffutils.FeatureDB(edb_file, keep_order=True)
            
            # command = [
            #     "python",
            #     "TxEnsDB103_layeredV6.py",
            #     "sorted_selected_events.csv",
            #     "principal_txs.csv",
            #     "temp_all_events_sashimi/FINAL_STATS_ALL_SASHIMIS.txt",
            #     f"{edb}"
            # ]
            
            # result = subprocess.run(command, capture_output=True, text=True)
            # print(result.stdout)

            command = [
                "R",
                "txens.R",
                "sorted_selected_events.csv",
                "principal_txs.csv",
                "temp_all_events_sashimi/FINAL_STATS_ALL_SASHIMIS.txt"
            ]
            
            result = subprocess.run(command, capture_output=True, text=True)
            print(result.stdout)
            
            with open("Summary_stats.txt","a") as fichier :
                fichier.write("************************ BACK FROM TxEnsDB103_layeredV6.R, CONTINUING pgp-a.sh \n")
             
        pd.set_option('display.max_columns', None)
            
        csv_data = pd.read_csv('all_tx_events.csv', header=None)
        csvi = 0
        samples = glob.glob('event_bedfiles/temp_*.bed')
        
        with open("Summary_stats.txt","a") as fichier :
            fichier.write("GENERATING BED FILES FOR EACH EVENT\n")
        
        for sample in samples :
            # Read csv entry
            csv_ln = csv_data.iloc[csvi]
            csvi = csvi + 1
        
            allexons = sample.split('/')[1].split('_')[1]
            gene_name1 = allexons.split('.')[0]
            gene_name = gene_name1.split('-')[0]
        
            # First sort the bed
            bedfile = pybedtools.BedTool(f'event_bedfiles/{allexons}')
            sorted_bed = bedfile.sort()
            sorted_bed.saveas(f'event_bedfiles/t{allexons}')
            
            # Also read Tx Files to retrieve selected Tx - should find better ways
            df = pd.read_csv(f'event_bedfiles/TxID{allexons}', sep='\\s+', header=None)
            TxID = df.iloc[0, 6]
            
            df = pd.read_csv(sample, sep='\\s+', header=None)
            strnd = df.iloc[0, 5]
        
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
        
            # Get up and down stream exon numbers
            upexon = us.iloc[:,10].tolist()[0]
            dnexon = ds.iloc[:,10].tolist()[0]
            
            # Events star and end
            event_st = us.iloc[:, 1].tolist()[0]
            event_end = us.iloc[:, 2].tolist()[0]
            diff_exon = int(upexon) - int(dnexon)
        
            # Take absolute value
            diff_exon_abs = abs(diff_exon)
        
            if diff_exon_abs>=1 : #ALL EVENTS THAT SPANS 2 OR MORE EXONS
            
                if strnd == '+' :
                    start = us.iloc[:, 7].tolist()[0]
                    end = ds.iloc[:, 8].tolist()[0]
                else :
                    start = ds.iloc[:, 7].tolist()[0]
                    end = us.iloc[:, 8].tolist()[0]
                
                # Also save
                # FORMAT IS: chr, start of up_ex,end of ds_ex, 1, 0, strand, gene_name, TxID
            	# First check if event lies between selected exons
                if (float(start)<=float(event_st) and float(end)>=float(event_end)) :
                    liste = [str(element) for element in us.values.tolist()[0]]
                    input_data = liste + [start, end, strnd, gene_name, TxID]
                    output_data = [input_data[0], str(input_data[13]), str(input_data[14]), '1', '0', input_data[15], input_data[16], input_data[17]]
                    
                    with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.bed", "a") as f:
                        f.write("\t".join(output_data) + "\n")
        
                    with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.bed", "a") as f:
                        f.write("\t".join(output_data) + "\n")
                    
                    with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.csv", "a", newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow(csv_ln)
                        
                    with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.csv", "a", newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow(csv_ln)
        
                else :
                    if strnd == '+' :
                        start = us.iloc[:, 7].tolist()[0]
                        end = ds.iloc[:, 8].tolist()[0] # Both are same
                        
                        # First check if star > event_start, then select upstream exon
                        if float(start)>=float(event_st) :
                            # Get exon (line number in bed file) to read
                            exon = ds.iloc[:, 10].tolist()[0]
                            exon = int(exon) - 2
                            
                            with open(f"event_bedfiles/{allexons}", 'r') as file:
                                bed_data = [line.strip() for line in file]
                            
                            bed_ln = bed_data[exon-1].split('\t')
                            
                            # Update start
                            start = bed_ln[1]
                            
                        # Now check if end <event_end
                        if float(end)<=float(event_end) :
                            # Get exon (line number in bed file) to read
                            exon = us.iloc[:, 10].tolist()[0]
        
                            with open(f"event_bedfiles/{allexons}", 'r') as file:
                                bed_data = [line.strip() for line in file]
                            
                            #print('bed_data ',bed_data)
                            bed_ln = bed_data[exon-1].split('\t')
                            #print('exon ', exon, ' bed_ln ', bed_ln)
        
                            # Update end
                            end = bed_ln[2]
                                
                        # Now one more time check if event lies between selected exons
                        if (float(start)<=float(event_st) and float(end)>=float(event_end)):
                            liste = [str(element) for element in us.values.tolist()[0]]
                            input_data = liste + [start, end, strnd, gene_name, TxID]
                            output_data = [input_data[0], str(input_data[13]), str(input_data[14]), '1', '0', input_data[15], input_data[16], input_data[17]]
                            
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.bed", "a") as f:
                                f.write("\t".join(output_data) + "\n")
        
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.bed", "a") as f:
                                f.write("\t".join(output_data) + "\n")
                            
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.csv", "a", newline='') as f:
                                writer = csv.writer(f)
                                writer.writerow(csv_ln)
                                
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.csv", "a", newline='') as f:
                                writer = csv.writer(f)
                                writer.writerow(csv_ln)
                                    
                        else :
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_progress2.txt", "a") as f:
                                f.write(f"ds 1 {ds}\n")
        
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_progress2.txt", "a") as f:
                                f.write(f"us 1 {us}\n")
                            
                            print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                                
                    else : # THIS IS FOR NEGATIVE STRAND
                        
                        start = ds.iloc[:, 7].tolist()[0]
                        end = us.iloc[:, 8].tolist()[0]
        
                        # First check if star > event_start, then select upstream exon
                        if int(start)>=int(event_st) :
                            
                            # Get exon (line number in bed file) to read
                            exon = ds.iloc[:, 10].tolist()[0]
                            exon = int(exon) # Tx on -ve strand has exons listed from bottom to top in increasing order
        
                            with open(f"event_bedfiles/{allexons}", 'r') as file:
                                bed_data = [line.strip() for line in file]
                            
                            bed_ln = bed_data[exon-1].split('\t')
        
                            # Update start
                            start = bed_ln[1]
                                
                        # Now check if end <event_end
                        if int(end)<=int(event_end) :
                            
                            # Get exon (line number in bed file) to read
                            exon = us.iloc[:, 10].tolist()[0]
                            exon = int(exon) - 2 # Tx on -ve strand has exons listed from bottom to top in increasing order
        
                            with open(f"event_bedfiles/{allexons}", 'r') as file:
                                bed_data = [line.strip() for line in file]
        
                            bed_ln = bed_data[exon-1].split('\t')
        
                            # Update end
                            end = bed_ln[2]
        
                        # Now one more time check if event lies between selected exons
                        if (int(start)<=int(event_st) and int(end)>=int(event_end)):
                            liste = [str(element) for element in us.values.tolist()[0]]
                            input_data = liste + [start, end, strnd, gene_name, TxID]
                            output_data = [input_data[0], str(input_data[13]), str(input_data[14]), '1', '0', input_data[15], input_data[16], input_data[17]]
                            
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.bed", "a") as f:
                                f.write("\t".join(output_data) + "\n")
        
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.bed", "a") as f:
                                f.write("\t".join(output_data) + "\n")
                            
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi2.csv", "a", newline='') as f:
                                writer = csv.writer(f)
                                writer.writerow(csv_ln)
                                
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.csv", "a", newline='') as f:
                                writer = csv.writer(f)
                                writer.writerow(csv_ln)
                                    
                        else :
                            with open(f"temp_all_events_sashimi/{splicing_events_file}_progress2.txt", "a") as f:
                                f.write(f"ds 2 {ds}\n")
                                f.write(f"us 2 {us}\n")
                            
                            print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                             
            elif diff_exon_abs==0 :
                if strnd == '+' :
                    start = us.iloc[:, 7].tolist()[0]
                    end = ds.iloc[:, 8].tolist()[0] # Both are same
        
            		 	# First check if star > event_start, then select upstream exon
                    if int(start)>=int(event_st) :
                        
                        # Get exon (line number in bed file) to read
                        exon = us.iloc[:, 10].tolist()[0]
                        exon = int(exon) - 2
        
                        with open(f"event_bedfiles/{allexons}", 'r') as file:
                            bed_data = [line.strip() for line in file]
        
                        bed_ln = bed_data[exon-1].split('\t')
        
                        # Update start
                        start = bed_ln[1]
        
                        # Now go on the other side
                        if int(start)>=int(event_st) : # Get the other exon
                            
                            # Get exon (line number in bed file) to read
                            exon = us.iloc[:, 10].tolist()[0]
            
                            with open(f"event_bedfiles/{allexons}", 'r') as file:
                                bed_data = [line.strip() for line in file]
            
                            bed_ln = bed_data[exon-1].split('\t')
            
                            # Update start
                            start = bed_ln[1]
        
                    # Now check if end <event_end
                    if int(end)<=int(event_end) :
                        
                        # Get exon (line number in bed file) to read
                        exon = us.iloc[:, 10].tolist()[0]
        
                        with open(f"event_bedfiles/{allexons}", 'r') as file:
                            bed_data = [line.strip() for line in file]
        
                        print(f"bed_data {bed_data}")
                        bed_ln = bed_data[exon-1].split('\t')
                        print(f"exon {exon} bed_ln {bed_ln}")
        
                        # Update end
                        end = bed_ln[2]
        
                        if int(end)<=int(event_end) :
                            # Get exon (line number in bed file) to read
                            exon = us.iloc[:, 10].tolist()[0]
                            exon = int(exon) - 2 # Reading line for readarray starts from 0
        
                            with open(f"event_bedfiles/{allexons}", 'r') as file:
                                bed_data = [line.strip() for line in file]
        
                            bed_ln = bed_data[exon-1].split('\t')
        
                            # Update end
                            end = bed_ln[2]
        
                    # Now one more time check if event lies between selected exons
                    if (int(start)<=int(event_st) and int(end)>=int(event_end)) :
                        liste = [str(element) for element in us.values.tolist()[0]]
                        input_data = liste + [start, end, strnd, gene_name, TxID]
                        output_data = [input_data[0], str(input_data[13]), str(input_data[14]), '1', '0', input_data[15], input_data[16], input_data[17]]
                        
                        with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi01.bed", "a") as f:
                            f.write("\t".join(output_data) + "\n")
        
                        with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.bed", "a") as f:
                            f.write("\t".join(output_data) + "\n")
        
                        with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi01.csv", "a", newline='') as f:
                            writer = csv.writer(f)
                            writer.writerow(csv_ln)
                            
                        with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.csv", "a", newline='') as f:
                            writer = csv.writer(f)
                            writer.writerow(csv_ln)
        
                    else :
                        with open(f"temp_all_events_sashimi/{splicing_events_file}_progress01.txt", "a") as f:
                            f.write(f"ds {ds}\n")
                            f.write(f"us {us}\n")
        
                        print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
        
                else : # THIS IS FOR NEGATIVE STRAND
                    start = ds.iloc[:, 7].tolist()[0]
                    end = us.iloc[:, 8].tolist()[0]
        
                    # First check if star > event_start, then select upstream exon
                    if int(start)>=int(event_st) :
                        # Get exon (line number in bed file) to read
                        exon = us.iloc[:, 10].tolist()[0]
        
                        with open(f"event_bedfiles/{allexons}", 'r') as file:
                            bed_data = [line.strip() for line in file]
        
                        bed_ln = bed_data[exon-1].split('\t')
        
                        # Update start
                        start = bed_ln[1]
                    
                    # Now check if end <event_end
                    if int(end)<=int(event_end) :
                        # Get exon (line number in bed file) to read
                        exon =us.iloc[:, 10].tolist()[0]
                        exon = int(exon) - 2 # Tx on -ve strand has exons listed from bottom to top in increasing order
                        
                        with open(f"event_bedfiles/{allexons}", 'r') as file:
                            bed_data = [line.strip() for line in file]
        
                        bed_ln = bed_data[exon-1].split('\t')
        
                        # Update end
                        end = bed_ln[2]
        
                    # Now one more time check if event lies between selected exons
                    if (int(start)<=int(event_st) and int(end)>=int(event_end)) :
                        liste = [str(element) for element in us.values.tolist()[0]]
                        input_data = liste + [start, end, strnd, gene_name, TxID]
                        output_data = [input_data[0], str(input_data[13]), str(input_data[14]), '1', '0', input_data[15], input_data[16], input_data[17]]
                        
                        with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi02.bed", "a") as f:
                            f.write("\t".join(output_data) + "\n")
        
                        with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.bed", "a") as f:
                            f.write("\t".join(output_data) + "\n")
        
                        with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi02.csv", "a", newline='') as f:
                            writer = csv.writer(f)
                            writer.writerow(csv_ln)
                            
                        with open(f"temp_all_events_sashimi/{splicing_events_file}_all_sashimi.csv", "a", newline='') as f:
                            writer = csv.writer(f)
                            writer.writerow(csv_ln)
        
                    else :
                        with open(f"temp_all_events_sashimi/{splicing_events_file}_progress02.txt", "a") as f:
                            f.write(f"ds {ds}\n")
                            f.write(f"us {us}\n")
        
                        print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                        
        
                
                # First check if star > event_start, then select upstream exon
        
            else :
                with open(f"temp_all_events_sashimi/{splicing_events_file}_progress_all.txt", "a") as f:
                    f.write(f"ds {ds}\n")
                    f.write(f"us {us}\n")   
        
        print("################################ DONE GENERATING BED AND OTHER RELATED FILES FOR ALL EVENTS TO CONTINUE FOR SASHIMI PLOTD")
        
        with open("Summary_stats.txt", "a") as f:
            f.write("################################ DONE GENERATING BED AND OTHER RELATED FILES FOR ALL EVENTS TO CONTINUE FOR SASHIMI PLOTD\n")
        
        if arg1==1 or arg1==2 or arg1==5 :
            print("################################ STARTED CREATING SASHIMI PLOTS FOR ALL EVENTS - MAY TAKE MANY HOURS ")
        
            with open("Summary_stats.txt", "a") as f:
                f.write("################################ STARTED CREATING SASHIMI PLOTS FOR ALL EVENTS - MAY TAKE MANY HOURS \n")
        
            # CALL ggsashimi for SASHIMI PLOTS
            command = [
                "python",
                "run_sashimiV1.py",
                f"temp_all_events_sashimi/{splicing_events_file}",
                f"temp_all_events_sashimi/{splicing_events_file}",
                "2"
            ]
            
            result = subprocess.run(command, capture_output=True, text=True)
            print(result.stdout)
    
            print("################################ DONE - RESULTANT pdf is in all_events_sashimi/all_events_sashimi.pdf ")
            
            with open("Summary_stats.txt", "a") as f :
                f.write("################################ DONE - RESULTANT pdf is in all_events_sashimi/all_events_sashimi.pdf \n")
    
        else :
            print("################################ SKIPPING SASHIMI PLOTS ")
            
            with open("Summary_stats.txt", "a") as f :
                f.write("################################ SKIPPING SASHIMI PLOTS \n") ##### Testé jusqu'ici et tout fonctionne #####
            
            
            
            
            
            
# # Step 2 - STARTING ce calculations
if arg1==3 or arg1==5 :
    
    with open ("Summary_stats.txt", "a") as f :
        f.write("############### NOW GENERATING BED FILES FOR ALL BAM SAMPLES, WILL TAKE A WHILE #########\n")
    
    print("############### NOW GENERATING BED FILES FOR ALL BAM SAMPLES, WILL TAKE A WHILE #########")
    
    # Step 1. Calculate TDP_SAMPLE_XXX.bam.bed files for each TDP43KD BAM file
 	# var to disable following
    file = open("all_bams.tsv")
    reader = csv.reader(file,delimiter='\t') 
    samples = []
    
    for row in reader :
        if "TDP43" in row[2] :
            samples.append(row[1])
    
    if "samples" in globals() and len(samples)!=0 :
        os.makedirs("bam_beds",exist_ok=True)
        # For file in bam_files/*.bam
        
        for file in samples :
            print(f"processing {file}")
            
            sample_bam = os.path.basename(file).split('.')[0]
            
            print(f"generating bam_beds/{sample_bam}.bam.bed")
            
            with open("Summary_stats.txt", "a") as f :
                f.write(f"generating bam_beds/{sample_bam}.bam.bed\n")
            
            bed = pybedtools.BedTool(file).bam_to_bed(split=True)
            bed.saveas(f"bam_beds/{sample_bam}.bam.bed")
            
            # Also create genome file containing chromosomes start and end for this bam file
            bam = pysam.AlignmentFile(file, "rb")
            chromosomes = []
            
            for ref in bam.references:
                length = bam.get_reference_length(ref)
                chromosomes.append((ref, length))
            
            df = pd.DataFrame(chromosomes, columns=['Chromosome', 'Length'])
            df.to_csv(f"bam_beds/{sample_bam}.chromosomes.sizes", sep='\t', index=False, header=False)
            
            # Now sort the BED file according to this genome file
            bed = pybedtools.BedTool(f"bam_beds/{sample_bam}.bam.bed")
            sorted_bed = bed.sort(g=f"bam_beds/{sample_bam}.chromosomes.sizes")
            sorted_bed.saveas(f"bam_beds/{sample_bam}-sorted.bam.bed")
            
        if arg1==3 : # Exit BAM to BED creation
            with open("Summary_stats.txt", "a") as f :
                f.write("DONE CREATING BED FIELS FOR BAM SAMPLES - NOW EXITING\n")
                
            sys.exit(1)
            
        else :
            with open("Summary_stats.txt", "a") as f :
                f.write("############### DONE CREATING BED FIELS FOR BAM SAMPLES #########\n")
            
            print("############### DONE CREATING BED FIELS FOR BAM SAMPLES #########")
        
    else :
        print("I did not get any BAM files from all_bams.tsv (from column 2), Please make sure that col3 contains string TDP43 for KD samples, ESITING")
        
        with open("Summary_stats", "a") as f :
            f.write("I did not get any BAM files from all_bams.tsv (from column 2), Please make sure that col3 contains string TDP43 for KD samples, ESITING")






if arg1==4 or arg1==5 :
#     if os.access(arg2,os.F_OK) :
#         print(arg2)
#         splicing_events_file = arg2.split('.')[0]
#         print(splicing_events_file)
#         selected_events = pd.read_csv(arg2, delimiter=',', header=None)
#         sorted_selected_events = selected_events.sort_values(by=selected_events.columns[4])
#         sorted_selected_events.to_csv("sorted_selected_events.csv", index=False, header=False)
        
#     else :
#         print("NO SPLICING_EVENTS csv file is provided, ABORTING!!!!")
#         print("PLEASE PROVIDE SPLICING_EVENTS csv file  'bash pgp_0.sh flag splicing_events.csv' and RERUN")
#         sys.exit(1)
        
#     # EVENTS CSV FILE CLEANING STARTS HERE
#     step00_flg = 1
    
#     if step00_flg==1 :
#         print("Now CLEANING INPUT FILE FOR DUPLICATES (EVENTS WITH SAME CHR#, START AND END) AND STAR > END OR START==END EVENTS- WILL TAKE AROUND 30 MINUTES")
#         if os.path.exists("file_clean_report.txt"):
#             os.remove("file_clean_report.txt")
        
#         # Get total records
#         with open(f"sorted_{splicing_events_file}.csv", 'r') as file :
#             t_recrds = sum(1 for line in file)
        
#         with open("file_clean_report.txt", "a") as f :
#             f.write(f"original file has {t_recrds} records")
        
#         # Delete if already exist
#         if os.path.exists(f"clean_{splicing_events_file}.csv"):
#             os.remove(f"clean_{splicing_events_file}.csv")
        
#         with open(f"sorted_{arg2}", 'r') as file:
#             all_splicing_events = file.read().splitlines() # For saving gene_ids as well to generate sashimi compatible csv file

#         i = 0
#         dups = 0
        
#         while i<t_recrds :
#             line = all_splicing_events[i]
#             chr = line.split(',')[0]
#             start = line.split(',')[1]
#             end = line.split(',')[2]
#             gene_id = line.split(',')[5]
#             lastline = i
#             i = i + 1
            
#             # Now go through each record
#             j = i
            
#             # Now go through rest of the data
#             flg = 0
            
#             while j<t_recrds :
#                 nextline = j + 1
#                 nline = all_splicing_events[j]
#                 nchr = nline.split(',')[0]
#                 nstart = nline.split(',')[1]
#                 nend = nline.split(',')[2]
#                 ngene_id = nline.split(',')[5]
                
#                 if (chr==nchr and start==nstart and end==nend and gene_id==ngene_id) :
#                     with open("file_clean_report.txt", "a") as f :
#                         f.write(f"Line {i} - {line} - has same start: {start} and end: {end} as line {nextline} - {nline} -, so ignoring\n")
                    
#                     dups = dups + 1
#                     flg = 1
                    
#                 j = j + 1
                
#             if (flg==0 and start<end) :
#                 # ALSO CHECK IF START < END
#                 if start<end :
#                     with open(f"clean_{splicing_events_file}.csv", "a") as f :
#                         f.write(f"{line}\n")
                    
#                 else :
#                     with open("file_clean_report.txt", "a") as f :
#                         f.write(f"Line {i} has problem with start: {start} and end: {end}")
                    
#                     dups = dups + 1
                    
#         with open(f"clean_{splicing_events_file}.csv", 'r') as file :
#             c_records = sum(1 for line in file)
        
#         print(" ")
#         print(" ")
#         print("DONE WITH CLEANING EVENTS FILE")
        
#         with open("file_clean_report.txt", "a") as f :
#             f.write(f"CLEANED FILE clean_{splicing_events_file}.csv HAS {c_records} RECORDS\n")
        
#         print(f"PLEASE CHECK file_clean_report.txt FILE FOR CLEANING PROCESSES AND clean_{splicing_events_file}.csv FOR RESULTING CLEAN FILE")
        
#         # Also add this information to Summart_stats.txt
#         with open("Summary_stats.txt", "a") as f :
#             f.write("******** START CLEANING REPORT ********\n")
#             f.write(f"ORIGINAL FILE HAS:  {t_recrds} RECORDS\n")
#             f.write(f"AFTER CLEANING, (REMOVING DUPLICATES (EVENTS WITH SAME CHR#, START AND END) AND STAR > END OR START==END EVENTS), FOUND {dups} events (Please see file_clean_report.txt file)\n")
#             f.write(f"Finally CLEANED FILE clean_{splicing_events_file}.csv HAS {c_records} RECORDS for further processing\n")
#             f.write("******** END CLEANING REPORT ********\n")
        
# # Step 2. Generate intronic_range and ce__all_scan_range_junctions.bed and ce_all_scan_range.bed files
# if arg1==4 or arg1==5 :
#     with open("Summary_stats.txt", "a") as f :
#         f.write("############### NOW STARTING EVENTS TYPE IDENTIFICATION, PLEASE MAKE SURE THAT bed files for all bams are present in bam_beds #########\n")
    
#     if arg1==5 :
#         if os.path.isfile(f"clean_{splicing_events_file}.csv") :
#             inpfile = f"clean_{splicing_events_file}.csv"
        
#         else :
#             print(f"FOUND EMPTY clean_{splicing_events_file}.csv SPLICING EVENTS FILE, PLEAES FIX THE PROBLEM AND RERUN")
#             with open("Summary_stats.txt", "a") as f :
#                 f.write(f"FOUND EMPTY clean_{splicing_events_file}.csv SPLICING EVENTS FILE, PLEAES FIX THE PROBLEM AND RERUN\n")
            
#             print("EXITING NOW")
            
#     else :
#         inpfile = f"clean_{splicing_events_file}.csv"
        
#     # Flag to run this section
#     cryptics_flg = 1
    
#     if cryptics_flg==1 :
#         with open("Summary_stats.txt", "a") as f :
#             f.write("############### NOW STARTING CE IDENTIFICATION #########\n")
        
#         print("############### NOW STARTING CE IDENTIFICATION #########")
        
#         os.makedirs("res_ce_all",exist_ok=True)
#         os.makedirs("event_bedfiles",exist_ok=True)
        
#         if os.access("events_to_tx_mapping_valid.csv",os.F_OK):
#             os.remove("events_to_tx_mapping_valid.csv")
        
#         if os.access("events_tx_mapping_invalid.csv",os.F_OK):
#             os.remove("events_tx_mapping_invalid.csv")
        
#         ######## IR SECTION
# 		###### END IR SECTION
        
#         print("CLEARING UP res_ce_all folder from previous run if any. SASHIMI FOLDERS WILL ONLY BE DELETED BEFORE CREATING NEW SASHIMI PLOTS.")
#         with open("Summary_stats.txt", "a") as f :
#             f.write("CLEARING UP res_ce_all folder from previous run if any. SASHIMI FOLDERS WILL ONLY BE DELETED BEFORE CREATING NEW SASHIMI PLOTS.\n")
        
#         if len(os.listdir("res_ce_all/"))!=0 :
#             for f in glob.glob("res_ce_all/*.*") :
#                 os.remove(f)
        
#         destination_path = os.path.join("res_ce_all", "Summary_stats.txt")
#         shutil.move("Summary_stats.txt", destination_path)
        
#         events_bed_create_flg = 1
        
#         if events_bed_create_flg==1 :
#             if len(os.listdir("event_bedfiles/"))!=0 :
#                 for f in glob.glob("event_bedfiles/*.*") :
#                     os.remove(f)
            
#             if os.access("EnsDB_tx_not_found.csv",os.F_OK):
#                 os.remove("EnsDB_tx_not_found.csv")
            
#             with open("res_ce_all/Summary_stats.txt", "a") as f :
#                 f.write("CALLING TxEnsDB103_layeredV6.R to generate bed files\n")

            # r_code = """
            # library(AnnotationHub)
            # library(AnnotationDbi)
            # library(ensembldb)
            # library(DBI)
            # library(RSQLite)

            # ah <- AnnotationHub()
            # edb <- query(ah, c("EnsDb", "Hsapiens", "103"))[[1]]
            
            # genes_df <- genes(edb)

            # conn <- dbConnect(RSQLite::SQLite(), "edb.db")
            # dbWriteTable(conn, "genes", genes_df)
            # dbDisconnect(conn)
            # """

            # ro.r(r_code)



            import rpy2.robjects as ro
            import rpy2.rinterface_lib.callbacks
            from rpy2.robjects.packages import importr

            # Configurer la gestion des messages de R
            rpy2.rinterface_lib.callbacks.logger.setLevel('INFO')

            # Importer les packages R nécessaires
            base = importr('base')
            utils = importr('utils')
            AnnotationHub = importr('AnnotationHub')
            GenomicFeatures = importr('GenomicFeatures')
            Repitools = importr('Repitools')
            ensembldb = importr('ensembldb')

            # Exécuter le code R pour obtenir l'objet 'edb'
            ro.r('''
            library(AnnotationHub)
            library(GenomicFeatures)
            library(Repitools)
            library(ensembldb)
            ah <- AnnotationHub()
            edb <- query(ah, c("EnsDb", "Hsapiens", "103"))[[1]]
            ''')

            # Importer l'objet 'edb' depuis R
            edb = ro.globalenv['edb']

            # Vous pouvez maintenant utiliser l'objet 'edb' en Python

            # Exemple: afficher les détails de l'objet 'edb'
            print(edb)




            print("CALLING TxEnsDB103_layeredV6.R to generate bed files")
            
            command = [
                "python",
                "TxEnsDB103_layeredV6.py",
                "sorted_selected_events.csv",
                "principal_txs.csv",
                "temp_all_events_sashimi/FINAL_STATS_ALL_SASHIMIS.txt"
            ]
            
            result = subprocess.run(command, capture_output=True, text=True)
            print(result.stdout)
            print(result.stderr)

            # command = [
            #     "R",
            #     "txens.R",
            #     "sorted_selected_events.csv",
            #     "principal_txs.csv",
            #     "temp_all_events_sashimi/FINAL_STATS_ALL_SASHIMIS.txt"
            # ]
            
            # result = subprocess.run(command, capture_output=True, text=True)
            # print(result.stdout)
            # print(result.stderr)

#             with open('votre_fichier.csv', 'r') as fichier:
#                 all_csv_data = [ligne.strip() for ligne in fichier] # For saving gene_ids as well to generate sashimi compatible csv file

#             event_i = 0
#             samples = glob.glob('event_bedfiles/temp_*.bed')
#             overlap_allowed = 5 # Overlap allowed (intron/exon) for ce events

#             for sample in samples :
#                 # Read line from csv file
#                 line_csv = all_csv_data[event_i]
#                 event_i = event_i + 1

#                 # Get  all exons bed
#                 allexons = sample.split('/')[1].split('_')[1]

#                 # Also get 5UTR and 3UTR bed files
#                 gene_name1 = allexons.split('.')[0]
#                 gene_name = gene_name1.split('-')[0]

#                 # First sort the bed
#                 bed = pybedtools.BedTool(f"event_bedfiles/{allexons}")
#                 sorted_bed = bed.sort()
#                 sorted_bed.saveas(f"event_bedfiles/t{allexons}")

#                 # Also read Tx Files to retrieve selected Tx - should find better ways
#                 with open(f"event_bedfiles/TxID{allexons}", 'r') as fichier:
#                     premiere_ligne = fichier.readline().strip()

#                 colonnes = premiere_ligne.split()
#                 TxID = colonnes[6]

#                 df = pd.read_csv(sample, sep='\\s+', header=None)
#                 strnd = df[5]

#                 print(sample)

#                 # Get distance to downstream exon (for ties, report first) from current reference and pick start, end and d
#                 a = pybedtools.BedTool(sample)
#                 b = pybedtools.BedTool(f'event_bedfiles/t{allexons}')
                
#                 closest = a.closest(b, s=True, D="a", iu=True, d=True, t="first")
#                 ds = closest.to_dataframe(names = [
#                     "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
#                     "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
#                 ])

#                 # Also get distance to upstream exon from current reference and pick start, end and d
#                 closest = a.closest(b, s=True, D="a", id=True, d=True, t="first")
#                 us = closest.to_dataframe(names = [
#                     "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
#                     "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
#                 ])

#                 # Get up/dn exon lengths
#                 upexonl = us.iloc[:,9].tolist()[0]
#                 dnexonl = ds.iloc[:,9].tolist()[0]

#                 # Get up and down stream exon numbers
#                 upexon = us.iloc[:,10].tolist()[0]
#                 dnexon = ds.iloc[:,10].tolist()[0]

#                 # Get overlap with up/dn exon - 0 means complete overlap which is assumed as exon skip event
#                 dsovlp1 = ds.iloc[:,12].tolist()[0]
#                 usovlp1 = us.iloc[:,12].tolist()[0]

#                 # Take absolute values
#                 dsovlp = abs(dsovlp1)
#                 usovlp = abs(usovlp1)

#                 diff_exon = int(upexon) - int(dnexon)

#                 # Take absolute value
#                 diff_exon_abs = abs(diff_exon)

#                 if (upexon==dnexon and upexonl!=dnexonl) : # condition [ $upexonl -ne $dnexonl ] REMOVES EXON_SKIP events
#                     with open("res_ce_all/ex1_ex2_ce.txt", "a") as f :
#                         f.write(f"{gene_name} is a CE event")
#                         f.write(f"ds {ds}")
#                         f.write(f"us {us}")
                    
#                     # Get distance to downstream exon (for ties, report first) from current reference and pick start, end and d
#                     closest = a.closest(b, s=True, D="a", iu=True, d=True, t="first", io=True)
#                     dsn = closest.to_dataframe(names = [
#                         "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
#                         "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
#                     ])

#                     # Also get ds overlap
#                     closest = a.closest(b, s=True, D="a", iu=True, d=True, t="first")
#                     dso = closest.to_dataframe(names = [
#                         "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
#                         "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
#                     ])

#                     # Get start and end of dso
#                     dso1 = dso.iloc[:,1].tolist()[0]
#                     dso2 = dso.iloc[:,8].tolist()[0]
#                     dsodiff = int(dso2) - int(dso1)

#                     # Also get distance to upstream exon from current reference and pick start, end and d
#                     closest = a.closest(b, s=True, D="a", id=True, d=True, t="last", io=True)
#                     usn = closest.to_dataframe(names = [
#                         "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
#                         "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
#                     ])

#                     # Also get us overlap
#                     closest = a.closest(b, s=True, D="a", id=True, d=True, t="last")
#                     uso = closest.to_dataframe(names = [
#                         "chr_a", "start_a", "end_a", "name_a", "row_line_a", "strand_a",
#                         "chr_b", "start_b", "end_b", "name_b", "row_line_b", "strand_b", "distance"
#                     ])

#                     # Get up/dn exon lengths
#                     usnexonl = usn.iloc[:,9].tolist()[0]
#                     dsnexonl = dsn.iloc[:,9].tolist()[0]

#                     with open("res_ce_all/ex1_ex2_ce.txt", "a") as f :
#                         f.write(f"now dsn {dsn}")
#                         f.write(f"now usn {usn}")
#                         f.write(f"dsodiff {dsodiff}")
                    
#                     if (usnexonl=="." or dsnexonl==".") :
#                         with open("res_ce_all/problematic_junctions.txt", "a") as f :
#                             f.write(f"{gene_name} is ce with up exon length {usnexonl} and dsnexonl {dsnexonl} , LEAVING it for now PLEASE MAKE SURE IT IS PROPERLY IDENTIFIED")
                        
#                         with open(sample, 'r') as infile, open("res_ce_all/IGV_problematic_junctionscsv", 'a', newline='') as outfile:
#                             writer = csv.writer(outfile)
                            
#                             for line in infile:
#                                 fields = line.strip().split()
#                                 combined_field = f"{fields[0]}:{fields[1]}-{fields[2]}"
#                                 output_fields = [combined_field] + fields + [gene_name]
#                                 writer.writerow(output_fields)
                            
#                     else : 
#                         if strnd=='-' :
#                             if dsodiff==1 :
#                                 if upexonl<60 :
#                                     comb = us.split('\t')                                
#                                     upex = f"{comb[0]}\t{comb[7]+1}\t{comb[8]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                
#                                 else :
#                                     comb = us.split('\t')                                
#                                     upex = f"{comb[0]}\t{comb[7]+1}\t{comb[7]+60+1}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                    
#                                 comb = us.split('\t')                                
#                                 ce = f"{comb[0]}\t{comb[1]+2}\t{comb[7]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                
#                                 if dnexonl<60 :
#                                     comb = ds.split('\t')                                
#                                     dsex = f"{comb[0]}\t{comb[1]-60+1}\t{comb[1]+1}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                    
#                                 else :
#                                     comb = ds.split('\t')                                
#                                     dsex = f"{comb[0]}\t{comb[1]-60+1}\t{comb[1]+1}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                    
#                             # Also get data for start and stop codons - only take -/+ 60 not the whole exon
#                             else :
#                                 if upexonl<60 :
#                                     comb = us.split('\t')                                
#                                     upex = f"{comb[0]}\t{comb[2]}\t{comb[2]+60}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                
#                                 else :
#                                     comb = us.split('\t')                                
#                                     upex = f"{comb[0]}\t{comb[2]}\t{comb[3]+60}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                    
#                                 comb = ds.split('\t')                                
#                                 ce = f"{comb[0]}\t{comb[8]+1}\t{comb[2]-1}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                
#                                 if dnexonl<60 :
#                                     comb = ds.split('\t')                                
#                                     dsex = f"{comb[0]}\t{comb[7]}\t{comb[8]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                    
#                                 else :
#                                     comb = ds.split('\t')                                
#                                     dsex = f"{comb[0]}\t{comb[8]-60}\t{comb[8]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                    
#                                 # Also get data for start and stop codons - only take -/+ 60 not the whole exon
#                                 junc1 = ds.iloc[:,1].tolist()[0]
#                                 junc2 = ds.iloc[:,2].tolist()[0]
#                                 upex1 = us.iloc[:,7].tolist()[0]
#                                 dsex2 = ds.iloc[:,8].tolist()[0]
                                
#                                 # This works better
#                                 if usd_abs<dsd_abs :
#                                     comb = ds.split('\t')
#                                     ce_scan1 = f"{comb[0]}\t{junc1}\t{comb[8]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                    
#                                 else :
#                                     comb = us.split('\t')
#                                     ce_scan1 = f"{comb[0]}\t{junc2}\t{comb[7]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                    
#                         else :
#                             if upexonl<60 :
#                                 comb = us.split('\t')                                
#                                 upex = f"{comb[0]}\t{comb[7]}\t{comb[8]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                            
#                             else :
#                                 comb = us.split('\t')                                
#                                 upex = f"{comb[0]}\t{comb[8]-60}\t{comb[8]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                
#                             use = us.iloc[:,8].tolist()[0]
#                             comb = dsn.split('\t')
#                             ce1 = f"{comb[0]}\t{use}\t{comb[7]-1}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t"
#                             comb1 = ce1.split('\t')                                
#                             ce = f"{comb[0]}\t{comb[1]}\t{comb[2]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                            
#                             if dsnexonl<60 :
#                                 comb = dsn.split('\t')                                
#                                 dsex = f"{comb[0]}\t{comb[7]-1}\t{comb[8]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                
#                             else :
#                                 comb = dsn.split('\t')                                
#                                 dsex = f"{comb[0]}\t{comb[7]-1}\t{comb[7]+60-1}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                
#                             # CE scan region
#                             junc1 = ds.iloc[:,1].tolist()[0]
#                             junc2 = ds.iloc[:,2].tolist()[0]
                            
#                             # Gave problem with single bp differences
#                             upex2 = us.iloc[:,8].tolist()[0]
#                             dsex1 = ds.iloc[:,7].tolist()[0]
                            
#                             # This works better
#                             if usd_abs<dsd_abs :
#                                 comb = ds.split('\t')
#                                 ce_scan1 = f"{comb[0]}\t{junc2}\t{comb[7]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                                
#                             else :
#                                 comb = us.split('\t')
#                                 ce_scan1 = f"{comb[0]}\t{junc1}\t{comb[8]}\t{comb[9]}\t{comb[4]}\t{comb[5]}\t{gene_name}"
                            
#                     with open(sample, 'r') as infile, open('res_ce_all/IGV_ce_inclusion.csv', 'a', newline='') as outfile:
#                         csv_writer = csv.writer(outfile)
                        
#                         for line in infile:
#                             fields = line.strip().split()
#                             first_column = f"{fields[0]}:{fields[1]}-{fields[2]}"
#                             row = [first_column, fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], gene_name, TxID]
#                             csv_writer.writerow(row)
                            
#                     # NOW REMOVE ANY EVENT FOR WHICH START>END
#                     upex_st = upex.iloc[:,1].tolist()[0]
#                     upex_en = upex.iloc[:,2].tolist()[0]
                    
#                     ce_st = ce.iloc[:,1].tolist()[0]
#                     ce_en = ce.iloc[:,2].tolist()[0]
                    
#                     dsex_st = dsex.iloc[:,1].tolist()[0]
#                     dsex_en = dsex.iloc[:,2].tolist()[0]
                    
#                     if (upex_st<upex_en and ce_st<ce_en and dsex_st<dsex_en) :
#                         # Finally paste three segments to a bed file
#                         # First remove if any such file exists
#                         os.remove(f"{gene_name}_nt.bed")
                        
#                         with open(f"{gene_name}_nt.bed", "a") as f :
#                             f.write(f"{upex}\n")
#                             f.write(f"{ce}\n")
#                             f.write(f"{dsex}\n")
                        
#                         os.remove(f"{gene_name}_nt.bed")
                        
#                     # Also coe ce boundary range
#                     with open("res_ce_all/ce_all_scan_range_junctions.bed", "a") as f :
#                         f.write(f"{upex}")
#                         f.write(f"{ce_scan1}")
#                         f.write(f"{dsex}")
                    
#                     with open("res_ce_all/ce_all_scan_range.bed", "a") as f :
#                         f.write(f"{ce_scan1}")
                    
#                 elif (diff_exon_abs==1 and dsovlp==5 and usovlp==5) : # Deals with exon_joining events few bp of
#                     with open(sample, 'r') as infile, open('res_ce_all/Exon_skip.csv', 'a', newline='') as outfile:
#                         csv_writer = csv.writer(outfile)
                        
#                         for line in infile:
#                             fields = line.strip().split()
#                             first_column = f"{fields[0]}:{fields[1]}-{fields[2]}"
#                             row = [first_column, fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], gene_name]
#                             csv_writer.writerow(row)
                            
#                     comb = line_csv.split(',')
#                     first = [comb[0], comb[1], comb[2], comb[3], comb[4], comb[5], comb[6]]
#                     combi = first.split(',')
#                     with open('res_ce_all/Annotated_junctions.csv', 'a', newline='') as outfile:
#                         csv_writer = csv.writer(outfile)
#                         csv_writer.writerow([combi[0], combi[1], combi[2], combi[3], combi[4], combi[5], combi[6]])
                    
#                 elif diff_exon_abs==1 :
#                     # THESE ARE CE INTRONIC events
#                     # CE event with difference of one exon
                    
#                     # Get start and end of dso
#                     usd = us.iloc[:,12].tolist()[0]
#                     usd_abs = abs(usd)
#                     dsd = ds.iloc[:,12].tolist()[0]
#                     dsd_abs = abs(dsd)
                    
#                     # Record 1 (downstream exon). chr# end_ce+d       end_ce+d_nuc_sz      strand  gene_name
#                     # ADDED TO REMOVE ALL EXONS THAT LIE INSIDE AN INTRON
#                     if (usd_abs<=5 or dsd_abs<=5) :
#                         if strnd=='+' :
#                             if upexonl<60 :
#                                 comb = us.split('\t')
#                                 upex = [comb[0], comb[7], comb[8], comb[9], comb[4], comb[5], gene_name]
                                
#                             else :
#                                 comb = us.split('\t')
#                                 upex = [comb[0], comb[8]-60, comb[8], comb[9], comb[4], comb[5], gene_name]
                                
#                             use = us.iloc[:,8].tolist()[0]
#                             comb = ds.split('\t')
#                             ce1 = [comb[0], use, comb[7]-1, comb[9], comb[4], comb[5]]
#                             comb1 = ce1.split('\t')
#                             ce = [comb1[0], comb1[1], comb1[2], comb1[9], comb1[4], comb1[5], gene_name]
                            
#                             if dnexonl<60 :
#                                 comb = ds.split('\t')
#                                 dsex = [comb[0], comb[7]-1, comb[8], comb[9], comb[4], comb[5], gene_name]
                                
#                             else :
#                                 comb = ds.split('\t')
#                                 dsex = [comb[0], comb[8]-1, comb[8]+60-1, comb[9], comb[4], comb[5], gene_name]
                            
#                             # CE scan region
#                             junc1 = ds.iloc[:,1].tolist()[0]
#                             junc2 = ds.iloc[:,2].tolist()[0]
                            
#                             # Gave problem with single bp differences
#                             upex2 = us.iloc[:,8].tolist()[0]
#                             dsex1 = ds.iloc[:,7].tolist()[0]
                            
#                             # This works better
#                             if (usd_abs < dsd_abs) :
#                                 comb = ds.split('\t')
#                                 ce_scan1 = [comb[0], junc2, comb[7], comb[9], comb[4], comb[5], gene_name]
                            
#                             else :
#                                 comb = us.split('\t')
#                                 ce_scan1 = [comb[0], junc1, comb[8], comb[9], comb[4], comb[5], gene_name]
                            
#                         if strnd=='-' :
#                             if upexonl<60 :
#                                 comb = us.split('\t')
#                                 upex = [comb[0], comb[7]-1, comb[8], comb[9], comb[4], comb[5], gene_name]
                                
#                             else :
#                                 comb = us.split('\t')
#                                 upex = [comb[0], comb[7]-1, comb[7]+60-1, comb[9], comb[4], comb[5], gene_name]
                            
#                             use = ds.iloc[:,8].tolist()[0]
#                             comb = ds.split('\t')
#                             ce1 = [comb[0], use, comb[7]-1, comb[9], comb[4], comb[5]]
#                             comb1 = ce1.split('\t')
#                             ce = [comb1[0], comb1[1], comb1[2], comb1[9], comb1[4], comb1[5], gene_name]
                            
#                             if dnexonl<60 :
#                                 comb = ds.split('\t')
#                                 dsex = [comb[0], comb[7]-1, comb[8], comb[9], comb[4], comb[5], gene_name]
                                
#                             else :
#                                 comb = ds.split('\t')
#                                 dsex = [comb[0], comb[8]-60, comb[8], comb[9], comb[4], comb[5], gene_name]
                            
#                             # CE scan region
#                             junc1 = ds.iloc[:,1].tolist()[0]
#                             junc2 = ds.iloc[:,2].tolist()[0]
#                             upex1 = us.iloc[:,7].tolist()[0]
#                             dsex2 = ds.iloc[:,8].tolist()[0]
                            
#                             # This works better
#                             if (usd_abs < dsd_abs) :
#                                 comb = ds.split('\t')
#                                 ce_scan1 = [comb[0], junc1, comb[8], comb[9], comb[4], comb[5], gene_name]
                            
#                             else :
#                                 comb = us.split('\t')
#                                 ce_scan1 = [comb[0], junc2, comb[7], comb[9], comb[4], comb[5], gene_name]
                            
#                         with open(sample, 'r') as infile, open('res_ce_all/IGV_ce_inclusion.csv', 'a', newline='') as outfile:
#                             csv_writer = csv.writer(outfile)
                            
#                             for line in infile:
#                                 fields = line.strip().split()
#                                 first_column = f"{fields[0]}:{fields[1]}-{fields[2]}"
#                                 row = [first_column, fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], gene_name, TxID]
#                                 csv_writer.writerow(row)
                        
#                         # NOW REMOVE ANY EVENT FOR WHICH START>END
#                         upex_st = upex.iloc[:,1].tolist()[0]
#                         upex_en = upex.iloc[:,2].tolist()[0]
                        
#                         ce_st = ce.iloc[:,1].tolist()[0]
#                         ce_en = ce.iloc[:,2].tolist()[0]
                        
#                         dsex_st = dsex.iloc[:,1].tolist()[0]
#                         dsex_en = dsex.iloc[:,2].tolist()[0]
                        
#                         # Also coe ce boundary range
#                         with open("res_ce_all/ce_all_scan_range_junctions.bed", "a") as f :
#                             f.write(f"{upex}\n")
#                             f.write(f"{ce_scan1}\n")
#                             f.write(f"{dsex}\n")
                        
#                         with open("res_ce_all/ce_all_scan_range.bed", "a") as f :
#                             f.write(f"{ce_scan1}")
                        
#                     else :
#                         # BOTH ENDS OF THESE EVENTS LIE INSIDE INTRON
#                         with open("res_ce_all/non_ce_events.txt", "a") as f :
#                             f.write(f"{gene_name} is unknown event type for gene {gene_name}")
                        
#                         comb = line_csv.split(',')
#                         add = [comb[0], comb[1], comb[2], comb[3], comb[4], comb[5], comb[6]].split(',')
                        
#                         with open(f"res_ce_all/Annotated_junctions.csv", "a", newline='') as f:
#                             writer = csv.writer(f)
#                             writer.writerow([add[0], add[1], add[2], add[3], add[4], add[5], add[6]])
                        
#                 else :
#                     with open("res_ce_all/non_ce_events.txt", "a") as f :
#                         f.write(f"{gene_name} is unknown event type for gene {gene_name}\n")
                    
#                     comb = line_csv.split(',')
#                     add = [comb[0], comb[1], comb[2], comb[3], comb[4], comb[5], comb[6]].split(',')
                    
#                     with open(f"res_ce_all/Annotated_junctions.csv", "a", newline='') as f:
#                         writer = csv.writer(f)
#                         writer.writerow([add[0], add[1], add[2], add[3], add[4], add[5], add[6]])
                    
#             ########## THIS SECTION WRITES NEW FILE res_ce_all/ce_all_scan_intron.bed to scan whole intron between two exons for probable ce events
#             with open("res_ce_all/ce_all_scan_range_junctions.bed", 'r') as file:
#                 nrecrds = sum(1 for line in file)
            
#             with open("res_ce_all/ce_all_scan_range_junctions.bed", 'r') as file:
#                 all_data = [line.strip() for line in file]
            
#             i = 0
#             while i<nrecrds :
#                 line11 = all_data[i]
#                 line12 = all_data[i+1]
#                 line13 = all_data[i+2]
#                 i = i + 3
#                 strand = line11.split('\t')[5]
                
#                 if strand=='+' :
#                     comb = [line11, line13].split('\t')
#                     with open("res_ce_all/ce_all_scan_intron.bed", "a") as f :
#                         f.write([comb[0], comb[2], comb[8], comb[10], comb[11], comb[12], comb[13]])
                    
#                 elif strand=='-' :
#                     comb = [line13, line11].split('\t')
#                     with open("res_ce_all/ce_all_scan_intron.bed", "a") as f :
#                         f.write([comb[0], comb[2], comb[8], comb[10], comb[11], comb[12], comb[13]])
                
#             ################## END SECTION WRITES NEW FILE res_ce_all/ce_all_scan_intron.bed
#             # Also copy for coverage calculations and avoid intronic_range calculations
#             with open("res_ce_all/ce_all_scan_range.bed", 'r') as infile, open("res_ce_all/ce_all_scan_unique_range.bed", 'a', newline='') as outfile:
#                 csv_writer = csv.writer(outfile)
#                 for line in infile:
#                     csv_writer.writerow(line)
        
# # CRYPTICS SECTION ENDS HERE

# # EVENT COVERAGE CALCULATIONS STARTS HERE
# # Step 3. Now calculate coverages for each ce-range from all TDP samples and aggregate them for each junction
# # Flag for step 3
# # FILES NEEDED FOR THIS STEP: res_ce_all/ce_all_scan_unique_range.bed, bam_files/*.bam, bam_beds/*.bed for all samples in bam_files folder
# coverages_bed_files_create_flg = 0

# if arg1==4 or arg1==5 :
#     print("NOW STARTED COVERAGE CALCULATIONS, TAKES LONG TIME OFTEN HOURS/DAYS FOR > 200 EVENTS")
    
#     with open("res_ce_all/Summary_stats.txt", "a") as f :
#         f.write("NOW STARTED COVERAGE CALCULATIONS, TAKES LONG TIME OFTEN HOURS/DAYS FOR > 200 EVENTS\n")
    
#     bamsamples_tdp = []
#     with open("all_bams.tsv", 'r') as file:
#         for line in file:
#             columns = line.strip().split('\t')
#             if len(columns) > 2 and re.search(r'TDP43', columns[2]) :
#                 filename = columns[1]
#                 basename = filename.split('/')[-1]
#                 name_without_ext = basename.split('.')[0]
#                 bamsamples_tdp.append(name_without_ext)
                
#                 if bamsamples_tdp :
#                     os.makedirs("coverages",exist_ok=True)
#                     i = 1 # Counter to use for calculations
#                     last_genename = ""
                    
#                     with open('res_ce_all/ce_all_scan_intron.bed', 'r') as file : # col2 col3 col4 col5 col6 col7
#                         for line in file:
#                             line = line.rstrip('\n')
                            
#                             # Now do calculations fot this ce across all replicates
#                             # Invert the range if it is not already
#                             genename = line[6]
                            
#                             if last_genename==genename :
#                                 i = i + 1
#                                 last_genename = genename
                                
#                             else :
#                                 i = 1
#                                 last_genename = genename
                            
#                             #########
#                             chrm_n = line[0]
#                             chrm_start = line[1]
#                             chrm_end = line[2]
#                             fnm_tdp = f"{chrm_n}_{chrm_start}_{chrm_end}_{genename}.cov.bed"
                            
#                             print(f"fnm_tdp {fnm_tdp}")
#                             if chrm_start<chrm_end :
#                                 comb = line.split('\t')
#                                 with open("temp_coord.bed", "w") as f :
#                                     f.write([comb[0], comb[1], comb[2], comb[3], comb[4], comb[5], comb[6]])
                                
#                             else :
#                                 comb = line.split('\t')
#                                 with open("temp_coord.bed", "w") as f :
#                                     f.write([comb[0], comb[2], comb[1], comb[3], comb[4], comb[5], comb[6]])
                                
#                             # And now for TDP43
#                             print("coverage calculations for tdp sampels for:")
#                             with open("temp_coord.bed", 'r') as file :
#                                 for line in file:
#                                     print(line, end='')
                            
#                             for bamsamp in bamsamples_tdp :
#                                 a = pybedtools.BedTool("temp_coord.bed")
#                                 b = pybedtools.BedTool(f"bam_beds/{bamsamp}-sorted.bam.bed")
#                                 result = a.coverage(b, g=f"bam_beds/{bamsamp}.chromosomes.txt", d=True, sorted=True)
#                                 result.saveas(f"temp_coord_{bamsamp}_tdp.bed.cov")
                            
#                             # Finally sum all coverages across replicates for this ce
#                             samples1 = ""
                            
#                             for bamsamp1 in bamsamples_tdp :
#                                 samples1 += f" temp_coord_{bamsamp1}_tdp.bed.cov"
                            
#                             print(f"got all coverage files for tdp {samples1}")
#                             samples2 = samples1.split()
                            
#                             dfs = [pd.read_csv(file, sep='\t', header=None) for file in samples2]
#                             combined_df = pd.concat(dfs, axis=1)
#                             numFiles = len(samples2)
#                             results = []
                            
#                             for index, row in combined_df.iterrows() :
#                                 row_data = row.iloc[:len(row) // numFiles]
#                                 sum_value = row[len(row) // numFiles:].sum()
                                
#                                 result_line = '\t'.join(row_data.astype(str)) + '\t' + str(sum_value)
#                                 results.append(result_line)
                            
#                             with open(f"coverages/{fnm_tdp}", 'w') as f:
#                                 f.write('\n'.join(results))
                            
#                             # Also remove temp_coord_ for each sample done
#                             for bamsamp in bamsamples_tdp :
#                                 if os.path.exists(f"temp_coord_{bamsamp}_tdp.bed.cov"):
#                                     os.remove(f"temp_coord_{bamsamp}_tdp.bed.cov")
                                
#                 else :
#                     print("I did not get any BAM files from all_bams.tsv (from column 2), Please make sure that col3 contains string TDP43 for KD samples, ESITING")
#                     with open("res_ce_all/Summary_stats.txt", "a") as f :
#                         f.write("I did not get any BAM files from all_bams.tsv (from column 2), Please make sure that col3 contains string TDP43 for KD samples, ESITING\n")
                    
#                     sys.exit(1)
                





# # CE IDENTIFICATION STARTS HERE

# # Step 4. using majiq_coverages_automateV1.R script identify final ce_extension_coord.bed and ce_extension_coord.bed files
# fact = 3

# with open(inpfile, 'r') as file:
#     total_splicing_events = sum(1 for line in file)

# unknow_events = 0

# if os.path.exists("res_ce_all/Annotated_junctions.csv") :
#     with open("res_ce_all/Annotated_junctions.csv", 'r') as file:
#         unknow_events = sum(1 for line in file)
    
# not_found = 0

# if os.path.exists("EnsDB_tx_not_found.csv") :
#     with open("EnsDB_tx_not_found.csv", 'r') as file:
#         not_found = sum(1 for line in file)

# problematic = 0

# if os.path.exists("res_ce_all/IGV_problematic_junctions.csv") :
#     with open("res_ce_all/IGV_problematic_junctions.csv", 'r') as file:
#         problematic = sum(1 for line in file)
        
# unseccesful_r1 = 0

# if os.path.exists("res_ce_all/skipped_ce.csv") :
#     with open("res_ce_all/skipped_ce.csv", 'r') as file:
#         unseccesful_r1 = sum(1 for line in file)

# unseccesful_r = unseccesful_r1/3

# with open("res_ce_all/Summary_stats.txt", "a") as f :
#     f.write(f"Total splicing events read are: {total_splicing_events}\n")
#     f.write(f"Out of these {total_splicing_events} total events\n")
#     f.write(f"Events not found in EnsDB are: {not_found}, please see EnsDB_tx_not_found.csv file\n")
#     f.write(f"Events that are not CE: {unknow_events} , please see res_ce_all/non_ce_events.txt and res_ce_all/Annotated_junctions.csv file\n")
#     f.write(f"Events that are somewhat problematic: {problematic} , please see res_ce_all/problematic_junctions.txt and res_ce_all/IGV_problematic_junctions.csv files\n")

# ce_boundary_events = total_splicing_events - not_found - unknow_events - problematic

# with open("res_ce_all/Summary_stats.txt", "a") as f :
#     f.write("#######################################\n")






# if arg1==4 or arg1==5 :
#     with open("res_ce_all/Summary_stats.txt", "a") as f :
#         f.write(f"Now starting ce_boundary calculations for a remaining total of: {ce_boundary_events} events, BY INVOKING Auto_CoverV4_layered_intronV3.R SCRIPT\n")
    
#     print("IDENTIFYING CE BOUNDARIES BY CALLING Auto_CoverV4_layered_intronV3.R script")
    
#     command = [
#         "python",
#         "Auto_CoverV4_layered_intronV3.py",
#         "res_ce_all/IGV_ce_inclusion.csv",
#         "res_ce_all/ce_all_scan_range_junctions.bed",
#         ".6",
#         "res_ce_all/ce_all_scan_intron.bed"
#     ]

#     result = subprocess.run(command, capture_output=True, text=True)
#     print(result.stdout)
    
#     print("BACK FROM CE_BOUNDARY CALCULATIONS")
    
#     with open("res_ce_all/Summary_stats.txt", "a") as f :
#         f.write("################################################\n")
    
#     # Step 4 ENDS
#     # Step 5. Finally identify nt and AA sequences
#     # Get nt and aa for ce boundaries
    
#     # First remove duplicates from ce' boundary coordinates
# 	# NEW CODE STARTS HERE
#     nrecrds = 0
#     nrecrdst = 0
    
#     if os.access("res_ce_all/ce_inclusion_coord.bed",os.F_OK) :
#         with open("res_ce_all/ce_inclusion_coord.bed", 'r') as file:
#             nrecrds = sum(1 for line in file)
#         nrecrdst = nrecrds/3
    
#     if nrecrds!=0 :
#         with open("res_ce_all/Summary_stats.txt", "a") as f :
#             f.write(f"Back From R Session, Now checking FOR DUPLICATES for a total of ce_inclusion events: {nrecrdst}\n")
        
#         with open("res_ce_all/ce_inclusion_coord.bed", 'r') as file:
#             all_data = [line.strip() for line in file]
        
#         # Also read ce_inclusion_coord_sashimi.bed to get TxID
#         with open("res_ce_all/ce_inclusion_coord_sashimi.bed", 'r') as file:
#             sashimi_data = [line.strip() for line in file]
        
#         # Also read csv file
#         with open("res_ce_all/IGV_R_returned_ce_inclusion.csv", 'r') as file:
#             csv_data = [line.strip() for line in file]
        
#         csvi = 0
#         i = 0
        
#         while i<nrecrds :
#             line11 = all_data[i]
#             line12 = all_data[i+1]
#             line13 = all_data[i+2]
            
#             # TxID
# 			# For sashimi
#             lines1 = sashimi_data[i]
#             lines2 = sashimi_data[i+1]
#             lines3 = sashimi_data[i+2]
            
#             csv_ln = csv_data[csvi]
#             csvi = csvi + 1
#             i = i + 3
#             j = i
            
#             # Now go through rest of the data
#             flg = 0
            
#             while j<nrecrds :
#                 line21 = all_data[j]
#                 line22 = all_data[j+1]
#                 line23 = all_data[j+2]
                
#                 if (line11==line21 and line12==line22 and line13==line23) :
#                     flg = 1
                
#                 j = j + 3
                
#             if flg==0 :
#                 comb11 = line11.split('\t')
#                 comb12 = line12.split('\t')
#                 comb13 = line13.split('\t')
                
#                 with open("res_ce_all/ce_inclusion_coord_uniq.bed", "a") as f :
#                     f.write(f"{comb11[0]}, {comb11[1]}, {comb11[2]}, {comb11[3]}, {comb11[4]}, {comb11[5]}, {comb11[6]}\n")
#                     f.write(f"{comb12[0]}, {comb12[1]}, {comb12[2]}, {comb12[3]}, {comb12[4]}, {comb12[5]}, {comb12[6]}\n")
#                     f.write(f"{comb13[0]}, {comb13[1]}, {comb13[2]}, {comb13[3]}, {comb13[4]}, {comb13[5]}, {comb13[6]}\n")
                    
#                 # Also save bed file for sashimi plots
#                 comb1 = lines1.split('\t')
#                 comb2 = lines2.split('\t')
#                 comb3 = lines3.split('\t')
                
#                 with open("res_ce_all/ce_inclusion_coord_uniq_sashimi.bed", "a") as f :
#                     f.write(f"{comb1[0]}, {comb1[1]}, {comb1[2]}, {comb1[3]}, {comb1[4]}, {comb1[5]}, {comb1[6]}, {comb1[7]}\n")
#                     f.write(f"{comb2[0]}, {comb2[1]}, {comb2[2]}, {comb2[3]}, {comb2[4]}, {comb2[5]}, {comb2[6]}, {comb2[7]}\n")
                
#                 # Also SAVING CE coordinates
#                 with open("res_ce_all/ce_inclusion_coord_only.bed", "a") as f :
#                     f.write(f"{comb2[0]}, {comb2[1]}, {comb2[2]}, {comb2[2]-comb2[1]}, {comb2[4]}, {comb2[5]}, {comb2[6]}, {comb2[7]}\n")
                
#                 with open("res_ce_all/ce_inclusion_coord_uniq_sashimi.bed", "a") as f :
#                     f.write(f"{comb3[0]}, {comb3[1]}, {comb3[2]}, {comb3[3]}, {comb3[4]}, {comb3[5]}, {comb3[6]}, {comb3[7]}\n")
                
#                 # Also save csv file
#                 with open("res_ce_all/CE_inclusion.csv", "a", newline='') as f:
#                     writer = csv.writer(f)
#                     writer.writerow(csv_ln)
                
#             else :
#                 with open("res_ce_all/ce_inclusion_coord_repeated.bed", "a") as f :
#                     f.write(f"{comb11[0]}, {comb11[1]}, {comb11[2]}, {comb11[3]}, {comb11[4]}, {comb11[5]}, {comb11[6]}\n")
#                     f.write(f"{comb12[0]}, {comb12[1]}, {comb12[2]}, {comb12[3]}, {comb12[4]}, {comb12[5]}, {comb12[6]}\n")
#                     f.write(f"{comb13[0]}, {comb13[1]}, {comb13[2]}, {comb13[3]}, {comb13[4]}, {comb13[5]}, {comb13[6]}\n")
                
#         with open("res_ce_all/ce_inclusion_coord_uniq.bed", 'r') as file:
#             nrecrds_uniqt = sum(1 for line in file)
        
#         nrecrds_uniq = nrecrds_uniqt/3
#         repeated_ce_boundary_events = nrecrdst - nrecrds_uniq
        
#         with open("res_ce_all/Summary_stats.txt", "a") as f :
#             f.write(f"Total unique ce_inclusion events are: {nrecrds_uniq}, please see res_ce_all/ce_inclusion_coord_uniq.bed and res_ce_all/CE_inclusion.csv\n")
#             f.write(f"Total repeated ce_inclusion events are: {repeated_ce_boundary_events}, please see res_ce_all/ce_inclusion_coord_repeated.bed\n")
        
#         print(f"total repeated ce_inclusion events were {repeated_ce_boundary_events}")
        
#     else :
#         with open("res_ce_all/Summary_stats.txt", "a") as f :
#             f.write(f"Back From R Session, total ce_inclusion events: {nrecrdst}\n")






# # First remove duplicates from IR events
# # NEW CODE STARTS HERE
# nrecrds_ir = 0
# nrecrds_irt = 0

# if os.access("res_ce_all/ce_extension_coord.bed",os.F_OK) :
#     with open("res_ce_all/ce_extension_coord.bed", 'r') as file:
#         nrecrds_ir = sum(1 for line in file)
#     nrecrds_irt = nrecrds_ir/3

# if nrecrds_ir!=0 :
#     with open("res_ce_all/Summary_stats.txt", "a") as f :
#         f.write(f"Total ce_extension events are: {nrecrds_irt}, Now checking for repeated ce_extension events\n")
    
#     with open("res_ce_all/ce_extension_coord.bed", 'r') as file:
#         all_data = [line.strip() for line in file]
    
#     # Also read ce_inclusion_coord_sashimi.bed to get TxID
#     with open("res_ce_all/ce_extension_coord_sashimi.bed", 'r') as file:
#         sashimi_data_ext = [line.strip() for line in file]
    
#     # Also read csv file
#     with open("res_ce_all/IGV_R_returned_ce_extension.csv", 'r') as file:
#         csv_data_ir = [line.strip() for line in file]
    
#     csvi = 0
#     i = 0
    
#     while i<nrecrds_ir :
#         line11 = all_data[i]
#         line12 = all_data[i+1]
#         line13 = all_data[i+2]

#         # For sashimi
#         lines1 = sashimi_data_ext[i]
#         lines2 = sashimi_data_ext[i+1]
#         lines3 = sashimi_data_ext[i+2]
        
#         # TxID
#         txid_recrd = sashimi_data_ext[i]
#         txid = txid_recrd.split('\t')[7]
        
#         csv_ln = csv_data_ir[csvi]
#         csvi = csvi + 1
#         i = i + 3
#         j = i
        
#         # Now go through rest of the data
#         flg = 0
        
#         while j<nrecrds_ir :
#             line21 = all_data[j]
#             line22 = all_data[j+1]
#             line23 = all_data[j+2]
            
#             if (line11==line21 and line12==line22 and line13==line23) :
#                 flg = 1
            
#             j = j + 3
            
#         if flg==0 :
#             comb11 = line11.split('\t')
#             comb12 = line12.split('\t')
#             comb13 = line13.split('\t')
            
#             with open("res_ce_all/ce_extension_coord_uniq.bed", "a") as f :
#                 f.write(f"{comb11[0]}, {comb11[1]}, {comb11[2]}, {comb11[3]}, {comb11[4]}, {comb11[5]}, {comb11[6]}\n")
#                 f.write(f"{comb12[0]}, {comb12[1]}, {comb12[2]}, {comb12[3]}, {comb12[4]}, {comb12[5]}, {comb12[6]}\n")
#                 f.write(f"{comb13[0]}, {comb13[1]}, {comb13[2]}, {comb13[3]}, {comb13[4]}, {comb13[5]}, {comb13[6]}\n")
                
#             # Also save bed file with TxID for sashimi plots
#             comb1 = lines1.split('\t')
#             comb2 = lines2.split('\t')
#             comb3 = lines3.split('\t')
            
#             with open("res_ce_all/ce_extension_coord_uniq_sashimi.bed", "a") as f :
#                 f.write(f"{comb1[0]}, {comb1[1]}, {comb1[2]}, {comb1[3]}, {comb1[4]}, {comb1[5]}, {comb1[6]}, {comb1[7]}\n")
#                 f.write(f"{comb2[0]}, {comb2[1]}, {comb2[2]}, {comb2[3]}, {comb2[4]}, {comb2[5]}, {comb2[6]}, {comb2[7]}\n")
            
#             # Also save ce coordinates only
#             with open("res_ce_all/ce_extension_coord_only.bed", "a") as f :
#                 f.write(f"{comb2[0]}, {comb2[1]}, {comb2[2]}, {comb2[2]-comb2[1]}, {comb2[4]}, {comb2[5]}, {comb2[6]}, {comb2[7]}\n")
            
#             with open("res_ce_all/ce_extension_coord_uniq_sashimi.bed", "a") as f :
#                 f.write(f"{comb3[0]}, {comb3[1]}, {comb3[2]}, {comb3[3]}, {comb3[4]}, {comb3[5]}, {comb3[6]}, {comb3[7]}\n")
            
#             # Also save csv file
#             with open("res_ce_all/CE_extension.csv", "a", newline='') as f:
#                 writer = csv.writer(f)
#                 writer.writerow(csv_ln)
            
#         else :
#             with open("res_ce_all/ce_extension_coord_repeated.bed", "a") as f :
#                 f.write(f"{comb11[0]}, {comb11[1]}, {comb11[2]}, {comb11[3]}, {comb11[4]}, {comb11[5]}, {comb11[6]}\n")
#                 f.write(f"{comb12[0]}, {comb12[1]}, {comb12[2]}, {comb12[3]}, {comb12[4]}, {comb12[5]}, {comb12[6]}\n")
#                 f.write(f"{comb13[0]}, {comb13[1]}, {comb13[2]}, {comb13[3]}, {comb13[4]}, {comb13[5]}, {comb13[6]}\n")
            
#     # Get total unique coordinates
#     with open("res_ce_all/ce_extension_coord_uniq.bed", 'r') as file:
#         nrecrds_uniq_irt = sum(1 for line in file)
    
#     nrecrds_uniq_ir = nrecrds_uniq_irt/3
#     repeated_IR_events = nrecrds_irt - nrecrds_uniq_ir
    
#     with open("res_ce_all/Summary_stats.txt", "a") as f :
#         f.write(f"Total unique ce_extension events are: {nrecrds_uniq_ir}, please see res_ce_all/ce_extension_coord_uniq.bed and res_ce_all/CE_extension.csv files\n")
#         f.write(f"Total repeated ce_extension events are: {repeated_IR_events}, please see res_ce_all/ce_extension_coord_repeated.bed\n")
    
#     print(f"Total repeated ce_extension events were {repeated_IR_events}")
    
#     # Now copy back data to ce_coord
#     # NEW CODE ENDS HERE - IR

# # Step 5 ENDS
# # NOW FOR IR EVENTS
# # First remove duplicates from IR events
# # NEW CODE STARTS HERE
# nrecrds_ir = 0
# nrecrds_irt = 0

# if os.access("res_ce_all/IR_coord.bed",os.F_OK) :
#     with open("res_ce_all/IR_coord.bed", 'r') as file:
#         nrecrds_ir = sum(1 for line in file)
    
# if nrecrds_ir!=0 :
#     with open("res_ce_all/Summary_stats.txt", "a") as f :
#         f.write(f"TOTAL IR EVENTS ARE: {nrecrds_ir}, NOW CHECKING FOR REPEATED IR EVENTS\n")
    
#     print(f"TOTAL IR EVENTS ARE: {nrecrds_ir}, NOW CHECKING FOR REPEATED IR EVENTS")
    
#     with open("res_ce_all/IR_coord.bed", 'r') as file:
#         all_data = [line.strip() for line in file]
    
#     # Also read ce_inclusion_coord_sashimi.bed to get TxID
#     with open("res_ce_all/IR_coord_sashimi.bed", 'r') as file:
#         sashimi_data_ext = [line.strip() for line in file]
    
#     # Also read csv file
#     with open("res_ce_all/IGV_R_returned_IR.csv", 'r') as file:
#         csv_data_ir = [line.strip() for line in file]
    
#     csvi = 0
#     i = 0
    
#     while i<nrecrds_ir :
#         line11 = all_data[i]
        
#         # For sashimi
#         lines1 = sashimi_data_ext[i]
        
# 		# TxID
#         txid_recrd = sashimi_data_ext[i]
#         txid = txid_recrd.split('\t')[7]
        
#         csv_ln = csv_data_ir[csvi]
#         csvi = csvi + 1
#         i = i + 1
#         j = i
        
#         # Now go through rest of the data
#         flg = 0
        
#         while j<nrecrds_ir :
#             line21 = all_data[j]
            
#             if line11==line21 :
#                 flg = 1
            
#             j = j + 1
            
#         if flg==0 :
#             comb11 = line11.split('\t')
#             comb1 = lines1.split('\t')
            
#             with open("res_ce_all/IR_coord_uniq.bed", "a") as f :
#                 f.write(f"{comb11[0]}, {comb11[1]}, {comb11[2]}, {comb11[3]}, {comb11[4]}, {comb11[5]}, {comb11[6]}\n")
            
#             with open("res_ce_all/IR_coord_uniq_sashimi.bed", "a") as f :
#                 f.write(f"{comb1[0]}, {comb1[1]}, {comb1[2]}, {comb1[3]}, {comb1[4]}, {comb1[5]}, {comb1[6]}, {comb1[7]}\n")
            
#             # Also save ce coordinates only
#             with open("res_ce_all/IR_coord_only.bed", "a") as f :
#                 f.write(f"{comb1[0]}, {comb1[1]}, {comb1[2]}, {comb1[2]-comb1[1]}, {comb1[4]}, {comb1[5]}, {comb1[6]}, {comb1[7]}\n")
            
#             # Also save csv file
#             comb = csv_ln.split(',')
#             combi = [comb[1], comb[2], comb[3], comb[6], comb[7], comb[8]].split(',')
            
#             with open("res_ce_all/IR.csv", "a") as f :
#                 f.write("{combi[0]}, {combi[1]}, {combi[2]}, {combi[3]}, {combi[4]}, {combi[5]}")
            
#         else :
#             with open("res_ce_all/IR_coord_repeated.bed", "a") as f :
#                 f.write(f"{comb11[0]}, {comb11[1]}, {comb11[2]}, {comb11[3]}, {comb11[4]}, {comb11[5]}, {comb11[6]}\n")
            
#     # Get total unique coordinates
#     with open("res_ce_all/IR_coord_uniq.bed", 'r') as file:
#         nrecrds_uniq_irt = sum(1 for line in file)
    
#     repeated_IR_events = nrecrds_ir - nrecrds_uniq_irt
    
#     with open("res_ce_all/Summary_stats.txt", "a") as f :
#         f.write(f"TOTAL UNIQUE IR EVENTS ARE: {nrecrds_uniq_irt}, please see res_ce_all/IR_coord_uniq.bed and res_ce_all/IR.csv files\n")
#         f.write(f"Total repeated IR events are: {repeated_IR_events}, please see res_ce_all/IR_coord_repeated.bed\n")
#         f.write("THOSE WERE ALL PERTINENT STAT - Please let us know if something is missing or some more stats can be useful!!!\n")
    
#     print(f"Total repeated IR events were {repeated_IR_events}")
    
#     # Now copy back data to ce_coord
#     # NEW CODE ENDS HERE - IR

# # SASHIMI PLOTS SECTION
# else :
#     print(f"CE_BOUNDARY_EVENTS ARE {ce_boundary_events} SO ABONDONED CE_BOUNDARY CALCULATIONS")

# ###### END IR EVENTS
# # FINAL STATISTICS
# # DELETE ALL sashimi files as we plot sashimi in part-b
# print("FINAL CLEANUP - Removing intermediate/unnecessary files")

# if len(os.listdir("res_skiptics"))!=0 :
#     for f in glob.glob("res_skiptics/*sashimi*") :
#         os.remove(f)
        
# if len(os.listdir("res_ce_all"))!=0 :
#     for f in glob.glob("res_ce_all/*sashimi*") :
#         os.remove(f)
        
# for f in glob.glob("*.bed") :
#     os.remove(f)

# if os.path.exists("all_events_bed_sashimi.tab"):
#     os.remove("all_events_bed_sashimi.tab")

# if os.path.exists("dump.fasta"):
#     os.remove("dump.fasta")

# if os.path.exists("res_ce_all/ce_all_scan_intron.bed"):
#     os.remove("res_ce_all/ce_all_scan_intron.bed")

# if os.path.exists("res_ce_all/ce_all_scan_range_junctions.bed"):
#     os.remove("res_ce_all/ce_all_scan_range_junctions.bed")

# if os.path.exists("res_ce_all/ce_all_scan_range.bed"):
#     os.remove("res_ce_all/ce_all_scan_range.bed")

# if os.path.exists("res_ce_all/ce_all_scan_unique_range.bed"):
#     os.remove("res_ce_all/ce_all_scan_unique_range.bed")

# if os.path.exists("res_ce_all/ce_extension_coord_only.bed"):
#     os.remove("res_ce_all/ce_extension_coord_only.bed")

# if os.path.exists("res_ce_all/ce_extension_coord.bed"):
#     os.remove("res_ce_all/ce_extension_coord.bed")

# if os.path.exists("res_ce_all/ce_inclusion_coord_only.bed"):
#     os.remove("res_ce_all/ce_inclusion_coord_only.bed")

# if os.path.exists("res_ce_all/ce_inclusion_coord_repeated.bed"):
#     os.remove("res_ce_all/ce_inclusion_coord_repeated.bed")

# if os.path.exists("res_ce_all/ce_inclusion_coord_uniq.bed"):
#     os.remove("res_ce_all/ce_inclusion_coord_uniq.bed")

# if os.path.exists("res_ce_all/ce_inclusion_coord.bed"):
#     os.remove("res_ce_all/ce_inclusion_coord.bed")

# if os.path.exists("res_ce_all/IR_coord_only.bed"):
#     os.remove("res_ce_all/IR_coord_only.bed")

# if os.path.exists("res_ce_all/IR_coord_repeated.bed"):
#     os.remove("res_ce_all/IR_coord_repeated.bed")

# if os.path.exists("res_ce_all/IR_coord_uniq.bed"):
#     os.remove("res_ce_all/IR_coord_uniq.bed")

# if os.path.exists("res_ce_all/IR_coord.bed"):
#     os.remove("res_ce_all/IR_coord.bed")

# if os.path.exists("res_ce_all/IGV_R_returned_ce_extension.csv"):
#     os.remove("res_ce_all/IGV_R_returned_ce_extension.csv")

# if os.path.exists("res_ce_all/IGV_R_returned_ce_inclusion.csv"):
#     os.remove("res_ce_all/IGV_R_returned_ce_inclusion.csv")

# if os.path.exists("res_ce_all/IGV_R_returned_IR.csv"):
#     os.remove("res_ce_all/IGV_R_returned_IR.csv")

# if os.path.exists("res_ce_all/IGV_ce_inclusion.csv"):
#     os.remove("res_ce_all/IGV_ce_inclusion.csv")

# print("For final statistics on ce events, please see res_ce_all/Summary_stats.txt file")
# print("ALL DONE - hopefully - Successfully")