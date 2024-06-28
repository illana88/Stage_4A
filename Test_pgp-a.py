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
            print("STARTED PROCESSING SAMPLE $sample")
            with open("Summary_stats.txt","a") as fichier :
                fichier.write("STARTED PROCESSING SAMPLE $sample\n")
            samp = os.path.splitext(os.path.basename(sample))[0]
            commande_shell = f"stringtie {sample} -p 8 -G gencode.v38.annotation.gtf -o iPSC_gtfs/{samp}.gtf"
            subprocess.run(commande_shell, shell=True, check=True) ##### Testé jusqu'ici et tout fonctionne #####
        if len(os.listdir("iPSC_gtfs"))!=0 :
            print("DONE WITH SringTie Calculations - NOW GENERATING LIST of ABUNDANT Txs")
            with open("Summary_stats.txt","a") as fichier :
                fichier.write("DONE WITH SringTie Calculations - NOW GENERATING LIST of ABUNDANT Txs\n")
            # Now get unique transcripts in csv files for each sample
            for i in glob.glob("iPSC_gtfs/*.gtf") :
                if (os.path.isfile(i)) :
                    samp = i.split('/')[1].split('.')[0]
                    print(samp)
                    print("PROCESSING SAMPLE $i")
                    with open("Summary_stats.txt","a") as fichier :
                        fichier.write("PROCESSING SAMPLE $i\n")
                #     # Step 1. Get all Txs having reference_id (ENST), ref_gene_id (ENSG) and ref_gene_name (HUGO SYMBOL), these are 24 column lines in stringtie's gtf file
                #     if len(i)==24 :
                #         col = [i[13], i[15], i[17], i[19], i[21], i[23]]
                #         for c in col :
                # else :
                #     break