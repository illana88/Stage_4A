######## Proteogenomic Pipeline for Biomarker Discovery in iPSC neurons
######### Inline Readme
# 1. Get List of most abundant Transcripts from KD iPSC samples
# This step assumes that all bam files are in current folder and follow have
# change directory to part-a and issue run part-a.sh with appropriate option and input file

import os
import glob
import csv
import subprocess

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
            samples.append(row[1]) # Testé jusqu'ici et tout fonctionne
    # Check if we got any samples
    if "samples" in globals and len(samples)!=0 :
        for sample in samples :
            print("STARTED PROCESSING SAMPLE $sample")
            with open("Summary_stats.txt","a") as fichier :
                fichier.write("STARTED PROCESSING SAMPLE $sample\n")
    #         samp = os.path.splitext(os.path.basename(sample))[0]
    #         commande_shell = f"stringtie {sample} -p 8 -G gencode.v38.annotation.gtf -o iPSC_gtfs/{samp}.gtf"
    #         subprocess.run(commande_shell, shell=True, check=True)
