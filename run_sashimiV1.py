########## run_sashimiV1.sh codé en Pyhton ##########

import os
import glob
import subprocess
import sys

# THIS IS FINAL SCRIPT FOR SASHIMI PLOTS FOR ALL PARTS OF PGP

# THIS SCRIPT CONTAINS SASHIMI PLOT CODE FOR
# 0. PLEASE NOTE THAT FOLLOWING 2 FILES (OR SOFT LINKS) SHOULD BE IN CURRENT FOLDER
    #01: ggsashimi_txV4.py
    #02: Homo_sapiens.GRCh38.103.chr.sorted_new.gtf
# 1. Skiptic Events
# 2. ALL MAJIQ EVENTS
# 3. CE (INCLUDING INCLUSION, EXTENSION AND IR) events

args = sys.argv[1:]

# First CHECK IF Called from pgp-a/b or pgp-c
# CHECK IF 3 ARGUMENTS ARE PROVIDED

if len(args)==3 :

    ######### NEW - Now read input csv and bed files and flag
    inp_csv = args[0]
    inp_bed = args[1]
    
    # Get folder
    inp_prefix = args[0].split('/')[0]
    
    if int(args[2])==1 :
        os.makedirs(f"{inp_prefix}/sashimi_plots",exist_ok=True)
        if len(os.listdir(f"{inp_prefix}/sashimi_plots/"))!=0 :
            for f in glob.glob(f"{inp_prefix}/sashimi_plots/*.*") :
                os.remove(f)
        
        bed = f"{inp_bed}"
        
        with open(bed, 'r') as file:
            all_bed_data = [line.strip() for line in file]
            
        with open(bed, 'r') as file:
            nrecrds = sum(1 for line in file)
            
        print("########### NRECRDS ###########")
        print(nrecrds)
        
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
            line = f"{field1}:{field2}-{field11}"
    
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
    
            result = subprocess.run(command)
            print(result.stdout)
    
        # Now merge all pdf's
        command = [
            "python", "merge_sashimis.py", f"{inp_prefix}/sashimi_plots/"
        ]
    
        result = subprocess.run(command)
        print(result.stdout)

    ###### THIS IS FOR ALL MAJIQ EVENTS
    if int(args[2])==2 :
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
            line1 = all_bed_data[i].split('\t')
            i = i + 1

            # Read strand
            strnd = '+'

            # Also read TxID
            TxID = line1[7]
            line = f"{line1[0]}:{int(line1[1])-50}-{int(line1[2])+50}" # This is the actual event

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
            event_identified = f"None -{line1[1]}-{line1[2]}"

            print(f"processing event num {eventn} and event {fn}")

            command = [
                "python",
                "ggsashimi_txV3.py",
                "-A", "median_j",
                "-b", "all_bams.tsv",
                "-c", line,
                "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
                "-GeneName", gene_name,
                "-MajiqStrnd", strnd,
                "-ORIG", "1",
                "-UEX", str(exon1),
                "-DEX", str(exon2),
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

            result = subprocess.run(command)
            print(result.stdout)

        # Now merge all pdf's
        command = ["python", "merge_sashimis.py", "all_events_sashimi/"]
        result = subprocess.run(command)
        print(result.stdout)
                
#                 # THIS SECTION IS FOR CE_INCLUSION EVENTS
#                 # WILL MERGE INCLUSION AND EXTENSION EVENTS
#                 if int(args[2])==3 :
#                     os.makedirs(f"{inp_prefix}/ce_incl_sashimi_plots",exist_ok=True)
#                     if len(os.listdir(f"{inp_prefix}/ce_incl_sashimi_plots/"))!=0 :
#                         for f in glob.glob(f"{inp_prefix}/ce_incl_sashimi_plots/*.*") :
#                             os.remove(f)
                
#                     bed = inp_bed

#                     with open(bed, 'r') as file:
#                         all_bed_data = [line.strip() for line in file]
                    
#                     with open(bed, 'r') as file:
#                         nrecrds = sum(1 for line in file)
                    
#                     nrecrdst = nrecrds/3

#                     print(f"read {nrecrdst} records")
#                     csv = inp_csv

#                     with open(csv, 'r') as file:
#                         all_csv_data = [line.strip() for line in file]
                    
#                     i = 0
#                     eventn = 0

#                     while i<nrecrds :
#                         # Construct string for ggsashimi
#                         line1 = all_bed_data[i]
#                         line2 = all_bed_data[i+1]
#                         line3 = all_bed_data[i+2]

#                         i = i + 3

#                         # Read strand
#                         strnd = line1.iloc[:, 5]

#                         # Also read TxID
#                         TxID = line1.iloc[7]

#                         if strnd=='+' :
#                             strndflg = "plus"

#                             # us exon length
#                             exon1 = line1.iloc[:, 3]

#                             # ds exon length
#                             exon2 = line3.iloc[:, 3]

#                             # Now modify to reflect whole up/dn exons
#                             up = line1.iloc[:, 1]
#                             upn = up - exon1
# 				             
#                             # ds exon
#                             dn = line3.iloc[:, 2]
#                             dnn = dn + exon2

#                             l123 = f"{line1} {line2} {line3}".split()
#                             fulltitle = f"{l123[0]}-{l123[1]}:{l123[2]}-{l123[9]}:{l123[10]}-{l123[17]}:{l123[18]}"

#                             comb = f"{line1} {upn} {dnn}".split()
#                             line = f"{comb[0]}:{comb[8]}-{comb[9]}"

#                         else :
#                             strndflg = "minus"

#                             # us exon length
#                             exon1 = line3.iloc[:, 3]

#                             # ds exon length
#                             exon2 = line1.iloc[:, 3]

#                             # Now modify to reflect whole up/dn exons
#                             dn = line3.iloc[:, 1]
#                             dnn = dn - exon1
# 				             
#                             # ds exon
#                             up = line1.iloc[:, 2]
#                             upn = up + exon2

#                             l123 = f"{line1} {line2} {line3}".split()
#                             fulltitle = f"{l123[0]}-{l123[17]}:{l123[18]}-{l123[9]}:{l123[10]}-{l123[1]}:{l123[2]}"

#                             comb = f"{line1} {dnn} {upn}".split()
#                             line = f"{comb[0]}:{comb[8]}-{comb[9]}"

#                         event = all_csv_data[eventn]

#                         fields = event.split(',')
#                         fn = f"{fields[7]}-{fields[1]}_{fields[2]}-{fields[3]}"

#                         # String for majiq event
#                         chr_name = event.split(',')[1]
#                         start = event.split(',')[2]
#                         end = event.split(',')[3]
#                         gene_name = event.split(',')[7]
#                         intron = end - start

#                         #gene_name-start-intron-end #for now using strand from ggsashimi
#                         comb_line = f"{chr_name} {start} {end}".split()
#                         majiq_event = f"{comb_line[0]}-{comb_line[1]}-{comb_line[2]}"

#                         # Also get actual event identified
#                         event_identified = f"{line2.iloc[:,0]}-{line2.iloc[:,1]}-{line2.iloc[:,2]}"
#                         eventn = eventn + 1

#                         print(f"processing event num {eventn} and event {fn}")

#                         command = [
#                             "./ggsashimi_txV3.py",
#                             "-A", "median_j",
#                             "-b", "all_bams.tsv",
#                             "-c", line,
#                             "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
#                             "-GeneName", gene_name,
#                             "-MajiqStrnd", strnd,
#                             "-ORIG", "1",
#                             "-UEX", exon1,
#                             "-DEX", exon2,
#                             "-FullTitle", fulltitle,
#                             "-MajiqTx", majiq_event,
#                             "-Majiq", fn,
#                             "-Tx", TxID,
#                             "-M", "1",
#                             "-C", "3",
#                             "-o", f"{inp_prefix}/ce_incl_sashimi_plots/{fn}",
#                             "-O", "3",
#                             "--alpha", "0.25",
#                             "--base-size=20",
#                             "--ann-height=2.5",
#                             "--height=2.5",
#                             "--width=18",
#                             "-P", "palette.txt"
#                         ]

#                         result = subprocess.run(command)
#                         print(result.stdout)

#                     # Now merge all pdf's
#                     command = ["python", "merge_sashimis.py", f"{inp_prefix}/ce_incl_sashimi_plots/"]
#                     result = subprocess.run(command)
#                     print(result.stdout)

#                 # CE_EXTENSION
#                 if int(args[2])==4 :
#                     os.makedirs(f"{inp_prefix}/ce_ext_sashimi_plots",exist_ok=True)
#                     if len(os.listdir(f"{inp_prefix}/ce_ext_sashimi_plots/"))!=0 :
#                         for f in glob.glob(f"{inp_prefix}/ce_ext_sashimi_plots/*.*") :
#                             os.remove(f)
                
#                     bed = inp_bed

#                     with open(bed, 'r') as file:
#                         all_bed_data = [line.strip() for line in file]
                    
#                     with open(bed, 'r') as file:
#                         nrecrds = sum(1 for line in file)
                    
#                     nrecrdst = nrecrds/3

#                     print(f"read {nrecrdst} records")
#                     csv = inp_csv

#                     with open(csv, 'r') as file:
#                         all_csv_data = [line.strip() for line in file]
                    
#                     i = 0
#                     eventn = 0

#                     while i<nrecrds :
#                         # Construct string for ggsashimi
#                         line1 = all_bed_data[i]
#                         line2 = all_bed_data[i+1]
#                         line3 = all_bed_data[i+2]

#                         i = i + 3

#                         # Read strand
#                         strnd = line1.iloc[:, 5]

#                         # Also read TxID
#                         TxID = line1.iloc[7]

#                         if strnd=='+' :
#                             strndflg = "plus"

#                             # us exon length
#                             exon1 = line1.iloc[:, 3]

#                             # ds exon length
#                             exon2 = line3.iloc[:, 3]

#                             # Now modify to reflect whole up/dn exons
#                             up = line1.iloc[:, 1]
#                             upn = up - exon1
# 				             
#                             # ds exon
#                             dn = line3.iloc[:, 2]
#                             dnn = dn + exon2

#                             l123 = f"{line1} {line2} {line3}".split()
#                             fulltitle = f"{l123[0]}-{l123[1]}:{l123[2]}-{l123[9]}:{l123[10]}-{l123[17]}:{l123[18]}"

#                             comb = f"{line1} {upn} {dnn}".split()
#                             line = f"{comb[0]}:{comb[8]}-{comb[9]}"

#                         else :
#                             # us exon length
#                             exon1 = line3.iloc[:, 3]

#                             # ds exon length
#                             exon2 = line1.iloc[:, 3]

#                             # Now modify to reflect whole up/dn exons
#                             dn = line3.iloc[:, 1]
#                             dnn = dn - exon1
# 				             
#                             # ds exon
#                             up = line1.iloc[:, 2]
#                             upn = up + exon2

#                             l123 = f"{line1} {line2} {line3}".split()
#                             fulltitle = f"{l123[0]}-{l123[17]}:{l123[18]}-{l123[9]}:{l123[10]}-{l123[1]}:{l123[2]}"

#                             comb = f"{line1} {dnn} {upn}".split()
#                             line = f"{comb[0]}:{comb[8]}-{comb[9]}"

#                         event = all_csv_data[eventn]

#                         fields = event.split(',')
#                         fn = f"{fields[7]}-{fields[1]}_{fields[2]}-{fields[3]}"

#                         # String for majiq event
#                         chr_name = event.split(',')[1]
#                         start = event.split(',')[2]
#                         end = event.split(',')[3]
#                         gene_name = event.split(',')[7]
#                         intron = end - start

#                         #gene_name-start-intron-end #for now using strand from ggsashimi
#                         comb_line = f"{chr_name} {start} {end}".split()
#                         majiq_event = f"{comb_line[0]}-{comb_line[1]}-{comb_line[2]}"

#                         # Also get actual event identified
#                         event_identified = f"{line2.iloc[:,0]}-{line2.iloc[:,1]}-{line2.iloc[:,2]}"
#                         eventn = eventn + 1

#                         print(f"processing event num {eventn} and event {fn}")

#                         command = [
#                             "./ggsashimi_txV3.py",
#                             "-A", "median_j",
#                             "-b", "all_bams.tsv",
#                             "-c", line,
#                             "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
#                             "-GeneName", gene_name,
#                             "-MajiqStrnd", strnd,
#                             "-ORIG", "1",
#                             "-UEX", exon1,
#                             "-DEX", exon2,
#                             "-FullTitle", fulltitle,
#                             "-MajiqTx", majiq_event,
#                             "-Majiq", fn,
#                             "-Tx", TxID,
#                             "-M", "1",
#                             "-C", "3",
#                             "-o", f"{inp_prefix}/ce_ext_sashimi_plots/{fn}",
#                             "-O", "3",
#                             "--alpha", "0.25",
#                             "--base-size=20",
#                             "--ann-height=2.5",
#                             "--height=2.5",
#                             "--width=18",
#                             "-P", "palette.txt"
#                         ]

#                         result = subprocess.run(command)
#                         print(result.stdout)

#                     # Now merge all pdf's
#                     command = ["python", "merge_sashimis.py", f"{inp_prefix}/ce_ext_sashimi_plots/"]
#                     result = subprocess.run(command)
#                     print(result.stdout)
                
#                 ######THIS IS FOR ALL IR EVENTS
#                 if int(args[2])==5 :
#                     os.makedirs(f"{inp_prefix}/ir_sashimi_plots",exist_ok=True)
#                     if len(os.listdir(f"{inp_prefix}/ir_sashimi_plots/"))!=0 :
#                         for f in glob.glob(f"{inp_prefix}/ir_sashimi_plots/*.*") :
#                             os.remove(f)
                
#                     bed = inp_bed

#                     with open(bed, 'r') as file:
#                         all_bed_data = [line.strip() for line in file]
                    
#                     with open(bed, 'r') as file:
#                         nrecrds = sum(1 for line in file)

#                     print(f"read {nrecrds} records")
#                     csv = inp_csv

#                     with open(csv, 'r') as file:
#                         all_csv_data = [line.strip() for line in file]
                    
#                     i = 0
#                     eventn = 0

#                     while i<nrecrds :
#                         # Construct string for ggsashimi
#                         line1 = all_bed_data[i]

#                         i = i + 1

#                         # Read strand
#                         strnd = line1.iloc[:, 5]

#                         # Also read TxID
#                         TxID = line1.iloc[7]

#                         line = f"{line1[0]}:{line1[1]-50}-{line1[2]+50}"

#                         # Get majiq event
#                         event = all_csv_data[eventn]
#                         eventn = eventn + 1

#                         # THIS IS TO KEEP DIFFERENT FILE NAMES FOR EVENTS WITH IDENTICAL COORDINATES AND GENE_NAMES
#                         gene_name = event.split(',')[4]

#                         if i==1 :
#                             temp_gene = event.split(',')[4]
#                             trackj = 1
                        
#                         elif temp_gene == gene_name :
#                             trackj = trackj + 1
                        
#                         else :
#                             temp_gene = event.split(',')[4]
#                             trackj = 1
                        
#                         # SHOULD NOT HAVE THIS BUT NEED TO KEEP AS PER REQUEST
#                         fields = event.split(',')
#                         fn = f"{fields[4]}-{fields[0]}_{fields[1]}-{fields[2]}-{trackj}" #also gene_id to avoid same file names for repeated events

#                         # String for majiq event
#                         chr_name = event.split(',')[0]
#                         start = event.split(',')[1]
#                         end = event.split(',')[2]

#                         comb_line = f"{chr_name} {start} {end}".split()
#                         majiq_event = f"{comb_line[0]}-{comb_line[1]}-{comb_line[2]}"
#                         exon1 = 0
#                         exon2 = 0

#                         # Also get actual event identified
#                         event_identified = f"{line1.iloc[:,0]}-{line1.iloc[:,1]}-{line1.iloc[:,2]}"

#                         print(f"processing event num {eventn} and event {fn}")

#                         command = [
#                             "./ggsashimi_txV3.py",
#                             "-A", "median_j",
#                             "-b", "all_bams.tsv",
#                             "-c", line,
#                             "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
#                             "-GeneName", gene_name,
#                             "-MajiqStrnd", strnd,
#                             "-ORIG", "1",
#                             "-UEX", exon1,
#                             "-DEX", exon2,
#                             "-MajiqTx", majiq_event,
#                             "-Majiq", fn,
#                             "-Tx", TxID,
#                             "-M", "1",
#                             "-C", "3",
#                             "-o", f"{inp_prefix}/ir_sashimi_plots/{fn}",
#                             "-O", "3",
#                             "--alpha", "0.25",
#                             "--base-size=20",
#                             "--ann-height=2.5",
#                             "--height=2.5",
#                             "--width=18",
#                             "-P", "palette.txt"
#                         ]

#                         result = subprocess.run(command)
#                         print(result.stdout)

#                     # Now merge all pdf's
#                     command = ["python", "merge_sashimis.py", f"{inp_prefix}/ir_sashimi_plots/"]
#                     result = subprocess.run(command)
#                     print(result.stdout)
                    
#             else :
#                 print("came in for peaksbackmapping")
                
#                 ###### THIS PART IS FOR BACKMAPPING
#                 flg_sashimi_files_only = 1
                
#                 if flg_sashimi_files_only==1 :
#                     step01_flg = 1
                    
#                     if step01_flg==1 :
#                         ################################ GENERATING BED FILE FOR SASHIMI PLOTS FOR ALL EVENTS
#                         # IMPORTANT - THIS CODE RELIES ON THE OUTPUT OF BEDTOOLS CLOSEST FUNCTION
#                         # GIVEN A MAJIQ EVENT (as a bed file) and BED FILE FOR THE TRANSCRIPT IT BELIEVES TO BE PART OF
#                         # BEDTOOLS CLOSEST FUNCTION USES INPUT RANGE (AMJIQ JUNTION) and finds closest exons (and their ranges) from the TRANSCRIPT BED FILE
#                         # BEDTOOLS CLOSEST RETURNS A BED FILE WITH FOLLOWING OUTPUT
#                         # INPUT: chr,start,end,1,0,+ (majiq event with strand)
#                         # OUTPUT: chr,start,end,1,0,+,chr,start,end,size,exon_num,+,distance
#                         # Here first 6 entries are the original majiq input and next 7 entries are the resulting closest exon coordinates, its size (in bp),exon_num,strand and distance from reference
#                         # SOULD ADD CHECKS ON EXON_NUM READING FROM FILE
                        
#                         inp = args[0].split('.')[0]
#                         print(f"run_sashimiV1.sh, got inp as {inp}")
                        
#                         # Also get folder
#                         folder = inp.split('/')[0]
#                         print(f"run_sashimiV1.sh, got folder as {folder}")
                        
#                         # Also get event TYPE
#                         eventyp = args[1].split('/')[1].split('.')[0]
#                         print(f"run_sashimiV1.sh, got eventype as {eventyp}")
                        
#                         if os.path.exists(f"{inp}_majiq.bed"):
#                             os.remove(f"{inp}_majiq.bed")
                        
#                         if os.path.exists(f"{inp}_majiq.csv"):
#                             os.remove(f"{inp}_majiq.csv")
                            
#                         if os.path.exists("all_tx_events.csv"):
#                             os.remove("all_tx_events.csv")
                            
#                         if os.path.exists("majiq_events_sashimi2.bed"):
#                             os.remove("majiq_events_sashimi2.bed")
                            
#                         if os.path.exists("majiq_events_sashimi2.csv"):
#                             os.remove("majiq_events_sashimi2.csv")
                            
#                         if os.path.exists("majiq_events_sashimi01.bed"):
#                             os.remove("majiq_events_sashimi01.bed")
                            
#                         if os.path.exists("majiq_events_sashimi01.csv"):
#                             os.remove("majiq_events_sashimi01.csv")
                            
#                         if os.path.exists("majiq_events_sashimi02.bed"):
#                             os.remove("majiq_events_sashimi02.bed")
                            
#                         if os.path.exists("majiq_events_sashimi02.csv"):
#                             os.remove("majiq_events_sashimi02.csv")
                            
#                         if os.path.exists("majiq_events_progress2.txt"):
#                             os.remove("majiq_events_progress2.txt")
                            
#                         if os.path.exists("majiq_events_progress01.txt"):
#                             os.remove("majiq_events_progress01.txt")
                            
#                         if os.path.exists("majiq_events_progress02.txt"):
#                             os.remove("majiq_events_progress02.txt")
                            
#                         if os.path.exists("majiq_events_progress_all.txt"):
#                             os.remove("majiq_events_progress_all.txt")
                            
#                         if os.path.exists("majiq_events_progress1.txt"):
#                             os.remove("majiq_events_progress1.txt")
                            
#                         flg = 1
                        
#                         if flg==1 :
#                             os.makedirs("event_bedfiles",exist_ok=True)
#                             if len(os.listdir("event_bedfiles/"))!=0 :
#                                 for f in glob.glob("event_bedfiles/*.*") :
#                                     os.remove(f)
                            
#                             command = [
#                                 "Rscript",
#                                 "txens.R",
#                                 args[0],
#                                 "principal_txs.csv"
#                             ]

#                             result = subprocess.run(command, capture_output=True, text=True)
#                             print(result.stdout)
                            
#                         with open("all_tx_events.csv", 'r') as file:
#                             csv_data = [line.strip() for line in file]
                        
#                         csvi = 0
#                         samples = glob.glob('event_bedfiles/temp_*.bed')
                        
#                         for sample in samples :
#                             # Read csv entry
#                             csv_ln = csv_data[csvi]
#                             csvi = csvi + 1
                            
#                             print(f"processing {sample}")
#                             allexons = sample.split('/')[1].split('_')[1]
#                             gene_name1 = allexons.split('.')[0]
#                             gene_name = gene_name1.split('-')[0]
                            
#                             # First sort the bed
#                             bed = pybedtools.BedTool(f"event_bedfiles/{allexons}")
#                             sorted_bed = bed.sort()
#                             sorted_bed.saveas(f"event_bedfiles/t{allexons}")
                            
#                             # Also read Tx Files to retrieve selected Tx - should find better ways
#                             with open(f"event_bedfiles/TxID{allexons}", 'r') as file:
#                                 first_line = file.readline().strip()
#                                 TxID = first_line.split()[6]
                                
#                             samp = pd.read_csv(sample, sep='\t', header=None)
#                             strnd = samp.iloc[:, 5]
                            
#                             # Get distance to downstream exon (for ties, report first) from current reference and pick start, end and d
#                             a = pybedtools.BedTool(sample)
#                             b = pybedtools.BedTool(f'event_bedfiles/t{allexons}')
#                             closest = a.closest(b, d=True, s=True, t='first')
#                             ds = [line.strip() for line in closest]
#                             ds = closest.to_dataframe()
                            
#                             # Also get distance to upstream exon from current reference and pick start, end and d
#                             closest = a.closest(b, d=True, s=True, t='last')
#                             us = [line.strip() for line in closest]
#                             us = closest.to_dataframe()
                            
#                             # Get up and down stream exon numbers
#                             upexon = us.iloc[:,10]
#                             dnexon = ds.iloc[:,10]
                            
#                             # Events star and end
#                             event_st = us.iloc[:,1]
#                             event_end = us.iloc[:,2]
#                             diff_exon = upexon - dnexon
                            
#                             # Take absolute value
#                             diff_exon_abs = abs(diff_exon)
                            
#                             if diff_exon_abs>=1 : # ALL EVENTS THAT SPANS 2 OR MORE EXONS
#                                 if strnd=="+" :
#                                     start = us.iloc[:,7]
#                                     end = ds.iloc[:,8]
                                    
#                                 else :
#                                     start = ds.iloc[:,7]
#                                     end = us.iloc[:,8]
                                
#                                 # Also save
#                                 # FORMAT IS: chr, start of up_ex,end of ds_ex, 1, 0, strand, gene_name, TxID
#                                 # First check if event lies between selected exons
#                                 if (start<=event_st and end>=event_end) :
#                                     input_data = us + [start, end, strnd] + gene_name + TxID
#                                     output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
#                                     with open(f"{inp}_majiq.bed", "a") as f:
#                                         f.write("\t".join(output_data) + "\n")
                                        
#                                     with open(f"{inp}_majiq.csv", "a") as f:
#                                         f.write(csv_ln)
                                    
#                                 else :
#                                     if strnd=="+" :
#                                         start = us.iloc[:,7]
#                                         end = ds.iloc[:,8] # Both are same
                                        
#                                         # First check if star > event_start, then select upstream exon
#                                         if start>=event_st :
#                                             # Get exon (line number in bed file) to read
#                                             exon = ds.iloc[:,10]
#                                             exon = exon - 2
                                            
#                                             with open(f"event_bedfiles/{allexons}", 'r') as file:
#                                                 bed_data = [line.strip() for line in file]
                                            
#                                             bed_ln = bed_data[exon]
                                            
#                                             # Update start
#                                             start = bed_ln.iloc[:,1]
                                            
#                                         # Now check if end <event_end
#                                         if end<=event_end :
#                                             # Get exon (line number in bed file) to read
#                                             exon = us.iloc[:,10]
                                            
#                                             with open(f"event_bedfiles/{allexons}", 'r') as file:
#                                                 bed_data = [line.strip() for line in file]
                                            
#                                             print(f"bed_data {bed_data}")
#                                             bed_ln = bed_data[exon]
#                                             print(f"exon {exon} bed_ln {bed_ln}")
                                            
#                                             # Update end
#                                             end = bed_ln.iloc[:,2]
                                            
#                                         # Now one more time check if event lies between selected exons
#                                         if (start<=event_st and end>=event_end):
#                                             input_data = us + [start, end, strnd] + gene_name + TxID
#                                             output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
#                                             with open(f"{inp}_majiq.bed", "a") as f:
#                                                 f.write("\t".join(output_data) + "\n")
                                                
#                                             with open(f"{inp}_majiq.csv", "a") as f:
#                                                 f.write(csv_ln)
                                            
#                                         else :
#                                             with open("majiq_events_progress2.txt", "a") as f:
#                                                 f.write(f"ds 1 {ds}\n")
#                                                 f.write(f"us 1 {us}\n")
                                                
#                                             print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                                                
#                                     else : # THIS IS FOR NEGATIVE STRAND
#                                         start = ds.iloc[:,7]
#                                         end = us.iloc[:,8]
                                        
#                                         # First check if star > event_start, then select upstream exon
#                                         if start>=event_st :
#                                             # Get exon (line number in bed file) to read
#                                             exon = ds.iloc[:,10]
#                                             with open(f"event_bedfiles/{allexons}", 'r') as file:
#                                                 bed_data = [line.strip() for line in file]
                                            
#                                             bed_ln = bed_data[exon]
                                            
#                                             # Update start
#                                             start = bed_ln.iloc[:,1]
                                            
#                                         # Now check if end <event_end
#                                         if end<=event_end :
#                                             # Get exon (line number in bed file) to read
#                                             exon = us.iloc[:,10]
#                                             exon = exon - 2
#                                             with open(f"event_bedfiles/{allexons}", 'r') as file:
#                                                 bed_data = [line.strip() for line in file]
                                            
#                                             bed_ln = bed_data[exon]
                                            
#                                             # Update end
#                                             end = bed_ln.iloc[:,2]
                                        
#                                         # Now one more time check if event lies between selected exons
#                                         if (start<=event_st and end>=event_end) :
#                                             input_data = us + [start, end, strnd] + gene_name + TxID
#                                             output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
#                                             with open(f"{inp}_majiq.bed", "a") as f:
#                                                 f.write("\t".join(output_data) + "\n")
                                                
#                                             with open(f"{inp}_majiq.csv", "a") as f:
#                                                 f.write(csv_ln)
                                            
#                                         else :
#                                             with open("majiq_events_progress2.txt", "a") as f:
#                                                 f.write(f"ds 2 {ds}\n")
#                                                 f.write(f"us 2 {us}\n")
                                                
#                                             print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                                        
#                             elif diff_exon_abs==0 :
#                                 if strnd=="+" :
#                                     start = us.iloc[:,7]
#                                     end = ds.iloc[:,8] # Both are same
                                    
#                                     # First check if star > event_start, then select upstream exon
#                                     if start>=event_st :
#                                         # Get exon (line number in bed file) to read
#                                         exon = us.iloc[:,10]
#                                         exon = exon - 2
                                        
#                                         with open(f"event_bedfiles/{allexons}", 'r') as file:
#                                             bed_data = [line.strip() for line in file]
                                        
#                                         bed_ln = bed_data[exon]
                                        
#                                         # Update start
#                                         start = bed_ln.iloc[:,1]
                                        
#                                         # Now go on the other side
#                                         if start>=event_st :
#                                             # Get exon (line number in bed file) to read
#                                             exon = us.iloc[:,10]
#                                             with open(f"event_bedfiles/{allexons}", 'r') as file:
#                                                 bed_data = [line.strip() for line in file]
                                            
#                                             bed_ln = bed_data[exon]
                                            
#                                             # Update start
#                                             start = bed_ln.iloc[:,1]
                                       
#                                     # Now check if end <event_end
#                                     if end<=event_end :
#                                         # Get exon (line number in bed file) to read
#                                         exon = us.iloc[:,10]
                                        
#                                         with open(f"event_bedfiles/{allexons}", 'r') as file:
#                                             bed_data = [line.strip() for line in file]
                                        
#                                         print(f"bed_data {bed_data}")
#                                         bed_ln = bed_data[exon]
#                                         print(f"exon {exon} bed_ln {bed_ln}")
                                        
#                                         # Update end
#                                         end = bed_ln.iloc[:,2]
                                        
#                                         if end<=event_end :
#                                             # Get exon (line number in bed file) to read
#                                             exon = us.iloc[:,10]
#                                             exon = exon - 2 # Reading line for readarray starts from 0
                                            
#                                             with open(f"event_bedfiles/{allexons}", 'r') as file:
#                                                 bed_data = [line.strip() for line in file]
                                            
#                                             bed_ln = bed_data[exon]
                                            
#                                             # Update end
#                                             end = bed_ln.iloc[:,2]
                                    
#                                     # Now one more time check if event lies between selected exons
#                                     if (start<=event_st and end>=event_end):
#                                         input_data = us + [start, end, strnd] + gene_name + TxID
#                                         output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
#                                         with open(f"{inp}_majiq.bed", "a") as f:
#                                             f.write("\t".join(output_data) + "\n")
                                            
#                                         with open(f"{inp}_majiq.csv", "a") as f:
#                                             f.write(csv_ln)
                                        
#                                     else :
#                                         with open("majiq_events_progress01.txt", "a") as f:
#                                             f.write(f"ds {ds}\n")
#                                             f.write(f"us {us}\n")
                                            
#                                         print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                                            
#                                 else : # THIS IS FOR NEGATIVE STRAND
#                                     start = ds.iloc[:,7]
#                                     end = us.iloc[:,8]
                                    
#                                     # First check if star > event_start, then select upstream exon
#                                     if start>=event_st :
#                                         # Get exon (line number in bed file) to read
#                                         exon = ds.iloc[:,10]
#                                         with open(f"event_bedfiles/{allexons}", 'r') as file:
#                                             bed_data = [line.strip() for line in file]
                                        
#                                         bed_ln = bed_data[exon]
                                        
#                                         # Update start
#                                         start = bed_ln.iloc[:,1]
                                        
#                                     # Now check if end <event_end
#                                     if end<=event_end :
#                                         # Get exon (line number in bed file) to read
#                                         exon = us.iloc[:,10]
#                                         exon = exon - 2
                                        
#                                         with open(f"event_bedfiles/{allexons}", 'r') as file:
#                                             bed_data = [line.strip() for line in file]
                                        
#                                         bed_ln = bed_data[exon]
                                        
#                                         # Update end
#                                         end = bed_ln.iloc[:,2]
                                    
#                                     # Now one more time check if event lies between selected exons
#                                     if (start<=event_st and end>=event_end) :
#                                         input_data = us + [start, end, strnd] + gene_name + TxID
#                                         output_data = [input_data[0], input_data[13], input_data[14], "1", "0", input_data[15], input_data[16], input_data[17]]
#                                         with open(f"{inp}_majiq.bed", "a") as f:
#                                             f.write("\t".join(output_data) + "\n")
                                            
#                                         with open(f"{inp}_majiq.csv", "a") as f:
#                                             f.write(csv_ln)
                                        
#                                     else :
#                                         with open("majiq_events_progress02.txt", "a") as f:
#                                             f.write(f"ds {ds}\n")
#                                             f.write(f"us {us}\n")
                                            
#                                         print(f"diff_exon_abs is {diff_exon_abs} selected event {sample} has event_st {event_st} selected start {start} event end {event_end} selected end {end} - please check")
                                        
#                                 # First check if start > event_st, then select upstream exon
                                
#                             else :
#                                 with open("majiq_events_progress_all.txt", "a") as f:
#                                     f.write(f"ds {ds}\n")
#                                     f.write(f"us {us}\n")
                                        
#                     flag = args[3]
                    
#                     if flag==1 :
#                         os.makedirs(f"{folder}/sashimi_plots/{eventyp}",exist_ok=True)
#                         if len(os.listdir(f"{folder}/sashimi_plots/{eventyp}/"))!=0 :
#                             for f in glob.glob(f"{folder}/sashimi_plots/{eventyp}/*.pdf") :
#                                 os.remove(f)
                            
#                             for f in glob.glob(f"{folder}/sashimi_plots/{eventyp}/*.svg") :
#                                 os.remove(f)
                            
#                         bed = f"{inp}_majiq.bed"
#                         with open(bed, 'r') as file:
#                             all_bed_data = [line.strip() for line in file]
                        
#                         with open(bed, 'r') as file:
#                             nrecrds = sum(1 for line in file)
                            
#                         print("read {nrecrds} records")
                        
#                         csv = f"{inp}_majiq.csv"
#                         with open(csv, 'r') as file:
#                             all_csv_data = [line.strip() for line in file]
                        
#                         # Read peaks bed file
#                         peak_bed = args[1]
#                         with open(peak_bed, 'r') as file:
#                             all_peak_data = [line.strip() for line in file]
                        
#                         # Also read peaks csv file (to get AA seq)
#                         with open(args[2], 'r') as file:
#                             all_peak_csv = [line.strip() for line in file]
                        
#                         i = 0
#                         eventn = 0
#                         if len(os.listdir("temp_pdfs/"))!=0 :
#                             for f in glob.glob("temp_pdfs/*.pdf") :
#                                 os.remove(f)
                        
#                         while i<nrecrds :
#                             # Construct string for ggsashimi
#                             line1 = all_bed_data[i]
#                             AA = all_peak_csv[i]
#                             i = i + 1
                            
#                             # Read strand
#                             strnd = line1.iloc[:,5]
                            
#                             # Also read TxID
#                             TxID = line1.iloc[:,7]
#                             line = f"{line1.iloc[:,0]}:{line1.iloc[:,1]-50}-{line1.iloc[:,2]+50}" # This is the actual event
                            
#                             # Get majiq event
#                             event = all_csv_data[eventn]
#                             peak_event = all_peak_data[eventn]
#                             eventn = eventn + 1
                            
#                             # THIS IS TO KEEP DIFFERENT FILE NAMES FOR EVENTS WITH IDENTICAL COORDINATES AND GENE_NAMES
#                             gene_name = event.split(',')[4]
                            
#                             if i==1 :
#                                 temp_gene = event.split(',')[4]
#                                 trackj = 1
                                
#                             elif temp_gene == gene_name :
#                                 trackj = trackj + 1
                                
#                             else :
#                                 temp_gene = event.split(',')[4]
#                                 trackj = 1
                                
#                             # SHOULD NOT HAVE THIS BUT NEED TO KEEP AS PER REQUEST
#                             fn = f"{event.split(',')[4]}-{event.split(',')[0]}-{event.split(',')[1]}-{event.split(',')[2]}-{trackj}" # Also gene_id to avoid same file names for repeated events
                            
#                             # String for majiq event
#                             chr_name = event.split(',')[0]
#                             start = event.split(',')[1]
#                             end = event.split(',')[2]
                            
#                             comb_line = f"{chr_name} {start} {end}".split()
#                             majiq_event = f"{comb_line[0]}-{comb_line[1]}-{comb_line[2]}"
                            
#                             # ALSO ADD PEAKS EVENT
#                             pchr_name = peak_event.split('\t')[0]
#                             pstart = peak_event.split('\t')[1]
#                             pend = peak_event.split('\t')[2]
                            
#                             pstart1 = peak_event.split('\t')[1]
#                             pend1 = peak_event.split('\t')[2]
                            
#                             # Also get length of event range
#                             sz1 = pend1 - pstart1
                            
#                             pcomb_line = f"{pchr_name} {pstart} {pend}".split()
#                             peaks_event = f"{pcomb_line[0]}-{pcomb_line[1]}-{pcomb_line[2]}"
                            
#                             # Also get AA seq
#                             AAseq1 = AA.split(',')[0]
                            
#                             # Now add new line if string is larger than 100 characters
#                             length_AA_flg = 0
                            
#                             if len(AAseq1)>100 :
#                                 AAseq2 = AAseq1[:100]
#                                 AAseq3 = AAseq1[100:]
                                
#                                 # Now concatenate
#                                 AAseq = f"{AAseq2}\n{AAseq3}"
#                                 length_AA_flg = 1
                                
#                             else :
#                                 AAseq = AAseq1
                            
#                             # Also get nt seq
#                             NTseq1 = AA.split(',')[6]
#                             length_nt_flg = 0
                            
#                             if len(NTseq1)-2>100 :
#                                 NTseq2 = NTseq1[:100]
#                                 NTseq3 = NTseq1[100:]
                                
#                                 # Now concatenate
#                                 NTseq = f"{NTseq2}\n{NTseq3}"
#                                 length_nt_flg = 1
                                
#                             else :
#                                 NTseq = NTseq1
                            
#                             # Now concatenate gene_name, NTseq and AA seq - here -4 is due to _num and two characters for new line
#                             if (length_nt_flg==1 and length_AA_flg==0) :
#                                 comb = f"{gene_name} {AAseq} {len(AAseq)} {NTseq} {len(NTseq)-4}".split()
#                                 gene_name1 = f"Gene : {comb[0]}"
                                
#                             elif (length_nt_flg==1 and length_AA_flg==1) :
#                                 comb = f"{gene_name} {AAseq} {len(AAseq)-2} {NTseq} {len(NTseq)-4}".split()
#                                 gene_name1 = f"Gene : {comb[0]}"
                                
#                             elif (length_nt_flg==0 and length_AA_flg==1) :
#                                 comb = f"{gene_name} {AAseq} {len(AAseq)-2} {NTseq} {len(NTseq)-2}".split()
#                                 gene_name1 = f"Gene : {comb[0]}"
                            
#                             else :
#                                 comb = f"{gene_name} {AAseq} {len(AAseq)} {NTseq} {len(NTseq)-2}".split()
#                                 gene_name1 = f"Gene : {comb[0]}"
                            
#                             ######## END PEAKS EVENT
#                             exon1 = 0
#                             exon2 = 0
                            
#                             # Also get actual event identified
#                             event_identified = f"None -{line1.iloc[:,1]}-{line1.iloc[:,2]}"
#                             print("processing event num {eventn} and event {fn}")
                            
#                             command = [
#                                 "python",
#                                 "ggsashimi_txV3.py",
#                                 "-A", "median_j",
#                                 "-b", "all_bams.tsv",
#                                 "-c", line,
#                                 "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
#                                 "PeaksAA", sz1,
#                                 "PeaksFlg", 1,
#                                 "PeaksTx", peaks_event,
#                                 "-GeneName", gene_name1,
#                                 "-MajiqStrnd", strnd,
#                                 "-ORIG", "1",
#                                 "-UEX", exon1,
#                                 "-DEX", exon2,
#                                 "-MajiqTx", majiq_event,
#                                 "PGPTx", event_identified,
#                                 "-Majiq", fn,
#                                 "-Tx", TxID,
#                                 "-M", "1",
#                                 "-C", "3",
#                                 "-o", f"{folder}/sashimi_plots/{eventyp}/{fn}",
#                                 "-O", "3",
#                                 "--alpha", "0.25",
#                                 "--base-size=20",
#                                 "--ann-height=4",
#                                 "--height=3",
#                                 "--width=18",
#                                 "-P", "palette.txt"
#                             ]

#                             result = subprocess.run(command)
#                             print(result.stdout)
                            
#                             if os.path.exists(f"{folder}/sashimi_plots/{eventyp}/{fn}.svg"):
#                                 os.remove(f"{folder}/sashimi_plots/{eventyp}/{fn}.svg")
                            
#                         # Now merge all pdf's
#                         command = ["python", "merge_sashimis.py", f"{folder}/sashimi_plots/{eventyp}/"]
#                         result = subprocess.run(command)
#                         print(result.stdout)
                        
#                     flag = args[3]
                    
#                     # THIS ADDS TWO EXONS IN PEAKS EVENTS SHARING EXON_CE OR CE_EXON
#             		# USES PeaksFlg=2 for ggsashimi_txV2.py
#             		# THIS SECTION DEALS WITH UEX_CE, CE_DEX for CE events (flag=0 or 1) and UEX_DEX for SKIPTICS (flag=2)
#                     if flag==2 :
#                         os.makedirs(f"{folder}/sashimi_plots/{eventyp}",exist_ok=True)
#                         if len(os.listdir(f"{folder}/sashimi_plots/{eventyp}/"))!=0 :
#                             for f in glob.glob(f"{folder}/sashimi_plots/{eventyp}/*.pdf") :
#                                 os.remove(f)
                            
#                             for f in glob.glob(f"{folder}/sashimi_plots/{eventyp}/*.svg") :
#                                 os.remove(f)
                        
#                         print(f"current wkd is {os.getcwd()}")
#                         print(f'came in "{flag} -eq 2" area with bed file {inp}.bed')
                        
#                         bed = f"{inp}_majiq.bed"
#                         with open(bed, 'r') as file:
#                             all_bed_data = [line.strip() for line in file]
                        
#                         with open(bed, 'r') as file:
#                             nrecrds = sum(1 for line in file)
                            
#                         print("read {nrecrds} records")
                        
#                         csv = f"{inp}_majiq.csv"
#                         with open(csv, 'r') as file:
#                             all_csv_data = [line.strip() for line in file]
                        
#                         # ALSO READ PEAKS EVENTS BED FILE, HERE IT HAS 2 EXONS PER EVENT
#                         peak_bed = args[1]
#                         with open(peak_bed, 'r') as file:
#                             all_peak_data = [line.strip() for line in file]
                        
#                         # Also read peaks csv file (to get AA seq)
#                         with open(args[2], 'r') as file:
#                             all_peak_csv = [line.strip() for line in file]
                        
#                         i = 0
#                         eventn = 0
#                         if len(os.listdir("temp_pdfs/"))!=0 :
#                             for f in glob.glob("temp_pdfs/*.pdf") :
#                                 os.remove(f)
                        
#                         while i<nrecrds :
#                             # Construct string for ggsashimi
#                             line1 = all_bed_data[i]
                            
#                             # Get majiq event
#                             event = all_csv_data[i]
#                             AA = all_peak_csv[i]
#                             i = i + 1
                            
#                             # Read strand
#                             strnd = line1.iloc[:,5]
                            
#                             # Also read TxID
#                             TxID = line1.iloc[:,7]
#                             line = f"{line1.iloc[:,0]}:{line1.iloc[:,1]-50}-{line1.iloc[:,2]+50}" # This is the actual event
                            
#                             peak_event = all_peak_data[eventn]
#                             eventn = eventn + 1
                            
#                             peak_event1 = all_peak_data[eventn]
#                             eventn = eventn + 1
                            
#                             # THIS IS TO KEEP DIFFERENT FILE NAMES FOR EVENTS WITH IDENTICAL COORDINATES AND GENE_NAMES
#                             gene_name = event.split(',')[4]
                            
#                             if i==1 :
#                                 temp_gene = event.split(',')[4]
#                                 trackj = 1
                                
#                             elif temp_gene == gene_name :
#                                 trackj = trackj + 1
                                
#                             else :
#                                 temp_gene = event.split(',')[4]
#                                 trackj = 1
                                
#                             # SHOULD NOT HAVE THIS BUT NEED TO KEEP AS PER REQUEST
#                             fn = f"{event.split(',')[4]}-{event.split(',')[0]}-{event.split(',')[1]}-{event.split(',')[2]}-{trackj}" # Also gene_id to avoid same file names for repeated events
                            
#                             # String for majiq event
#                             chr_name = event.split(',')[0]
#                             start = event.split(',')[1]
#                             end = event.split(',')[2]
                            
#                             comb_line = f"{chr_name} {start} {end}".split()
#                             majiq_event = f"{comb_line[0]}-{comb_line[1]}-{comb_line[2]}"
                            
#                             # ALSO ADD 2 PEAKS EVENT - one for us/ds ans other for ce
#                             pchr_name = peak_event.split('\t')[0]
#                             pstart = peak_event.split('\t')[1]
#                             pend = peak_event.split('\t')[2]
                            
#                             pchr_name1 = peak_event1.split('\t')[0]
#                             pstart1 = peak_event1.split('\t')[1]
#                             pend1 = peak_event1.split('\t')[2]
                            
#                             pcomb_line = f"{pchr_name} {pstart} {pend} {pstart1} {pend1}".split()
#                             peaks_event = f"{pcomb_line[0]}-{pcomb_line[1]}-{pcomb_line[2]}-{pcomb_line[3]}-{pcomb_line[4]}"
                            
#                             # Also get AA seq
#                             # Also get nt seq
#                             NTseq1 = AA.split(',')[6]
#                             length_nt_flg = 0
                            
#                             if len(NTseq1)-2>100 :
#                                 NTseq2 = NTseq1[:100]
#                                 NTseq3 = NTseq1[100:]
                                
#                                 # Now concatenate
#                                 NTseq = f"{NTseq2}\n{NTseq3}"
#                                 length_nt_flg = 1
                                
#                             else :
#                                 NTseq = NTseq1
                            
#                             AAseq1 = AA.split(',')[0]
#                             length_AA_flg = 0
                            
#                             if len(AAseq1)>100 :
#                                 AAseq2 = AAseq1[:100]
#                                 AAseq3 = AAseq1[100:]
                                
#                                 # Now concatenate
#                                 AAseq = f"{AAseq2}\n{AAseq3}"
#                                 length_AA_flg = 1
                                
#                             else :
#                                 AAseq = AAseq1
                            
#                             # Now concatenate gene_name and NTseq
#                             ###########
#                             if (length_nt_flg==1 and length_AA_flg==0) :
#                                 comb = f"{gene_name} {AAseq} {len(AAseq)} {NTseq} {len(NTseq)-4}".split()
#                                 gene_name1 = f"Gene : {comb[0]}"
                                
#                             elif (length_nt_flg==1 and length_AA_flg==1) :
#                                 comb = f"{gene_name} {AAseq} {len(AAseq)-2} {NTseq} {len(NTseq)-4}".split()
#                                 gene_name1 = f"Gene : {comb[0]}"
                                
#                             elif (length_nt_flg==0 and length_AA_flg==1) :
#                                 comb = f"{gene_name} {AAseq} {len(AAseq)-2} {NTseq} {len(NTseq)-2}".split()
#                                 gene_name1 = f"Gene : {comb[0]}"
                            
#                             else :
#                                 comb = f"{gene_name} {AAseq} {len(AAseq)} {NTseq} {len(NTseq)-2}".split()
#                                 gene_name1 = f"Gene : {comb[0]}"
                            
#                             ############
#                             # Now adding upexon and ce part or ce and dnexon part separately
#                             if args[4]==0 :
#                                 NTseque1 = AA.split(',')[7]
#                                 if len(NTseque1)-2>100 :
#                                     NTseque2 = NTseque1[:100]
#                                     NTseque3 = NTseque1[100:]
                                    
#                                     # Now concatenate
#                                     NTseque = f"{NTseque2}\n{NTseque3}"
                                    
#                                 else :
#                                     NTseque = NTseque1
                                
#                                 NTseqce1 = AA.split(',')[8]
#                                 if len(NTseqce1)-2>100 :
#                                     NTseqce2 = NTseqce1[:100]
#                                     NTseqce3 = NTseqce1[100:]
                                    
#                                     # Now concatenate
#                                     NTseqce = f"{NTseqce2}\n{NTseqce3}"
                                    
#                                 else :
#                                     NTseqce = NTseqce1
                                    
#                                 AAseque1 = AA.split(',')[9]
#                                 if len(AAseque1)-2>100 :
#                                     AAseque2 = AAseque1[:100]
#                                     AAseque3 = AAseque1[100:]
                                    
#                                     # Now concatenate
#                                     AAseque = f"{AAseque2}\n{AAseque3}"
                                    
#                                 else :
#                                     AAseque = AAseque1
                                
#                                 AAseqce1 = AA.split(',')[9]
#                                 if len(AAseqce1)-2>100 :
#                                     AAseqce2 = AAseqce1[:100]
#                                     AAseqce3 = AAseqce1[100:]
                                    
#                                     # Now concatenate
#                                     AAseqce = f"{AAseqce2}\n{AAseqce3}"
                                    
#                                 else :
#                                     AAseqce = AAseqce1
                                
#                                 comb = f"{gene_name} {AAseque} {len(AAseque)} {NTseque} {len(NTseque)-2} {AAseqce} {len(AAseqce)} {NTseqce} {len(NTseqce)-2} {AAseq} {len(AAseq)} {NTseq} {len(NTseq)-2}".split()
#                                 gene_name2 = f"Gene : {comb[0]}"
                                
#                                 # Also get length of each range
#                                 len1 = AA.split(',')[12]
#                                 len2 = AA.split(',')[13]
#                                 sz1 = len2 - len1
#                                 len1 = AA.split(',')[14]
#                                 len2 = AA.split(',')[15]
#                                 sz2 = len2 - len1
                                
#                                 # Also add UEX_CE coordinates
#                                 comb = f"{AA} {sz1} {sz2}".split()
#                                 JuncSpanning_Coord = f"{comb[11]}:{comb[12]}-{comb[13]}:{comb[16]}\n{comb[11]}:{comb[14]}-{comb[15]}:{comb[17]}"
                                
#                             elif args[4]==1 :
#                                 NTseque = AA.split(',')[7]
#                                 NTseqce = AA.split(',')[8]
#                                 AAseque = AA.split(',')[9]
#                                 NTseqce = AA.split(',')[10]
                                
#                                 comb = f"{gene_name} {AAseque} {len(AAseque)} {NTseque} {len(NTseque)-2} {AAseqce} {len(AAseqce)} {NTseqce} {len(NTseqce)-2} {AAseq} {len(AAseq)} {NTseq} {len(NTseq)-2}".split()
#                                 gene_name2 = f"Gene : {comb[0]}"
                                
#                                 # Also add UEX_CE coordinates
#                                 # Also get length of each range
#                                 len1 = AA.split(',')[12]
#                                 len2 = AA.split(',')[13]
#                                 sz1 = len2 - len1
#                                 len1 = AA.split(',')[14]
#                                 len2 = AA.split(',')[15]
#                                 sz2 = len2 - len1
                                
#                                 comb = f"{AA} {sz1} {sz2}".split()
#                                 JuncSpanning_Coord = f"{comb[11]}:{comb[12]}-{comb[13]}:{comb[16]}\n{comb[11]}:{comb[14]}-{comb[15]}:{comb[17]}"
                                
#                             else :
#                                 NTseque = AA.split(',')[7]
#                                 NTseqce = AA.split(',')[8]
#                                 AAseque = AA.split(',')[9]
#                                 NTseqce = AA.split(',')[10]
                                
#                                 comb = f"{gene_name} {AAseque} {len(AAseque)} {NTseque} {len(NTseque)-2} {AAseqce} {len(AAseqce)} {NTseqce} {len(NTseqce)-2} {AAseq} {len(AAseq)} {NTseq} {len(NTseq)-2}".split()
#                                 gene_name2 = f"Gene : {comb[0]}"
                                
#                                 # Also get length of each range
#                                 len1 = AA.split(',')[12]
#                                 len2 = AA.split(',')[13]
#                                 sz1 = len2 - len1
#                                 len1 = AA.split(',')[14]
#                                 len2 = AA.split(',')[15]
#                                 sz2 = len2 - len1
                                
#                                 # Also add UEX_CE coordinates
#                                 comb = f"{AA} {sz1} {sz2}".split()
#                                 JuncSpanning_Coord = f"{comb[11]}:{comb[12]}-{comb[13]}:{comb[16]}\n{comb[11]}:{comb[14]}-{comb[15]}:{comb[17]}"
                                
#                             ##### END PEAKS EVENT
#                             exon1 = 0
#                             exon2 = 0
                            
#                             # Also get actual event identified
#                             event_identified = f"None -{line1.iloc[:,1]}-{line1.iloc[:,2]}"
#                             print("processing event num {i} and event {fn}")
                            
#                             ### Should remove this
#                             command = [
#                                 "./ggsashimi_txV3.py",
#                                 "-A", "median_j",
#                                 "-b", "all_bams.tsv",
#                                 "-c", line,
#                                 "-g", "Homo_sapiens.GRCh38.103.chr.sorted_new.gtf",
#                                 "PeaksJuncSpanningCoord", JuncSpanning_Coord,
#                                 "PeaksAA", AAseq,
#                                 "PeaksFlg", 2,
#                                 "PeaksTx", peaks_event,
#                                 "-GeneName", gene_name2,
#                                 "-MajiqStrnd", strnd,
#                                 "-ORIG", "1",
#                                 "-UEX", exon1,
#                                 "-DEX", exon2,
#                                 "-MajiqTx", majiq_event,
#                                 "PGPTx", event_identified,
#                                 "-Majiq", fn,
#                                 "-Tx", TxID,
#                                 "-M", "1",
#                                 "-C", "3",
#                                 "-o", f"{folder}/sashimi_plots/{eventyp}/{fn}",
#                                 "-O", "3",
#                                 "--alpha", "0.25",
#                                 "--base-size=20",
#                                 "--ann-height=4",
#                                 "--height=3",
#                                 "--width=18",
#                                 "-P", "palette.txt"
#                             ]

#                             result = subprocess.run(command)
#                             print(result.stdout)
                            
#                             if os.path.exists(f"{folder}/sashimi_plots/{eventyp}/{fn}.svg") :
#                                 os.remove(f"{folder}/sashimi_plots/{eventyp}/{fn}.svg")
                            
#                         # Now merge all pdf's
#                         command = [
#                             "python", "merge_sashimis.py", f"{folder}/sashimi_plots/{eventyp}/"
#                         ]

#                         result = subprocess.run(command)
#                         print(result.stdout)
                        
#             ########## FIN run_sashimiV1.sh codé en Pyhton ##########