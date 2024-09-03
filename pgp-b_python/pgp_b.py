import os
import pandas as pd
import sys
import subprocess

# HOW TO RUN: bash pgp_b_ce_inclusion.sh ce_extension_pgp1.csv COV_CUTOFF NUM_EVENTS_PROCESSED
# EXAMPLE: bash pgp_b_ce_inclusion.sh ce_extension_pgp1.csv .6 1

# PLEASE MAKE SURE TO CONCATENATE
# 1. CE_INCLUSION_FUSED_AA.fasta FOR EACH CUTOFF (WITH 40% FOLLOWED BY NEXT CUTOFF and so on untill all hand_curated events are accounted for)
# 2. ce_inclusion_fused.transeq_in.fasta
# 3. FINAL_CE_INCLUSION_FUSED_AA.fasta
# 4. IGV_unique_ce_inclusion.csv
###########################

# CHECK IF 3 ARGUMENTS ARE PROVIDED
args = ["sorted_selected_events.csv"]

if len(args)==1 :
    if os.path.exists("res_skiptics/Summary.Skiptics.txt"):
        os.remove("res_skiptics/Summary.Skiptics.txt")
    
    os.makedirs("res_skiptics",exist_ok=True)
    
    if os.path.isfile(args[0]):
        splicing_events_file = args[0].split('.')[0]
    
    else :
        print("PLEASE PROVIDE SPLICING_EVENTS csv file  'bash pgp_b.sh splicing_events.csv' and RERUN")
        
        with open("res_skiptics/Summary.Skiptics.txt", "w") as f :
            f.write("PLEASE PROVIDE SPLICING_EVENTS csv file  'bash pgp_b.sh splicing_events.csv' and RERUN\n")
        
        sys.exit(1)
    
    print("################################ STARTED SKIPTICS IDENTIFICATION ")
    
    with open("res_skiptics/Summary.Skiptics.txt", "a") as f :
        f.write("################################ STARTED SKIPTICS IDENTIFICATION \n")
    
    if os.path.exists("EnsDB_tx_not_found.csv"):
        os.remove("EnsDB_tx_not_found.csv")
    
    skiptics_flg = 1
    
    if skiptics_flg==1 :
        print("#############################")
        print(f"INVOKING SKIPTICS Script esV5_layered_CDSV3.sh FOR FILE sorted_{splicing_events_file}.csv from pgp_b.sh")
        
        with open("res_skiptics/Summary.Skiptics.txt", "a") as f :
            f.write(f"INVOKING SKIPTICS Script esV5_layered_CDSV3.sh FOR FILE sorted_{splicing_events_file}.csv from pgp_b.sh\n")
        
        print("#############################") ##### Teste jusqu'ici et tout fonctionne #####
        
        command = [
            "python",
            "esV5_layered_CDSV3.py",
            f"sorted_{splicing_events_file}.csv",
            "principal_txs.csv"
        ]
        
        result = subprocess.run(command, capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)
        
        