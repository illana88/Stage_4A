import pandas as pd
import numpy as np
import os
import sys

args = sys.argv[1:]

csv_record = pd.read_csv(args[0], sep='\t', header=None, dtype=str)
ce_coord = pd.read_csv(args[1], sep='\t', header=None, dtype=str)
intron_coord = pd.read_csv(args[3], sep='\t', header=None, dtype=str)
folder = args[0].split('/')[0]
pcnt = -1*float(args[2])

def pct(x) :
    return np.diff(x)/x[:-1]

def pct1(x) :
    return np.diff(x)/x[1:] # For - strands scanning reverse

if os.path.exists(os.path.join(folder, "ce_inclusion_coord.bed")) :
    os.remove(os.path.join(folder, "ce_inclusion_coord.bed"))
    
if os.path.exists(os.path.join(folder, "/ce_extension_coord.bed")) :
    os.remove(os.path.join(folder, "/ce_extension_coord.bed"))

if os.path.exists(os.path.join(folder, "/coverage_file_not_found.txt")) :
    os.remove(os.path.join(folder, "/coverage_file_not_found.txt"))

if os.path.exists(os.path.join(folder, "/coverage_file_not_found.bed")) :
    os.remove(os.path.join(folder, "/coverage_file_not_found.bed"))

if os.path.exists(os.path.join(folder, "/skipped_ce.csv")) :
    os.remove(os.path.join(folder, "/skipped_ce.csv"))

if os.path.exists(os.path.join(folder, "/IGV_R_returned_ce_inclusion.csv")) :
    os.remove(os.path.join(folder, "/IGV_R_returned_ce_inclusion.csv"))

if os.path.exists(os.path.join(folder, "/IGV_coverage_file_notfound_ce.csv")) :
    os.remove(os.path.join(folder, "/IGV_coverage_file_notfound_ce.csv"))

if os.path.exists(os.path.join(folder, "/IGV_R_returned_ce_extension.csv")) :
    os.remove(os.path.join(folder, "/IGV_R_returned_ce_extension.csv"))

if os.path.exists(os.path.join(folder, "/IGV_R_ce1_ce2.csv")) :
    os.remove(os.path.join(folder, "/IGV_R_ce1_ce2.csv"))

if os.path.exists(os.path.join(folder, "/ce_inclusion_coord_sashimi.bed")) :
    os.remove(os.path.join(folder, "/ce_inclusion_coord_sashimi.bed"))

if os.path.exists(os.path.join(folder, "/ce_extension_coord_sashimi.bed")) :
    os.remove(os.path.join(folder, "/ce_extension_coord_sashimi.bed"))

if os.path.exists(os.path.join(folder, "/IGV_R_returned_IR.csv")) :
    os.remove(os.path.join(folder, "/IGV_R_returned_IR.csv"))

if os.path.exists(os.path.join(folder, "/IR_coord_sashimi.bed")) :
    os.remove(os.path.join(folder, "/IR_coord_sashimi.bed"))

if os.path.exists(os.path.join(folder, "/IR_coord.bed")) :
    os.remove(os.path.join(folder, "/IR_coord.bed"))

if os.path.exists(os.path.join(folder, "/skipped_ce_sashimi.csv")) :
    os.remove(os.path.join(folder, "/skipped_ce_sashimi.csv"))

last_genename = ""
k = 0
IRFlag = 0
scan_flag = 0 # 0 for start<end, 1 for start>end
start_cov = 0
end_cov = 0

ce_coord[1] = pd.to_numeric(ce_coord[1])
ce_coord[2] = pd.to_numeric(ce_coord[2])

# Read every 3rd entry starting from 2
cnt_ce = 0
cnt_ir = 0
cnt_ext = 0
cnt_f = 0
cnt_cf = 0
ii = 0 # For csv_record

for i in range(2, len(ce_coord)+1, 3) :
    ii = ii + 1
    gene_name = ce_coord.iloc[i,6].tolist()[0]
    chr_n = ce_coord.iloc[i,0].tolist()[0]
    start = ce_coord.iloc[i,1].tolist()[0]
    end = ce_coord.iloc[i,2].tolist()[0]
    strnd = ce_coord.iloc[i,5].tolist()[0]
    
    star_int = intron_coord.iloc[ii,1].tolist()[0]
    end_int = intron_coord.iloc[ii,2].tolist()[0]
    
    # Also get intron start and end
    if strnd=='+' :
        cstart = ce_coord.iloc[i-1,2].tolist()[0]
        cend = ce_coord.iloc[i+1,1].tolist()[0]
        
    else :
        cstart = ce_coord.iloc[i+1,2].tolist()[0]
        cend = ce_coord.iloc[i-1,1].tolist()[0]
    
    if start<end :
        scan_flag = 0
    
    else :
        scan_flag = 1
    
    k = 1
    
    # NOW COVERAGE FILE NAME FOR EACH EVENT IS CHANGED TO chr_start_end_gene_name to make each cov file unique
    fname = f"coverages/{chr_n}_{star_int}_{end_int}_{gene_name}.cov.bed"
    
    # Now read coverage file
    if os.path.exists(fname) :
        CovData = pd.read_csv(fname, sep='\t', header=None)
        
        # Now get ce coverage data for the two ends
        if scan_flag==0 :
            data_val = CovData.iloc[:,8].tolist()[0]
            
            # Change start_cov position as we are scanning whole intron now
            begt = pd.to_numeric(CovData[1])
            beg = begt.iloc[0].tolist()[0]
            offv = float(start) - float(beg)
            start_cov = data_val[offv]
            data_val1 = data_val[offv-1:]
            
            end_cov_index = np.min(np.where(pct(data_val1) <= pcnt)[0]) # Same as befor - just used -sign with pcnt
            
            ce1 = start
            if end_cov_index==np.inf :
                # NOW CHECK FOR WHOLE INTRON COVERAGE
                # FIRST CHECK FOR THE COVERAGES ON THE RIGHT TILL NEXT EXON
                data_val1 = data_val[offv-1::-1] # Get data on the right
                end_cov_index1 = np.min(np.where(pct(data_val1) <= pcnt)[0]) # Same as befor - just used -sign with pcnt
                
                if end_cov_index1==np.inf : # THIS IS IR EVENT
                    if strnd=='+' :
                        ce1 = ce_coord.iloc[i-1,2].tolist()[0] + 1
                        ce2 = ce_coord.iloc[i+1,1].tolist()[0] - 1
                    
                    else :
                        ce1 = ce_coord.iloc[i-1,1].tolist()[0] - 1
                        ce2 = ce_coord.iloc[i+1,2].tolist()[0] + 1
                    
                    IRFlag = 2
                    
                else : # CE_EXTENSION
                    if strnd=='+' :
                        ce2 = ce_coord.iloc[i+1,1].tolist()[0] - 1
                    
                    else :
                        ce2 = start
                        ce1 = ce_coord.iloc[i-1,1].tolist()[0] - 1
                    
                    IRFlag = 1
            else :
                IRFlag = 0
                
                # Changed due to whole intron
                ce2 = ce1 + end_cov_index - 1
            
            if ce1>ce2 : # Avoid start>end
                tempt = ce1
                ce1 = ce2
                ce2 = tempt
            
        else :
            # Here we first reverse the dataframe w.r.to rows
            # Transpose of dataframe
            transpose = CovData.T
            
            # Converting the result to dataframe
            transpose = pd.DataFrame(transpose)
            
            # Calculating reverse of dataframe
            rev_data_frame = transpose.iloc[::-1]
            
            # Transpose of reverse dataframe
            rev_data_frame = rev_data_frame.T
            
            # Converting the result to dataframe
            rev_data_frame = pd.DataFrame(rev_data_frame)
            rev_data_frame[7] = rev_data_frame[7].astype(float)
            rev_data_frame[8] = rev_data_frame[8].astype(float)
            
            data_val = rev_data_frame.iloc[:,8].tolist()[0]
            
            # Change start_cov position as we are scanning whole intron now
            endt = float(rev_data_frame.iloc[:,2]) # Get end position
            et = endt[0]
            offv = et - float(start)
            start_cov = data_val[offv]
            data_val1 = data_val[offv-1:]
            
            end_cov_index = np.min(np.where(pct(data_val1) <= pcnt)[0]) # New criterion
            ce2 = start
            
            if end_cov_index==np.inf :
                # NOW CHECK FOR WHOLE INTRON COVERAGE
                # NOTE HERE WE ARE SCANNING IN THE FORWARD DIRECTION, SO USING pct instead of pct1
                data_val1 = data_val[offv-1:]
                end_cov_index1 = np.min(np.where(pct(data_val1) <= pcnt)[0]) # Same as befor - just used -sign with pcnt
                
                if end_cov_index1==np.inf : # THIS IS IR EVENT
                    IRFlag = 2
                    
                    if strnd=='+' :
                        ce1 = ce_coord.iloc[i-1,2].tolist()[0] + 1
                        ce2 = ce_coord.iloc[i+1,1].tolist()[0] - 1
                    
                    else :
                        ce1 = ce_coord.iloc[i+1,2].tolist()[0] + 1
                        ce2 = ce_coord.iloc[i-1,1].tolist()[0] - 1
                    
                else : # CE_EXTENSION
                    if strnd=='+' :
                        ce2 = start
                        ce1 = ce_coord.iloc[i-1,2].tolist()[0] + 1
                    
                    else :
                        ce1 = ce_coord.iloc[i+1,2].tolist()[0] + 1
                        ce2 = start
                    
                    IRFlag = 1
                    
            else :
                IRFlag = 0
                
                if strnd=='-' :
                    ce1 = float(ce2) - end_cov_index + 2
                
                else :
                    ce1 = float(ce2) - end_cov_index + 4
                
            if ce1>ce2 : # Avoid star>end
                tempt = ce1
                ce1 = ce2
                ce2 = tempt
            
        # Save bedfile for current junction
        if (IRFlag==0 and ce1!=ce2) : # Avoid ce1=ce2
            cnt_ce = cnt_ce + 1
            
            SJdat = pd.DataFrame({
                'seqnames': [ce_coord.iloc[i,0]],
                'start': [ce1-1], # No idea why -1
                'end': [ce2],
                'name': [1], # Reverted back ce_coord[i,]$V4, #1, #changed to exon size
                'score': [0],
                'strand': [ce_coord.iloc[i,5]],
                'gene_name': [ce_coord.iloc[i,6]]
            })
            
            # Also new file with 5',ce and 3'
            recrd_1 = ce_coord.iloc[i-1,:].tolist()[0]
            recrd_1[3]=1 # Change exon length to 1 - causes trouble
            
            recrd_1.to_csv(f"{folder}/ce_inclusion_coord.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            SJdat.to_csv(f"{folder}/ce_inclusion_coord.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            
            recrd_2 = ce_coord.iloc[i+1,:].tolist()[0]
            recrd_2[3] = 1
            
            recrd_2.to_csv(f"{folder}/ce_inclusion_coord.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            
            # Also write bed file for sashimi plots including TxID
            c1 = ce_coord.iloc[i-1,:].tolist()[0]
            c1[7] = csv_record.iloc[ii,0].split(',')[8]
            SJdat["TxID"] = csv_record.iloc[ii,0].split(',')[8]
            c2 = ce_coord.iloc[i+1,:].tolist()[0]
            c2[7] = csv_record.iloc[ii,0].split(',')[8]
            
            c1.to_csv(f"{folder}/ce_inclusion_coord_sashimi.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            SJdat.to_csv(f"{folder}/ce_inclusion_coord_sashimi.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            c2.to_csv(f"{folder}/ce_inclusion_coord_sashimi.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            
            # Also write csv file
            csv_record.iloc[ii,:].to_csv(f"{folder}/IGV_R_returned_ce_inclusion.csv", sep='\t', mode='a', header=False, index=False, quoting=3)
        
        elif (IRFlag==1 and ce1!=ce2) : # Extension event
            cnt_ext = cnt_ext + 1
            
            SJdat = pd.DataFrame({
                'seqnames': [ce_coord.iloc[i,0]],
                'start': [ce1],
                'end': [ce2-1], # To match with DNM1 - for exton extension events
                'name': [1], #ce_coord[i,]$V4, #1,
                'score': [0],
                'strand': [ce_coord.iloc[i,5]],
                'gene_name': [ce_coord.iloc[i,6]]
            })
            
            # Also new file with 5',ce and 3'
            recrd_1 = ce_coord.iloc[i-1,:].tolist()[0]
            recrd_1[3]=1
            
            recrd_1.to_csv(f"{folder}/ce_extension_coord.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            SJdat.to_csv(f"{folder}/ce_extension_coord.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            
            recrd_2 = ce_coord.iloc[i+1,:].tolist()[0]
            recrd_2[3] = 1
            
            recrd_2.to_csv(f"{folder}/ce_extension_coord.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            
            # Also write bed file for sashimi plots including TxID
            c1 = ce_coord.iloc[i-1,:].tolist()[0]
            c1[7] = csv_record.iloc[ii,0].split(',')[8]
            SJdat["TxID"] = csv_record.iloc[ii,0].split(',')[8]
            c2 = ce_coord.iloc[i+1,:].tolist()[0]
            c2[7] = csv_record.iloc[ii,0].split(',')[8]
            
            c1.to_csv(f"{folder}/ce_extension_coord_sashimi.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            SJdat.to_csv(f"{folder}/ce_extension_coord_sashimi.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            c2.to_csv(f"{folder}/ce_extension_coord_sashimi.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            
            # Also write csv file
            csv_record.iloc[ii,:].to_csv(f"{folder}/IGV_R_returned_ce_extension.csv", sep='\t', mode='a', header=False, index=False, quoting=3)
            
        elif (IRFlag==2 and ce1!=ce2) : # IR events
            cnt_ir = cnt_ir + 1
            
            SJdat = pd.DataFrame({
                'seqnames': [ce_coord.iloc[i,0]],
                'start': [ce1],
                'end': [ce2-1], # To match with DNM1 - for exton extension events
                'name': [1], #ce_coord[i,]$V4, #1,
                'score': [0],
                'strand': [ce_coord.iloc[i,5]],
                'gene_name': [ce_coord.iloc[i,6]]
            })
            
            # Also new file with 5',ce and 3'
            SJdat.to_csv(f"{folder}/IR_coord.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            
            # Also write bed file for sashimi plots including TxID
            SJdat["TxID"] = csv_record.iloc[ii,0].split(',')[8]
            SJdat.to_csv(f"{folder}/IR_coord_sashimi.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
            
            # Also write csv file
            csv_record.iloc[ii,:].to_csv(f"{folder}/IGV_R_returned_IR.csv", sep='\t', mode='a', header=False, index=False, quoting=3)
            
        else :
            cnt_f = cnt_f + 1
            print(f"ce1==ce2 {ce1}{ce2} for gene: {gene_name} so skipping for ce and events are in folder/skipped_ce.csv")
            
            SJdat = pd.DataFrame({
                'seqnames': [ce_coord.iloc[i,0]],
                'start': [ce1],
                'end': [ce2],
                'name': [1], #ce_coord[i,]$V4, #1,
                'score': [0],
                'strand': [ce_coord.iloc[i,5]],
                'gene_name': [ce_coord.iloc[i,6]]
            })
            
            # Also new file with 5',ce and 3'
            recrd_1 = ce_coord.iloc[i-1,:].tolist()[0]
            recrd_1[3]=1
            
            recrd_1.to_csv(f"{folder}/skipped_ce.csv", sep=',', mode='a', header=False, index=False, quoting=3)
            SJdat.to_csv(f"{folder}/skipped_ce.csv", sep=',', mode='a', header=False, index=False, quoting=3)
            
            recrd_2 = ce_coord.iloc[i+1,:].tolist()[0]
            recrd_2[3] = 1
            
            recrd_2.to_csv(f"{folder}/skipped_ce.csv", sep=',', mode='a', header=False, index=False, quoting=3)
            
            # Also write bed file for sashimi plots including TxID
            c1 = ce_coord.iloc[i-1,:].tolist()[0]
            c1[7] = csv_record.iloc[ii,0].split(',')[8]
            SJdat["TxID"] = csv_record.iloc[ii,0].split(',')[8]
            c2 = ce_coord.iloc[i+1,:].tolist()[0]
            c2[7] = csv_record.iloc[ii,0].split(',')[8]
            
            c1.to_csv(f"{folder}/skipped_ce_sashimi.csv", sep='\t', mode='a', header=False, index=False, quoting=3)
            SJdat.to_csv(f"{folder}/skipped_ce_sashimi.csv", sep='\t', mode='a', header=False, index=False, quoting=3)
            c2.to_csv(f"{folder}/skipped_ce_sashimi.csv", sep='\t', mode='a', header=False, index=False, quoting=3)
            
            csv_record.iloc[ii,:].to_csv(f"{folder}/IGV_R_ce1_ce2.csv", sep=',', mode='a', header=False, index=False, quoting=3)
            
    else :
        cnt_cf = cnt_cf + 1
        print(f"coverage file: {fname} does not exists, so skipping")
        
        with open(f"{folder}/coverage_file_not_found.txt", "a") as f :
            f.write(f"coverage file: {fname} does not exists, so skipping")
        
        # Also write bed file to generate coverage files for missing junctions
        # Also new file with 5',ce and 3'
        ce_coord.iloc[i,:].to_csv(f"{folder}/coverage_file_not_found.bed", sep='\t', mode='a', header=False, index=False, quoting=3)
        
        # Also write csv file
        csv_record.iloc[ii,:].to_csv(f"{folder}/IGV_coverage_file_notfound_ce.csv", sep=',', mode='a', header=False, index=False, quoting=3)
        
print(f"total events read are: {ce_coord.shape[0]/3}")

with open(f"{folder}/Summary_stats.txt", "a") as f :
    f.write(f"Finished R session for CE_boundary calculations of a total of: {ce_coord.shape[0]/3} successful events, please see folder/IGV_ce_inclusion.csv\n")

print(f"Total events for which ce_coverage files not found are : {cnt_cf} - please see folder/coverage_file_not_found.bed and folder/IGV_coverage_file_notfound_ce.csv")

with open(f"{folder}/Summary_stats.txt", "a") as f :
    f.write(f"Total events for which ce_coverage files not found are : {cnt_cf} (please see folder/coverage_file_not_found.bed)\n")

print(f"Total events unsuccessful (due to ce1==ce2) are: {cnt_f} - Please see file folder/skipped_ce.csv file")

with open(f"{folder}/Summary_stats.txt", "a") as f :
    f.write(f"Total events unsuccessful (due to ce1==ce2) are: {cnt_f} Please see file folder/skipped_ce.csv file. It contains us, ce and ds exon coordinates for each event\n")

print(f"Total events successfully processed are: {cnt_ce+cnt_ir+cnt_ext}")

with open(f"{folder}/Summary_stats.txt", "a") as f :
    f.write(f"Total events successfully processed are: {cnt_ce+cnt_ir+cnt_ext}\n")

print(f"Out of total {cnt_ce+cnt_ir+cnt_ext} successfull events processed, ce_inclusion events are: {cnt_ce} ce_extension events are: {cnt_ext} and IR events are: {cnt_ir}")

with open(f"{folder}/Summary_stats.txt", "a") as f :
    f.write(f"Out of total: {cnt_ce+cnt_ir+cnt_ext} successfull events processed\n")
    f.write(f"            ce_inclusion events are: {cnt_ce}\n")
    f.write(f"            ce_extension events are: {cnt_ext}\n")
    f.write(f"            IR events are: {cnt_ir}\n")