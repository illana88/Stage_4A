import os
import glob

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
    if len(os.listdir("iPSC_gtfs"))!=0 :
        for f in glob.glob("*.*") :
            os.remove(f)
        os.makedirs("iPSC_gtfs",exist_ok=True)
#     file = open("all_bams.tsv")
#     reader = csv.reader(file,delimiter='\t')