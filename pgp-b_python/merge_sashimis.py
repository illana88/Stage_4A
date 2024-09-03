from PyPDF2 import PdfFileMerger
import os
import sys

#Create and instance of PdfFileMerger() class
merger = PdfFileMerger()

#Create a list with PDF file names
path_to_files = sys.argv[1] #r'pdf_files/'

#Get the file names in the directory
for root, dirs, file_names in os.walk(path_to_files):
    #Iterate over the list of file names
    for file_name in file_names:
        #Append PDF files
        file_path = os.path.join(root, file_name)
        # Append PDF files using a context manager to ensure files are closed
        with open(file_path, 'rb') as pdf_file:
            merger.append(pdf_file)

#Write out the merged PDF
with open(os.path.join(path_to_files, "all_sashimis.pdf"), 'wb') as output_pdf:
    merger.write(output_pdf)
    
merger.close()