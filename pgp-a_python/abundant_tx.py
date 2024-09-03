import os
import pandas as pd

if os.access("principal_tx.csv",os.F_OK):
    # Delete file if it exists
    os.remove(("principal_tx.csv"))
    
# Create a list of the files from your target directory
file_list = [f for f in os.listdir("iPSC_gtfs/") if f.endswith('.csv')]

ddf = pd.DataFrame({
    'TxID': pd.Series(dtype='str'),
    'GeneID': pd.Series(dtype='str'),
    'Gene_Name': pd.Series(dtype='str'),
    'cov': pd.Series(dtype='float'),
    'FPKM': pd.Series(dtype='float'),
    'TPM': pd.Series(dtype='float')
})

all_ddf = pd.DataFrame({
    'TxID': pd.Series(dtype='str'),
    'GeneID': pd.Series(dtype='str'),
    'Gene_Name': pd.Series(dtype='str'),
    'cov': pd.Series(dtype='float'),
    'FPKM': pd.Series(dtype='float'),
    'TPM': pd.Series(dtype='float')
})

for i in range(len(file_list)):
    print(f"reading file : {file_list[i]}")
    rec = pd.read_csv(f"iPSC_gtfs/{file_list[i]}", sep="\\s+", header=None, dtype=str)
    rec.columns = ["TxID","GeneID","Gene_Name","cov","FPKM","TPM"]
    
    rec = rec.astype({
        'TxID': 'str',
        'GeneID': 'str',
        'Gene_Name': 'str',
        'cov': 'float',
        'FPKM': 'float',
        'TPM': 'float'
    })
    
    all_ddf = pd.concat([all_ddf, rec], ignore_index=True)

# And sort for fast retrieval
all_ddf = all_ddf.sort_values(by='Gene_Name')

# Get unique gene names
unique_genes = pd.unique(all_ddf['Gene_Name'])

Tx_ddf = pd.DataFrame({
    'TxID': pd.Series(dtype='str'),
    'GeneID': pd.Series(dtype='str'),
    'Gene_Name': pd.Series(dtype='str')
})

# Select row with max cov for each gene
print('Now Generating Txs Table, Will take a while !!!!!!!')

for i in range (len(unique_genes)):
    TxSubset = all_ddf[all_ddf['Gene_Name'] == unique_genes[i]]
    tx_gene = TxSubset[TxSubset['cov'] == TxSubset['cov'].max()]
    tx_gene = tx_gene.iloc[:, :3]
    Tx_ddf = pd.concat([Tx_ddf, tx_gene], ignore_index=True)

# Write in file
# Also remove trailing version numbers
Tx_ddf['GeneID'] = Tx_ddf['GeneID'].str.split('.', n=1, expand=True)[0]
Tx_ddf['TxID'] = Tx_ddf['TxID'].str.split('.', n=1, expand=True)[0]
Tx_ddf.to_csv("principal_txs.csv", index=False, sep=',', header=False)
print('Done With Txs Table: principal_txs.csv')