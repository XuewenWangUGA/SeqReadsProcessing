#!/usr/bin/python
import sys
import gffutils

# pip install gffutils

'''
FUNCTION:
This script will take the gff3 annotation fiel and filter/remove out short CDS of any gene based on given length in bp, report to a new gff3 file with longer (>=) CDS
Parameters: outfileName in.gff3 cdslength_cutoff
Version 2024-May 1st
# use python version 3
'''
__author__ ="Xuewen Wang"

usage=""" python3 gff3_cds_filter.py OUT_gff3.txt IN.gff3 CdsMinLength
          e.g.: python3 gff3_cds_filter.py filtered.gff3 Urhy.gene.gff3 300
          """
print(usage)

outf=sys.argv[1] #'filtered.gff3'
inf=sys.argv[2] #'Urhy.gene.gff3'
minlengthCds=int(sys.argv[3])  #300
db = gffutils.create_db(inf, dbfn='temp.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

# Create a dictionary to store the IDs of the genes with CDS length >= 300
valid_gene_ids = {}
valid_mRNA_ids = {}

# Iterate over all CDS features
for feature in db.all_features(featuretype='CDS'):
    # If the CDS length is >= 300, add the parent gene ID to the dictionary
    if len(feature) >= minlengthCds:
        mRNA = db[feature.attributes['Parent'][0]]
        gene_id = mRNA.attributes['Parent'][0]
        valid_gene_ids[gene_id] = valid_gene_ids.get(gene_id, 0) + len(feature)
        valid_mRNA_ids[mRNA.id] = valid_mRNA_ids.get(mRNA.id, 0) + len(feature)

# Create a list of valid feature types
valid_feature_types = ['gene', 'CDS', 'mRNA', 'exon', 'intron', 'start_codon', 'stop_codon']

with open(outf, 'w') as outfile:
    # Iterate over all features
    for feature in db.all_features(order_by=('seqid', 'start', 'end')):
        # If the feature is a gene and its ID is in the dictionary, write it to the output file
        if feature.featuretype == 'gene' and feature.id in valid_gene_ids:
            outfile.write(str(feature) + "\n")
            # After writing the gene feature, find and write the corresponding mRNA feature
            for mRNA_feature in db.children(feature, featuretype='mRNA'):
                if mRNA_feature.id in valid_mRNA_ids:
                    outfile.write(str(mRNA_feature) + "\n")
                    # After writing the mRNA feature, find and write the corresponding other features
                    for child_feature in db.children(mRNA_feature):
                        if child_feature.featuretype in valid_feature_types:
                            outfile.write(str(child_feature) + "\n")
