#!/usr/bin/python3
#Author Gaurav Sablok
# Universitat Potsdam
# Date 2024-4-12
import pandas as pd
def generatingAlignments(inputgff, fasta = None, file = None):
    """
    docstring: given a miniprot alignment and the fasta
    it will extract the aligned regions and will make a fasta
    of the same. Only the aligned regions are extracted as a 
    part of the substring.
    usage: generatingAlignments("/home/gaurav/Desktop/sample.gff")
    it optionally takes a fasta and a file name to write the corresponding
    patterns
    @inputgff = gff aligned from the miniprot
    @fasta = fasta sequences from which you want to extract the pattern
    @file = file to which you want to write the sequences of the extracted pattern
    generatingAlignments("/home/gaurav/Desktop/final_code_push/multi.gff", 
                        "/home/gaurav/Desktop/final_code_push/multi.fasta", 
                               "/home/gaurav/Desktop/final_code_push/multiout.fasta")
    """
    readfile = [i for i in open(inputgff, "r").readlines() if "#" not in i]
    with open(inputgff + ".clipped.gff", "w") as writegff:
        writegff.write("col0 \t col1 \t col2 \t col3 \t col4 \t col5 \t col6 \t col7 \t col8 \t col9\n")
        for line in readfile:
            writegff.write(line)
        writegff.close()
    with open(inputgff + ".clipped.gff", "r") as readgff: 
        readdataframe = pd.read_csv(readgff, sep = "\t")
        mRNAannotations = []
        mRNAstart_coordinate = []
        mRNAend_coordinate = []
        for i in range(len(readdataframe.iloc[:,2])):
            if readdataframe.iloc[:,2][i] == "mRNA":
                mRNAannotations.append(readdataframe.iloc[:,2][i])
        for i in range(len(readdataframe.iloc[:,3])):
            if readdataframe.iloc[:,2][i] == "mRNA":
               mRNAstart_coordinate.append(readdataframe.iloc[:,3][i])
        for i in range(len(readdataframe.iloc[:,4])):
            if readdataframe.iloc[:,2][i] == "mRNA":
               mRNAend_coordinate.append(readdataframe.iloc[:,4][i])
        extract_pattern = {}
        for i in range(len(fasta_sequences)):
            extract_pattern[fasta_names[i]] = fasta_sequences[mRNAstart_coordinate[i]:mRNAend_coordinate[i]]
        extractkeys = list(extract_pattern.keys())
        extractvalues = list(extract_pattern.values())
        finalextractvalues = [str(''.join(extractvalues[i])) for i in range(len(list(extractvalues)))]
        with open(file, "w") as fastawrite:
            for i in range(len(extractkeys)):
                fastawrite.write(f">{extractkeys[i]}\n{finalextractvalues[i]}\n")
            fastawrite.close()
