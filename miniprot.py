#!/usr/bin/python3
# Author Gaurav Sablok
# Universitat Potsdam
# Date 2024-3-21
def miniprotalignment(alignment, splicetype=None):
    """
    an minprot to the alignment for the 
    paf to allow for the extraction of the
    annotations coming from the mapping of
    the paf file and the gff files. 
    """
    # /home/gaurav/Desktop/week/miniprot/test/aln.gff
    # miniprotalignment("/home/gaurav/Desktop/aln.gff", splicetype="mRNA")
    with open(alignment, "r") as readminprot:
        with open(alignment + "modified.txt", "w") as outminprot:
            for line in readminprot.readlines():
                if "Rank" in line:
                    outminprot.write(line) 
                else:
                    pass 
            outminprot.close()
        readminprot.close()
    with open(alignment + "modified.txt", "r") as alignment:
        reading_mRNA = []
        reading_CDS = []
        reading_stop_codon = []
        for line in alignment.readlines():
            if line.strip().split()[2] == "mRNA":
                reading_mRNA.append([line.strip().split()[0],[line.strip().split()[3],line.strip().split()[4]]])
            elif line.strip().split()[2] == "CDS":
                reading_CDS.append([line.strip().split()[0],[line.strip().split()[3],line.strip().split()[4]]])
            elif line.strip().split()[2] == "stop_codon":
                reading_stop_codon.append([line.strip().split()[0],[line.strip().split()[3],line.strip().split()[4]]])
    if alignment and splicetype == "mRNA":
        return reading_mRNA
    elif alignment and splicetype == "CDS":
        return reading_CDS 
    elif alignment and splicetype == "stop_codon":
        return reading_stop_codon 
    if __name__ == "__main__":
        miniprotalignment(alignmentfile, splicetype="CDS")
