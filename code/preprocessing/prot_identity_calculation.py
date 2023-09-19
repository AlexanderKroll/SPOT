from Bio.Emboss.Applications import NeedleCommandline
import os
from os.path import join
import pandas as pd
import sys
import time
import numpy as np


arg = int(sys.argv[1])

CURRENT_DIR = "/gpfs/scratch/alkro105/Fasta_files/"
    
def calculate_identity(fasta_file_1, fasta_file_2):
    needle_cline = NeedleCommandline(asequence = fasta_file_1, bsequence = fasta_file_2,
                                     gapopen=10, gapextend=0.5,  filter = True)

    out = needle_cline()[0]
    out = out[out.find("Identity"):]
    out = out[:out.find("\n")]
    percent = float(out[out.find("(")+1 :out.find(")")-1].replace(" ", ""))
    return(percent)

finished_files = os.listdir(join("/gpfs/scratch/alkro105/", "Transporter_classes_ident"))
if not "test_seq" + str(arg) + ".txt" in finished_files:

    identities = []
    for i in range(13710):
        print(arg, i)
        ident = calculate_identity(fasta_file_1 = join(CURRENT_DIR, "test_seq_" + str(arg) + ".fasta"),
                   fasta_file_2 = join(CURRENT_DIR, "train_seq_" + str(i) + ".fasta"))
        identities.append(ident)


    ofile = open(join("/gpfs/scratch/alkro105/", "Transporter_classes_ident", "test_seq" + str(arg) + ".txt"), "w")
    ofile.write(str(np.argmax(identities)) + " " + str(max(identities)))
    ofile.close()
    


