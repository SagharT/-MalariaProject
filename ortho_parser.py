#!/usr/bin/env python3
"""
Script name: ortho_parser.py
Author: Defne YANARTAS 
Created: 2023-03-01

Usage: python ortho_parser.py 
        Please edit the file list manually in the script file!!!!!
        Also make sure that the ortho file name in the script matches 
        your tsv output.
    
   
Description: The script takes the proteinortho output, and parses it to fetch 
        the gene_ids that are orthologoues and one copy for every species 
        and /or double or none for outgroup Tg. And those selected orthologue
        gene sequences are fetcged from the correspinding species files and 
        written to a new file.
        
Input: A list of fasta files with genes
Output: fasta files as many as orthologues, they are numerically named 
        with ".fasta" extension. 
        Orthologues that had double genes for the Tg is labeled "double_tg"
        Orthologues that had no genes for the Tg is labeled "no_tg"
Procedure:
    1. Make a list of species files.
    2. Use that list to loop in the files, get the sequences and gene ids, 
    and make one dictionary containing all the gene ids and their sequences. 
    The gene ids are renamed as the species abbreviation is appended to it.
    3. Loop in the proteinortho output file and fetch the orthologoues of 
    interest, use that information to fetch the corresponding sequences from 
    sequence dictionary, generate a fasta file for each orthologue.

              
Quality Control:
    - Some random orthologues have been checked if the output file contents 
    correctly match the original files.
    
"""
#%% 1. Make a list of species files.
faa_files=["Ht_new.faa","Pb_new.faa","Pc_new.faa","Pf_new.faa",
           "Pk_new.faa","Pv_new.faa","Py_new.faa","Tg_new.faa"]

                
#%% 2. Use that list to loop in the files, get the sequences and gene ids, 
#and make one dictionary containing all the gene ids and their sequences. 
#The gene ids are renamed as the species abbreviation is appended to it.

seq_dict={}
for i in faa_files:
    with open (i,"r") as read:
        for line in read:
            if line.startswith(">"):
                id=line.strip(">").rstrip()
                seq=next(read).rstrip()
                species=i.replace("_new.faa","")
                seq_dict[f'>{id}_{species}']=seq

#%%  3. Loop in the proteinortho output file and fetch the orthologoues of 
#interest, use that information to fetch the corresponding sequences from 
#sequence dictionary, generate a fasta file for each orthologue.

ortho_file="myproject.proteinortho.tsv"             

with open(ortho_file ,"r") as ortho:                                            # get the output file                                                                             
    header=ortho.readline().split("\t")                                         # the header is made into a list so we can fetch the species name each time we are at a orthologue
    counter=1                                                                   # counter will serve as the name for output fasta files
    for line in ortho:                                                          # looping in the ortho output starting after the header
        l=line.rstrip().split("\t")                                             # make a list of the items in the line spearated by tab.
        if l[0]=="8" and l[1]=="8":                                             # if the number of species and genes match ( every orthologue is represented once in every species)
            with open(f'{counter}.fasta',"w") as write:                         # then we want that orthologue to be written to a file, the naming follows the latest number not taken
                for i in range(3,11):                                           # I use numbers to loop in the line because then I can use the same index to refer to the species name in the header
                    id=l[i].rstrip()                                            # get the current index if the line (which is a gene id)  get rid of newlines 
                    species=header[i].replace("_new.faa","").rstrip()           # get the corresponding column name from header (whihc is species name)
                    print(f">{id}_{species}\n{seq_dict[f'>{id}_{species}']}",   
                          file=write)                                           # write the gene id and species name in the desired format as fasta header and fetch the sequence from the sequence dictionary and write that to the next line
            counter+=1                                                          # if we have written a file, we need to change the counter so the next fasta file will have a different name
        elif l[0]=="8" and l[1]=="9":                                           # this bit includes the orthologues that the outgroup Tg gave double hits, in that case if we decide to include those, we can fetch only the first hit.
            with open(f'{counter}.fasta',"w") as write:               # rest is almost the same with the first condition of this loop.
                for i in range(3,11):                                           
                    if i!=10:                                                   # operation is the same until the outgroup column
                        id=l[i].rstrip()
                    elif i==10:                                                 # in the outgroup, we only pick the first gene id (genes were comma separated)
                        id=l[i].split(",")[0]
                    species=header[i].replace("_new.faa","").rstrip()           # rest is the same
                    print(f">{id}_{species}\n{seq_dict[f'>{id}_{species}']}",
                          file=write)
            counter+=1
'''
        elif l[0]=="7" and l[1]=="7" and "*" in l[10]:                          # we can also choose to get the orthologues that the out group gave no hits, even though we may not use it for phylogeny, it might be interesting to have them
            with open(f'{counter}_no_tg.fasta',"w") as write:
                for i in range(3,10):                                           # we loop until before the outgroup column
                    id=l[i].rstrip()
                    species=header[i].replace("_new.faa","").rstrip()
                    print(f">{id}_{species}\n{seq_dict[f'>{id}_{species}']}",
                          file=write) 
            counter+=1
                    
'''
            


