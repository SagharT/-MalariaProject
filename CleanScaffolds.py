#!/user/bin/python3
#make a list of all contigs in the file scaffolds.txt

contigs=[]
with open("scaffolds2.txt") as f:
    for line in f:
        contigs.append(line.strip())
#find this contig in the file gffParse.fna
#make a new file without this contig
with open("Haemoproteus_tartakovskyi.fasta", "r") as f:
    with open("Tartakovskyi_filtered.fasta", "w") as f1:
      #check if all contigs are in the file
      #check if the line starts with >
      #check the third column of the line if it is scaffold or contig 
      #>  length=     scaffold=contig
      for line in f:
            if line.startswith(">"):
                if line[1:].split()[0] not in contigs:
                    f1.write(line) 
                    line=next(f)
                    #print(line) 
                    f1.write(line)     
