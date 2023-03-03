# Malaria study
## Plasmodium data Collection
The files for this study were downloaded from the BINP29 course server under /resources/binp29/Data/malaria and saved in the Malaria folder.  
Because prediction takes time, students of the course did the gene prediction program GeneMark on particular genome files and were stored in a shared folder on the server. 

```bash
mkdir Malaria
cd Malaria
cp -r /resources/binp29/Data/malaria .
cp -r /tmp//Prediction/*.gtf .
```

# Processing of Haemoproteus tartakovskyi data
The H. tartakovskyi genome was found on the server. Because malaria is a parasite, the file contain reads from parasites and birds too. We need to remove bird scaffolds.
This was done by difference in GC content between parasite reads and bird reads.
Haemoproteus tartakovskyi has a lower GC content compared to the host.

An appropriate GC-content 20% and remove scaffolds above
this threshold. Also remove scaffolds that are less than 3000 nucleotides. 
write a  python-script.
Furthermore, all reads shorter than 3000 nucleotides are sorted out as well.
This is done with the python script Scaffold.py. The minimum scaffold length of 3000 was a given, the GC-content was chosen based on the lecture material. Run with


```bash
gunzip Haemoproteus_tartakovskyi.raw.genome.gz
#filter by GC content
python3 removeScaffold.py Haemoproteus_tartakovskyi.raw.genome 26 Haemoproteus_tartakovskyi.fasta 3000
#---->Haemoproteus_tartakovskyi.fasta
```


## Make a gene prediction
You will probably still have some scaffolds that derive from the bird. These should be short. Why?

```bash
#prediction
gmes_petap.pl --min_contig 3000 --Es --sequence Haemoproteus_tartakovskyi.fasta 
#--->genemark.gtf
#remove from firts colum to not get error 
cat genemark.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht2.gff
#--->Ht2.gff
#
cd backup_results
gunzip Ht.blastout.gz
#Create fasta sequences from the gff file 
chmod +x gffParse.pl 
#This get sequence of each gene that predicted
perl gffParse.pl -i Haemoproteus_tartakovskyi.fasta -g Ht2.gff -c -p
#--->gffParse.faa (protein)
#--->gffParse.fna (nt)
#Run blast
#If you don’t want to wait download the Ht.blastp file from the course server.
#blast on protein gffParse.faa 
#nohup will send it to the background
nohup blastp -db SwissProt -query gffParse.faa & 2>&1 
#see the progress
jobs
[1]+  Stopped                 top
[2]-  Done                    nohup blastp -db SwissProt -query gffParse.faa SwissProt -query gffParse.faa &
#---> nohup.out

```
Use the taxonomy for the top hit in the BLAST output to decide if the query sequence orginates from the host he five characters after the underscore is a species abbreviation (what is before the underscore?).

  sp|P14223|ALF_PLAFA Fructose-bisphosphate aldolase OS=Plasmodiu...  560     0.0   
  sp|Q7KQL9|ALF_PLAF7 Fructose-bisphosphate aldolase OS=Plasmodiu...  560     0.0   
  sp|P49577|ALF2_PLABA Fructose-bisphosphate aldolase 2 OS=Plasmo...  545     0.0   
  sp|Q86A67|ALF_DICDI Fructose-bisphosphate aldolase OS=Dictyoste...  422     2e-146
  sp|O65581|ALFC5_ARATH Fructose-bisphosphate aldolase 5, cytosol...  417     6e-145
  
In SwissProt, the five characters before the underscore in a protein identifier correspond to the organism abbreviation. The top hit in the BLAST output (sp|P14223|ALF_PLAFA) corresponds to the fructose-bisphosphate aldolase protein from the organism Plasmodium falciparum.

The script is datParser.py and is found on the server at /resources/binp29/Data/malaria/
```bash

#on my data:
python datParser.py nohup.out gffParse.fna taxonomy.dat uniprot_sprot.dat > scaffolds2.txt

```
I got different result than the one from the server but it makes sence due to different cut off

```bash
cat scaffolds2.txt

> contig00029
> contig00065
> contig00067
> contig00122
> contig00142
> contig00329
> contig00352
> contig00356
> contig00569
> contig00661
> contig01393
> contig01894
> contig01918
> contig02015
```
Now we run python on this file
```bash
python3 CleanScaffolds.py 
#input: Haemoproteus_tartakovskyi.fasta
#---->Tartakovskyi_filtered.fasta

cat Haemoproteus_tartakovskyi.fasta | wc -l
2740
cat Tartakovskyi_filtered.fasta | wc -l
2712
```

## Redo gene Prediction on filtered file
I ran this on filtered file
```bash
gmes_petap.pl --min_contig 3000 --Es --sequence Tartakovskyi_filtered.fasta
#----> genemark.gtf
```
### Fill the table:  
Genome size:  
```bash
gunzip plasmodiumGenomes.tgz
tar -xf plasmodiumGenomes.tar
for file in *.genome; do echo $file; cat $file | grep -v '>' | tr -d '\n\s' | wc -c; done
for file in *.fasta; do echo $file; cat $file | grep -v '>' | tr -d '\n\s' | wc -c; done
```
Genes:
```bash
for file in *.gtf; do echo $file; cat $file | grep 'CDS' | cut -f 9 | cut -d ';' -f 1 | sort -u | wc -l; done
for file in *.gff; do echo $file; cat $file | grep 'CDS' | cut -f 9 | cut -d ';' -f 1 | sort -u | wc -l; done

```
Genomic GC:
```bash
for file in *.genome; do cat $file | grep -v "^>" | tr -d "\n" | awk '{s+=gsub(/[GC]/,"")}END{printf "%.2f\n",s*100/(s+gsub(/[ATCG]/,""))}'; done

cat Tartakovskyi_filtered.fasta | grep -v "^>" | tr -d "\n" | \
awk '{s+=gsub(/[GC]/,"")}END{printf "%.2f\n",s*100/(s+gsub(/[ATCG]/,""))}'
```
## Result table

| # |Species | Host | Genome size | Genes |Genomic GC
| :---   | :---: | :---: | :---: | :---:   | :---: | 
1 | Plasmodium berghei  | rodents| 17954629 |7235 | 23.72| 
2 | Plasmodium cynomolgi  | macaques| 26181343| 5787|40.38 | 
3 | Plasmodium falciparum  | humans | 23270305| 5207|19.36 | 
4 | Plasmodium knowlesi  | lemures | 23462346|4953 |38.83 | 
5 | Plasmodium vivax  |  humans| 27007701| 5682| 42.28| 
6 | Plasmodium yoelii  | rodents| 22222369| 4919| 21.77| 
7 | Haemoproteus tartakovskyi  | birds| 8947331 |2060 | 24.16| 
8 | Toxoplasma gondii  | humans| 128105889|15892 | 52.35| 

# Phylogenetic trees
Print the protein sequences in fasta format for each genome with the aid of the
gffParse.pl program downloaded earlier. Use the -b option to give the protein
file a shorter name as well as the -c option.   

Pb, Plasmodium berghei  
Pc, Plasmodium chabaudi  
Ht, Haemoproteus tartakovskyi  
Pf, Plasmodium falciparum  
Pk, Plasmodium knowlesi   
Pv, Plasmodium vivax  
Py, Plasmodium yoelii  
Tg, Toxoplasma gondii

```bash
cat genemark.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > ht_filtered.gff

perl gffParse.pl -c -p -i Haemoproteus_tartakovskyi.fasta  -g ht_filtered.gff -b Ht #---> Ht.faa and Ht.fna
perl gffParse.pl -c -p -i PlasmodiumGenom/Plasmodium_berghei.genome -g GenePrediction_Plasmodium/P_berghei.gtf -b Pb #---> Pb.faa and Pb.fna
perl gffParse.pl -c -p -i PlasmodiumGenom/Toxoplasma_gondii.genome -g GenePrediction_Plasmodium/Tg.gff -b Tg #---> Tg.faa and Tg.fna
perl gffParse.pl -c -p -i PlasmodiumGenom/Plasmodium_knowlesi.genome -g GenePrediction_Plasmodium/Plasmodium_knowlesi_corrected.gtf -b Pk #---> Pk.faa and Pk.fna
perl gffParse.pl -c -p -i PlasmodiumGenom/Plasmodium_vivax.genome -g GenePrediction_Plasmodium/P_vivax.gtf -b Pv #---> Pv.faa and Pv.fna
perl gffParse.pl -c -p -i PlasmodiumGenom/Plasmodium_yoelii.genome -g GenePrediction_Plasmodium/P_yoelii.gtf -b Py #---> Py.faa and Py.fna
perl gffParse.pl -c -p -i PlasmodiumGenom/Plasmodium_faciparum.genome -g GenePrediction_Plasmodium/P_falciparum.gtf -b Pf #---> Pf.faa and Pf.fna
perl gffParse.pl -c -p -i PlasmodiumGenom/Plasmodium_cynomolgi.genome -g GenePrediction_Plasmodium/P_cynomolgi.gtf -b Pc #---> Pc.faa and Pc.fna
```

## Identify orthologs using Proteinortho
In Proteinortho you will get extra one here that is not in the Busco set. It tries to find matches between genomes compare to other genomes.  
So by evalues each gene get connections to other genes. In the end you get allignment between these.  
proteinortho can be installed over conda. The version we are using is 6.1.7

```bash
conda create -n proteinortho -c bioconda -c conda-forge proteinortho=6.0.33 -y

for f in *.faa; do cat ${f} | sed s'/\t.*//g' > ${f/.faa/_new.faa}; done
#proteinortho
proteinortho6.pl gffparse_All/{Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}_new.faa 
#----> myproject.proteinortho.tsv
```
## Make the Tree
Add the heather and all that starts with 8 to a new file.  
I use just columns for 8 - 8 and 8 - 9
```bash
grep "^#" myproject.proteinortho.tsv > myproject.proteinortho8.tsv
grep "^8" myproject.proteinortho.tsv >> myproject.proteinortho8.tsv
```
Get the python parser from server, to find sequence   
To run this program these files should be availabe:  
1. faa_files=[Ht_new.faa,Pb_new.faa,Pc_new.faa,Pf_new.faa,
           Pk_new.faa,Pv_new.faa,Py_new.faa,Tg_new.faa]  
2. myproject.proteinortho.tsv

```bash
cp /tmp/Programs/ortho_parser.py .
#Edit the script for my data
nano ortho_parser.py
python ortho_parser.py
```
Download clustal from conda

```bash
conda create -n Clustal
conda activate clustal
conda install -c bioconda clustalo
```
Run clustalo on the proteinortho protein fasta files 
```bash
for f in *.fasta; do clustalo -i ${f} -o ${f/.fasta/_aligned.fasta} -v;done

for f in *_aligned.fasta; do cat ${f} | sed 's/.*_/>/g' > ${f/_aligned/renamed_aligned} ; done
#Run for raxml
for f in *renamed_aligned.fasta; do raxmlHPC -s ${f} -n ${f/renamed_aligned.fasta/} tre -o Tg -m PROTGAMMABLOSUM62 -p 12345; done
``` 
phylip can be installed by conda
```bash
conda install -c bioconda phylip
#merge all files to one which is input.tre
cat RAxML_bestTree* > input.tre
#Run Phylip to make tree
#Use consense from the phylip package to merge all individual trees
phylip consense
#it askes for input witch will be input.tre
``` 
Download the tree file from server and uppload it in this website: http://itol.embl.de/
The trees can be visualized at: https://itol.embl.de/tree/1302358190436061677844647   


![image](https://user-images.githubusercontent.com/112621611/222738867-6dd9ea93-b7eb-4057-96d2-0fb45912e153.png)

