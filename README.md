# Vc_Paraguay2024 <!-- omit in toc -->
Here you'll find all the commands necessary to solve the *Vibrio cholerae* command line exercise. Let's start!

## Objectives <!-- omit in toc -->
1) Determine the quality of the sequenced reads
2) Assemble and annotate your genomes
3) Look for AMR and virulence determinants
4) Determine serogroup and serotype
5) Determine Vc lineage and sublineage

```
You will be working in a user of the server named "nw07", so all the directories mentioned in this manual will refer to directories present in this user.
```
## Exercise <!-- omit in toc -->
### Quality control
First, enter the working directory where we will be working throughout the whole course

Type:
```
cd Vibrio
```

Let's see what's inside the folder, type:

```
ls
```

You can see you have different fastq files. We are interested in checking the quality of these reads, type:

```
mkdir fastqc_raw

fastqc -t 8 -o fastqc_raw *.fastq.gz
```

To see the results of these commands `cd` to the "fastq_raw" folder you've created previously. You'll probably need to download the .html files to analyze them in your computer.

> **Remember**: to do this, drop the .html files to the left part of your WinSCP window.

**Question 1: Analysing the FastQC reports, how many reads are there? What's the mean coverage of depth for each sample?**

**Question 2: With respect to quality scores, which of the two files (for each sample) has better-quality data?**

**Question 3: Are these datasets contaminated with any Illumina sequencing adapter oligonucleotides?**

> **Remember**: once you've enter a new directory, you can go to the previous one by typing: `cd ..`

Type:

```
cd ..

```

Now, we want to trim our fastq files by quality as well as checking/removing Illumina adapter sequences, so type:

```
mkdir trimmed

for f in *_1.fastq.gz; do trim_galore -j 8 -q 30 --paired --fastqc --illumina -o trimmed $f ${f%_1.fastq.gz}_2.fastq.gz; done

```

Once trim_galore has finished, check the outputs. You should see that the new files are named something like this:

ARIMX_1_val_1.fq.gz

ARIMX_2_val_2.fq.gz

Let's change the names of these files, type:

```
cd trimmed

for f in *_1_val_1.fq.gz; do mv $f ${f%_1_val_1.fq.gz}_R1.fastq.gz; done

for f in *_2_val_2.fq.gz; do mv $f ${f%_2_val_2.fq.gz}_R2.fastq.gz; done
```

Now, we'll run a software named Kraken2 which matches each read to the lowest common ancestor (LCA) of all genomes, asigning specific taxa. Being in the "trimmed" directory, type: 

```
mkdir kraken2

for f in *_R1.fastq.gz; do kraken2 -db /mnt/Netapp/KRKDB/KRKDB_st8 --threads 8 --gzip-compressed --paired --report kraken2/${f%_R1.fastq.gz}_kraken2.txt --use-names $f ${f%_R1.fastq.gz}_R2.fastq.gz; done
```

To visualize all the Kraken reports at once you can use Pavian: [https://fbreitwieser.shinyapps.io/pavian/](https://fbreitwieser.shinyapps.io/pavian/).
Upload them by selecting all the reports you want to analyze and afterwards click on "Results overview".

### Assembly

We now want to obtain de novo assemblies for our genomes. This process is computationally demanding, so it will take a while to finish for all the samples. Being in the "trimmed" directory, type: 

```
for f in *_R1.fastq.gz; do unicycler -1 $f -2 ${f%_R1.fastq.gz}_R2.fastq.gz -o assemblies/${f%_R1.fastq.gz}.uni.out --verbosity 2 -t 12
```

Once it is over, type:

```
cd assemblies

for f in *.out/assembly.fasta; do cp $f ${f%_uni.out/assembly.fasta}.fasta; done
```

Let's check the quality of our newly obtained assemblies, being in the assembly folder type:

```
quast.py -t 4 *.fasta
```

### Annotation

We will annotate the de novo assemblies created previously. Being in the "assemblies" folder, type:

```
for f in *.fasta; do prokka --prefix ${f%.fasta} --outdir prokka/${f%.fasta}.annotation --addgenes --mincontiglen 300 --cpus 8 $f
```

To explore the files produced, enter the "prokka" folder.

### Looking for AMR and virulence determinants

There are lots of softwares for this purpose but we will be using ariba (with resfinder and vfdb databases) and AMRFinderPlus. To start, being in the "trimmed" folder, type:

```
mkdir ariba

for f in *_R1.fastq.gz; do ariba run --threads 6 /home/inei/secuencias/prepareref/resfinder_db.out/ $f ${f%_R1.fastq.gz}_R2.fastq.gz ariba/${f%_R1.fastq.gz}.res.out.dir; done

for f in *_R1.fastq.gz; do ariba run --threads 6 /home/inei/secuencias/prepareref/vfdb_core_db.out/ $f ${f%_R1.fastq.gz}_R2.fastq.gz ariba/${f%_R1.fastq.gz}.vfdb.out.dir; done

ariba summary ariba/out_res ariba/*res.out.dir/report.tsv

ariba summary ariba/out_vfdb ariba/*vfdb.out.dir/report.tsv
```

AMRFinderPlus requires as input a fasta file, so we will try this software using the de novo assemblies. Being in the "assemblies" folder, type:

```
for f in *.fasta; do amrfinder -O Vibrio_cholerae --plus -n $f -o ${f%.fasta}_amrfinderplus.tsv
```

### Determining serogroup and serotype

To determine the serotype we will check the reports obtained with ariba and vfdb database:

> We assume that an Inaba phenotype would be conferred on isolates in which ariba was unable to detect or assemble *wbeT* in its
totality, and if a mutation in *wbeT* was detected that was predicted to frameshift or truncate translated wbeT (N62fs, N165fs, F244fs, Q274trunc), was associated with Inaba phenotypes (I206K), or was otherwise known to confer an Inaba phenotype (S158P). [(Dorman et al 2020 PMID: 33004800)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7530988/)

To determine the serogroup, being in the "assemblies" folder, type:

```
mkdir serogroup/

cd serogroup/

cp /home/sh12/Analisis/VC/AnalisisVC/serogroupDB/DB_Vc_Oserogroup.fasta .

makeblastdb -in DB_Vc_Oserogroup.fasta -dbtype nucl -out serogroup_db

cd ..

blastn -query T_VC4.fasta -db serogroup/serogroup_db -out T_VC4_blast_O.txt -outfmt 6 -max_target_seqs 3
```

### Bonus! <!-- omit in toc -->

We would like to determine the *ctxB* variant, for that we will use ariba with a custome database. Being in the trimmed folder, type:

```
for f in *_R1.fastq.gz; do ariba run --threads 4 /home/sh12/Analisis/VC/AnalisisVC/ctxB_out_prepareref/ $f ${f%_R1.fastq.gz}_R2.fastq.gz ariba/${f%_R1.fastq.gz}.ctx.out.dir; done 
```

### Determining the lineage of a *Vibrio cholerae* genome

We will use a small collection of annotated genomes in order to build a phylogenetic tree from core genome SNPs. You can find the metadata here: [https://docs.google.com/spreadsheets/d/1nyGB5RA3eu7xlJ_BBXrF5trsQyJmqKetKuBfPbN4WvE/edit?usp=sharing](https://docs.google.com/spreadsheets/d/1nyGB5RA3eu7xlJ_BBXrF5trsQyJmqKetKuBfPbN4WvE/edit?usp=sharing)

To obtain the core genome alignment we will use [Roary](https://github.com/sanger-pathogens/Roary). Being in the "Vc_training2023" folder, type: 

```
roary -e --mafft -f roary070723 -p 12 roary_dataset/*.gff
```

Now, we want to keep the variant sites from the alignment, for that we will use [snp-sites](https://github.com/sanger-pathogens/snp-sites). Being in the "Vc_training2023" folder, type:

```
snp-sites -o roary070723/core_gene_alignment_snps.aln roary070723/core_gene_alignment.aln
```

For building up the tree we will use [IQ-TREE](https://github.com/iqtree/iqtree2) which is an efficient and versatile phylogenomic software by maximum likelihood. Being in "roary070723" folder, type:

```
iqtree -s core_gene_alignment_snps.aln -m TEST -pre vc_roary -bb 1000 -nt 12
```

### Determining the sublineage of a *Vibrio cholerae* 7PET genome

For determining the sublineage we will build up a tree from the SNPs obtained from the mapping of the reads against a 7PET reference genome (N16961). We will use [snippy](https://github.com/tseemann/snippy) to do the reference mapping. The metadata of the context dataset used in this task is in: [https://docs.google.com/spreadsheets/d/17j9f_cYUxsbDTbTWhxfRB1oPt_aXP1GjwrnmjtDdecA/edit#gid=0](https://docs.google.com/spreadsheets/d/17j9f_cYUxsbDTbTWhxfRB1oPt_aXP1GjwrnmjtDdecA/edit#gid=0)

Being in the "Vc_training2023", type:

```
snippy --cpus 4 --outdir snippy_dataset/T_VC4_snippy --ref N16961.fna --R1 trimmed/T_VC4_R1.fastq.gz --R2 trimmed/T_VC4_R2.fastq.gz
```

> Repeat this command for every 7PET genome.

Now lets use 'snippy-core' to summarise all these genomes and create a multiple sequence alignment. Type: 

```
cd snippy_dataset
snippy-core --ref ../N16961.fna *_snippy
```

As we want everything masked in the same way, let’s change that so anything uncertain is marked as ’N’ using the 'snippy-clean_full_aln' script that comes with 'snippy'. Type:

```
snippy-clean_full_aln core.full.aln > clean.full.aln
```

We can use [gubbins](https://github.com/nickjcroucher/gubbins) to infer recombining sites by looking for increased SNP density that occurs in specific ancestral nodes. Type:

```
conda activate
run_gubbins.py --threads 10 -p gubbins clean.full.aln
conda deactivate
```
Again, to remove all the invariant sites and create a SNP-only multiple sequence alignment we will use [snp-sites](https://github.com/sanger-pathogens/snp-sites). Type: 

```
snp-sites -c gubbins.filtered_polymorphic_sites.fasta -o clean.core.aln
```

Now, let's build our tree. Type: 

```
iqtree -s clean.core.aln -m TEST -pre vc_snippy -bb 1000 -nt 12
```

