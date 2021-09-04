
# coding: utf-8

# # Metaviromics
# 
# In this exercise we will use a very basic workflow for the search of viral metagenomes from RNA sequencing data.
# We are going to use a small, dummy, RNA sequencing dataset that is derived from a larger, real one. In this example we are going to be virologists who have returned from a field expedition in a warm, tropical country where the african tiger mosquito, Aedes aegypti, is present. Aedes mosquitoes are vectors for pathogens like Dengue, Zika and Chikungunya and our aim is to discover which (if any) RNA viruses are infecting the mosquitoes in this region. We have collected mosquito samples, extracted RNA from them and sent them to a company for RNA sequencing with an Illumina instrument. Our starting data are:
#  - **Fastq files**. We have sequenced in paired-end so we have files ending with R1 and R2, which are pairs (sequences from the same fragment, facing each other). In this example we also have more lanes (the L1, L2... suffixes) indicating that our data is split between more lanes.
#  - **Reference genome**. For this exercise we'll use a small portion of the real *Ae. aegypti* genome.
# 
# 
# ## Before starting
# 
# Let's enter the folder that contains the data for this exercise (metaviromics_data) and check the content with ls.
# We start by creating the folders in which we'll store the data we produce. There are two folders already: Fastq and ReferenceGenome whose content is self-explicatory! Let's create folders for the alignment, the unmapped reads, Clark output table and the viral reads we'll select. SPAdes will create is own output folder. 
# We also have to extract the pre-made clark database I conveniently created for you (to save you from downloading the entire NCBI taxonomy). 

# In[ ]:


cd /home/student/DATA/metaviromics_data/
ls
mkdir alignment unmapped clark viral
cd /home/student/DATA/metaviromics_data/CLARKSCV1.2.6.1/
tar -zxvf database.tar.gz
tar -xvf DBdata.tar
ls


# Step 1: There are 8 fastq files, 4 for each pair, corresponding to 4 different lanes used for the sequencing. When working with large files it may be useful to process the lanes in parallel, but our fastq files are easy to handle (small) so we'll merge all the lanes together and produce only two fastq files, one for R1 and one for R2 reads. Be careful to keep the same merging order (from L1 to L4) to avoid errors. Keep the files compressed, all our programs can use them.

# In[ ]:


cd /home/student/DATA/metaviromics_data/Fastq/
cat *R1_L*.fastq.gz > RNAseq_R1.fastq.gz
cat *R2_L*.fastq.gz > RNAseq_R2.fastq.gz


# Step 2: Check if the fastq files are good. Let's use fastqc. The program creates two html files, that can be opened in different ways.
# 1) you can just download the html file on your computer and open it with any browser (if you are using a drag-and drop system) and check if the parameters are ok.
# 2) If you are using mobaxterm, navigate to the files using the sftp browser on the left. Click on the html files with the RMB and click on "open with". Choose your favourite browser, or edge if you have no other choice or want to make microsoft happy.
# 3) On Ubuntu ...
# 4) If you can't download it on your computer, use xdg-open to read the html files.

# In[ ]:


fastqc RNAseq_R1.fastq.gz RNAseq_R2.fastq.gz
xdg-open RNAseq_R1_fastqc.html
xdg-open RNAseq_R2_fastqc.html


# Step 3: 
# Move into the Reference genome folder, look at the genome file and build the index for Hisat2 (it will be used to map RNA-sequencing reads). The two parameters we must provide to the hisat2-build commands are the fasta reference and the name of the new index we will create.

# In[ ]:


cd ../ReferenceGenome
ls
hisat2-build AaegL5_Chr1.fa AaegL5_index
ls


# Step 4: 
# You can see that our command produced 8 files with the .ht2 suffix. These are the index files.
# Map the RNA-seq raw reads to the host reference genome. In our case a fragment of the African Tiger Mosquito (Aedes albopictus) Chromosome 1. 
# For hisat2 to work, we must provide some parameters:
# -x : the reference genome
# -1 and -2 : the R1 and R2 fastq files (can be .gz)
# -S : the sam output file
# -p : the number of threads use by the program

# In[ ]:


cd ..
hisat2 -x ReferenceGenome/AaegL5_index -1 Fastq/RNAseq_R1.fastq.gz -2 Fastq/RNAseq_R2.fastq.gz -p 2 -S alignment/RNAseq_AlignedRefGenome.sam


# Step 5: The output of step 4 is a SAM file (which is uncompressed, so it takes a lot of space) that contains all the reads from the starting fastq files. Some of these reads will be mapped as concordant pair (R1 and R2 facing each other in a meaningful way), others will be mapped as singletons (alone, without a mate) and some will not be mapped at all. 
# Now, extract the unmapped reads. We are interested in reads that are NOT derived from the host genome because they *may* be derived from viruses. Keep in mind that these unmapped reads may also come from bacteria, parasites, contamination or, simply, from the mosquito genome.

# In[ ]:


samtools fastq -@ 2 -f 4 alignment/RNAseq_AlignedRefGenome.sam -1 unmapped/Unmapped_R1.fastq -2 unmapped/Unmapped_R2.fastq -s unmapped/Singletons.fastq


# Step 6: The command should have produced three files, stored in the "unmapped" folder. There is on fastq for each group of reads: R1, R2 and singletons. 
# Time to start searching for viruses! Assemble contigs using SPAdes in viral RNA mode.

# In[ ]:


rnaviralspades.py -t 2 -m 2 -1 unmapped/Unmapped_R1.fastq -2 unmapped/Unmapped_R2.fastq -s unmapped/Singletons.fastq -o SPAdes_output


# Step 7: the output of rnaSPAdes will be a folder with several files and subfolders. For our exercise we are only interested in the contigs.fasta file, which contains the contigs assembled by the programs. Each contig is named as NODE_contig# and the name also reports the length and average RNA coverage of the node.
# Now, classify all the contigs produced by SPAdes using a custom database of viruses (including Flaviviruses, Alphaviruses and Vesiculoviruses) and see if we find any contig from these viral genera! Of course clark can be used with a complete database of viruses, bacteria, or anything the user wants.

# In[ ]:


#DO NOT RUN THIS: $ bash CLARKSCV1.2.6.1/set_targets.sh CLARKSCV1.2.6.1/database/ custom
chmod +x CLARKSCV1.2.6.1/*.sh
chmod +x CLARKSCV1.2.6.1/exe/*
bash CLARKSCV1.2.6.1/classify_metagenome.sh -m 1 -n 2 --light -O SPAdes_output/contigs.fasta -R clark_output/all_contigs_clark


# Step 8: The output of clark is a table with the name of the contig (NODE_#), the length and the taxonomy ID (or NA). 
# Create a new table containing only the names and taxIDS of the contigs classified as viral by Clark. This may be useful if you want to create a table of these contigs on Excel, R or any other software. While you're there, also create a text file of the contigs names, which is required by the ExtractFasta python script to extract the fasta sequences of the viral assemblies.

# In[ ]:


awk -F, '$3 !~ /NA/' clark_output/all_contigs_clark.csv > viral/viral_contigs_clark.csv
cut -f1 -d, viral/viral_contigs_clark.csv | tail -n+2 > viral/virIDs.txt
python ExtractFasta.py -f SPAdes_output/contigs.fasta -k viral/virIDs.txt > viral/Viral_contigs.fasta


# ## Now that we have clark output and the viral contigs fasta...
# 
# Now there is a set of (putative) viral sequences. What can you do with them? I suggest to download on your computer the viral_contigs_clark.csv and import it on Excel or a similar software. You have the fasta of these sequences so you can check what they are using BLAST, NCBI viral taxonomy etc. Also, in "real life" you would have a lot more of results that you can visualize using metagenomics analysis tools (e.g. MEGAN), online tools, Krona... it all depends on your skills and scientific questions.
