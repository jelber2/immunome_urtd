# Immunome_URTD
=====
## Get Fastq files from BaseSpace
### Rename fastq files
    # make run and combined folders
    cd /media/immunome_2014/work/jelber2/immunome_urtd/
    mkdir run1
    mkdir run2
    mkdir combined
    # rename the files to remove spaces
    cd ~/Desktop/
    # run1
    mv Bioo\ NEXTflex-23135384.zip Bioo-NEXTflex-23135384.zip
    ## copy the file to external drive
    cp ~/Desktop/Bioo-NEXTflex-23135384.zip /media/immunome_2014/work/jelber2/immunome_urtd/run1/.
    ## go to the directory
    cd /media/immunome_2014/work/jelber2/immunome_urtd/run1/
    ## unzip the archive
    unzip Bioo-NEXTflex-23135384.zip
    ## rename the unzipped folder
    mv Bioo\ NEXTflex-23135384 Bioo-NEXTflex-23135384
    ## grab all of the file names and save them as a document
    find /media/immunome_2014/work/jelber2/immunome_urtd/run1/Bioo-NEXTflex-23135384/ \
    -maxdepth 5 -type f -print > rename-fastq-run1
    ## make a directory called fastq
    mkdir fastq
    ## use regular expressions to create rename-fastq-run1.sh
    perl -pe s"/(.+\/BaseCalls\/)(\w+)(_\w+_\w+_)(R\d)(_\d+.fastq.gz)/mv \1\2\3\4\5 \/media\/immunome_2014\/work\/jelber2\/immunome_urtd\/run1\/fastq\/\2-\4.fastq.gz/" rename-fastq-run1 > rename-fastq-run1.sh
    # add #!/bin/bash to first line
    sed -i '1 i\#!/bin/bash' rename-fastq-run1.sh
    bash rename-fastq-run1.sh
    # run2
    ## copy fastq files to /media/immunome_2014/work/jelber2/immunome_urtd/run2/
    cd /media/immunome_2014/work/jelber2/immunome_urtd/run2/
    zip immunome_urtd_run2.zip *.gz
    ## make a directory called fastq
    mkdir fastq
    ## move data to fastq
    mv *.gz ./fastq/.
    ## use regular expressions to create rename-fastq-run2.sh
    cd fastq
    ls *.gz > rename-fastq-run2
    perl -pe s"/(\w+)(_\w+_\w+_)(R\d)(_\d+.fastq.gz)/mv \1\2\3\4 \1-\3.fastq.gz/" rename-fastq-run2 > rename-fastq-run2.sh
    # add #!/bin/bash to first line
    sed -i '1 i\#!/bin/bash' rename-fastq-run2.sh
    bash rename-fastq-run2.sh
## Put Data on SuperMikeII
    # run1
    rsync --stats --progress --archive \
    /media/immunome_2014/work/jelber2/immunome_urtd/run1/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_urtd/run1/ -n
    # run2
    rsync --stats --progress --archive \
    /media/immunome_2014/work/jelber2/immunome_urtd/run2/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_urtd/run2/ -n
## Install programs and get reference genome
### trimmomatic-0.32
    cd /home/jelber2/bin/
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
    unzip Trimmomatic-0.32.zip
    mv Trimmomatic-0.32.zip Trimmomatic-0.32
    #PATH=~/home/jelber2/bin/Trimmomatic-0.32/trimmomatic-0.32.jar
### bbmerge-5.4 (part of bbmap-34.33)
    cd /home/jelber2/bin/
    mkdir bbmap-34.33
    cd bbmap-34.33/
    wget http://downloads.sourceforge.net/project/bbmap/BBMap_34.33.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2F%3Fsource%3Ddlp&ts=1421955805&use_mirror=iweb
    mv BBMap_34.33.tar.gz\?r\=http\:%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2F\?source\=dlp BBMap_34.33.tar.gz
    tar xzf BBMap_34.33.tar.gz
    cd bbmap/
    mv * ..
    cd ..
    rm -r bbmap
    #PATH=~/bin/bbmap-34.33/bbmerge.sh
### bwa-0.7.12
    cd /home/jelber2/bin/
    wget https://github.com/lh3/bwa/archive/0.7.12.tar.gz
    mv 0.7.12 bwa-0.7.12.tar.gz
    tar xzf bwa-0.7.12.tar.gz
    mv bwa-0.7.12.tar.gz bwa-0.7.12
    cd bwa-0.7.12/
    make
    #PATH=~/bin/bwa-0.7.12/bwa
### stampy-1.0.23
    cd /home/jelber2/bin/
    wget http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz
    tar xzf Stampy-latest.tgz
    cd stampy-1.0.23
    make
    #PATH=~/bin/stampy-1.0.23/stampy.py
### java jre1.7.0
    #had to download using firefox on my Centos machine
    #saved in /home/jelber2/bin/
    rsync --stats --archive --progress /home/jelber2/bin/jre-7-linux-x64.tar.gz jelber2@mike.hpc.lsu.edu:/home/jelber2/bin/ -n
    #switched to SuperMikeII
    cd /home/jelber2/bin/
    tar xzf jre-7-linux-x64.tar.gz
    mv jre-7-linux-x64.tar.gz jre1.7.0
    #add  PATH += /home/jelber2/bin/jre1.7.0/bin  to .soft file
    nano ~/.soft
    #then resoft
    #resoft
    #PATH=~/bin/jre1.7.0/bin
### picard-1.128
    #on my Centos machine
    cd /home/jelber2/bin/
    wget https://github.com/broadinstitute/picard/releases/download/1.128/picard-tools-1.128.zip
    rsync --stats --archive --progress /home/jelber2/bin/picard-tools-1.128.zip jelber2@mike.hpc.lsu.edu:/home/jelber2/bin/ -n
    #switched to SuperMikeII
    cd /home/jelber2/bin/
    unzip picard-tools-1.128.zip
    mv picard-tools-1.128.zip picard-tools-1.128
    #PATH=~/bin/picard-tools-1.128/picard.jar
### GATK-3.3.0
    #had to download using firefox on my Centos machine
    #saved in /home/jelber2/bin/GATK-3.3.0
    rsync --stats --archive --progress /home/jelber2/bin/GATK-3.3.0/ jelber2@mike.hpc.lsu.edu:/home/jelber2/bin/GATK-3.3.0/ -n
    #switched to SuperMikeII
    cd /home/jelber2/bin/
    cd GATK-3.3.0
    tar xjf GenomeAnalysisTK-3.3-0.tar.bz2
    #PATH=~/bin/GATK-3.3.0/GenomeAnalysisTK.jar
### samtools-1.1
    cd /home/jelber2/bin/
    wget http://downloads.sourceforge.net/project/samtools/samtools/1.1/samtools-1.1.tar.bz2?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2Ffiles%2Fsamtools%2F1.1%2F&ts=1421967581&use_mirror=softlayer-dal
    tar xjf samtools-1.1.tar.bz2 
    mv samtools-1.1.tar.bz2 samtools-1.1
    cd samtools-1.1
    make
    nano ~/.soft #add the following line to .soft file using nano
    PATH += /home/jelber2/bin/samtools-1.1/
### parallel-20150122
    cd /home/jelber2/bin/
    wget ftp://ftp.gnu.org/gnu/parallel/parallel-20150122.tar.bz2
    tar xjf parallel-20150122.tar.bz2
    mv parallel-20150122.tar.bz2 parallel-20150122
    #PATH=~/bin/parallel-20150122/src/parallel
#### Get bedtools2.22.1
    cd ~/bin/
    wget https://github.com/arq5x/bedtools2/releases/download/v2.22.1/bedtools-2.22.1.tar.gz
    tar xzf bedtools-2.22.1.tar.gz
    mv bedtools2 bedtools-2.22.1
    mv bedtools-2.22.1.tar.gz bedtools-2.22.1
    cd bedtools-2.22.1
    make
### Got painted turtle reference genome (on SuperMikeII)
    cd /work/jelber2/reference/
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
    gunzip GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
## STEPS FOR QUALITY CONTROL, MAPPING, & SNP CALLING
### 1.Make indexes (on SuperMikeII)
    qsub /home/jelber2/scripts/immunome_urtd/01-make_indexes.sh
### 2.Quality control and adapter trimming (on SuperMikeII)
    # run1
    cd /work/jelber2/immunome_urtd/run1/fastq/
    ~/scripts/immunome_urtd/02-trimmomatic.py *.fastq.gz
    # run2
    cd /work/jelber2/immunome_urtd/run2/fastq/
    ~/scripts/immunome_urtd/02-trimmomatic-run2.py *.fastq.gz
### 3.BWA alignment (on SuperMikeII)
    # run1
    cd /work/jelber2/immunome_urtd/run1/trimmed-data/
    ~/scripts/immunome_urtd/03-bwa.py *.trim.fastq.gz
    # run2
    cd /work/jelber2/immunome_urtd/run2/trimmed-data/
    ~/scripts/immunome_urtd/03-bwa-run2.py *.trim.fastq.gz
### 4.STAMPY alignment (on SuperMikeII)
    # run1
    cd /work/jelber2/immunome_urtd/run1/bwa-alignment/
    ~/scripts/immunome_urtd/04-stampy.py *.bwa.sam
    # run2
    cd /work/jelber2/immunome_urtd/run2/bwa-alignment/
    ~/scripts/immunome_urtd/04-stampy-run2.py *.bwa.sam
### 5a.Clean,Sort,Add Read Groups (on SuperMikeII)
    # run1
    cd /work/jelber2/immunome_urtd/run1/stampy-alignment/
    ~/scripts/immunome_urtd/05a-clean_sort_addRG.py *.stampy.bam
    # run2
    cd /work/jelber2/immunome_urtd/run2/stampy-alignment/
    ~/scripts/immunome_urtd/05a-clean_sort_addRG-run2.py *.stampy.bam
### 5b.Clean,Sort,Add Read Groups, DeDup,Realign Around Indels(on SuperMikeII)
    cd /work/jelber2/immunome_urtd/run1/clean-sort-addRG/
    ~/scripts/immunome_urtd/05b-clean_sort_addRG_markdup_realign.py *-CL-RG.bam
### 6.Merge BAM files, Call SNPs initially (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/realign-around-indels/
    ~/scripts/immunome_urtd/06-mergeBAM_callSNPs_initial.py *-realigned.bam
### 7.Quality score recalibration 1 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/realign-around-indels/
    find . -name '*-realigned.bam' -not -name 'ALL-samples-*' \
     -exec ~/scripts/immunome_urtd/07-qual_score_recal01.py {} \;
### 8.Merge BAM files, Call SNPs recalibrated 1 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal01/
    ~/scripts/immunome_urtd/08-mergeBAM_callSNPs_recal01.py *-recal01.bam
### 9.Quality score recalibration 2 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal01/
    find . -name '*-recal01.bam' -not -name 'ALL-samples-*' \
    -exec ~/scripts/immunome_urtd/09-qual_score_recal02.py {} \;
### 10.Merge BAM files, Call SNPs recalibrated 2 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal02/
    ~/scripts/immunome_urtd/10-mergeBAM_callSNPs_recal02.py *-recal02.bam
### 11.Quality score recalibration 3 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal02/
    find . -name '*-recal02.bam' -not -name 'ALL-samples-*' \
    -exec ~/scripts/immunome_urtd/11-qual_score_recal03.py {} \;
### 12.Merge BAM files, Call SNPs recalibrated 3 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal03/
    ~/scripts/immunome_urtd/12-mergeBAM_callSNPs_recal03.py *-recal03.bam
### 13.Sequencing metrics
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal03/
    ~/scripts/immunome_2014/13-seq_metrics.py ALL-samples-recal03.bam
### 14.plot coverage for each sample
#### run bedtools coverage on all bam files, then keep only lines with 'all' on them
    #FROM http://gettinggeneticsdone.blogspot.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
    #ALSO FROM https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md
#####make samplelist
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/call-SNPs-recal03/
    ls *.bam | grep -v "ALL"| perl -pe "s/-recal03.bam//g" > samplelist
#####calculate coverage
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir plot_coverage
    cd plot_coverage/
    #command below takes 5-10 minutes
    while read i;do
    ~/bin/bedtools-2.22.1/bin/bedtools coverage -abam ../call-SNPs-recal03/$i-recal03.bam \
    -b /media/immunome_2014/work/jelber2/reference/immunome_baits_C_picta-3.0.3.bed \
    -hist | grep ^all > $i.baitcoverage.all.txt
    done < ../call-SNPs-recal03/samplelist
    #now use modified R scripts from links above to plot coverage
### 15.Need to use featureCounts to summarize number of genes, reads per gene, etc
#### Get Subread
    #featureCounts is part of the Subread package http://bioinf.wehi.edu.au/featureCounts/
    cd ~/bin/
    wget http://downloads.sourceforge.net/project/subread/subread-1.4.6/subread-1.4.6-Linux-x86_64.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsubread%2Ffiles%2Fsubread-1.4.6%2F&ts=1423446986&use_mirror=iweb
    mv subread-1.4.6-Linux-x86_64.tar.gz?r=http:%2F%2Fsourceforge.net%2Fprojects%2Fsubread%2Ffiles%2Fsubread-1.4.6%2F subread-1.4.6-Linux-x86_64.tar.gz
    tar xzf subread-1.4.6-Linux-x86_64.tar.gz
    mv subread-1.4.6-Linux-x86_64.tar.gz subread-1.4.6-Linux-x86_64
#### Get genometools-1.5.4 to annotate introns in gff file
    cd ~/bin/
    wget http://genometools.org/pub/genometools-1.5.4.tar.gz
    tar xzf genometools-1.5.4.tar.gz
    mv genometools-1.5.4.tar.gz genometools-1.5.4
    cd genometools-1.5.4
    #on MacOSX
    make
    #on CentOS
    #install ruby first
    #become superuser
    su
    yum install ruby.x86_64
    #stop being a super user
    exit
    #make the executable using 64bit mode but without cairo
    make 64bit=yes cairo=no
    #test the install - will take a long time (>30min?)
    make 64bit=yes cairo=no test
#####Use genometools to get introns
    ~/bin/genometools-1.5.4/bin/gt gff3 -addintrons yes \
    -o /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.introns \
    /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff
#### Need to prefilter the GFF file for immune genes
    cd /media/immunome_2014/work/jelber2/reference/
    #intersect the gff file
    ~/bin/bedtools-2.22.1/bin/bedtools intersect \
    -a GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.introns \
    -b immunome_baits_C_picta-3.0.3.bed \
    > GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gff
#### Convert GFF file of gene annotations to GTF
    cd /media/immunome_2014/work/jelber2/reference/
    perl -pe "s/\S+=GeneID:(\d+).+/gene_id \"\1\";/g" \
    GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gff \
    > GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf
#### run subread
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir featureCounts
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/featureCounts/
#####all samples at once
    #on CentOS machine
    #at the gene level
    ~/bin/subread-1.4.6-Linux-x86_64/bin/featureCounts \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf \
    -o ALL.gene -F GTF -T 2 --ignoreDup \
    ../call-SNPs-recal03/ALL-samples-recal03.bam
    #at exon level
    ~/bin/subread-1.4.6-Linux-x86_64/bin/featureCounts \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf \
    -o ALL.exon -F GTF -T 2 -f --ignoreDup \
    ../call-SNPs-recal03/ALL-samples-recal03.bam
#####each sample separately
    #at the gene level
    while read i;do
    ~/bin/subread-1.4.6-Linux-x86_64/bin/featureCounts \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf \
    -o $i.gene -F GTF -T 2 --ignoreDup ../call-SNPs-recal03/$i-recal03.bam
    done < ../call-SNPs-recal03/samplelist
    #at the exon level
    while read i;do
    ~/bin/subread-1.4.6-Linux-x86_64/bin/featureCounts \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf \
    -o $i.exon -F GTF -T 2 -f --ignoreDup ../call-SNPs-recal03/$i-recal03.bam
    done < ../call-SNPs-recal03/samplelist
#### write an R function to do the following
#####count number of possible immune genes
    wc -l CF53.gene
    #total genes = 632 (after subtracting 2 header lines)
#####count number of possible immune gene exons
    wc -l CF53.exon
    #total exons = 37275 (after subtracting 2 header lines)
#####how many different immune genes were captured
    grep -Pv "\t0$" ALL.gene | wc -l
    #615 (after subtracting 2 header lines)
#####how many different immune gene exons were captured
    grep -Pv "\t0$" ALL.exon | wc -l
    #4860 (after subtracting 2 header lines)
#####count number of genes per sample
    while read i;do
    test=$(grep -Pv "\t0$" $i.gene | wc -l)
    echo -e $i'\t'$test > $i.gene.count
    done < ../call-SNPs-recal03/samplelist
    cat *.gene.count > gene.counts.per.sample
#####count number of exons per sample
    while read i;do
    grep -Pv "\t0$" $i.exon | wc -l > $i.exon.count
    done < ../call-SNPs-recal03/samplelist
    cat *.exon.count > exon.counts.per.sample
### 16.Haplotype Caller (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal03/
    #excludes file ALL-samples-recal03.bam
    find . -name '*-recal03.bam' -not -name 'ALL-samples-*' \
    -exec ~/scripts/immunome_urtd/14-haplotypecaller.py {} \;
### 17.Ran GenotypeGVCFs to perform joint genotyping
    #on Cenots machine
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/hc/
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    --max_alternate_alleles 32 \
    --variant CF219-raw-snps-indels.vcf \
    --variant CF53-raw-snps-indels.vcf \
    --variant CF69-raw-snps-indels.vcf \
    --variant CF72-raw-snps-indels.vcf \
    --variant CF80-raw-snps-indels.vcf \
    --variant CF90-raw-snps-indels.vcf \
    --variant FC13-raw-snps-indels.vcf \
    --variant FC15-raw-snps-indels.vcf \
    --variant FC19-raw-snps-indels.vcf \
    --variant FC47-raw-snps-indels.vcf \
    --variant FC58-raw-snps-indels.vcf \
    --variant OLD106-raw-snps-indels.vcf \
    --variant OLD107-raw-snps-indels.vcf \
    --variant OLD65-raw-snps-indels.vcf \
    --variant OLD77-raw-snps-indels.vcf \
    --variant OLD92-raw-snps-indels.vcf \
    -o ALL-samples-raw-snps-indels.vcf
### 18.Added expressions to filter variants
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -L /media/immunome_2014/work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -V ALL-samples-raw-snps-indels.vcf \
    --clusterWindowSize 10 \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
    --filterName "Bad_Validation" \
    --filterExpression "QUAL < 30.0" \
    --filterName "LowQual" \
    --genotypeFilterExpression "DP < 10.0" \
    --genotypeFilterName "Low_Read_Depth_Over_Sample" \
    --genotypeFilterExpression "GQ < 20.0" \
    --genotypeFilterName "Low_GenotypeQuality" \
    -o ALL-samples-Q30-snps-indels.vcf
### 19.Got only Indel variants
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -T SelectVariants \
    -V ALL-samples-Q30-snps-indels.vcf \
    -o ALL-samples-Q30-indels.vcf \
    -selectType INDEL
### 20.Got only SNP variants
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -T SelectVariants \
    -V ALL-samples-Q30-snps-indels.vcf \
    -o ALL-samples-Q30-snps.vcf \
    -selectType SNP
### 21.Ran Variant Recalibrator
#### Note: used SNPS and Indels from previous immunome_2014 data for "truthing" and filtration
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir vqsr
    cd vqsr
#####Recalibrated snps
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -input ../hc/ALL-samples-Q30-snps.vcf \
    -resource:concordantSet,VCF,known=true,training=true,truth=true,prior=10.0 /media/immunome_2014/work/jelber2/immunome_2014/combined/beagle/ALL-samples-Q30-snps-recal-beagle.vcf \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
    -recalFile VQSR-snps.recal \
    -mode SNP \
    -tranchesFile VQSR-snps.tranches \
    -rscriptFile VQSR-snps.plots.R \
    --maxGaussians 4
#####Recalibrated indels
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -input ../hc/ALL-samples-Q30-indels.vcf \
    -resource:concordantSet,VCF,known=true,training=true,truth=true,prior=10.0 /media/immunome_2014/work/jelber2/immunome_2014/combined/beagle/ALL-samples-Q30-indels-recal-beagle.vcf \
    -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
    -recalFile VQSR-indels.recal \
    -mode INDEL \
    -tranchesFile VQSR-indels.tranches \
    -rscriptFile VQSR-indels.plots.R \
    --maxGaussians 4
#### Applied the recalibration on snps
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -input ../hc/ALL-samples-Q30-snps.vcf \
    --ts_filter_level 99.5 \
    -tranchesFile VQSR-snps.tranches \
    -recalFile VQSR-snps.recal \
    -o ALL-samples-Q30-snps-recal.vcf
#### Applied the recalibration on indels
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -input ../hc/ALL-samples-Q30-indels.vcf \
    --ts_filter_level 99.0 \
    -tranchesFile VQSR-indels.tranches \
    -recalFile VQSR-indels.recal \
    -o ALL-samples-Q30-indels-recal.vcf
### 22.Needed to use beagle to improve SNPs (using Linkage Disequilibrium) called by Haplotype Caller
#### Downloaded beagle
    cd ~/bin
    wget http://faculty.washington.edu/browning/beagle/beagle.r1398.jar
#### Ran beagle on snps and indels separately
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir beagle
    cd beagle
    #snps
    java -Xmx4000m -jar ~/bin/beagle.r1398.jar \
    gtgl=/media/immunome_2014/work/jelber2/immunome_urtd/combined/vqsr/ALL-samples-Q30-snps-recal.vcf\
    nthreads=2 \
    out=/media/immunome_2014/work/jelber2/immunome_urtd/combined/beagle/ALL-samples-Q30-snps-recal-beagle
    #indels
    java -Xmx4000m -jar ~/bin/beagle.r1398.jar \
    gtgl=/media/immunome_2014/work/jelber2/immunome_urtd/combined/vqsr/ALL-samples-Q30-indels-recal.vcf\
    nthreads=2 \
    out=/media/immunome_2014/work/jelber2/immunome_urtd/combined/beagle/ALL-samples-Q30-indels-recal-beagle
### 23.Get only polymorphic loci
    # Because many SNP/indel loci will occur because "fixed" differeces
    # between Gopher tortoise and Western Painted turtle
#### Downloaded vcftools
    cd /home/jelber2/bin/
    wget http://downloads.sourceforge.net/project/vcftools/vcftools_0.1.12b.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fvcftools%2Ffiles%2F&ts=1411515317&use_mirror=superb-dca2
    tar -xzf vcftools_0.1.12b.tar.gz 
    mv vcftools_0.1.12b.tar.gz vcftools_0.1.12b
    cd vcftools_0.1.12b/
    nano ~/.soft #add the following two lines to using nano.soft file
    PATH+=/home/jelber2/bin/tabix-0.2.6
    PERL5LIB = /home/jelber2/bin/vcftools_0.1.12b/perl
    resoft #to refresh soft file
    cd /home/jelber2/bin/vcftools_0.1.12b/
    make #compile vcftools
    # Path to vcftools executable
    /home/jelber2/bin/vcftools_0.1.12b/bin/vcftools
#### Remove loci with AF=1
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/beagle/
#####snps
    #removes loci with AF=1
    zcat ../beagle/ALL-samples-Q30-snps-recal-beagle.vcf.gz | grep -v "AF=1" \
    > ALL-samples-Q30-snps-recal-beagle2.vcf
    #still leaves some unwanted non-polymorphic snps?
    #try calculating allele frequencies then
    ~/bin/vcftools_0.1.12b/bin/vcftools \
    --vcf ALL-samples-Q30-snps-recal-beagle2.vcf \
    --freq --out ALL-samples-Q30-snps-recal-beagle
    #get only non-polymorphic loci
    grep ":1" ALL-samples-Q30-snps-recal-beagle.frq | cut -f 1-2 > nonpolymorphicsnps
    grep -v "\s2\s" ALL-samples-Q30-snps-recal-beagle.frq | cut -f 1-2 > multiallelicloci
    cat nonpolymorphicloci multiallelicloci > multiallelic_or_nonpolymorphicloci
    #calculate number of di-,tri-,tetra-allelic loci
    cut -f 3 ALL-samples-Q30-snps-recal-beagle.frq | sort | uniq -c
    #    di = 19355 (includes non-polymorphic loci = 3186)
    #   tri = 651
    # tetra = 1
    #
    #filter out nonpolymorphicsnps
    #might take a few minutes
    while read i;do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/remove\tlocus\t\n/" ALL-samples-Q30-snps-recal-beagle2.vcf
    done < multiallelic_or_nonpolymorphicloci
    grep -v "remove" ALL-samples-Q30-snps-recal-beagle2.vcf \
    > ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf
#####indels
    #removes loci with AF=1
    zcat ../beagle/ALL-samples-Q30-indels-recal-beagle.vcf.gz | grep -v "AF=1" \
    > ALL-samples-Q30-indels-recal-beagle2.vcf
    #still leaves some unwanted non-polymorphic indels?
    #try calculating allele frequencies then
    ~/bin/vcftools_0.1.12b/bin/vcftools \
    --vcf ALL-samples-Q30-indels-recal-beagle2.vcf \
    --freq --out ALL-samples-Q30-indels-recal-beagle
    #get only non-polymorphic loci
    grep ":1" ALL-samples-Q30-indels-recal-beagle.frq | cut -f 1-2 > nonpolymorphicindels
    #filter out nonpolymorphicindels
    #might take a few minutes
    while read i;do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/remove\tlocus\t\n/" ALL-samples-Q30-indels-recal-beagle2.vcf
    done < nonpolymorphicindels
    grep -v "remove" ALL-samples-Q30-indels-recal-beagle2.vcf \
    > ALL-samples-Q30-indels-recal-beagle-polymorphic.vcf
=====
## STEPS FOR VARIANT PREDICTION
### 1.Download Tools First
#### Downloaded snpEff version 4.0e
    #ideally want to know if variants will affect protein structure and possibly immune gene function
    cd /work/jelber2/reference
    wget http://iweb.dl.sourceforge.net/project/snpeff/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
#### Added Chrysemys_picta_bellii-3.0.3 to snpEff.config using nano
    cd /media/immunome_2014/work/jelber2/reference/snpEff
    nano snpEff.config # added the following four lines after the Capsella_rubella_v1.0 entry (remove 4 spaces on left if cut and pasting)
    # Chrysemys_picta_bellii-3.0.3
    Chrysemys_picta_bellii-3.0.3.genome : western painted turtle
    	Chrysemys_picta_bellii-3.0.3.reference : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/
    	Chrysemys_picta_bellii-3.0.3.M.codonTable : Standard
#### Created data directory for Chrysemys_picta_bellii-3.0.3 genome
    cd /media/immunome_2014/work/jelber2/reference/snpEff
    mkdir data
    cd data
    mkdir Chrysemys_picta_bellii-3.0.3
    cd Chrysemys_picta_bellii-3.0.3
    # downloaded FASTA file
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
    # snpEff requires genome.fa file to be called "sequences.fa"
    mv GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz sequences.fa.gz
    # have to unzip sequences.fa.gz
    gunzip sequences.fa.gz
    # downloaded gff3 file (i.e., gene annotation file)
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.gz
    # snpEff requires gene annotation file be called "genes.gff"
    mv GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.gz genes.gff.gz
    # unzipped genes.gff.gz
    gunzip genes.gff.gz
    # download protein sequences
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_protein.faa.gz
    mv GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_protein.faa.gz protein.fa.gz
    gunzip protein.fa.gz
#### Built snpEff database for Chrysemys_picta_bellii-3.0.3
    cd /media/immunome_2014/work/jelber2/reference/snpEff/
    # used snpEff_build.py script to implement command below, which took < 30 minutes
    java -jar -Xmx4g /media/immunome_2014/work/jelber2/reference/snpEff/snpEff.jar build -gff3 -v Chrysemys_picta_bellii-3.0.3 2>&1 | tee Chrysemys_picta_bellii-3.0.3.build
#### Downloaded bcftools
    cd ~/bin/
    git clone --branch=develop git://github.com/samtools/htslib.git
    git clone --branch=develop git://github.com/samtools/bcftools.git
    cd bcftools; make
### 2.Need to look for protein altering variants shared by samples in the same phenotype
#### a.Split vcf file for snpEff
    #snpEff needs ALL-samples*.vcf file split by sample (i.e., into Sample1.vcf, Sample2.vcf)
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir split-vcfs
    cd split-vcfs/
    #compress vcf files
    ~/bin/samtools-1.1/htslib-1.1/bgzip -f ../beagle/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf
    ~/bin/samtools-1.1/htslib-1.1/bgzip -f ../beagle/ALL-samples-Q30-indels-recal-beagle-polymorphic.vcf
    #index vcf.gz with tabix
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ../beagle/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ../beagle/ALL-samples-Q30-indels-recal-beagle-polymorphic.vcf.gz
    #split files
    #code to split each vcf file
    #snps
    while read i;do
    ~/bin/bcftools/bcftools view -s $i ../beagle/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz > $i-snps.vcf
    done < ../call-SNPs-recal03/samplelist
    #indels
    while read i;do
    ~/bin/bcftools/bcftools view -s $i ../beagle/ALL-samples-Q30-indels-recal-beagle-polymorphic.vcf.gz > $i-indels.vcf
    done < ../call-SNPs-recal03/samplelist
#### b.Ran snpEff on each split vcf file
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/split-vcfs/
    # command below to run snpEff on all samples in samplelist
    # not implemented on SuperMikeII b/c process was < 15 min
    #snps
    while read i;do
    java -Xmx4g -jar /media/immunome_2014/work/jelber2/reference/snpEff/snpEff.jar \
    -v -i vcf -o gatk \
    Chrysemys_picta_bellii-3.0.3 \
    $i-snps.vcf > $i-snps-snpeff.vcf
    mv snpEff_genes.txt $i-snps-snpeff-genes.txt
    mv snpEff_summary.html $i-snps-snpeff-summary.html
    done < ../call-SNPs-recal03/samplelist
    #indels
    while read i;do
    java -Xmx4g -jar /media/immunome_2014/work/jelber2/reference/snpEff/snpEff.jar \
    -v -i vcf -o gatk \
    Chrysemys_picta_bellii-3.0.3 \
    $i-indels.vcf > $i-indels-snpeff.vcf
    mv snpEff_genes.txt $i-indels-snpeff-genes.txt
    mv snpEff_summary.html $i-indels-snpeff-summary.html
    done < ../call-SNPs-recal03/samplelist
#### c.Ran VariantAnnotator on each snpeff file
    #snps
    while read i;do
    rm $i-snps.vcf.idx
    rm $i-snps-snpeff.vcf.idx
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -A SnpEff \
    --variant $i-snps.vcf \
    --snpEffFile $i-snps-snpeff.vcf \
    -L $i-snps.vcf \
    -o $i-snps-annotated.vcf
    done < ../call-SNPs-recal03/samplelist
    #indels
    while read i;do
    rm $i-indels.vcf.idx
    rm $i-indels-snpeff.vcf.idx
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -A SnpEff \
    --variant $i-indels.vcf \
    --snpEffFile $i-indels-snpeff.vcf \
    -L $i-indels.vcf \
    -o $i-indels-annotated.vcf
    done < ../call-SNPs-recal03/samplelist
#### d.Merge split, annotated vcfs
    #compress then index split files
    #snps
    while read i;do
    ~/bin/samtools-1.1/htslib-1.1/bgzip -f $i-snps-annotated.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf $i-snps-annotated.vcf.gz
    done < ../call-SNPs-recal03/samplelist
    #indels
    while read i;do
    ~/bin/samtools-1.1/htslib-1.1/bgzip -f $i-indels-annotated.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf $i-indels-annotated.vcf.gz
    done < ../call-SNPs-recal03/samplelist
#### e.Merge vcf files then index
    #snps
    ~/bin/bcftools/bcftools merge \
    -o ALL-samples-snps-annotated.vcf.gz \
    -O z -m none \
     ../split-vcfs/*-snps-annotated.vcf.gz
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-snps-annotated.vcf.gz
    #indels
    ~/bin/bcftools/bcftools merge \
    -o ALL-samples-indels-annotated.vcf.gz \
    -O z -m none \
     ../split-vcfs/*-indels-annotated.vcf.gz
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-indels-annotated.vcf.gz
#### f.Get only high quality non-synonymous alleles
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir shared-variants
    cd shared-variants
    #script to automate unzipping bgzipped files
    #and get only  with non-synonymous variants
    #snps
    ## what the proceeding complicated code does:
    # 1.Reads in each line of samplelist passes that value to "i"
    # 2.Unzips the file ../split-vcfs/i-snps-annotated.vcf.gz
    # 3.Looks through the unzipped file stream, ignoring lines with "#" 
    #   but keeping lines with "SNPEFF_AMINO_ACID_CHANGE=CapitalLetterNumbersCapitalLetter"
    # 4.Grabs only the first, second, and tenth columns of text
    # 5.Converts haplotypes from format 1|1 to 1\t1 (where \t is a tab)
    # 6.Sets the delimiter as a tab,
    #   sets the 3rd column text as index 1 of array "a",
    #   sets the 4th col text as index 2 of array "a",
    #   sorts array "a" numerically, using gnu awk (gawk)
    #   prints a new tab-delimited line in the form:
    #   col 1, col 2, array a value 1, array a value 2
    #   note that col 3 and 4 are sorted numerically
    # 7.Converts columns 3 and 4 into a single column, and
    #   saves the file as i-nonsyn-snps-genotype.txt
    # 8.Adds the sample name to the first line
    # 9.Repeats steps 1-7 for all lines of samplelist
    while read i;do
    zcat ../split-vcfs/$i-snps-annotated.vcf.gz | \
    grep -v '#' | grep -P 'SNPEFF_AMINO_ACID_CHANGE=\w*[A-Z]\d+\w*[A-Z]' | \
    cut -f 1-2,10 | \
    perl -pe "s/(\w+\.\d)\t(\d+)\t(\d)\|(\d).+\n/\1\t\2\t\3\t\4\n/" | \
    gawk -v OFS='\t' '{a[1]=$3;a[2]=$4;asort(a);print $1,$2,a[1],a[2];}' - | \
    awk -v OFS=' ' '{b=$3$4;print $1,$2,b;}' - | \
    echo -e "$i\n$(cat - )" > $i-nonsyn-snps-genotype.txt
    done < ../call-SNPs-recal03/samplelist
    # then combine the files into one file with each file being a separate column
    paste *snps* > ALL-samples-nonsyn-snps-genotype.txt
    # how many SNP loci have non-synonymous variants
    wc -l CF219-nonsyn-snps-genotype.txt
    #3946 (after subtracting 1 for header)
    #indels
    while read i;do
    zcat ../split-vcfs/$i-indels-annotated.vcf.gz | \
    grep -v '#' | grep -P 'SNPEFF_AMINO_ACID_CHANGE=\w*[A-Z]\d+\w*[A-Z]' | \
    cut -f 1-2,10 | \
    perl -pe "s/(\w+\.\d)\t(\d+)\t(\d)\|(\d).+\n/\1\t\2\t\3\t\4\n/" | \
    gawk -v OFS='\t' '{a[1]=$3;a[2]=$4;asort(a);print $1,$2,a[1],a[2];}' - | \
    awk -v OFS=' ' '{b=$3$4;print $1,$2,b;}' - | \
    echo -e "$i\n$(cat - )" > $i-nonsyn-indels-genotype.txt
    done < ../call-SNPs-recal03/samplelist
    # then combine the files into one file with each file being a separate column
    paste *indels* > ALL-samples-nonsyn-indels-genotype.txt
    # how many indel loci have non-synonymous variants
    wc -l CF219-nonsyn-indels-genotype.txt
    #231 (after subtracting 1 for header)
    # phenotypes
    # Subclinical/Recovered (n=4): CF090, CF072, FC015, OLD065
    # Sick at last observation (n=6): CF080, CF219, FC013, FC047, OLD092, OLD107
    # Healthy (n=6): CF053, CF069, FC019, FC058, OLD077, OLD106
#### For getting snp alleles shared amongst phenotypes
#####http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#TOC-Venn-Diagrams
=====
## STEPS FOR LOOKING FOR SNPs UNDER SELECTION
### 1.Download tools
#### a.Download BayeScan
    cd ~/bin/
    wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
    unzip BayeScan2.1.zip
    mv BayeScan2.1.zip BayeScan2.1
    cd BayeScan2.1/
    cd binaries/
    chmod u+x BayeScan2.1_linux64bits # makes the file executable
    # Path to BayeScan
    /home/jelber2/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits
#### b.Download Simple Fool's Guide (SFG) to RNA-seq scripts to convert vcf file to BayeScan input format
    cd ~/scripts/immunome_2014/
    mkdir fromSFG
    cd fromSFG
    wget http://sfg.stanford.edu/Scripts.zip
    unzip Scripts.zip 
    mv Scripts\ for\ SFG/ Scripts_for_SFG
#### c. make directory
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir bayescan
    cd bayescan
### 2.Run for BayeScan for snps
#### a.Add Genotype Qualities to ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    zcat ../beagle/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz | perl -pe "s/(GT:DS:GP)/\1:GQ/" \
    > ALL-samples-Q30-snps-recal-beagle-polymorphic-fixed.vcf
    perl -pe "s/(\d\|\d:\d:\d,\d,\d)/\1:30/g" \
    ALL-samples-Q30-snps-recal-beagle-polymorphic-fixed.vcf \
    > ALL-samples-Q30-snps-recal-beagle-polymorphic-fixed2.vcf
#### b.Make populations.txt file
    text file with samplename\tpopulation (samplename tab population)
    CF219   CF
    CF53    CF
    CF69    CF
    CF72    CF
    CF80    CF
    CF90    CF
    FC13    FC
    FC15    FC
    FC19    FC
    FC47    FC
    FC58    FC
    OLD106  OLD
    OLD107  OLD
    OLD65   OLD
    OLD77   OLD
    OLD92   OLD
#### c.Make phenotypes.txt file
    text file with samplename\tphenotype (samplename tab phenotype)
    CF219   sick
    CF53    healthy
    CF69    healthy
    CF72    subclinical
    CF80    sick
    CF90    subclinical
    FC13    sick
    FC15    subclinical
    FC19    healthy
    FC47    sick
    FC58    healthy
    OLD106  healthy
    OLD107  sick
    OLD65   subclinical
    OLD77   healthy
    OLD92   sick
#### c2.Make phenotypes.txt file
    text file with samplename\tphenotype (samplename tab phenotype)
    CF219   sick
    CF53    healthy
    CF69    healthy
    CF72    sick
    CF80    sick
    CF90    sick
    FC13    sick
    FC15    sick
    FC19    healthy
    FC47    sick
    FC58    healthy
    OLD106  healthy
    OLD107  sick
    OLD65   sick
    OLD77   healthy
    OLD92   sick
#### d.Ran make_bayescan_input.py on snps and populations.txt
    #30 = min genotype quality
    #4 = min number of good quality genotype required from each population in order for a given SNP to be included in the analysis
    #1 = min number of copies of the minor allele that are necc. for a locus to be considered trustworthy enough to be used in BayeScan
    #1 = make outfile file (used_snp_genos.txt) showing what snp genotype were used
    #> = creates a file so you know the values for each population
    #output = bayes_input.tx, snpkey.txt, low_freq_snps.txt, used_snp_genos.txt
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    python ~/scripts/immunome_2014/fromSFG/Scripts_for_SFG/make_bayescan_input_using_phased_data.py \
    ../bayescan/ALL-samples-Q30-snps-recal-beagle-polymorphic-fixed2.vcf \
    populations.txt 30 4 1 1 > population-info.txt
    mv bayes_input.txt bayes_input.txt.snps.pop
    mv low_freq_snps.txt low_freq_snps.txt.snps.pop
    mv population-info.txt population-info.txt.snps.pop
    mv snpkey.txt snpkey.txt.snps.pop
    mv used_snp_genos.txt used_snp_genos.txt.snps.pop
#### e.Ran make_bayescan_input.py on snps and phenotypes2.txt
    #30 = min genotype quality
    #4 = min number of good quality genotype required from each population in order for a given SNP to be included in the analysis
    #1 = min number of copies of the minor allele that are necc. for a locus to be considered trustworthy enough to be used in BayeScan
    #1 = make outfile file (used_snp_genos.txt) showing what snp genotype were used
    #> = creates a file so you know the values for each population
    #output = bayes_input.tx, snpkey.txt, low_freq_snps.txt, used_snp_genos.txt
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    python ~/scripts/immunome_2014/fromSFG/Scripts_for_SFG/make_bayescan_input_using_phased_data.py \
    ../bayescan/ALL-samples-Q30-snps-recal-beagle-polymorphic-fixed2.vcf \
    phenotypes2.txt 30 4 1 1 > population-info.txt
    mv bayes_input.txt bayes_input.txt.snps.pheno
    mv low_freq_snps.txt low_freq_snps.txt.snps.pheno
    mv population-info.txt population-info.txt.snps.pheno
    mv snpkey.txt snpkey.txt.snps.pheno
    mv used_snp_genos.txt used_snp_genos.txt.snps.pheno
#### f.Run BayeScan on snps and pops (on SuperMikeII)
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /work/jelber2/immunome_urtd/combined/bayescan/bayes_input.txt.snps.pop \
    -snp \
    -d low_freq_snps.txt.snps.pop \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles.snps.pop \
    -threads 16
#####i.View Bayescan results for pop
    R
    #source the plot_R.r script from Bayescan
    source("/home/jelber2/bin/BayeScan2.1/R functions/plot_R_no_plot.r")
    #plot fst values without minor alleles below minor allele frequency of 1 copy
    noMAF_snps_results <- plot_bayescan("bayescan_no_loci_with_low_freq_minor_alleles.snps.pop_fst.txt", FDR=0.1)
    #save the candidate loci to a text file
    write(noMAF_snps_results$outliers, file= "noMAF_loci_FDR_0.1_outlier_snps.pop.txt", ncolumns= 1,append= FALSE)
    q()
#####ii.View Bayescan results in IGV
    #create a copy of snpkey.txt, so it can be modified
    cp snpkey.txt.snps.pop snpkey.txt.snps2.pop
    #code to create IGV batch file for noMAF loci
    while read i;do
    perl -pi -e "s/^$i\t(.+)_(.+)\n/goto \1:\2\n/" snpkey.txt.snps2.pop
    done < noMAF_loci_FDR_0.1_outlier_snps.pop.txt
    grep 'goto' snpkey.txt.snps2.pop > noMAF_loci_FDR_0.1_outlier_snps_pop_igv.txt
    #view in IGV
    ~/bin/IGV_2.3.40/igv.sh /media/immunome_2014/work/jelber2/immunome_urtd/split-vcfs/ALL-samples-snps-annotated.vcf.gz
    #open noMAF_loci_FDR_0.1_outlier_snps__pop_igv.txt
#####iii.Filter annotated VCF file by outlier snps
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    zcat ../split-vcfs/ALL-samples-snps-annotated.vcf.gz > ALL-samples-snps-pop-annotated2.vcf
    perl -pe "s/goto (\w+\.\d):(\d+)\n/\1\t\2\n/" noMAF_loci_FDR_0.1_outlier_snps_pop_igv.txt > noMAF_loci_FDR_0.1_outlier_snps_pop_vcf.txt
    while read i;do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/\1\tOUTLIER_SNP\t\2\n/" ALL-samples-snps-pop-annotated2.vcf
    done < noMAF_loci_FDR_0.1_outlier_snps_pop_vcf.txt
    grep 'OUTLIER_SNP\|^#' ALL-samples-snps-pop-annotated2.vcf | grep -v "contig" > ALL-samples-outlier-snps-pop.vcf
    # get only gene names, note that some SNPs are intergenic
    grep -v "#" ALL-samples-outlier-snps-pop.vcf | \
    perl -pe "s/.+SNPEFF_GENE_NAME=(\w+);.+\n/\1\n/" | \
    perl -pe "s/NW_.+\n/intergenic\n/g" | \
    sort | uniq -c | \
    perl -pe "s/( )+/\t/g" > ALL-samples-outlier-snps-pop-gene-names.txt
    # how many SNPs are under selection?
    perl -ane '$sum += $F[0]; END {print $sum; print "\n"}' ALL-samples-outlier-snps-pop-gene-names.txt
    # 2
    # how many genes have SNPs under selection
    grep -v "intergenic" ALL-samples-outlier-snps-pop-gene-names.txt | wc -l
    # note: the gff3 files says that both genic snps while snpeffect says only
    #       ones is genic?
    # 2
    # what genes?
    # interferon-induced protein with tetratricopeptide repeats 1-like (IFIT1-like)
    # interferon-induced protein 44-like (IFI44-like)
    # how many SNPs are intergenic
    grep "intergenic" ALL-samples-outlier-snps-pop-gene-names.txt | perl -ane '$sum += $F[0]; END {print $sum; print "\n"}'
    # 0 - based off of the gff3
    # how many SNPs are contained in genes
    grep -v "intergenic" ALL-samples-outlier-snps-pop-gene-names.txt | perl -ane '$sum += $F[0]; END {print $sum; print "\n"}'
    # 2 - based off of gff3 file
#### g.Run BayeScan on snps and pheno (on SuperMikeII)
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /work/jelber2/immunome_urtd/combined/bayescan/bayes_input.txt.snps.pheno \
    -snp \
    -d low_freq_snps.txt.snps.pheno \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles.snps.pheno \
    -threads 16
#####i.View Bayescan results for pheno
    #initiate R in the terminal
    R
    setwd("/media/immunome_2014/work/jelber2/immunome_urtd/bayescan/")
    #source the plot_R.r script from Bayescan
    source("/home/jelber2/bin/BayeScan2.1/R functions/plot_R_no_plot.r")
    #plot fst values without minor alleles below minor allele frequency of 1 copy
    noMAF_snps_results <- plot_bayescan("bayescan_no_loci_with_low_freq_minor_alleles.snps.pheno_fst.txt", FDR=0.1)
    #save the candidate loci to a text file
    write(noMAF_snps_results$outliers, file= "noMAF_loci_FDR_0.05_outlier_snps.pheno.txt", ncolumns= 1,append= FALSE)
    q()
    # no outliers! even using FDR=0.5!
### 3.Run BayeScan for indels
#### a.Add Genotype Qualities to ALL-samples-Q30-indels-recal-beagle-polymorphic.vcf.gz
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    zcat ../beagle/ALL-samples-Q30-indels-recal-beagle-polymorphic.vcf.gz | perl -pe "s/(GT:DS:GP)/\1:GQ/" \
    > ALL-samples-Q30-indels-recal-beagle-polymorphic-fixed.vcf
    perl -pe "s/(\d\|\d:\d:\d,\d,\d)/\1:30/g" \
    ALL-samples-Q30-indels-recal-beagle-polymorphic-fixed.vcf \
    > ALL-samples-Q30-indels-recal-beagle-polymorphic-fixed2.vcf
#### b.Ran make_bayescan_input.py on indels and populations.txt
    #30 = min genotype quality
    #4 = min number of good quality genotype required from each population in order for a given SNP to be included in the analysis
    #1 = min number of copies of the minor allele that are necc. for a locus to be considered trustworthy enough to be used in BayeScan
    #1 = make outfile file (used_snp_genos.txt) showing what snp genotype were used
    #> = creates a file so you know the values for each population
    #output = bayes_input.tx, snpkey.txt, low_freq_indels.txt, used_snp_genos.txt
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    python ~/scripts/immunome_2014/fromSFG/Scripts_for_SFG/make_bayescan_input_using_phased_data.py \
    ../bayescan/ALL-samples-Q30-indels-recal-beagle-polymorphic-fixed2.vcf \
    populations.txt 30 4 1 1 > population-info.txt
    mv bayes_input.txt bayes_input.txt.indels.pop
    mv low_freq_snps.txt low_freq_indels.txt.indels.pop
    mv population-info.txt population-info.txt.indels.pop
    mv snpkey.txt snpkey.txt.indels.pop
    mv used_snp_genos.txt used_snp_genos.txt.indels.pop
#### c.Ran make_bayescan_input.py on indels and phenotypes2.txt
    #30 = min genotype quality
    #4 = min number of good quality genotype required from each population in order for a given SNP to be included in the analysis
    #1 = min number of copies of the minor allele that are necc. for a locus to be considered trustworthy enough to be used in BayeScan
    #1 = make outfile file (used_snp_genos.txt) showing what snp genotype were used
    #> = creates a file so you know the values for each population
    #output = bayes_input.tx, snpkey.txt, low_freq_indels.txt, used_snp_genos.txt
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    python ~/scripts/immunome_2014/fromSFG/Scripts_for_SFG/make_bayescan_input_using_phased_data.py \
    ../bayescan/ALL-samples-Q30-indels-recal-beagle-polymorphic-fixed2.vcf \
    phenotypes2.txt 30 4 1 1 > population-info.txt
    mv bayes_input.txt bayes_input.txt.indels.pheno
    mv low_freq_snps.txt low_freq_indels.txt.indels.pheno
    mv population-info.txt population-info.txt.indels.pheno
    mv snpkey.txt snpkey.txt.indels.pheno
    mv used_snp_genos.txt used_snp_genos.txt.indels.pheno
#### d.Run BayeScan on indels and pops (on SuperMikeII)
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /work/jelber2/immunome_urtd/combined/bayescan/bayes_input.txt.indels.pop \
    -snp \
    -d low_freq_indels.txt.indels.pop \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles.indels.pop \
    -threads 16
#####i.View Bayescan results for pop
    #initiate R in the terminal
    R
    source("/home/jelber2/bin/BayeScan2.1/R functions/plot_R_no_plot.r")
    #plot fst values without minor alleles below minor allele frequency of 1 copy
    noMAF_indels_results <- plot_bayescan("bayescan_no_loci_with_low_freq_minor_alleles.indels.pop_fst.txt", FDR=0.1)
    #save the candidate loci to a text file
    write(noMAF_indels_results$outliers, file= "noMAF_loci_FDR_0.05_outlier_indels.pop.txt", ncolumns= 1,append= FALSE)
    q()
    # no outlier loci even using a FDR=0.2!
#### e.Run BayeScan on indels and pheno (on SuperMikeII)
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /work/jelber2/immunome_urtd/combined/bayescan/bayes_input.txt.indels.pheno \
    -snp \
    -d low_freq_indels.txt.indels.pheno \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles.indels.pheno \
    -threads 16
#####i.View Bayescan results for pheno
    #initiate R in the terminal
    R
    #source the plot_R.r script from Bayescan
    source("/home/jelber2/bin/BayeScan2.1/R functions/plot_R_no_plot.r")
    #plot fst values without minor alleles below minor allele frequency of 1 copy
    noMAF_indels_results <- plot_bayescan("bayescan_no_loci_with_low_freq_minor_alleles.indels.pheno_fst.txt", FDR=0.1)
    #save the candidate loci to a text file
    write(noMAF_indels_results$outliers, file= "noMAF_loci_FDR_0.05_outlier_indels.pheno.txt", ncolumns= 1,append= FALSE)
    q()
    # no outlier loci even using a FDR=0.4!
=====
## STEPS FOR LOOKING FOR Genes UNDER SELECTION
### Get FASTA read depth 20, min length 60
### Get intervals at least read depth 20
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir fasta-intervals
    mkdir fasta-seqs
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-intervals/
    #Use GATK Callableloci
    while read i;do
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T CallableLoci \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I ../call-SNPs-recal03/$i-recal03.bam \
    -L /media/immunome_2014/work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    --summary $i.callableloci.summary \
    --minBaseQuality 20 \
    --minMappingQuality 10 \
    --minDepth 20 \
    --minDepthForLowMAPQ 10 \
    --format BED \
    --out $i.callableloci
    #Use grep to get only "CALLABLE" intervals
    #Use bedtools merge to get contiguous regions
    grep "CALLABLE" $i.callableloci | ~/bin/bedtools-2.22.1/bin/bedtools merge > $i.callableloci.cont.bed
    #Use awk to get interval lengths
    awk -v OFS='\t' '{a=$3-$2;print $1,$2,$3,a;}' $i.callableloci.cont.bed > $i.callableloci.bylength.bed
    done < ../call-SNPs-recal03/samplelist
### Get only regions at least 60bp long using R
    R
    #set the working directory
    setwd("/media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-intervals/")
    #get the desired files
    print(files <- list.files(pattern=".bylength.bed$"))
    #convert files to samples (i.e., AL102 instead of AL102.callableloci.bylength.bed)
    print(samples <- gsub("prefixToTrash-0|\\.callableloci\\.bylength\\.bed",
          "", files, perl=TRUE), sep="")
    #for each sample:
    #1.Create the string to name the input file
    #2.Read in the file
    #3.Get only intervals of length 60 or greater
    #4.Create the string to name the output file
    #5.Write the output file
    for (i in samples){
      filein <- paste(i, ".callableloci.bylength.bed", sep="")
      i.callable <- read.table(filein)
      i.callablesixty = i.callable[i.callable$V4 > 59, ]
      fileout <- paste(i, ".callableloci.depth20.len60.bed", sep="")
      write.table(i.callablesixty,
                  file= fileout,
                  append = FALSE, quote = FALSE, sep = "\t",
                  eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                  col.names = FALSE)
    }
    quit()
### Use bedtools then grep to get regions shared among all samples
    ~/bin/bedtools-2.22.1/bin/bedtools multiinter -header -i *.callableloci.depth20.len60.bed > ALLsamples
    #use grep to get intervals shared by all samples (n total = 16)
    grep -P "\t16\t" ALLsamples | cut -f 1-3 > ALLsamples.callableloci.bed
    #get only intervals that are the same length for all samples
    ~/bin/bedtools-2.22.1/bin/bedtools intersect \
    -a ALLsamples.callableloci.bed \
    -b *.callableloci.depth20.len60.bed \
    -f 1.0 -r -u > ALLsamples.callableloci.samelength.bed
    #use awk to calculate lengths
    awk -v OFS='\t' '{a=$3-$2;print $1,$2,$3,a;}' ALLsamples.callableloci.samelength.bed \
    > ALLsamples.callableloci.bylength.bed
### Use R to get intervals >= 60 bp
    R
    setwd("/media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-intervals/")
    callable <- read.table("ALLsamples.callableloci.bylength.bed")
    callablesixty = callable[callable$V4 > 59, ]
    write.table(callablesixty,
                file= "ALLsamples.callableloci.depth20.len60.bed",
                append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE)
    quit()
### use awk to convert from 0-based to 1-based positions then use perl to convert format for samtools region
    awk -v OFS='\t' '{a=$2+1;print $1,a,$3,$4;}' ALLsamples.callableloci.depth20.len60.bed | 
    perl -pe "s/(\w+_\w+\.\d)\t(\d+)\t(\d+)\t\d+/\1:\2-\3/" > ../fasta-seqs/loci2filterbams.txt
### use GATK's FastaAlternateReferenceMaker to output consensus.fa with degenerate bases
    #first combine Q30 snps and indels after beagle genotype imputation
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/beagle/
    #concatentate snps and indels with GATK
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    --variant ALL-samples-Q30-snps-recal-beagle.vcf \
    --variant ALL-samples-Q30-indels-recal-beagle.vcf \
    -o ALL-samples-Q30-snps-indels-recal-beagle.vcf \
    --assumeIdenticalSamples \
    -genotypeMergeOptions UNSORTED
    #updated files with rsync
    rsync --stats --progress --archive /media/immunome_2014/work/jelber2/immunome_urtd/combined/beagle/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_urtd/combined/beagle/ -n
    rsync --stats --progress --archive /media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-seqs/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_urtd/combined/fasta-seqs/ -n
    #ran the following on SuperMikeII using create_fasta.sh
    #using 16 cores for at most 32 hours - actually took 6.5 hrs
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-seqs/
    ~/bin/parallel-20150122/src/parallel \
    'while read i;do
    java -Xmx2g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -T FastaAlternateReferenceMaker \
    -o $i.{}.fa \
    -L {} \
    --variant ../beagle/ALL-samples-Q30-snps-indels-recal-beagle.vcf \
    -IUPAC $i \
    --lineWidth 10000
    perl -pi -e "s/>.+\\n/>$i.{}\\n/" $i.{}.fa
    done < ../call-SNPs-recal03/samplelist
    cat *.{}.fa > {}.fa
    rm *.{}.fa' :::: loci2filterbams.txt
    #get files from SuperMikeII
    rsync --stats --progress --archive jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_urtd/combined/fasta-seqs/ \
    /media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-seqs/ -n
### add referenc/outgroup sequence
    ~/bin/parallel-20150122/src/parallel \
    'java -Xmx2g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T FastaReferenceMaker \
    -R /work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -o ref.{}.fa \
    -L {} \
    --lineWidth 10000
    perl -pi -e "s/>.+\\n/>ref.{}\\n/" ref.{}.fa
    cat ref.{}.fa >> {}.fa
    rm ref.{}.fa' :::: loci2filterbams.txt
### get PHASE
    cd ~/bin/
    #note had to manually download and place in ~/bin/ because had to enter a password for download
    wget http://stephenslab.uchicago.edu/phase/phasecode/phase.2.1.1.linux.tar.gz
    tar xzf phase.2.1.1.linux.tar.gz 
    mv phase.2.1.1.linux.tar.gz phase.2.1.1.linux
### get SeqPHASE to generate PHASE input file and process PHASE output
    cd ~/bin/
    mkdir seqphase
    cd seqphase/
    wget http://seqphase.mpg.de/seqphase/seqphase2014.zip
    unzip seqphase2014.zip
### get muscle
    cd ~/bin/
    mkdir muscle-3.8.31
    cd muscle-3.8.31/
    wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
    tar xzf muscle3.8.31_i86linux64.tar.gz
### run muscle,SeqPHASE,PHASE,then SeqPHASE using gnu parallel
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-seqs/
    ~/bin/parallel-20150122/src/parallel \
    '~/bin/muscle-3.8.31/muscle3.8.31_i86linux64 -in {}.fa -out {}.fa.aln
    ~/bin/seqphase/seqphase1.pl -1 {}.fa.aln -p {}
    ~/bin/phase.2.1.1.linux/PHASE {}.inp {}.out
    ~/bin/seqphase/seqphase2.pl -c {}.const -i {}.out_pairs -o {}.fa.phased' :::: loci2filterbams.txt
#### Rename fasta headers to CF219a and CF219b
    perl -pi -e "s/>(\w+).\w+_\d+.\d:\d+-\d+(\w)_.+/>\1\2/" *.fa.phased
#### created file phased for file conversion below
    ls *.phased | perl -pe "s/(.+).fa.phased/\1/" > phased
### Results
    #How many regions are there among the 16 samples that are at least 60bp and
    #have a read depth of 20 reads per individual?
    ls *.fa | wc -l
    #1680
    #How many of these regions are polymorphic, and thus phaseable?
    #Also includes how many regions have a ref seq the same length as the tortoises
    ls *.fa.phased | wc -l
    #1556
    #How many genes and exons do these 1764 regions represent
    ls *.fa.phased | perl -pe "s/(\w+_\d+\.\d):(\d+)-(\d+)\.fa\.phased\n/\1\t\2\t\3\n/g" |
    awk -v OFS='\t' '{a=$2-1;print $1,a,$3;}' - | sort -k 1,1 -k2,2n > ALLsamples.callableloci.depth20.len60.poly.bed
    #use bedtools to intersect the gff and bed file
    ~/bin/bedtools-2.22.1/bin/bedtools intersect \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.introns \
    -b ALLsamples.callableloci.depth20.len60.poly.bed -wao > ALLsamples.callableloci.depth20.len60.poly2.bed
    grep -v "RefSeq" ALLsamples.callableloci.depth20.len60.poly2.bed | grep -Pv ".\t-1\t-1\t0" | \
    grep -Pv "Gnomon\tmRNA" | grep -Pv "Gnomon\tCDS" | \
    grep -Pv "Gnomon\ttranscript" > ALLsamples.callableloci.depth20.len60.poly3.bed
### Run PGDSpider to convert FASTA to PHYLIP
#### make selection folder
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir selection
    cd selection
#### now make the spid template
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-seqs/
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-seqs/NW_007359913.1:1047555-1047630.fa.phased \
    -inputformat FASTA \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/selection/NW_007359913.1:1047555-1047630.fa.phased.phy \
    -outputformat PHYLIP
#####Manually edit ../selection/template_FASTA_PHYLIP.spid
    # spid-file generated: Tue Jul 21 21:31:38 CDT 2015
    # FASTA Parser questions
    PARSER_FORMAT=FASTA
    # Select the type of the data:
    FASTA_PARSER_DATA_TYPE_QUESTION=DNA
    # PHYLIP (RAxML) Writer questions
    WRITER_FORMAT=PHYLIP
    # Numeric SNP data: enter the integer that codes for nucleotide T:
    PHYLIP_WRITER_CODE_T_QUESTION=
    # Numeric SNP data: enter the integer that codes for nucleotide C:
    PHYLIP_WRITER_CODE_C_QUESTION=
    # Numeric SNP data: enter the integer that codes for nucleotide A:
    PHYLIP_WRITER_CODE_A_QUESTION=
    # Save relaxed PHYLIP format (e.g. for RAxML)?
    PHYLIP_WRITER_RELAXED_QUESTION=
    # Select the kind of file you want to write:
    PHYLIP_WRITER_DATA_KIND_QUESTION=
    # Numeric SNP data: enter the integer that codes for nucleotide G:
    PHYLIP_WRITER_CODE_G_QUESTION=
    # Specify the DNA locus you want to write to the PHYLIP (RAxML) file or write "CONCAT" for concatenation:
    ARLEQUIN_WRITER_ONCATENATE_QUESTION=
#### Save as fasta2phylip.spid
### Convert all *.fa.phased to *.phy
    #note used newer version of PGDSPider (i.e., 2.0.8.3)
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/selection/
    while read i;do
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-seqs/$i.fa.phased \
    -inputformat FASTA \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/selection/$i.fa.phased.phy \
    -outputformat PHYLIP \
    -spid /media/immunome_2014/work/jelber2/immunome_urtd/combined/selection/fasta2phylip.spid
    done < /media/immunome_2014/work/jelber2/immunome_urtd/combined/fasta-seqs/phased
### Get and install VariScan
    cd ~/bin/
    wget http://www.ub.es/softevol/variscan/variscan-2.0.3.tar.gz
    tar xzf variscan-2.0.3.tar.gz
    mv variscan-2.0.3.tar.gz variscan-2.0.3
    cd variscan-2.0.3/
   ./configure
    make distclean
    ./configure
    make
#### Make VariScan config files for all populations and combinations
    #note actual files do not have leading four spaces on each line below
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/selection
    #CF.conf
    StartPos = 1
    EndPos = 0
    RefPos = 0
    BlockDataFile = none
    SeqChoice = 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    Outgroup = none
    RunMode = 12
    UseMuts = 1
    UseLDSinglets = 0
    CompleteDeletion = 1
    FixNum = 1
    NumNuc = 4
    SlidingWindow = 0
    WidthSW = 10
    JumpSW = 10
    WindowType = 0
    IndivNames = 
    RefSeq = 1
    #CF.outgroup.conf
    StartPos = 1
    EndPos = 0
    RefPos = 0
    BlockDataFile = none
    SeqChoice = 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
    Outgroup = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
    RunMode = 22
    UseMuts = 1
    UseLDSinglets = 0
    CompleteDeletion = 1
    FixNum = 1
    NumNuc = 4
    SlidingWindow = 0
    WidthSW = 10
    JumpSW = 10
    WindowType = 0
    IndivNames = 
    RefSeq = 1
    #FC.conf
    StartPos = 1
    EndPos = 0
    RefPos = 0
    BlockDataFile = none
    SeqChoice = 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0
    Outgroup = none
    RunMode = 12
    UseMuts = 1
    UseLDSinglets = 0
    CompleteDeletion = 1
    FixNum = 1
    NumNuc = 4
    SlidingWindow = 0
    WidthSW = 10
    JumpSW = 10
    WindowType = 0
    IndivNames = 
    RefSeq = 1
    #FC.outgroup.conf
    StartPos = 1
    EndPos = 0
    RefPos = 0
    BlockDataFile = none
    SeqChoice = 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1
    Outgroup = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
    RunMode = 22
    UseMuts = 1
    UseLDSinglets = 0
    CompleteDeletion = 1
    FixNum = 1
    NumNuc = 4
    SlidingWindow = 0
    WidthSW = 10
    JumpSW = 10
    WindowType = 0
    IndivNames = 
    RefSeq = 1
    #OLD.conf
    StartPos = 1
    EndPos = 0
    RefPos = 0
    BlockDataFile = none
    SeqChoice = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0
    Outgroup = none
    RunMode = 12
    UseMuts = 1
    UseLDSinglets = 0
    CompleteDeletion = 1
    FixNum = 1
    NumNuc = 4
    SlidingWindow = 0
    WidthSW = 10
    JumpSW = 10
    WindowType = 0
    IndivNames = 
    RefSeq = 1
    #OLD.outgroup.conf
    StartPos = 1
    EndPos = 0
    RefPos = 0
    BlockDataFile = none
    SeqChoice = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1
    Outgroup = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
    RunMode = 22
    UseMuts = 1
    UseLDSinglets = 0
    CompleteDeletion = 1
    FixNum = 1
    NumNuc = 4
    SlidingWindow = 0
    WidthSW = 10
    JumpSW = 10
    WindowType = 0
    IndivNames = 
    RefSeq = 1
#### Run VariScan
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/selection/
    #first make file to read file names
    ls *.phy | perl -pe "s/(.+).fa.phased.phy/\1/" > phylip
    #get headers
    ~/bin/variscan-2.0.3/src/variscan NC_024218.1:27914240-27914483.fa.phased.phy CF.conf | \
    grep "Eta" | perl -pe "s/#/chr/" > header.no.og
    ~/bin/variscan-2.0.3/src/variscan NC_024218.1:27914240-27914483.fa.phased.phy CF.outgroup.conf | \
    grep "Eta" | perl -pe "s/#/chr/" > header.og
    #CF
    while read i;do
    ~/bin/variscan-2.0.3/src/variscan $i.fa.phased.phy CF.conf | \
    grep -P "^#( )+" | grep -v "Eta" | \
    perl -pe "s/( )+/\t/g" > $i.CF.vs
    echo $i | cat - $i.CF.vs > temp && mv temp $i.CF.vs
    perl -pi -e "s/(\w+_\d+.\d:\d+-\d+)\n/\1/" $i.CF.vs
    perl -pi -e "s/#//" $i.CF.vs
    done < phylip
    cat header.no.og *.CF.vs > CF.out
    rm *.vs
    #CF.outgroup
    while read i;do
    ~/bin/variscan-2.0.3/src/variscan $i.fa.phased.phy CF.outgroup.conf | \
    grep -P "^#( )+" | grep -v "Eta" | \
    perl -pe "s/( )+/\t/g" > $i.CF.outgroup.vs
    echo $i | cat - $i.CF.outgroup.vs > temp && mv temp $i.CF.outgroup.vs
    perl -pi -e "s/(\w+_\d+.\d:\d+-\d+)\n/\1/" $i.CF.outgroup.vs
    perl -pi -e "s/#//" $i.CF.outgroup.vs
    done < phylip
    cat header.og *.CF.outgroup.vs > CF.outgroup.out
    rm *.outgroup.vs
    #FC
    while read i;do
    ~/bin/variscan-2.0.3/src/variscan $i.fa.phased.phy FC.conf | \
    grep -P "^#( )+" | grep -v "Eta" | \
    perl -pe "s/( )+/\t/g" > $i.FC.vs
    echo $i | cat - $i.FC.vs > temp && mv temp $i.FC.vs
    perl -pi -e "s/(\w+_\d+.\d:\d+-\d+)\n/\1/" $i.FC.vs
    perl -pi -e "s/#//" $i.FC.vs
    done < phylip
    cat header.no.og *.FC.vs > FC.out
    rm *.vs
    #FC.outgroup
    while read i;do
    ~/bin/variscan-2.0.3/src/variscan $i.fa.phased.phy FC.outgroup.conf | \
    grep -P "^#( )+" | grep -v "Eta" | \
    perl -pe "s/( )+/\t/g" > $i.FC.outgroup.vs
    echo $i | cat - $i.FC.outgroup.vs > temp && mv temp $i.FC.outgroup.vs
    perl -pi -e "s/(\w+_\d+.\d:\d+-\d+)\n/\1/" $i.FC.outgroup.vs
    perl -pi -e "s/#//" $i.FC.outgroup.vs
    done < phylip
    cat header.og *.FC.outgroup.vs > FC.outgroup.out
    rm *.outgroup.vs
    #OLD
    while read i;do
    ~/bin/variscan-2.0.3/src/variscan $i.fa.phased.phy OLD.conf | \
    grep -P "^#( )+" | grep -v "Eta" | \
    perl -pe "s/( )+/\t/g" > $i.OLD.vs
    echo $i | cat - $i.OLD.vs > temp && mv temp $i.OLD.vs
    perl -pi -e "s/(\w+_\d+.\d:\d+-\d+)\n/\1/" $i.OLD.vs
    perl -pi -e "s/#//" $i.OLD.vs
    done < phylip
    cat header.no.og *.OLD.vs > OLD.out
    rm *.vs
    #OLD.outgroup
    while read i;do
    ~/bin/variscan-2.0.3/src/variscan $i.fa.phased.phy OLD.outgroup.conf | \
    grep -P "^#( )+" | grep -v "Eta" | \
    perl -pe "s/( )+/\t/g" > $i.OLD.outgroup.vs
    echo $i | cat - $i.OLD.outgroup.vs > temp && mv temp $i.OLD.outgroup.vs
    perl -pi -e "s/(\w+_\d+.\d:\d+-\d+)\n/\1/" $i.OLD.outgroup.vs
    perl -pi -e "s/#//" $i.OLD.outgroup.vs
    done < phylip
    cat header.og *.OLD.outgroup.vs > OLD.outgroup.out
    rm *.outgroup.vs
#### Use R to combine CF.out CF.outgroup, etc.
    R
    #read in the neutrality stats into R
    #stats calculated without outgroup
    CF <-  read.table("CF.out",header=TRUE)
    #stats calculated with outgroup
    CF2 <- read.table("CF.outgroup.out",header=TRUE)
    #combine the two into 1 data.frame
    CF <- cbind(CF,CF2$FuLi_D,CF2$FuLi_F,CF2$FayWu_H)
    write.table(x=CF,
                file="CF.neutrality.stats.txt",
                sep = "\t",
                eol= "\n",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE,
                na = "NA",
                append = FALSE,
                dec = ".")
    #for FC samples
    FC <-  read.table("FC.out",header=TRUE)
    FC2 <- read.table("FC.outgroup.out",header=TRUE)
    FC <- cbind(FC,FC2$FuLi_D,FC2$FuLi_F,FC2$FayWu_H)
    write.table(x=FC,
                file="FC.neutrality.stats.txt",
                sep = "\t",
                eol= "\n",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE,
                na = "NA",
                append = FALSE,
                dec = ".")
    #for OLD samples
    OLD <-  read.table("OLD.out",header=TRUE)
    OLD2 <- read.table("OLD.outgroup.out",header=TRUE)
    OLD <- cbind(OLD,OLD2$FuLi_D,OLD2$FuLi_F,OLD2$FayWu_H)
    write.table(x=OLD,
                file="OLD.neutrality.stats.txt",
                sep = "\t",
                eol= "\n",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE,
                na = "NA",
                append = FALSE,
                dec = ".")
    quit()
#### Use cat to combine neutrality.stats.txt files
    cat *neutrality* > neutrality_stats.txt
    # remove extra column headers
    grep -v "chr" neutrality_stats.txt > neutrality_stats.txt2
    # get the header
    grep "chr" CF.neutrality.stats.txt > header
    # manually remove CF2$
    nano header
    # combine header and neutrality_stats.txt2
    cat header neutrality_stats.txt2 > neutrality_stats.txt
## What Genes have extreme values for Tajima's D?
### Use ms to simulate the coalescent and extreme values
    # see https://www.biostars.org/p/12227/#12233
    # for program http://home.uchicago.edu/rhudson1/source/mksamples.html
    # download manually by clicking the link
    tar xzf ms.tar.gz
    #saved tar and pdf file to ~/bin/msdir/
    mv ms.tar.gz msdir/.
    mv msdoc.pdf msdir/.
    # compile ms with random number generator #1
    cd msdir/
    gcc -O3 -o ms ms.c streec.c rand1.c -lm
    # compile sample_stats
    gcc -o sample_stats sample_stats.c tajd.c -lm
    #path to ms
    ~/bin/msdir/ms
    #path to sample_stats
    ~/bin/msdir/sample_stats
    # simulate 1-138 segragating sites from population of 5 individuals
    # with 1000 simulations
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir ms
    cd ms/
    for i in `seq 1 138`;do
    ~/bin/msdir/ms 5 1000 -s $i | ~/bin/msdir/sample_stats > $i.tsv
    done
    cat *.tsv > ms_simulations.txt
    rm *.tsv
### use R script to calculate p values for Tajima's D
    #neutrality_tests.R
    setwd("C:/Users/jelber2/Dropbox/LSU/Dissertation/Manuscripts/immunome_URTD/")
    ms.sims <- read.table("ms_simulations.txt") #read in the simulated data
    AL <- read.table("neutrality_stats.txt", header = TRUE) #read in the empirical data
    AL <- AL[AL$S != 0,] #get rid of empirical data that have 0 segragating sites
    AL <- AL[with(AL, order(S, chr)), ] #sort the empirical data first by increasing segragating sites then by increasing chromosome
    #
    AL2 <- "" #make an empty list AL2
    for (i in 1:138) { #go through the values of S, 1-138, S=# of segragating sites
      a <- ms.sims[ms.sims$V4 == i,] #get rows where the fourth column (S) is equal to "i" for the simulated data
      b <- a$V6 #get the values for Tajima's D for the ith value of S
      den.fun <- approxfun(density(b)) #get the density function for the ith value of S
      den.fun.approx <- approx(density(b)) #get the density values for the ith value of S
      den.fun.max <- max(den.fun.approx$x) #get the max value for the ith value of S
      den.fun.min <- min(den.fun.approx$x) #get the min value for the ith value of S
      pFmax <- function(q) integrate(den.fun, q, den.fun.max)$value #calculate the one-sided p value for values > 0
      pFmin <- function(q) integrate(den.fun, den.fun.min, q)$value #calculate the one-sided p value for values < 0
      y <- "" #make and empty list or clear y
      c <- AL[AL$S == i,] #get only rows where the values of S are "i" for the empirical data
      d <- c$Tajima_D #get only the 
      for (j in d) { #for each value of Tajima's D,
        if (is.na(j)) y <- rbind(y,"NA") #if the value is NA, then put an NA
        else { #else 
          if (j < 0) { #if the value of Tajima's D is negative,
            res <- try(y <- rbind(y, pFmin(j))) #determine if you can calculate the integral without an error
            if(inherits(res, "try-error")) { #if there is an error, 
              y <- rbind(y, 0) #then put p-value as "0"
            }
          }
          else { #if the value of Tajima's D is positive,
            res <- try(y <- rbind(y, pFmax(j))) #determine if you can calculate the integral without an error
            if(inherits(res, "try-error")) { #if there is an error,
              y <- rbind(y, 0) #then put the p-value as "0"
            }
          }
        }
      }
      y <- y[-1] #get rid of the empty first value
      AL2 <- c(AL2,y) #combine AL2 and y
    }
    AL2 = AL2[-1] #get rid of the first empty value of AL2
    AL3 <- p.adjust(AL2,method = "fdr",n = length(AL2)) #adjust the p-values using false discovery rate method
    AL <- cbind(AL, "prob_Tajima_D"= AL3) #add the AL3 data as a column to the AL data frame as "prob_Tajima_D"
    #
    #
    BL2 <- "" #make an empty list BL2
    for (i in 1:138) { #go through the values of S, 1-138, S=# of segragating sites
      a <- ms.sims[ms.sims$V4 == i,] #get rows where the fourth column (S) is equal to "i" for the simulated data
      b <- a$V10 #get the values for Way&Fu'sH for the ith value of S
      den.fun <- approxfun(density(b)) #get the density function for the ith value of S
      den.fun.approx <- approx(density(b)) #get the density values for the ith value of S
      den.fun.max <- max(den.fun.approx$x) #get the max value for the ith value of S
      den.fun.min <- min(den.fun.approx$x) #get the min value for the ith value of S
      pFmax <- function(q) integrate(den.fun, q, den.fun.max)$value #calculate the one-sided p value for values > 0
      pFmin <- function(q) integrate(den.fun, den.fun.min, q)$value #calculate the one-sided p value for values < 0
      y <- "" #make and empty list or clear y
      c <- AL[AL$S == i,] #get only rows where the values of S are "i" for the empirical data
      d <- c$FayWu_H #get only the 
      for (j in d) { #for each value of Way&Fu'sH,
        if (is.na(j)) y <- rbind(y,"NA") #if the value is NA, then put an NA
        else { #else 
          if (j < 0) { #if the value of Way&Fu'sH is negative,
            res <- try(y <- rbind(y, pFmin(j))) #determine if you can calculate the integral without an error
            if(inherits(res, "try-error")) { #if there is an error, 
              y <- rbind(y, 0) #then put p-value as "0"
            }
          }
          else { #if the value of Way&Fu'sH is positive,
            res <- try(y <- rbind(y, pFmax(j))) #determine if you can calculate the integral without an error
            if(inherits(res, "try-error")) { #if there is an error,
              y <- rbind(y, 0) #then put the p-value as "0"
            }
          }
        }
      }
      y <- y[-1] #get rid of the empty first value
      BL2 <- c(BL2,y) #combine BL2 and y
    }
    BL2 = BL2[-1] #get rid of the first empty value of BL2
    BL3 <- p.adjust(BL2,method = "fdr",n = length(BL2)) #adjust the p-values using false discovery rate method
    AL <- cbind(AL, "prob_Way_Fu_H"= BL3) #add the BL3 data as a column to the AL data frame as "prob_Way_Fu_H"
    #
    AL4 <- AL[AL$prob_Tajima_D != 0,] #get rid of values of prob ==0 because you couldn't properly calculate the integral
    #
    sig.TD <- AL4[AL4$prob_Tajima_D < 0.05,] #how many regions have Tajima's D values that are significant
    nrow(sig.TD)
    #70 gene regions but including duplicates because 3 populations analyzed separately
    #
    nrow(sig.TD[sig.TD$Tajima_D < 0,])
    #27 gene regions with negative Tajima's D
    # these are regions with rare alleles present at low frequencies indicating
    # Recent selective sweep, population expansion after a recent bottleneck, linkage to a swept gene
    #
    nrow(sig.TD[sig.TD$Tajima_D > 0,])
    #43 gene regions with positive Tajima's D
    # these are regions with multiple alleles present, some at low, others at high frequencies
    # Balancing selection, sudden population contraction
    #
    #
    unique.sig.TD <- unique(sig.TD$chr) #how many sig. gene regions are unique
    length(unique.sig.TD)
    #50 unique regions deviating from neutral expectations
    #
    #but how many genes? #use bedtools intersect
    #
    write.table(x=sig.TD,
                file="TajimaD.sig.txt",
                sep = "\t",
                eol= "\n",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE,
                na = "NA",
                append = FALSE,
                dec = ".")
    quit()
#### convert to BED format for bedtools
    cut -f 1 TajimaD.sig.txt | grep -v "chr" | sort -u | \
    perl -pe "s/(\w+_\d+.\d):(\d+)-(\d+)\n/\1\t\2\t\3\n/" > TajimaD.extr.regions.bed
#### bedtools intersect
    ~/bin/bedtools-2.22.1/bin/bedtools intersect \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.introns \
    -b TajimaD.extr.regions.bed -wao > TajimaD.extr.regions.bed.overlap
    grep -v "RefSeq" TajimaD.extr.regions.bed.overlap | grep -Pv ".\t-1\t-1\t0" | \
    grep -Pv "Gnomon\tmRNA" | grep -Pv "Gnomon\tCDS" | \
    grep -Pv "Gnomon\ttranscript" > TajimaD.extr.regions.bed.overlap.genes_exons_introns
    # get only genes
    grep "gene" TajimaD.extr.regions.bed.overlap.genes_exons_introns | \
    grep -v "exon" > TajimaD.extr.regions.bed.overlap.genes
    #get only unique genes
    perl -pe "s/.+Name=(\w+);.+\n/\1\n/" TajimaD.extr.regions.bed.overlap.genes | \
    perl -pe "s/.+gene=(\w+);.+\n/\1\n/" | sort -u > TajimaD.extr.regions.bed.overlap.genes.unique
    #how many regions have extreme values for Tajima's D?
    wc -l TajimaD.extr.regions.bed
    #50
    #how many genes have extreme values for Tajima's D?
    wc -l TajimaD.extr.regions.bed.overlap.genes.unique
    #35
#### Get protein accessions
    # get protein accessions
    grep -Pv ".\t-1\t-1\t0" TajimaD.extr.regions.bed.overlap | \
    grep "XP_" | \
    cut -f 9 | \
    perl -pe "s/.+Name=(XP_\d+\.\d).+/\1/" | \
    sort | uniq > TajimaD.extr.regions.bed.overlap.proteins.txt
=====
## STEPS FOR POPGEN
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir popgen
### 1.Make copy of ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/beagle/
    cp ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz ../popgen/.
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/
    gunzip ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz
### 2.Get only non-synonymous SNPs
    zcat ../split-vcfs/ALL-samples-snps-annotated.vcf.gz| grep '#' > header
    zcat ../split-vcfs/ALL-samples-snps-annotated.vcf.gz| grep -P 'SNPEFF_AMINO_ACID_CHANGE=\w*[A-Z]\d+\w*[A-Z]' > non-synonymous
    cat header non-synonymous > ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf
### 3.Check loci for linkage disequilibrium and Hardy-Weinberg Equilibrium
#### Had to make populations.txt file and population-specific files for vcf filtering
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/
    grep "CF" ../bayescan/populations.txt > CF
    grep "OLD" ../bayescan/populations.txt > OLD
    grep "FC" ../bayescan/populations.txt > FC
    cp ../bayescan/populations.txt .
#### Hardy-Weinberg Equilibrium test
    #all snps
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.CF --keep CF
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.OLD --keep OLD
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.FC --keep FC
    #only non-synonymous snps
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf --hardy --out hwe.CF-nonsyn --keep CF
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf --hardy --out hwe.OLD-nonsyn --keep OLD
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf --hardy --out hwe.FC-nonsyn --keep FC
#####R function to count number of sites out of HWE (i.e., p_HWE < 0.05 after fdr correction)
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/
    #all snps
    R
    CFhwe <-read.table(file ="hwe.CF.hwe", header = TRUE)
    CFhwe.fdr <- p.adjust(p = CFhwe$P_HWE, method = "fdr", n = length(CFhwe$P_HWE))
    summary(CFhwe.fdr)
    FChwe <-read.table(file ="hwe.FC.hwe", header = TRUE)
    FChwe.fdr <- p.adjust(p = FChwe$P_HWE, method = "fdr", n = length(FChwe$P_HWE))
    summary(FChwe.fdr)
    OLDhwe <-read.table(file ="hwe.OLD.hwe", header = TRUE)
    OLDhwe.fdr <- p.adjust(p = OLDhwe$P_HWE, method = "fdr", n = length(OLDhwe$P_HWE))
    summary(OLDhwe.fdr)
    #only nonsyn snps
    CFhwe <-read.table(file ="hwe.CF-nonsyn.hwe", header = TRUE)
    CFhwe.fdr <- p.adjust(p = CFhwe$P_HWE, method = "fdr", n = length(CFhwe$P_HWE))
    summary(CFhwe.fdr)
    FChwe <-read.table(file ="hwe.FC-nonsyn.hwe", header = TRUE)
    FChwe.fdr <- p.adjust(p = FChwe$P_HWE, method = "fdr", n = length(FChwe$P_HWE))
    summary(FChwe.fdr)
    OLDhwe <-read.table(file ="hwe.OLD-nonsyn.hwe", header = TRUE)
    OLDhwe.fdr <- p.adjust(p = OLDhwe$P_HWE, method = "fdr", n = length(OLDhwe$P_HWE))
    summary(OLDhwe.fdr)
    # no snps or nonsyn snps out of Hardy-Weinberg Equilibrium
#### outputs linkage disequilibrium pvalues
    #all snps
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --geno-chisq --out geno.CF --keep CF
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --geno-chisq --out geno.FC --keep FC
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --geno-chisq --out geno.OLD --keep OLD
    #only nonsyn snps
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf --geno-chisq --out geno.CF-nonsyn --keep CF
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf --geno-chisq --out geno.FC-nonsyn --keep FC
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf --geno-chisq --out geno.OLD-nonsyn --keep OLD
#####R function to count number of sites out of linkage equi (i.e., PVAL < 0.05 after fdr correction)
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/
    #all snps
    R
    FCld = na.omit(read.table(file ="geno.FC.geno.chisq", header = TRUE))
    FCld.fdr=p.adjust(p = FCld$PVAL, method = "fdr", n = length(FCld$PVAL))
    summary(FCld.fdr)
    CFld = na.omit(read.table(file ="geno.CF.geno.chisq", header = TRUE))
    CFld.fdr=p.adjust(p = CFld$PVAL, method = "fdr", n = length(CFld$PVAL))
    summary(CFld.fdr)
    OLDld = na.omit(read.table(file ="geno.OLD.geno.chisq", header = TRUE))
    OLDld.fdr=p.adjust(p = OLDld$PVAL, method = "fdr", n = length(OLDld$PVAL))
    summary(OLDld.fdr)
    #only nonsyn snps
    R
    FCld = na.omit(read.table(file ="geno.FC-nonsyn.geno.chisq", header = TRUE))
    FCld.fdr=p.adjust(p = FCld$PVAL, method = "fdr", n = length(FCld$PVAL))
    summary(FCld.fdr)
    CFld = na.omit(read.table(file ="geno.CF-nonsyn.geno.chisq", header = TRUE))
    CFld.fdr=p.adjust(p = CFld$PVAL, method = "fdr", n = length(CFld$PVAL))
    summary(CFld.fdr)
    OLDld = na.omit(read.table(file ="geno.OLD-nonsyn.geno.chisq", header = TRUE))
    OLDld.fdr=p.adjust(p = OLDld$PVAL, method = "fdr", n = length(OLDld$PVAL))
    summary(OLDld.fdr)
    # no snps or nonsyn snps out of Linkage Equilibrium
### 3.Convert VCF file to structure, fstat
#### VCF to Structure
#####use PGDSpider
    #get the program
    cd ~/bin/
    wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.8.3.zip
    unzip PGDSpider_2.0.8.3.zip
    mv PGDSpider_2.0.8.3.zip PGDSpider_2.0.8.3
    #use following command to generate spider.conf.xml and spid file in ../popgen directory
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/structure-input.txt \
    -outputformat STRUCTURE
    #edit spider.conf.xml
    nano ~/bin/PGDSpider_2.0.8.3/spider.conf.xml #to add path to samtools
    #change <entry key="PathBcftools"></entry>
    #to <entry key="PathBcftools">/home/jelber2/bin/samtools-0.1.19/bcftools/bcftools</entry>
    #change <entry key="PathSamtools"></entry>
    #to <entry key="PathSamtools">/home/jelber2/bin/samtools-0.1.19/samtools</entry>
    #save and exit
    #edit the spid file
    nano /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/template_VCF_STRUCTURE.spid
    #contents after editing, minus the leading spaces
        # spid-file generated: Wed Jan 14 18:47:06 CST 2015
        # VCF Parser questions
        PARSER_FORMAT=VCF
        # Do you want to include a file with population definitions?
        VCF_PARSER_POP_QUESTION=true
        # Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
        VCF_PARSER_REGION_QUESTION=
        # What is the ploidy of the data?
        VCF_PARSER_PLOIDY_QUESTION=DIPLOID
        # Only output following individuals (ind1, ind2, ind4, ...):
        VCF_PARSER_IND_QUESTION=
        # Output genotypes as missing if the read depth of a position for the sample is below:
        VCF_PARSER_READ_QUESTION=
        # Take most likely genotype if "PL" or "GL" is given in the genotype field?
        VCF_PARSER_PL_QUESTION=
        # Do you want to exclude loci with only missing data?
        VCF_PARSER_EXC_MISSING_LOCI_QUESTION=
        # Select population definition file:
        VCF_PARSER_POP_FILE_QUESTION=/media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/populations.txt
        # Only output SNPs with a phred-scaled quality of at least:
        VCF_PARSER_QUAL_QUESTION=
        # Do you want to include non-polymorphic SNPs?
        VCF_PARSER_MONOMORPHIC_QUESTION=false
        # Output genotypes as missing if the phred-scale genotype quality is below:
        VCF_PARSER_GTQUAL_QUESTION=
        #
        # STRUCTURE Writer questions
        WRITER_FORMAT=STRUCTURE
        # Save more specific fastSTRUCTURE format?
        STRUCTURE_WRITER_FAST_FORMAT_QUESTION=false
        # Specify the locus/locus combination you want to write to the STRUCTURE file:
        STRUCTURE_WRITER_LOCUS_COMBINATION_QUESTION=
        # Specify which data type should be included in the STRUCTURE file  (STRUCTURE can only analyze one data type per file):
        STRUCTURE_WRITER_DATA_TYPE_QUESTION=SNP
        # Do you want to include inter-marker distances?
        STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=false
    #saved as vcf2structure.spid
    #edit vcf2structure.spid
    nano vcf2structure.spid
    #change the following line from populations.txt to phenotypes.txt
        VCF_PARSER_POP_FILE_QUESTION=/media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/phenotypes.txt
#####Do VCF to STRUCTURE file conversion
    #all snps using populations.txt
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/structure-input-allsnps.txt \
    -outputformat STRUCTURE \
    -spid /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/vcf2structure.spid > structure-input-allsnps.log
    #all snps using phenotypes.txt
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/structure-input-pheno.txt \
    -outputformat STRUCTURE \
    -spid /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/vcf2structurepheno.spid > structure-input-pheno.log
    #all snps using phenotypes2.txt
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/structure-input-pheno2.txt \
    -outputformat STRUCTURE \
    -spid /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/vcf2structurepheno2.spid > structure-input-pheno2.log
#### VCF to FSTAT
#####create the spid file
    #replace the STRUCTURE section of vcf2structure.spid with the following for FSTAT
    #minus the leading spaces
    nano /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/vcf2structure.spid
        # FSTAT Writer questions
        WRITER_FORMAT=FSTAT
        # Specify which data type should be included in the FSTAT file  (FSTAT can only analyze one data type per file):
        FSTAT_WRITER_DATA_TYPE_QUESTION=SNP
        # Save label file
        FSTAT_WRITER_LABEL_FILE_QUESTION=
        # Do you want to save an additional file with labels (population names)?
        FSTAT_WRITER_INCLUDE_LABEL_QUESTION=false
        # Specify the locus/locus combination you want to write to the FSTAT file:
        FSTAT_WRITER_LOCUS_COMBINATION_QUESTION=
    #saved as vcf2fstat.spid
#####Do VCF to FSTAT file conversion
    #all snps
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/fstat-input-allsnps.txt \
    -outputformat FSTAT \
    -spid /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/vcf2fstat.spid > fstat-input-allsnps.log
    #only nonsyn snps
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/fstat-input-nonsyn.txt \
    -outputformat FSTAT \
    -spid /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/vcf2fstat.spid > fstat-input-nonsyn.log
### 4.Look at population structure with structure
#### Installed structure from source
#####Get source
    cd ~/bin/
    mkdir structure
    cd structure
    wget http://pritchardlab.stanford.edu/structure_software/release_versions/v2.3.4/structure_kernel_source.tar.gz
    tar xzf structure_kernel_source.tar.gz
    cd structure_kernel_src/
#####compile
    make
#####used the following setting for mainparams file created using Windows front-end version
        #define OUTFILE /work/jelber2/immunome_urtd/combined/popgen2/structure-results-001-k1
        #define INFILE /work/jelber2/immunome_urtd/combined/popgen2/structure-input-allsnps.txt
        #define NUMINDS 16
        #define NUMLOCI 16821
        #define LABEL 1
        #define POPDATA 1
        #define POPFLAG 0
        #define LOCDATA 0
        #define PHENOTYPE 0
        #define MARKERNAMES 1
        #define MAPDISTANCES 0
        #define ONEROWPERIND 0
        #define PHASEINFO 0
        #define PHASED 0
        #define RECESSIVEALLELES 0
        #define EXTRACOLS 0
        #define MISSING 
        #define PLOIDY 2
        #define MAXPOPS 1
        #define BURNIN 100000
        #define NUMREPS 1000000
        #define NOADMIX 0
        #define LINKAGE 0
        #define USEPOPINFO 0
        #define LOCPRIOR 0
        #define INFERALPHA 1
        #define ALPHA 1.0
        #define POPALPHAS 0
        #define UNIFPRIORALPHA 1
        #define ALPHAMAX 10.0
        #define ALPHAPROPSD 0.025
        #define FREQSCORR 1
        #define ONEFST 0
        #define FPRIORMEAN 0.01
        #define FPRIORSD 0.05
        #define INFERLAMBDA 0
        #define LAMBDA 1.0
        #define COMPUTEPROB 1
        #define PFROMPOPFLAGONLY 0
        #define ANCESTDIST 0 
        #define STARTATPOPINFO 0
        #define METROFREQ 10
        #define UPDATEFREQ 1
#### Ran structure using the following command for k1,k2,k3,k4 for 20 reps each
    # note had to make 80 mainparams.test.0* files
    # ex: mainparams.test.001.k1, mainparams.test.002.k1, etc.
    # run structure for populations.txt
    cd /work/jelber2/immunome_urtd/combined/popgen2/
    ~/bin/structure/structure_kernel_src/structure \
    -m mainparams.test.001.k1 \
    -e ~/bin/structure/structure_kernel_src/extraparams
    # etc.
    # note we used default settings for extraparams
    # (i.e., the correlated allele frequency and the admixture ancestry models)
    # implemented on SuperMike II using /home/jelber2/scripts/immunome_urtd/16-structure.py
#### Ran structure using phenotypes.txt for k1,k2,k3,k4 for 20 reps each
    # copy param files for phenotypes
    cd /work/jelber2/immunome_urtd/combined/popgen2/
    mkdir phenotypes
    cp mainparams.test.0* phenotypes/.
    cd phenotypes/
    # edit the INFILE
    perl -pi -e "s/structure-input-allsnps.txt/structure-input-pheno.txt/g" mainparams.test.0*
    # edit the OUTFILE
    perl -pi -e "s/(structure-results-0\w+.k\w)/\phenotypes\/\1/g" mainparams.test.0*
    # run structure
    cd /work/jelber2/immunome_urtd/combined/popgen/phenotypes/
    ~/bin/structure/structure_kernel_src/structure \
    -m mainparams.test.001.k1 \
    -e ~/bin/structure/structure_kernel_src/extraparams
    # etc.
    # note we used default settings for extraparams
    # (i.e., the correlated allele frequency and the admixture ancestry models)
    # implemented on SuperMike II using /home/jelber2/scripts/immunome_urtd/16-structure-pheno.py
#### Ran structure using phenotypes2.txt (i.e., for 2 phenotypes) for k1,k2,k3,k4 for 20 reps each
    # copy param files for phenotypes
    cd /work/jelber2/immunome_urtd/combined/popgen2/
    mkdir phenotypes2
    cp mainparams.test.0* phenotypes2/.
    cd phenotypes2/
    # edit the INFILE
    perl -pi -e "s/structure-input-allsnps.txt/structure-input-pheno2.txt/g" mainparams.test.0*
    # edit the OUTFILE
    perl -pi -e "s/(structure-results-0\w+.k\w)/\phenotypes2\/\1/g" mainparams.test.0*
    #sync files
    rsync --stats --progress --archive /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen2/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_urtd/combined/popgen2/ -n
    #sync script
    rsync --stats --progress --archive /home/jelber2/scripts/immunome_urtd/16-structure-pheno2.py \
    jelber2@mike.hpc.lsu.edu:/home/jelber2/scripts/immunome_urtd/ -n
    # run structure on SuperMikeII
    cd /work/jelber2/immunome_urtd/combined/popgen2/phenotypes2/
    ~/bin/structure/structure_kernel_src/structure \
    -m mainparams.test.001.k1 \
    -e ~/bin/structure/structure_kernel_src/extraparams
    # etc.
    # note we used default settings for extraparams
    # (i.e., the correlated allele frequency and the admixture ancestry models)
    # implemented on SuperMike II using /home/jelber2/scripts/immunome_urtd/16-structure-pheno2.py
    #sync files back to my computer
    rsync --stats --progress --archive jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_urtd/combined/popgen/phenotypes2/ \
    /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/phenotypes2/ -n
    # zip structure-results on my computer
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/phenotypes2/
    zip structure-results-* structure-results-phenotype2.zip
#### Used STRUCTURE HARVESTER web v0.6.94 to select best K values
#### Used CLUMPAK web to visualize population assignments
=====
# FUNCTIONAL ENRICHMENT ANALYSIS of gene deviating from neutrality
## 1.Get BLAST
    cd ~/bin/
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.31+-x64-linux.tar.gz
    tar xzf ncbi-blast-2.2.31+-x64-linux.tar.gz 
    mv ncbi-blast-2.2.31+-x64-linux.tar.gz ncbi-blast-2.2.31+
## 2.Get proteins in FASTA
    cd /media/immunome_2014/work/jelber2/blast/
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Chrysemys_picta/protein/protein.fa.gz
    gunzip protein.fa.gz
    #change FASTA header from gi|#######|ref|XP_#######.#| to XP_######.#
    perl -pi -e "s/gi\|\w+\|ref\|(XP_\d+.\d)\| /\1 /g" protein.fa
## 3.Make BLAST database
    ~/bin/ncbi-blast-2.2.31+/bin/makeblastdb \
    -in protein.fa \
    -input_type fasta \
    -dbtype prot \
    -parse_seqids \
    -out C_picta_protein.db \
    -title C_picta_protein
## 4.Perform BLAST to make input for BLAST2GO in xml format (on SuperMikeII)
    ~/bin/ncbi-blast-2.2.31+/bin/blastp \
    -db C_picta_protein.db \
    -query protein.fa \
    -num_threads 4 \
    -out prot.blast2go.input.xml \
    -outfmt 5
## 5.Use BLAST2GO version 3.1
    # import BLAST XML reults using default settings
    File->Load->Load Blast Results->Xml Files
    /media/immunome_2014/work/jelber2/blast/prot.blast2go.input.xml
    # wait for xml file to load in table
    # then import FASTA sequences using default settings
    File->Load->Load Sequences (e.g.:. fasta)
    Select "Add to the existing project"
    Select "Protein Sequences"
    /media/immunome_2014/work/jelber2/blast/protein.fa
## 6.BLAST2GO InterproScan
    # use default settings
    InterProScan->InterProScan
    InterProScan->Merge InterProScan GOs to Annotation
## 7.Enrichment Analysis
    Analysis->Enrichment Analysis (Fisher's Exact Test)
    # use default options
    # proteins with gene regions that deviate from neutrality
    TajimaD.extr.regions.bed.overlap.proteins.txt
=====
# GWAS with plink
# 1 install plink.107 and make dir for gwas
    #install plink
    cd ~/bin/
    wget http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-1.07-x86_64.zip
    unzip plink-1.07-x86_64.zip 
    mv plink-1.07-x86_64.zip plink-1.07-x86_64
    mv plink-1.07-x86_64/ plink-1.07
    #mkdir
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir gwas
    cd gwas
# 2 Use vcftools to create plink files
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/gwas/
    ~/bin/vcftools_0.1.12b/bin/vcftools \
    --vcf ../popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --plink
# 3 Rename plink files
    mv out.map ALL-samples-Q30-snps-recal-beagle-polymorphic.map
    mv out.ped ALL-samples-Q30-snps-recal-beagle-polymorphic.ped
    mv out.log ALL-samples-Q30-snps-recal-beagle-polymorphic.log
# 4 Coded phenotypes into .ped file
    nano ALL-samples-Q30-snps-recal-beagle-polymorphic.ped
    #edited 6th column (1=unaffected/asymptomatic, 2=affected/symptomatic)
# 5 Run plink association all samples
    ~/bin/plink-1.07/plink \
    --file ALL-samples-Q30-snps-recal-beagle-polymorphic \
    --assoc \
    --out ALL-samples-Q30-snps-recal-beagle-polymorphic \
    --allow-no-sex --adjust
# 6 Look for FDR_BH less than 0.05
    less -S ALL-samples-Q30-snps-recal-beagle-polymorphic.assoc.adjusted
    # no p values < 0.05
# 7 Repeat steps 5and6 with removing CF72 (4th line/sample)
    sed '4d' ALL-samples-Q30-snps-recal-beagle-polymorphic.ped \
    > ALL-samples-Q30-snps-recal-beagle-polymorphic-noCF72.ped
    #make a map file copy but name it with -noCF72 filename
    cp ALL-samples-Q30-snps-recal-beagle-polymorphic.map \
    ALL-samples-Q30-snps-recal-beagle-polymorphic-noCF72.map
    #run plink without CF72
    ~/bin/plink-1.07/plink \
    --file ALL-samples-Q30-snps-recal-beagle-polymorphic-noCF72 \
    --assoc \
    --out ALL-samples-Q30-snps-recal-beagle-polymorphic-noCF72 \
    --allow-no-sex --adjust
# 8 Look for FDR_BH less than 0.05
    less -S ALL-samples-Q30-snps-recal-beagle-polymorphic-noCF72.assoc.adjusted
    # no p values < 0.05
# 9 Get non-syn snps and indels
    #get header
    zcat ../split-vcfs/ALL-samples-snps-annotated.vcf.gz | grep "#" > header
    #get non-syn snps
    zcat ../split-vcfs/ALL-samples-snps-annotated.vcf.gz | \
    grep -v '#' | \
    grep -P 'SNPEFF_AMINO_ACID_CHANGE=\w*[A-Z]\d+\w*[A-Z]' \
    > ALL-samples-snps-annotated-nonsyn.txt
    #combine the header and file
    cat header ALL-samples-snps-annotated-nonsyn.txt \
    > ALL-samples-snps-annotated-nonsyn.vcf
    #get non-syn indels
    zcat ../split-vcfs/ALL-samples-indels-annotated.vcf.gz | \
    grep -v '#' | \
    grep -P 'SNPEFF_AMINO_ACID_CHANGE=\w*[A-Z]\d+\w*[A-Z]' \
    > ALL-samples-indels-annotated-nonsyn.txt
    #combine the header and file
    cat header ALL-samples-indels-annotated-nonsyn.txt \
    > ALL-samples-indels-annotated-nonsyn.vcf
# 10 Repeat steps #2-6 on non-syn snps
    #make plink input files
    ~/bin/vcftools_0.1.12b/bin/vcftools \
    --vcf ALL-samples-snps-annotated-nonsyn.vcf --plink \
    --out ALL-samples-snps-annotated-nonsyn
    #edit ped file to add phenotypes
    nano ALL-samples-snps-annotated-nonsyn.ped
    #edited 6th column (1=unaffected/asymptomatic, 2=affected/symptomatic)
    #run plink association on all samples
    ~/bin/plink-1.07/plink \
    --file ALL-samples-snps-annotated-nonsyn \
    --assoc \
    --out ALL-samples-snps-annotated-nonsyn \
    --allow-no-sex --adjust
    #look for FDR_BH < 0.05
    less -S ALL-samples-snps-annotated-nonsyn.assoc.adjusted
    #no significant SNPs after correcting for multiple tests
# 11 Repeat steps #2-6 on indels
    #make plink input files
    ~/bin/vcftools_0.1.12b/bin/vcftools \
    --gzvcf ../beagle/ALL-samples-Q30-indels-recal-beagle-polymorphic.vcf.gz \
    --plink \
    --out ALL-samples-Q30-indels-recal-beagle-polymorphic
    #edit ped file to put phenotypes
    cut -f 1-6 ALL-samples-Q30-snps-recal-beagle-polymorphic.ped \
    > first6rows
    cut -f 7- ALL-samples-Q30-indels-recal-beagle-polymorphic.ped \
    > ALL-samples-Q30-indels-recal-beagle-polymorphic2.ped
    paste first6rows ALL-samples-Q30-indels-recal-beagle-polymorphic2.ped \
    > ALL-samples-Q30-indels-recal-beagle-polymorphic.ped
    #run plink association on all samples
    ~/bin/plink-1.07/plink \
    --file ALL-samples-Q30-indels-recal-beagle-polymorphic \
    --assoc \
    --out ALL-samples-Q30-indels-recal-beagle-polymorphic \
    --allow-no-sex --adjust
    #look for FDR_BH < 0.05
    less -S ALL-samples-Q30-indels-recal-beagle-polymorphic.assoc.adjusted
    #no significant SNPs after correcting for multiple tests
# 12 Repeat steps #2-6 on non-imputated snps
    #make plink input files
    ~/bin/vcftools_0.1.12b/bin/vcftools \
    --vcf ../vqsr/ALL-samples-Q30-snps-recal.vcf \
    --plink \
    --out ALL-samples-Q30-snps-recal
    #edit ped file to put phenotypes
    cut -f 7- ALL-samples-Q30-snps-recal.ped \
    > ALL-samples-Q30-snps-recal2.ped
    paste first6rows ALL-samples-Q30-snps-recal2.ped \
    > ALL-samples-Q30-snps-recal.ped
    #run plink association on all samples
    ~/bin/plink-1.07/plink \
    --file ALL-samples-Q30-snps-recal \
    --assoc \
    --out ALL-samples-Q30-snps-recal \
    --allow-no-sex --adjust
    #look for FDR_BH < 0.05
    less -S ALL-samples-Q30-snps-recal.assoc.adjusted
    #no significant SNPs after correcting for multiple tests
# 13 Repeat steps #2-6 on non-imputated indels
    #make plink input files
    ~/bin/vcftools_0.1.12b/bin/vcftools \
    --vcf ../vqsr/ALL-samples-Q30-indels-recal.vcf \
    --plink \
    --out ALL-samples-Q30-indels-recal
    #edit ped file to put phenotypes
    cut -f 7- ALL-samples-Q30-indels-recal.ped \
    > ALL-samples-Q30-indels-recal2.ped
    paste first6rows ALL-samples-Q30-indels-recal2.ped \
    > ALL-samples-Q30-indels-recal.ped
    #run plink association on all samples
    ~/bin/plink-1.07/plink \
    --file ALL-samples-Q30-indels-recal \
    --assoc \
    --out ALL-samples-Q30-indels-recal \
    --allow-no-sex --adjust
    #look for FDR_BH < 0.05
    less -S ALL-samples-Q30-indels-recal.assoc.adjusted
    #no significant indels after correcting for multiple tests
=====
# GWAS with plink except use SNP-sets
# 1 Make snp.list and indel.list files
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/gwas/
    #snps
    while read i;do
        echo $i | \
        perl -pe "s/0\s(\w+_\d+\.\d):(\d+)\s0\s\d+/\1\t\2/" | \
        awk -v OFS='\t' '{a=$2-1;print $1,a,$2;}' - | \
        ~/bin/bedtools-2.22.1/bin/bedtools intersect \
        -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff \
        -b stdin | \
        grep -P "Gnomon\tgene" | \
        sed -n '1p' | \
        cut -f 1,4,5,9 | \
        perl -pe "s/(\w+_\d+\.\d)\t(\d+)\t(\d+)\t.+Name=(\w+);.+/\1\t\2\t\3\t\4/" >> snp.list
    done < ALL-samples-Q30-snps-recal-beagle-polymorphic.map
    #indels
    while read i;do
        echo $i | \
        perl -pe "s/0\s(\w+_\d+\.\d):(\d+)\s0\s\d+/\1\t\2/" | \
        awk -v OFS='\t' '{a=$2-1;print $1,a,$2;}' - | \
        ~/bin/bedtools-2.22.1/bin/bedtools intersect \
        -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff \
        -b stdin | \
        grep -P "Gnomon\tgene" | \
        sed -n '1p' | \
        cut -f 1,4,5,9 | \
        perl -pe "s/(\w+_\d+\.\d)\t(\d+)\t(\d+)\t.+Name=(\w+);.+/\1\t\2\t\3\t\4/" >> indel.list
    done < ALL-samples-Q30-indels-recal-beagle-polymorphic.map
# 2 Make setlists with plink
    #snps
    ~/bin/plink-1.07/plink \
    --file ALL-samples-Q30-snps-recal-beagle-polymorphic \
    --make-set snp.list \
    --write-set \
    --out ALL-samples-Q30-snps-recal-beagle-polymorphic
    #indels
    ~/bin/plink-1.07/plink \
    --file ALL-samples-Q30-indels-recal-beagle-polymorphic \
    --make-set indel.list \
    --write-set \
    --out ALL-samples-Q30-indels-recal-beagle-polymorphic
# 3 Do association tests with SNP- and Indel-sets
    #snps
    ~/bin/plink-1.07/plink \
    --file ALL-samples-Q30-snps-recal-beagle-polymorphic \
    --set-test \
    --set-max 1000 \
    --set ALL-samples-Q30-snps-recal-beagle-polymorphic.set \
    --mperm 10000 \
    --assoc \
    --allow-no-sex \
    --out ALL-samples-Q30-snps-recal-beagle-polymorphic
    #indels
    ~/bin/plink-1.07/plink \
    --file ALL-samples-Q30-indels-recal-beagle-polymorphic \
    --set-test \
    --set-max 1000 \
    --set ALL-samples-Q30-indels-recal-beagle-polymorphic.set \
    --mperm 10000 \
    --assoc \
    --allow-no-sex \
    --out ALL-samples-Q30-indels-recal-beagle-polymorphic


# GWAS with ROADTRIPS
## Controls for unknown population structure and unknown relatedness
# 1 install ROADTRIPS make dir roadtrips
    #install
    cd ~/bin/
    wget http://faculty.washington.edu/tathornt/software/ROADTRIPS2/ROADTRIPS2.0.tar.gz
    tar xzf ROADTRIPS2.0.tar.gz 
    mv ROADTRIPS2.0.tar.gz ROADTRIPS/
    #roadtrips
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir roadtrips
# 2 install FORMAT_PED_PHENO
    cd ~/bin/
    wget http://www.stat.uchicago.edu/~mcpeek/software/FORMAT_PED_PHENO/FORMAT_PED_PHENO1.0.tar.gz
    tar xzf FORMAT_PED_PHENO1.0.tar.gz 
    mv FORMAT_PED_PHENO1.0.tar.gz FORMAT
# 3 install KinINbcoef
    cd ~/bin/
    wget http://www.stat.uchicago.edu/~mcpeek/software/KinInbcoef/v1.1/KinInbcoef.tar.gz
    tar xzf KinInbcoef.tar.gz
    mv KinInbcoef.tar.gz KinInbcoef
# 4 SNPs convert PLINK ped file to tped format and tfam format
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/roadtrips/
    # make transposed ped file
    ~/bin/plink-1.07/plink \
    --file ../gwas/ALL-samples-Q30-snps-recal-beagle-polymorphic \
    --recode12 \
    --output-missing-genotype 0 \
    --transpose \
    --out ALL-samples-Q30-snps-recal-beagle-polymorphic
# 5 Snps make .pedpheno, kinpedigree, and .kinlist files
    ~/bin/FORMAT/FORMAT \
    -f ALL-samples-Q30-snps-recal-beagle-polymorphic.tfam \
    -o ALL-samples-Q30-snps-recal-beagle-polymorphic
# 6 Snps make kinfile to run ROADTRIPS
    ~/bin/KinInbcoef/KinInbcoef \
    ALL-samples-Q30-snps-recal-beagle-polymorphic.kinpedigree \
    ALL-samples-Q30-snps-recal-beagle-polymorphic.kinlist \
    ALL-samples-Q30-snps-recal-beagle-polymorphic.kinfile
# 7 Snps make prevalence file - accurate prevalence increase statistical power
    touch ALL-samples-Q30-snps-recal-beagle-polymorphic.prevalence
    nano ALL-samples-Q30-snps-recal-beagle-polymorphic.prevalence
    #type 0enter0 then save file
# 8 Snps run ROADTRIPS
    ~/bin/ROADTRIPS/ROADTRIPS \
    -g ALL-samples-Q30-snps-recal-beagle-polymorphic.tped \
    -p ALL-samples-Q30-snps-recal-beagle-polymorphic.pedpheno \
    -k ALL-samples-Q30-snps-recal-beagle-polymorphic.kinfile \
    -r ALL-samples-Q30-snps-recal-beagle-polymorphic.prevalence
    #rename output files
    mv ROADTRIPS_Software.err ROADTRIPS_Software.snps.err
    mv ROADTRIPStest.out ROADTRIPStest.snps.out
    mv ROADTRIPStest.pvalues ROADTRIPStest.snps.pvalues
    mv ROADTRIPStest.testvalues ROADTRIPStest.snps.testvalues
    mv ROADTRIPStest.top ROADTRIPStest.snps.top
# 10 Indels convert PLINK ped file to tped format and tfam format
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/roadtrips/
    # make transposed ped file
    ~/bin/plink-1.07/plink \
    --file ../gwas/ALL-samples-Q30-indels-recal-beagle-polymorphic \
    --recode12 \
    --output-missing-genotype 0 \
    --transpose \
    --out ALL-samples-Q30-indels-recal-beagle-polymorphic
# 11 Indels make .pedpheno, kinpedigree, and .kinlist files
    ~/bin/FORMAT/FORMAT \
    -f ALL-samples-Q30-indels-recal-beagle-polymorphic.tfam \
    -o ALL-samples-Q30-indels-recal-beagle-polymorphic
# 12 Indels make kinfile to run ROADTRIPS
    ~/bin/KinInbcoef/KinInbcoef \
    ALL-samples-Q30-indels-recal-beagle-polymorphic.kinpedigree \
    ALL-samples-Q30-indels-recal-beagle-polymorphic.kinlist \
    ALL-samples-Q30-indels-recal-beagle-polymorphic.kinfile
# 13 Indels make prevalence file - accurate prevalence increase statistical power
    touch ALL-samples-Q30-indels-recal-beagle-polymorphic.prevalence
    nano ALL-samples-Q30-indels-recal-beagle-polymorphic.prevalence
    #type 0enter0 then save file
# 14 Indels run ROADTRIPS
    ~/bin/ROADTRIPS/ROADTRIPS \
    -g ALL-samples-Q30-indels-recal-beagle-polymorphic.tped \
    -p ALL-samples-Q30-indels-recal-beagle-polymorphic.pedpheno \
    -k ALL-samples-Q30-indels-recal-beagle-polymorphic.kinfile \
    -r ALL-samples-Q30-indels-recal-beagle-polymorphic.prevalence
    #rename output files
    mv ROADTRIPS_Software.err ROADTRIPS_Software.indels.err
    mv ROADTRIPStest.out ROADTRIPStest.indels.out
    mv ROADTRIPStest.pvalues ROADTRIPStest.indels.pvalues
    mv ROADTRIPStest.testvalues ROADTRIPStest.indels.testvalues
    mv ROADTRIPStest.top ROADTRIPStest.indels.top
### GOT THE SAME Top SNPs and Indels as PLINK analysis
    #see ROADTRIPStest.indels.top and ROADTRIPStest.snps.top
