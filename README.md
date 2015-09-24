#Immunome_URTD
=====
##Get Fastq files from BaseSpace
###Rename fastq files
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
    perl -pe s"/(.+\/BaseCalls\/)(\w+)(_\w+_\w+_)(R\d)(_\d+.fastq.gz)/mv \1\2\3\4\5 \/media\/immunome_2014\/work\/jelber2\/immunome_urtd\/run1/fastq\/\2-\4.fastq.gz/" rename-fastq-run1 > rename-fastq-run1.sh
    # add #!/bin/bash to first line
    sed -i '1 i\#!/bin/bash' rename-fastq-run1.sh
    bash rename-fastq-run1.sh
    # run2
    mv Bioo\ NEXTflex-23135384.zip Bioo-NEXTflex-23135384.zip
    ## copy the file to external drive
    cp ~/Desktop/Bioo-NEXTflex-23135384.zip /media/immunome_2014/work/jelber2/immunome_urtd/run2/.
    ## go to the directory
    cd /media/immunome_2014/work/jelber2/immunome_urtd/run2/
    ## unzip the archive
    unzip Bioo-NEXTflex-23135384.zip
    ## rename the unzipped folder
    mv Bioo\ NEXTflex-23135384 Bioo-NEXTflex-23135384
    ## grab all of the file names and save them as a document
    find /media/immunome_2014/work/jelber2/immunome_urtd/run2/Bioo-NEXTflex-23135384/ \
    -maxdepth 5 -type f -print > rename-fastq-run2
    ## make a directory called fastq
    mkdir fastq
    ## use regular expressions to create rename-fastq-run2.sh
    perl -pe s"/(.+\/BaseCalls\/)(\w+)(_\w+_\w+_)(R\d)(_\d+.fastq.gz)/mv \1\2\3\4\5 \/media\/immunome_2014\/work\/jelber2\/immunome_urtd\/run1/fastq\/\2-\4.fastq.gz/" rename-fastq-run2 > rename-fastq-run2.sh
    # add #!/bin/bash to first line
    sed -i '1 i\#!/bin/bash' rename-fastq-run2.sh
    bash rename-fastq-run2.sh
##Put Data on SuperMikeII
    # run1
    rsync --stats --progress --archive \
    /media/immunome_2014/work/jelber2/immunome_urtd/run1/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_urtd/run1/ -n
    # run2
    rsync --stats --progress --archive \
    /media/immunome_2014/work/jelber2/immunome_urtd/run2/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_urtd/run2/ -n
##Install programs and get reference genome
###trimmomatic-0.32
    cd /home/jelber2/bin/
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
    unzip Trimmomatic-0.32.zip
    mv Trimmomatic-0.32.zip Trimmomatic-0.32
    #PATH=~/home/jelber2/bin/Trimmomatic-0.32/trimmomatic-0.32.jar
###bbmerge-5.4 (part of bbmap-34.33)
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
###bwa-0.7.12
    cd /home/jelber2/bin/
    wget https://github.com/lh3/bwa/archive/0.7.12.tar.gz
    mv 0.7.12 bwa-0.7.12.tar.gz
    tar xzf bwa-0.7.12.tar.gz
    mv bwa-0.7.12.tar.gz bwa-0.7.12
    cd bwa-0.7.12/
    make
    #PATH=~/bin/bwa-0.7.12/bwa
###stampy-1.0.23
    cd /home/jelber2/bin/
    wget http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz
    tar xzf Stampy-latest.tgz
    cd stampy-1.0.23
    make
    #PATH=~/bin/stampy-1.0.23/stampy.py
###java jre1.7.0
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
###picard-1.128
    #on my Centos machine
    cd /home/jelber2/bin/
    wget https://github.com/broadinstitute/picard/releases/download/1.128/picard-tools-1.128.zip
    rsync --stats --archive --progress /home/jelber2/bin/picard-tools-1.128.zip jelber2@mike.hpc.lsu.edu:/home/jelber2/bin/ -n
    #switched to SuperMikeII
    cd /home/jelber2/bin/
    unzip picard-tools-1.128.zip
    mv picard-tools-1.128.zip picard-tools-1.128
    #PATH=~/bin/picard-tools-1.128/picard.jar
###GATK-3.3.0
    #had to download using firefox on my Centos machine
    #saved in /home/jelber2/bin/GATK-3.3.0
    rsync --stats --archive --progress /home/jelber2/bin/GATK-3.3.0/ jelber2@mike.hpc.lsu.edu:/home/jelber2/bin/GATK-3.3.0/ -n
    #switched to SuperMikeII
    cd /home/jelber2/bin/
    cd GATK-3.3.0
    tar xjf GenomeAnalysisTK-3.3-0.tar.bz2
    #PATH=~/bin/GATK-3.3.0/GenomeAnalysisTK.jar
###samtools-1.1
    cd /home/jelber2/bin/
    wget http://downloads.sourceforge.net/project/samtools/samtools/1.1/samtools-1.1.tar.bz2?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2Ffiles%2Fsamtools%2F1.1%2F&ts=1421967581&use_mirror=softlayer-dal
    tar xjf samtools-1.1.tar.bz2 
    mv samtools-1.1.tar.bz2 samtools-1.1
    cd samtools-1.1
    make
    nano ~/.soft #add the following line to .soft file using nano
    PATH += /home/jelber2/bin/samtools-1.1/
###parallel-20150122
    cd /home/jelber2/bin/
    wget ftp://ftp.gnu.org/gnu/parallel/parallel-20150122.tar.bz2
    tar xjf parallel-20150122.tar.bz2
    mv parallel-20150122.tar.bz2 parallel-20150122
    #PATH=~/bin/parallel-20150122/src/parallel
####Get bedtools2.22.1
    cd ~/bin/
    wget https://github.com/arq5x/bedtools2/releases/download/v2.22.1/bedtools-2.22.1.tar.gz
    tar xzf bedtools-2.22.1.tar.gz
    mv bedtools2 bedtools-2.22.1
    mv bedtools-2.22.1.tar.gz bedtools-2.22.1
    cd bedtools-2.22.1
    make
###Got painted turtle reference genome (on SuperMikeII)
    cd /work/jelber2/reference/
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
    gunzip GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
##STEPS FOR QUALITY CONTROL, MAPPING, & SNP CALLING
###1.Make indexes (on SuperMikeII)
    qsub /home/jelber2/scripts/immunome_urtd/01-make_indexes.sh
###2.Quality control and adapter trimming (on SuperMikeII)
    # run1
    cd /work/jelber2/immunome_urtd/run1/fastq
    ~/scripts/immunome_urtd/02-trimmomatic.py *.fastq.gz
    # run2
    cd /work/jelber2/immunome_urtd/run2/fastq
    ~/scripts/immunome_urtd/02-trimmomatic-run2.py *.fastq.gz
###3.BWA alignment (on SuperMikeII)
    # run1
    cd /work/jelber2/immunome_urtd/run1/trimmed-data/
    ~/scripts/immunome_urtd/03-bwa.py *.trim.fastq.gz
    # run2
    cd /work/jelber2/immunome_urtd/run2/trimmed-data/
    ~/scripts/immunome_urtd/03-bwa-run2.py *.trim.fastq.gz
###4.STAMPY alignment (on SuperMikeII)
    # run1
    cd /work/jelber2/immunome_urtd/run1/bwa-alignment/
    ~/scripts/immunome_urtd/04-stampy.py *.bwa.sam
    # run2
    cd /work/jelber2/immunome_urtd/run2/bwa-alignment/
    ~/scripts/immunome_urtd/04-stampy-run2.py *.bwa.sam
###5a.Clean,Sort,Add Read Groups (on SuperMikeII)
    # run1
    cd /work/jelber2/immunome_urtd/run1/stampy-alignment/
    ~/scripts/immunome_urtd/05a-clean_sort_addRG.py *.stampy.bam
    # run2
    cd /work/jelber2/immunome_urtd/run2/stampy-alignment/
    ~/scripts/immunome_urtd/05a-clean_sort_addRG-run2.py *.stampy.bam
###5b.Clean,Sort,Add Read Groups, DeDup,Realign Around Indels(on SuperMikeII)
    cd /work/jelber2/immunome_urtd/run1/clean-sort-addRG/
    ~/scripts/immunome_urtd/05b-clean_sort_addRG_markdup_realign.py *-CL-RG.bam
###6.Merge BAM files, Call SNPs initially (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/realign-around-indels/
    ~/scripts/immunome_urtd/06-mergeBAM_callSNPs_initial.py *-realigned.bam
###7.Quality score recalibration 1 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/realign-around-indels/
    find . -name '*-realigned.bam' -not -name 'ALL-samples-*' \
     -exec ~/scripts/immunome_urtd/07-qual_score_recal01.py {} \;
###8.Merge BAM files, Call SNPs recalibrated 1 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal01/
    ~/scripts/immunome_urtd/08-mergeBAM_callSNPs_recal01.py *-recal01.bam
###9.Quality score recalibration 2 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal01/
    find . -name '*-recal01.bam' -not -name 'ALL-samples-*' \
    -exec ~/scripts/immunome_urtd/09-qual_score_recal02.py {} \;
###10.Merge BAM files, Call SNPs recalibrated 2 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal02/
    ~/scripts/immunome_urtd/10-mergeBAM_callSNPs_recal02.py *-recal02.bam
###11.Quality score recalibration 3 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal02/
    find . -name '*-recal02.bam' -not -name 'ALL-samples-*' \
    -exec ~/scripts/immunome_urtd/11-qual_score_recal03.py {} \;
###12.Merge BAM files, Call SNPs recalibrated 3 (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal03/
    ~/scripts/immunome_urtd/12-mergeBAM_callSNPs_recal03.py *-recal03.bam
###13.Sequencing metrics
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal03/
    ~/scripts/immunome_2014/13-seq_metrics.py ALL-samples-recal03.bam
###14.plot coverage for each sample
####run bedtools coverage on all bam files, then keep only lines with 'all' on them
    #FROM http://gettinggeneticsdone.blogspot.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
    #ALSO FROM https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md
#####make samplelist
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/call-SNPs-recal03/
    ls *.bam | grep -v "ALL"| perl -pe "s/-recal03.bam//g" > samplelist
#####calculate coverage
    cd /media/immunome_2014/work/jelber2/immunome_urtd/
    mkdir plot_coverage
    cd plot_coverage/
    #command below takes 5-10 minutes
    while read i;do
    ~/bin/bedtools-2.22.1/bin/bedtools coverage -abam ../call-SNPs-recal03/$i-recal03.bam \
    -b /media/immunome_2014/work/jelber2/reference/immunome_baits_C_picta-3.0.3.bed \
    -hist | grep ^all > $i.baitcoverage.all.txt
    done < ../call-SNPs-recal03/samplelist
    #now use modified R scripts from links above to plot coverage
###15.Need to use featureCounts to summarize number of genes, reads per gene, etc
####Get Subread
    #featureCounts is part of the Subread package http://bioinf.wehi.edu.au/featureCounts/
    cd ~/bin/
    wget http://downloads.sourceforge.net/project/subread/subread-1.4.6/subread-1.4.6-Linux-x86_64.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsubread%2Ffiles%2Fsubread-1.4.6%2F&ts=1423446986&use_mirror=iweb
    mv subread-1.4.6-Linux-x86_64.tar.gz?r=http:%2F%2Fsourceforge.net%2Fprojects%2Fsubread%2Ffiles%2Fsubread-1.4.6%2F subread-1.4.6-Linux-x86_64.tar.gz
    tar xzf subread-1.4.6-Linux-x86_64.tar.gz
    mv subread-1.4.6-Linux-x86_64.tar.gz subread-1.4.6-Linux-x86_64
####Get genometools-1.5.4 to annotate introns in gff file
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
####Need to prefilter the GFF file for immune genes
    cd /media/immunome_2014/work/jelber2/reference/
    #intersect the gff file
    ~/bin/bedtools-2.22.1/bin/bedtools intersect \
    -a GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.introns \
    -b immunome_baits_C_picta-3.0.3.bed \
    > GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gff
####Convert GFF file of gene annotations to GTF
    cd /media/immunome_2014/work/jelber2/reference/
    perl -pe "s/\S+=GeneID:(\d+).+/gene_id \"\1\";/g" \
    GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gff \
    > GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf
####run subread
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
####write an R function to do the following
#####count number of possible immune genes
    wc -l CF53.gene
    #total genes = 632 (after subtracting 2 header lines)
#####count number of possible immune gene exons
    wc -l AL102.exon
    #total exons = 37275 (after subtracting 2 header lines)
#####how many different immune genes were captured
    grep -Pv "\t0$" ALL.gene | wc -l
    #609 (after subtracting 2 header lines)
#####how many different immune gene exons were captured
    grep -Pv "\t0$" ALL.exon | wc -l
    #4676 (after subtracting 2 header lines)
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
###16.Haplotype Caller (on SuperMikeII)
    cd /work/jelber2/immunome_urtd/combined/call-SNPs-recal03/
    #excludes file ALL-samples-recal03.bam
    find . -name '*-recal03.bam' -not -name 'ALL-samples-*' \
    -exec ~/scripts/immunome_urtd/14-haplotypecaller.py {} \;
###17.Ran GenotypeGVCFs to perform joint genotyping
    #on Cenots machine
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/hc/
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
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
###18.Added expressions to filter variants
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
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
###19.Got only Indel variants
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -T SelectVariants \
    -V ALL-samples-Q30-snps-indels.vcf \
    -o ALL-samples-Q30-indels.vcf \
    -selectType INDEL
###20.Got only SNP variants
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -T SelectVariants \
    -V ALL-samples-Q30-snps-indels.vcf \
    -o ALL-samples-Q30-snps.vcf \
    -selectType SNP
###21.Ran Variant Recalibrator
####Note: used SNPS and Indels from previous immunome_2014 data for "truthing" and filtration
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir vqsr
    cd vqsr
#####Recalibrated snps
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
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
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
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
####Applied the recalibration on snps
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -input ../hc/ALL-samples-Q30-snps.vcf \
    --ts_filter_level 99.5 \
    -tranchesFile VQSR-snps.tranches \
    -recalFile VQSR-snps.recal \
    -o ALL-samples-Q30-snps-recal.vcf
####Applied the recalibration on indels
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -input ../hc/ALL-samples-Q30-indels.vcf \
    --ts_filter_level 99.0 \
    -tranchesFile VQSR-indels.tranches \
    -recalFile VQSR-indels.recal \
    -o ALL-samples-Q30-indels-recal.vcf
###22.Needed to use beagle to improve SNPs (using Linkage Disequilibrium) called by Haplotype Caller
####Downloaded beagle
    cd ~/bin
    wget http://faculty.washington.edu/browning/beagle/beagle.r1398.jar
####Ran beagle on snps and indels separately
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir beagle
    cd beagle
    #snps
    java -Xmx8000m -jar ~/bin/beagle.r1398.jar \
    gtgl=/media/immunome_2014/work/jelber2/immunome_urtd/vqsr/ALL-samples-Q30-snps-recal.vcf\
    nthreads=2 \
    out=/media/immunome_2014/work/jelber2/immunome_urtd/beagle/ALL-samples-Q30-snps-recal-beagle
    #indels
    java -Xmx8000m -jar ~/bin/beagle.r1398.jar \
    gtgl=/media/immunome_2014/work/jelber2/immunome_urtd/vqsr/ALL-samples-Q30-indels-recal.vcf\
    nthreads=2 \
    out=/media/immunome_2014/work/jelber2/immunome_urtd/beagle/ALL-samples-Q30-indels-recal-beagle
###23.Get only polymorphic loci
    # Because many SNP/indel loci will occur because "fixed" differeces
    # between Gopher tortoise and Western Painted turtle
####Downloaded vcftools
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
####Remove loci with AF=1
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
    #filter out nonpolymorphicsnps
    #might take a few minutes
    while read i;do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/remove\tlocus\t\n/" ALL-samples-Q30-snps-recal-beagle2.vcf
    done < nonpolymorphicsnps
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
##STEPS FOR VARIANT PREDICTION
###1.Download Tools First
####Downloaded snpEff version 4.0e
    #ideally want to know if variants will affect protein structure and possibly immune gene function
    cd /work/jelber2/reference
    wget http://iweb.dl.sourceforge.net/project/snpeff/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
####Added Chrysemys_picta_bellii-3.0.3 to snpEff.config using nano
    cd /media/immunome_2014/work/jelber2/reference/snpEff
    nano snpEff.config # added the following four lines after the Capsella_rubella_v1.0 entry (remove 4 spaces on left if cut and pasting)
    # Chrysemys_picta_bellii-3.0.3
    Chrysemys_picta_bellii-3.0.3.genome : western painted turtle
    	Chrysemys_picta_bellii-3.0.3.reference : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/
    	Chrysemys_picta_bellii-3.0.3.M.codonTable : Standard
####Created data directory for Chrysemys_picta_bellii-3.0.3 genome
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
####Built snpEff database for Chrysemys_picta_bellii-3.0.3
    cd /media/immunome_2014/work/jelber2/reference/snpEff/
    # used snpEff_build.py script to implement command below, which took < 30 minutes
    java -jar -Xmx8g /media/immunome_2014/work/jelber2/reference/snpEff/snpEff.jar build -gff3 -v Chrysemys_picta_bellii-3.0.3 2>&1 | tee Chrysemys_picta_bellii-3.0.3.build
####Downloaded bcftools
    cd ~/bin/
    git clone --branch=develop git://github.com/samtools/htslib.git
    git clone --branch=develop git://github.com/samtools/bcftools.git
    cd bcftools; make
###2.Need to look for protein altering variants shared by samples in the same phenotype
####a.Split vcf file for snpEff
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
####b.Ran snpEff on each split vcf file
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/split-vcfs/
    # command below to run snpEff on all samples in samplelist
    # not implemented on SuperMikeII b/c process was < 15 min
    #snps
    while read i;do
    java -Xmx8g -jar /media/immunome_2014/work/jelber2/reference/snpEff/snpEff.jar \
    -v -i vcf -o gatk \
    Chrysemys_picta_bellii-3.0.3 \
    $i-snps.vcf > $i-snps-snpeff.vcf
    mv snpEff_genes.txt $i-snps-snpeff-genes.txt
    mv snpEff_summary.html $i-snps-snpeff-summary.html
    done < ../call-SNPs-recal03/samplelist
    #indels
    while read i;do
    java -Xmx8g -jar /media/immunome_2014/work/jelber2/reference/snpEff/snpEff.jar \
    -v -i vcf -o gatk \
    Chrysemys_picta_bellii-3.0.3 \
    $i-indels.vcf > $i-indels-snpeff.vcf
    mv snpEff_genes.txt $i-indels-snpeff-genes.txt
    mv snpEff_summary.html $i-indels-snpeff-summary.html
    done < ../call-SNPs-recal03/samplelist
####c.Ran VariantAnnotator on each snpeff file
    #snps
    while read i;do
    rm $i-snps.vcf.idx
    rm $i-snps-snpeff.vcf.idx
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
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
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -A SnpEff \
    --variant $i-indels.vcf \
    --snpEffFile $i-indels-snpeff.vcf \
    -L $i-indels.vcf \
    -o $i-indels-annotated.vcf
    done < ../call-SNPs-recal03/samplelist
####d.Merge split, annotated vcfs
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
####e.Merge vcf files then index
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
####f.Get only high quality non-synonymous alleles
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
    #3764 (after subtracting 1 for header)
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
    #220 (after subtracting 1 for header)
    # phenotypes
    # Subclinical/Recovered (n=4): CF090, CF072, FC015, OLD065
    # Sick at last observation (n=6): CF080, CF219, FC013, FC047, OLD092, OLD107
    # Healthy (n=6): CF053, CF069, FC019, FC058, OLD077, OLD106
####For getting snp alleles shared amongst phenotypes
#####http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#TOC-Venn-Diagrams
=====
##STEPS FOR LOOKING FOR SNPs UNDER SELECTION
###1.Download tools
####a.Download BayeScan
    cd ~/bin/
    wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
    unzip BayeScan2.1.zip
    mv BayeScan2.1.zip BayeScan2.1
    cd BayeScan2.1/
    cd binaries/
    chmod u+x BayeScan2.1_linux64bits # makes the file executable
    # Path to BayeScan
    /home/jelber2/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits
####b.Download Simple Fool's Guide (SFG) to RNA-seq scripts to convert vcf file to BayeScan input format
    cd ~/scripts/immunome_2014/
    mkdir fromSFG
    cd fromSFG
    wget http://sfg.stanford.edu/Scripts.zip
    unzip Scripts.zip 
    mv Scripts\ for\ SFG/ Scripts_for_SFG
####c. make directory
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir bayescan
    cd bayescan
###2.Run for BayeScan for snps
####a.Add Genotype Qualities to ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    zcat ../beagle/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz | perl -pe "s/(GT:DS:GP)/\1:GQ/" \
    > ALL-samples-Q30-snps-recal-beagle-polymorphic-fixed.vcf
    perl -pe "s/(\d\|\d:\d:\d,\d,\d)/\1:30/g" \
    ALL-samples-Q30-snps-recal-beagle-polymorphic-fixed.vcf \
    > ALL-samples-Q30-snps-recal-beagle-polymorphic-fixed2.vcf
####b.Make populations.txt file
    text file with samplename\tpopulation (samplename tab population)
####c.Ran make_bayescan_input.py
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
    mv bayes_input.txt bayes_input.txt.snps
    mv low_freq_snps.txt low_freq_snps.txt.snps
    mv population-info.txt population-info.txt.snps
    mv snpkey.txt snpkey.txt.snps
    mv used_snp_genos.txt used_snp_genos.txt.snps
####d.Run BayeScan
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /media/immunome_2014/work/jelber2/immunome_urtd/bayescan/bayes_input.txt.snps \
    -snp \
    -d low_freq_snps.txt.snps \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles.snps \
    -threads 2
####e.View bayescan results
    #initiate R in the terminal
    R
    setwd("/media/immunome_2014/work/jelber2/immunome_urtd/bayescan/")
    #source the plot_R.r script from Bayescan
    source("/home/jelber2/bin/BayeScan2.1/R functions/plot_R.r")
    #plot fst values without minor alleles below minor allele frequency of 1 copy
    noMAF_snps_results <- plot_bayescan("bayescan_no_loci_with_low_freq_minor_alleles.snps_fst.txt", FDR=0.05)
    #save the candidate loci to a text file
    write(noMAF_snps_results$outliers, file= "noMAF_loci_FDR_0.05_outlier_snps.txt", ncolumns= 1,append= FALSE)
    q()
####f.View bayescan results in IGV
    #create a copy of snpkey.txt, so it can be modified
    cp snpkey.txt.snps snpkey.txt.snps2
    #code to create IGV batch file for noMAF loci
    while read i;do
    perl -pi -e "s/^$i\t(.+)_(.+)\n/goto \1:\2\n/" snpkey.txt.snps2
    done < noMAF_loci_FDR_0.05_outlier_snps.txt
    grep 'goto' snpkey.txt.snps2 > noMAF_loci_FDR_0.05_outlier_snps_igv.txt
    #view in IGV
    ~/bin/IGV_2.3.40/igv.sh /media/immunome_2014/work/jelber2/immunome_urtd/split-vcfs/ALL-samples-snps-annotated.vcf.gz
    #open noMAF_loci_FDR_0.05_outlier_snps_igv.txt
####g.Filter annotated VCF file by outlier snps
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    zcat ../split-vcfs/ALL-samples-snps-annotated.vcf.gz > ALL-samples-snps-annotated2.vcf
    perl -pe "s/goto (\w+\.\d):(\d+)\n/\1\t\2\n/" noMAF_loci_FDR_0.05_outlier_snps_igv.txt > noMAF_loci_FDR_0.05_outlier_snps_vcf.txt
    while read i;do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/\1\tOUTLIER_SNP\t\2\n/" ALL-samples-snps-annotated2.vcf
    done < noMAF_loci_FDR_0.05_outlier_snps_vcf.txt
    grep 'OUTLIER_SNP\|^#' ALL-samples-snps-annotated2.vcf | grep -v "contig" > ALL-samples-outlier-snps.vcf
    # get only gene names, note that some SNPs are intergenic
    grep -v "#" ALL-samples-outlier-snps.vcf | \
    perl -pe "s/.+SNPEFF_GENE_NAME=(\w+);.+\n/\1\n/" | \
    perl -pe "s/NW_.+\n/intergenic\n/g" | \
    sort | uniq -c | \
    perl -pe "s/( )+/\t/g" > ALL-samples-outlier-snps-gene-names.txt
    # how many SNPs are under selection?
    perl -ane '$sum += $F[0]; END {print $sum; print "\n"}' ALL-samples-outlier-snps-gene-names.txt
    # 2
    # how many genes have SNPs under selection
    grep -v "intergenic" ALL-samples-outlier-snps-gene-names.txt | wc -l
    # 2
    # CD74 and C-type lectin domain family 2 member D
    # how many SNPs are intergenic
    grep "intergenic" ALL-samples-outlier-snps-gene-names.txt | perl -ane '$sum += $F[0]; END {print $sum; print "\n"}'
    # 0
    # how many SNPs are contained in genes
    grep -v "intergenic" ALL-samples-outlier-snps-gene-names.txt | perl -ane '$sum += $F[0]; END {print $sum; print "\n"}'
    # 2
###3.Run BayeScan for indels
####a.Add Genotype Qualities to ALL-samples-Q30-indels-recal-beagle-polymorphic.vcf.gz
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    zcat ../beagle/ALL-samples-Q30-indels-recal-beagle-polymorphic.vcf.gz | perl -pe "s/(GT:DS:GP)/\1:GQ/" \
    > ALL-samples-Q30-indels-recal-beagle-polymorphic-fixed.vcf
    perl -pe "s/(\d\|\d:\d:\d,\d,\d)/\1:30/g" \
    ALL-samples-Q30-indels-recal-beagle-polymorphic-fixed.vcf \
    > ALL-samples-Q30-indels-recal-beagle-polymorphic-fixed2.vcf
####b.Ran make_bayescan_input.py
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
    mv bayes_input.txt bayes_input.txt.indels
    mv low_freq_snps.txt low_freq_indels.txt.indels
    mv population-info.txt population-info.txt.indels
    mv snpkey.txt snpkey.txt.indels
    mv used_snp_genos.txt used_snp_genos.txt.indels
####c.Run BayeScan
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /media/immunome_2014/work/jelber2/immunome_urtd/bayescan/bayes_input.txt.indels \
    -snp \
    -d low_freq_indels.txt.indels \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles.indels \
    -threads 2
####d.View bayescan results
    #initiate R in the terminal
    R
    setwd("/media/immunome_2014/work/jelber2/immunome_urtd/bayescan/")
    #source the plot_R.r script from Bayescan
    source("/home/jelber2/bin/BayeScan2.1/R functions/plot_R.r")
    #plot fst values without minor alleles below minor allele frequency of 1 copy
    noMAF_indels_results <- plot_bayescan("bayescan_no_loci_with_low_freq_minor_alleles.indels_fst.txt", FDR=0.05)
    #save the candidate loci to a text file
    write(noMAF_indels_results$outliers, file= "noMAF_loci_FDR_0.05_outlier_indels.txt", ncolumns= 1,append= FALSE)
    q()
####e.View bayescan results in IGV
    #create a copy of snpkey.txt, so it can be modified
    cp snpkey.txt.indels snpkey.txt.indels2
    #code to create IGV batch file for noMAF loci
    while read i;do
    perl -pi -e "s/^$i\t(.+)_(.+)\n/goto \1:\2\n/" snpkey.txt.indels2
    done < noMAF_loci_FDR_0.05_outlier_indels.txt
    grep 'goto' snpkey.txt.indels2 > noMAF_loci_FDR_0.05_outlier_indels_igv.txt
    #view in IGV
    ~/bin/IGV_2.3.40/igv.sh /media/immunome_2014/work/jelber2/immunome_urtd/split-vcfs/ALL-samples-indels-annotated.vcf.gz
    #open noMAF_loci_FDR_0.05_outlier_indels_igv.txt
####f.Filter annotated VCF file by outlier indels
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/bayescan/
    zcat ../split-vcfs/ALL-samples-indels-annotated.vcf.gz > ALL-samples-indels-annotated2.vcf
    perl -pe "s/goto (\w+\.\d):(\d+)\n/\1\t\2\n/" noMAF_loci_FDR_0.05_outlier_indels_igv.txt > noMAF_loci_FDR_0.05_outlier_indels_vcf.txt
    while read i;do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/\1\tOUTLIER_SNP\t\2\n/" ALL-samples-indels-annotated2.vcf
    done < noMAF_loci_FDR_0.05_outlier_indels_vcf.txt
    grep 'OUTLIER_SNP\|^#' ALL-samples-indels-annotated2.vcf | grep -v "contig" > ALL-samples-outlier-indels.vcf
    # get only gene names, note that some indels are intergenic
    grep -v "#" ALL-samples-outlier-indels.vcf | \
    perl -pe "s/.+SNPEFF_GENE_NAME=(\w+);.+\n/\1\n/" | \
    perl -pe "s/NW_.+\n/intergenic\n/g" | \
    sort | uniq -c | \
    perl -pe "s/( )+/\t/g" > ALL-samples-outlier-indels-gene-names.txt
    # how many indels are under selection?
    perl -ane '$sum += $F[0]; END {print $sum; print "\n"}' ALL-samples-outlier-indels-gene-names.txt
    # 0
    # how many genes have indels under selection
    grep -v "intergenic" ALL-samples-outlier-indels-gene-names.txt | wc -l
    # 0
    # how many indels are intergenic
    grep "intergenic" ALL-samples-outlier-indels-gene-names.txt | perl -ane '$sum += $F[0]; END {print $sum; print "\n"}'
    # 0
    # how many indels are contained in genes
    grep -v "intergenic" ALL-samples-outlier-indels-gene-names.txt | perl -ane '$sum += $F[0]; END {print $sum; print "\n"}'
    # 0
=====
##STEPS FOR POPGEN
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/
    mkdir popgen
###1.Make copy of ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/beagle/
    cp ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz ../popgen/.
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/
    gunzip ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz
###2.Get only non-synonymous SNPs
    zcat ../split-vcfs/ALL-samples-snps-annotated.vcf.gz| grep '#' > header
    zcat ../split-vcfs/ALL-samples-snps-annotated.vcf.gz| grep -P 'SNPEFF_AMINO_ACID_CHANGE=\w*[A-Z]\d+\w*[A-Z]' > non-synonymous
    cat header non-synonymous > ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf
###3.Check loci for linkage disequilibrium and Hardy-Weinberg Equilibrium
####Had to make populations.txt file and population-specific files for vcf filtering
    cd /media/immunome_2014/work/jelber2/immunome_urtd/combined/popgen/
    grep "CF" ../bayescan/populations.txt > CF
    grep "OLD" ../bayescan/populations.txt > OLD
    grep "FC" ../bayescan/populations.txt > FC
    cp ../bayescan/populations.txt .
####Hardy-Weinberg Equilibrium test
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
####outputs linkage disequilibrium pvalues
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
###3.Convert VCF file to structure, fstat
####VCF to Structure
#####use PGDSpider
    #get the program
    cd ~/bin/
    wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.8.3.zip
    unzip PGDSpider_2.0.8.3.zip
    mv PGDSpider_2.0.8.3.zip PGDSpider_2.0.8.3
    #use following command to generate spider.conf.xml and spid file in ../popgen directory
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/popgen/structure-input.txt \
    -outputformat STRUCTURE
    #edit spider.conf.xml
    nano ~/bin/PGDSpider_2.0.8.3/spider.conf.xml #to add path to samtools
    #change <entry key="PathBcftools"></entry>
    #to <entry key="PathBcftools">/home/jelber2/bin/samtools-0.1.19/bcftools/bcftools</entry>
    #change <entry key="PathSamtools"></entry>
    #to <entry key="PathSamtools">/home/jelber2/bin/samtools-0.1.19/samtools</entry>
    #save and exit
    #edit the spid file
    nano /media/immunome_2014/work/jelber2/immunome_urtd/popgen/template_VCF_STRUCTURE.spid
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
        VCF_PARSER_POP_FILE_QUESTION=/media/immunome_2014/work/jelber2/immunome_urtd/popgen/populations.txt
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
#####Do VCF to STRUCTURE file conversion
    #all snps
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/popgen/structure-input-allsnps.txt \
    -outputformat STRUCTURE \
    -spid /media/immunome_2014/work/jelber2/immunome_urtd/popgen/vcf2structure.spid > structure-input-allsnps.log
    #only nonsyn snps
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/popgen/structure-input-nonsyn.txt \
    -outputformat STRUCTURE \
    -spid /media/immunome_2014/work/jelber2/immunome_urtd/popgen/vcf2structure.spid > structure-input-nonsyn.log
####VCF to FSTAT
#####create the spid file
    #replace the STRUCTURE section of vcf2structure.spid with the following for FSTAT
    #minus the leading spaces
    nano /media/immunome_2014/work/jelber2/immunome_urtd/popgen/vcf2structure.spid
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
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/popgen/fstat-input-allsnps.txt \
    -outputformat FSTAT \
    -spid /media/immunome_2014/work/jelber2/immunome_urtd/popgen/vcf2fstat.spid > fstat-input-allsnps.log
    #only nonsyn snps
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.8.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_urtd/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic-nonsyn.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_urtd/popgen/fstat-input-nonsyn.txt \
    -outputformat FSTAT \
    -spid /media/immunome_2014/work/jelber2/immunome_urtd/popgen/vcf2fstat.spid > fstat-input-nonsyn.log
###4.Look at population structure with structure
####Installed structure from source
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
        #define OUTFILE /work/jelber2/immunome_2014/combined/popgen/structure-results
        #define INFILE /work/jelber2/immunome_2014/combined/popgen/structure-input.txt
        #define NUMINDS 16
        #define NUMLOCI 15891
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
        #define BURNIN 50000
        #define NUMREPS 100000
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
####Ran structure using the following command for k1,k2,k3,k4
    cd /work/jelber2/immunome_2014/combined/popgen/
    ~/bin/structure/structure_kernel_src/structure \
    -m mainparams.test.k1 \
    -e ~/bin/structure/structure_kernel_src/extraparams
    #etc.
    #note we used default settings for extraparams
    #(i.e., the correlated allele frequency and the admixture ancestry models)
    #implemented on SuperMike II using /home/jelber2/scripts/immunome_2014/16-structure.py
####Used STRUCTURE HARVESTER weeb v0.6.94 to select best K values
####Used CLUMPP v1.1.2b to average data from multiple runs
    #get CLUMPP
    wget http://web.stanford.edu/group/rosenberglab/software/CLUMPP_Linux64.1.1.2.tar.gz
    tar xzf CLUMPP_Linux64.1.1.2.tar.gz 
    mv CLUMPP_Linux64.1.1.2.tar.gz CLUMPP_Linux64.1.1.2/.
    #used CLUMPP
    #
    #
####Used DISTRUCT v1.1 to visualize population assignments.
    #get DISTRUCT
    cd ~/bin/
    wget http://web.stanford.edu/group/rosenberglab/software/distruct1.1.tar.gz
    tar xzf distruct1.1.tar.gz 
    mv distruct1.1.tar.gz distruct1.1/.
    #used DISTRUCT
    #
    #

