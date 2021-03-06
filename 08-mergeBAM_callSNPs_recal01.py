#! /usr/bin/env python

# PBS cluster job submission in Python
# Use picard to merge recal01 BAM files then makes index
# Then call SNPs with GATK-3.3.0
# By Jean P. Elbers
# jean.elbers@gmail.com
# Last modified 24 Sep 2015
###############################################################################
Usage = """

08-mergeBAM_callSNPs_recal01.py - version 1.0

Command:
cd InDir = /work/jelber2/immunome_urtd/combined/call-SNPs-recal01
1.Merge Bam files
    java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar MergeSamFiles \
    SO=coordinate \
    AS=true \
    CREATE_INDEX=true \
    I=Sample1-recal01.bam \
    I=Sample2-recal01.bam \
    O=ALL-samples-recal01.bam

2.Call SNPs
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I ALL-samples-recal01.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -maxAltAlleles 32 \
    -gt_mode DISCOVERY \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -o ALL-samples-recal01-Q30-rawSNPS.vcf

3.Call Indels
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I ALL-samples-recal01.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -maxAltAlleles 32 \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -o ALL-samples-recal01-Q30-indels.vcf
    
4.Filter SNP calls around indels
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -V ALL-samples-recal01-Q30-rawSNPS.vcf \
    --mask ALL-samples-recal01-Q30-indels.vcf \
    --maskExtension 5 \
    --maskName InDel \
    --clusterWindowSize 10 \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
    --filterName "Bad Validation" \
    --filterExpression "QUAL < 30.0" \
    --filterName "LowQual" \
    --filterExpression "QD < 2.0" \
    --filterName "Low Variant Confidence" \
    --genotypeFilterExpression "DP < 10.0" \
    --genotypeFilterName "Low Read Depth Over Sample" \
    --genotypeFilterExpression "GQ < 20.0" \
    --genotypeFilterName "Low GenotypeQuality" \
    -o ALL-samples-recal01-Q30-SNPs.vcf


File Info:
InDir = /work/jelber2/immunome_urtd/combined/call-SNPs-recal01
Input Files = *-recal01.bam
OutDir = InDir
Output Files = ALL-samples-recal01.bam
               ALL-samples-recal01-Q30-SNPs.vcf


Usage (execute following code in InDir):

~/scripts/immunome_urtd/08-mergeBAM_callSNPs_recal01.py *-recal01.bam

"""
###############################################################################
import os, sys, subprocess #imports os, sys, subprocess modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    IFileList =[]
    for File in FileList:
        IFile = "I="+File+" \\"
        IFileList.append(IFile)
    IFileListString = '\n'.join(IFileList)
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome_urtd/combined/call-SNPs-recal01"
    os.chdir(InDir)
    # Customize your options here
    Queue = "single"
    Allocation = "hpc_gopo04"
    Processors = "nodes=1:ppn=4"
    WallTime = "04:00:00"
    LogOut = "/work/jelber2/immunome_urtd/combined/call-SNPs-recal01"
    LogMerge = "oe"
    JobName = "mergeBAM_callSNPs_recal01"
    Command ="""
    java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar MergeSamFiles \
    SO=coordinate \
    AS=true \
    CREATE_INDEX=true \
    %s
    O=ALL-samples-recal01.bam

    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I ALL-samples-recal01.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -maxAltAlleles 32 \
    -gt_mode DISCOVERY \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -o ALL-samples-recal01-Q30-rawSNPS.vcf

    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I ALL-samples-recal01.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -maxAltAlleles 32 \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -o ALL-samples-recal01-Q30-indels.vcf

    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -V ALL-samples-recal01-Q30-rawSNPS.vcf \
    --mask ALL-samples-recal01-Q30-indels.vcf \
    --maskExtension 5 \
    --maskName InDel \
    --clusterWindowSize 10 \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
    --filterName "Bad Validation" \
    --filterExpression "QUAL < 30.0" \
    --filterName "LowQual" \
    --filterExpression "QD < 2.0" \
    --filterName "Low Variant Confidence" \
    --genotypeFilterExpression "DP < 10.0" \
    --genotypeFilterName "Low Read Depth Over Sample" \
    --genotypeFilterExpression "GQ < 20.0" \
    --genotypeFilterName "Low GenotypeQuality" \
    -o ALL-samples-recal01-Q30-SNPs.vcf""" % \
    (IFileListString,
    RefDir,
    RefDir,
    RefDir)

    JobString = """
    #!/bin/bash
    #PBS -q %s
    #PBS -A %s
    #PBS -l %s
    #PBS -l walltime=%s
    #PBS -o %s
    #PBS -j %s
    #PBS -N %s

    cd %s
    %s\n""" % (Queue, Allocation, Processors, WallTime, LogOut, LogMerge, JobName, InDir, Command)

    #Create pipe to qsub
    proc = subprocess.Popen(['qsub'], shell=True,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
    (child_stdout, child_stdin) = (proc.stdout, proc.stdin)

    #Print JobString
    JobName = proc.communicate(JobString)[0]
    print JobString
    print JobName
