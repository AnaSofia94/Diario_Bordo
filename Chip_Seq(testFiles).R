#################################################################################
#Genome- wide data####

#cell-line: Jurkat ==> single end 
#CHIPseq_Pol1_: (GSM722164)==> (SRR443847)
#CHIPseq_Pol2_: (GSM722165)==> (SRR443848)
#CHIPseq_Input: (GSM722173)==> (SRR443856)

#-------------SBATCH file liftover-STEP1-------------------------
#!/bin/bash
#SBATCH --job-name=IniAnali
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/


srun shifter fastq-dump SRR443847.sra
##srun shifter fastq-dump SRR443848.sra
##srun shifter fastq-dump SRR443856.sra

##srun shifter fastqc SRR443847.fastq
##srun shifter fastqc SRR443848.fastq
##srun shifter fastqc SRR443856.fastq

##utilizado para retirar os ultimos nucleotidos 
##srun shifter trim_galore SRR3114088.fastq

##echo "Statistics for job $SLURM_JOB_ID:"
##sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#----------------------SBATCH file liftover- STEP2------------------
#!/bin/bash
#SBATCH --job-name=bowtie2
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:ummidock/bowtie2:latest
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

srun shifter bowtie2 -p 10 -x $GENOME -U SRR443847.fastq -S SRR443847.sam
srun shifter bowtie2 -p 10 -x $GENOME -U SRR443848.fastq -S SRR443848.sam
srun shifter bowtie2 -p 10 -x $GENOME -U SRR443856.fastq -S SRR443856.sam


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


#----------------------SBATCH Unique & peak- STEP3------------------
#!/bin/bash
#SBATCH --job-name=uni_peak
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/peakcalling
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

#converter de sam para bam 
srun shifter samtools view -bt SRR443847.sam > SRR443847.bam
srun shifter samtools view -bt SRR443848.sam > SRR443848.bam
srun shifter samtools view -bt SRR443856.sam > SRR443856.bam


##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

srun shifter macs2 callpeak -q 0.01 -t SRR443847.bam -c SRR443856.bam -f BAM -g hs -n SRR443847_peakCal.bed 
srun shifter macs2 callpeak -q 0.01 -t SRR443848.bam -c SRR443856.bam -f BAM -g hs -n SRR443848_peakCal.bed 

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

###################################################################################
#################################14/09/2017#######################################
#Genome Wide- Reconstitution-####

#cell-line: Jurkat ==> single end 
#CHIPseq_Pol1_: (GSM722164)==> (SRR443847)
#CHIPseq_Pol2_: (GSM722165)==> (SRR443848)
#CHIPseq_Input: (GSM722173)==> (SRR443856)

#-------------SBATCH file liftover-STEP1-------------------------
#!/bin/bash
#SBATCH --job-name=RepIA
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/


srun shifter fastq-dump SRR443847.sra
srun shifter fastq-dump SRR443848.sra
srun shifter fastq-dump SRR443856.sra

srun shifter fastqc SRR443847.fastq
srun shifter fastqc SRR443848.fastq
srun shifter fastqc SRR443856.fastq

##utilizado para retirar os ultimos nucleotidos 
##srun shifter trim_galore SRR3114088.fastq

##echo "Statistics for job $SLURM_JOB_ID:"
##sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#-------------SBATCH file liftover-MdStep2_Fastx-------------------------
#!/bin/bash
#SBATCH --job-name=RepFastx
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/


srun shifter fastq_quality_filter -q 20 -i SRR443847.fastq -o SRR443847CUT.fastq
srun shifter fastq_quality_filter -q 20 -i SRR443848.fastq -o SRR443848CUT.fastq
srun shifter fastq_quality_filter -q 20 -i SRR443856.fastq -o SRR443856CUT.fastq


##utilizado para retirar os ultimos nucleotidos 
##srun shifter trim_galore SRR3114088.fastq

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


#----------------------SBATCH Bowtie2--STEP2------------------
#!/bin/bash
#SBATCH --job-name=bowtie2Rep
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:ummidock/bowtie2:latest
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

srun shifter bowtie2 -p 10 -x $GENOME -U SRR443847CUT.fastq -S SRR443847.sam
srun shifter bowtie2 -p 10 -x $GENOME -U SRR443848CUT.fastq -S SRR443848.sam
srun shifter bowtie2 -p 10 -x $GENOME -U SRR443856CUT.fastq -S SRR443856.sam


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#----------------------See_Uniqe file liftover- STEP3------------------
#nao foi usado, tem este codigo so para identificar
#!/bin/bash
#SBATCH --job-name=SeeUni_Te
#SBATCH --time=2:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/peakcalling
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome


srun shifter grep -E "@|NM:" SRR443847.sam > SRR443847_unique.sam
srun shifter grep -E "@|NM:" SRR443848.sam > SRR443848_unique.sam
srun shifter grep -E "@|NM:" SRR443856.sam > SRR443856_unique.sam

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#----------------------SBATCH Unique & peak- STEP4------------------
#!/bin/bash
#SBATCH --job-name=DelUniq
#SBATCH --time=2:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/peakcalling
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

#ja nao precisa de srun shifter porque o grep ja é uma linha de comandos
grep -v "XS:i" SRR443847.sam > uniq_firstSampl.sam
grep -v "XS:i" SRR443848.sam > uniq_secSampl.sam
grep -v "XS:i" SRR443856.sam > uniq_input.sam

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


#----------------------SBATCH Unique & peak- STEP5------------------
#!/bin/bash
#SBATCH --job-name=uni_peak
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/peakcalling
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

#converter de sam para bam 
srun shifter samtools view -bt uniq_firstSampl.sam > uniq_firstSampl.bam
srun shifter samtools view -bt uniq_secSampl.sam > uniq_secSampl.bam
srun shifter samtools view -bt uniq_input.sam > uniq_input.bam
#eleminar os sam 

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

srun shifter macs2 callpeak -q 0.01 -t uniq_firstSampl.bam -c uniq_input.bam -f BAM -g hs -n Sample1_peakCal.bed 
srun shifter macs2 callpeak -q 0.01 -t uniq_secSampl.bam -c uniq_input.bam -f BAM -g hs -n Sample2_peakCal.bed 

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

###################################################################################
#################################15/09/2017#######################################
#Genome Wide- Reconstitution-####

#cell-line: Jurkat ==> single end 
#CHIPseq_Pol1_: (GSM722164)==> (SRR443847)
#CHIPseq_Pol2_: (GSM722165)==> (SRR443848)
#CHIPseq_Input: (GSM722173)==> (SRR443856)

#-------------SBATCH file liftover-STEP1-------------------------
#!/bin/bash
#SBATCH --job-name=RepIA
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/


srun shifter fastq-dump SRR443847.sra
srun shifter fastq-dump SRR443848.sra
srun shifter fastq-dump SRR443856.sra

srun shifter fastqc SRR443847.fastq
srun shifter fastqc SRR443848.fastq
srun shifter fastqc SRR443856.fastq

##utilizado para retirar os ultimos nucleotidos 
##srun shifter trim_galore SRR3114088.fastq

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#-------------SBATCH file liftover-MdStep2_Fastx-------------------------
#!/bin/bash
#SBATCH --job-name=RepFastx
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##second time=RepFastx2
srun shifter fastq_quality_filter -q 20 -Q 33 -i SRR443847.fastq -o SRR443847CUT.fastq
srun shifter fastq_quality_filter -q 20 -Q 33 -i SRR443848.fastq -o SRR443848CUT.fastq
srun shifter fastq_quality_filter -q 20 -Q 33 -i SRR443856.fastq -o SRR443856CUT.fastq

##srun shifter fastx_trimmer -l 20 -i SRR443847.fastq -o SRR443847CUT.fastq
##srun shifter fastx_trimmer -l 20 -i SRR443848.fastq -o SRR443848CUT.fastq
##srun shifter fastx_trimmer -l 20 -i SRR4438456.fastq -o SRR443856CUT.fastq


srun shifter fastqc SRR443847CUT.fastq
srun shifter fastqc SRR443848CUT.fastq
srun shifter fastqc SRR443856CUT.fastq

##srun shifter fastx_quality_stats -i SRR443847CUT.fastq -o SRR443847CUT_stats.txt
##srun shifter fastx_quality_stats -i SRR443848CUT.fastq -o SRR443848CUT_stats.txt
##srun shifter fastx_quality_stats -i SRR443856CUT.fastq -o SRR443856CUT_stats.txt

##utilizado para retirar os ultimos nucleotidos 
##srun shifter trim_galore SRR3114088.fastq

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#----------------------SBATCH Bowtie2--STEP2------------------
#!/bin/bash
#SBATCH --job-name=bowtie2Rep
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:ummidock/bowtie2:latest
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

srun shifter bowtie2 -p 10 -x $GENOME -U SRR443847CUT.fastq -S SRR443847.sam
srun shifter bowtie2 -p 10 -x $GENOME -U SRR443848CUT.fastq -S SRR443848.sam
srun shifter bowtie2 -p 10 -x $GENOME -U SRR443856CUT.fastq -S SRR443856.sam


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#----------------------SBATCH Unique & peak- STEP4------------------
#!/bin/bash
#SBATCH --job-name=DelUniq
#SBATCH --time=2:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/peakcalling
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

#ja nao precisa de srun shifter porque o grep ja é uma linha de comandos
grep -v "XS:i" SRR443847.sam > uniq_firstSampl.sam
grep -v "XS:i" SRR443848.sam > uniq_secSampl.sam
grep -v "XS:i" SRR443856.sam > uniq_input.sam

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#----------------------SBATCH Unique & peak- STEP4------------------
#!/bin/bash
#SBATCH --job-name=uniq_Flag
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/peakcalling
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/
##job ID-6196

#get the uniques and convert SAM-->BAM 
srun shifter samtools view -c -F 4 uniq_firstSampl.sam > uniq_firstSampl.bam
srun shifter samtools view -c -F 4 uniq_secSampl.sam > uniq_secSampl.bam
srun shifter samtools view -c -F 4 uniq_input.sam > uniq_input.bam

#eleminar os sam 

##final unique regarding the initial aligned reads;
##Does a full pass through the input file to calculate and print statistics to stdout.
srun shifter samtools flagstat uniq_firstSampl.bam
srun shifter samtools flagstat uniq_secSampl.bam
srun shifter samtools flagstat uniq_input.bam


#----------------------SBATCH Unique & peak- STEP4------------------
#!/bin/bash
#SBATCH --job-name=uni_peak
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/peakcalling
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

#see the peaks
srun shifter macs2 callpeak -q 0.05 -t uniq_firstSampl.bam -c uniq_input.bam -f BAM -g hs -n Sample1_peakCal
srun shifter macs2 callpeak -q 0.05 -t uniq_secSampl.bam -c uniq_input.bam -f BAM -g hs -n Sample2_peakCal 

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

###################################################################################
#################################18/09/2017#######################################
#Genome Wide- Reconstitution-####

#cell-line: Jurkat ==> single end 
#CHIPseq_Pol1_: (GSM722164)==> (SRR443847)
#CHIPseq_Pol2_: (GSM722165)==> (SRR443848)
#CHIPseq_Input: (GSM722173)==> (SRR443856)

#-------------SBATCH file liftover-STEP1-------------------------
#!/bin/bash
#SBATCH --job-name=RepIA
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/


srun shifter fastq-dump SRR443847.sra
srun shifter fastq-dump SRR443848.sra
srun shifter fastq-dump SRR443856.sra

srun shifter fastqc SRR443847.fastq
srun shifter fastqc SRR443848.fastq
srun shifter fastqc SRR443856.fastq

#run fastx in order to use a quality filter
srun shifter fastq_quality_filter -p 80 -q 20 -Q33  -i SRR443847.fastq -o SRR443847CUT.fastq 
srun shifter fastq_quality_filter -p 80 -q 20 -Q33  -i SRR443848.fastq -o SRR443848CUT.fastq 
srun shifter fastq_quality_filter -p 80 -q 20 -Q33  -i SRR443856.fastq -o SRR443856CUT.fastq 

#volar a correr o fastq
srun shifter fastqc SRR443847CUT.fastq
srun shifter fastqc SRR443848CUT.fastq
srun shifter fastqc SRR443856CUT.fastq 


##utilizado para retirar os ultimos nucleotidos 
##srun shifter trim_galore SRR3114088.fastq

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#-------------SBATCH file alignment-STEP2-------------------------
#!/bin/bash
#SBATCH --job-name=btie2Rep
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:ummidock/bowtie2:latest
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

#neste job como nao fiz o fastx, as reads nao estao mais curtas
srun shifter bowtie2 -p 10 -x $GENOME -U SRR443847CUT.fastq -S SRR443847CUT.sam
srun shifter bowtie2 -p 10 -x $GENOME -U SRR443848CUT.fastq -S SRR443848CUT.sam
srun shifter bowtie2 -p 10 -x $GENOME -U SRR443856CUT.fastq -S SRR443856CUT.sam

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#----------------------SBATCH Unique & peak- STEP4------------------
#!/bin/bash
#SBATCH --job-name=DelUniq
#SBATCH --time=2:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/peakcalling
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

#ja nao precisa de srun shifter porque o grep ja é uma linha de comandos
grep -v "XS:i" SRR443847CUT.sam > uniq_firstSampl.sam
grep -v "XS:i" SRR443848CUT.sam > uniq_secSampl.sam
grep -v "XS:i" SRR443856CUT.sam > uniq_input.sam

#get the uniques and convert SAM-->BAM 
srun shifter samtools view -b -F 4 uniq_firstSampl.sam > uniq_firstSampl.bam
srun shifter samtools view -b -F 4 uniq_secSampl.sam > uniq_secSampl.bam
srun shifter samtools view -b -F 4 uniq_input.sam > uniq_input.bam

#eleminar os sam 

##final unique regarding the initial aligned reads;
##Does a full pass through the input file to calculate and print statistics to stdout.
srun shifter samtools flagstat uniq_firstSampl.bam
srun shifter samtools flagstat uniq_secSampl.bam
srun shifter samtools flagstat uniq_input.bam

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#----------------------SBATCH Unique & peak- STEP4------------------
#!/bin/bash
#SBATCH --job-name=uni_peak
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/peakcalling
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

##Indice 
export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

#see the peaks

srun shifter macs2 callpeak -t uniq_firstSampl.bam -c uniq_input.bam -f BAM -g hs -n Sample1_peakCal -q 0.05 -s 40
srun shifter macs2 callpeak -t uniq_secSampl.bam   -c uniq_input.bam -f BAM -g hs -n Sample2_peakCal -q 0.05 -s 40

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#----------------------SBATCH Enriched regions- STEP5--------------####
#!/bin/bash
#SBATCH --job-name=Ftools
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/peakcalling
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/

export BLACKLIST=hg38_mergedBlacklists.bed

##Script to prepare bigWig for genomebrowser (you just need the bam from bowtie and the bed from MACS2)

##Filter only reads overlapping enriched regions
srun shifter intersectBed -wa -ubam -abam Sample1_peakCal.bam -b MACS2_output_Sample1_peakCal_highSignPeaks.bed > Sample1_peakCal_enrichedReads.bam
srun shifter intersectBed -wa -ubam -abam Sample2_peakCal.bam -b MACS2_output_Sample1_peakCal_highSignPeaks.bed > Sample2_peakCal_enrichedReads.bam

##Produce index for bam file
srun shifter samtools sort -o Sample1_peakCal_enrichedReads_sorted.bam Sample1_peakCal_enrichedReads.bam
srun shifter samtools index Sample1_peakCal_enrichedReads_sorted.bam

srun shifter samtools sort -o Sample2_peakCal_enrichedReads_sorted.bam Sample2_peakCal_enrichedReads.bam
srun shifter samtools index Sample2_peakCal_enrichedReads_sorted.bam

## Produce bigwig file for USCS Genome Browser
srun shifter bamCoverage -b Sample1_peakCal_enrichedReads_sorted.bam -o Sample1_peakCal_deepTools.bw -of bigwig --binSize 10 --numberOfProcessors max/2 --normalizeUsingRPKM --extendReads 200 --ignoreDuplicates --hg38_mergedBlacklists.bed $BLACKLIST
srun shifter bamCoverage -b Sample2_peakCal_enrichedReads_sorted.bam -o Sample2_peakCal_deepTools.bw -of bigwig --binSize 10 --numberOfProcessors max/2 --normalizeUsingRPKM --extendReads 200 --ignoreDuplicates --hg38_mergedBlacklists.bed $BLACKLIST


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

###################################################################################
#################################18/09/2017#######################################
#Genome Wide- deeptools-####

#!/bin/bash
#SBATCH --job-name=temp1
#SBATCH --time=3:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS
#SBATCH --exclusive

##----------------------------------------------------------- INPUTS --------------------------------------------- ###
export wk=$SLURM_SUBMIT_DIR ## working directory

## ----------------------------------------------------------- INPUTS --------------------------------------------- ###
##criar ficheiros para por nas pastas
mkdir sra_files fastq_files fastqc_files 

##export directory into different folders
export sra_dir=$wk/sra_files
export fastq_dir=$wk/fastq_files
export fastqc_dir=$wk/fastqc_files

##ficheiros a correr em paralelo;
srun="srun -N1 -n1"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog /home/anasofiamoreira/scratch/MasterAnaS/$SLURM_JOB_ID.log"

time ls *.sra | $parallel 'mv {} $sra_dir'

## Transform .sra to .fastq
time ls $sra_dir/*.sra | $parallel '$srun shifter --image=argrosso/htspreprocessing:0.1.1 fastq-dump {} --outdir $fastq_dir'

## Check the quality of the reads with fastqc
time ls $fastq_dir/*.fastq | $parallel '$srun shifter --image=argrosso/htspreprocessing:0.1.1 fastqc {} -o $fastqc_dir'

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


###################################################################################
#################################18/09/2017#######################################

#----------------------SBATCH Enriched regions- STEP5--------------####
#!/bin/bash
#SBATCH --job-name=cut_xls
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/htstools:latest
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/


##export BLACKLIST=/mnt/nfs/lobo/SALMEIDA-NFS/GeneArchive/salmeida_data/Annotations/mm10/mm10_mergedBlacklists.bed ##MOUSE:GRCm38/mm10
##export BLACKLIST=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome/hg38_mergedBlacklists.bed

export BLACKLIST=hg38_mergedBlacklists.bed


##cut the files in macs in order to get a bed completed
srun shifter sed '1,28d' Sample1_peakCal_peaks.xls > Sample1_peakCal_cut.xls
cut -f 1,2,3 Sample1_peakCal_cut.xls > Sample1_peakCal_cut.bed
srun shifter sed '1,28d' Sample2_peakCal_peaks.xls > Sample2_peakCal_cut.xls
cut -f 1,2,3 Sample2_peakCal_cut.xls > Sample2_peakCal_cut.bed


##Script to prepare bigWig for genomebrowser (you just need the bam from bowtie and the bed from MACS2)

##Filter only reads overlapping enriched regions
srun shifter intersectBed -wa -ubam -abam uniq_firstSampl.bam -b Sample1_peakCal_cut.bed > Sample1_peakCal_enrichedReads.bam
srun shifter intersectBed -wa -ubam -abam uniq_secSampl.bam -b Sample2_peakCal_cut.bed > Sample2_peakCal_enrichedReads.bam

##Produce index for bam file
srun shifter samtools sort -o Sample1_peakCal_enrichedReads_sorted.bam Sample1_peakCal_enrichedReads.bam
srun shifter samtools index Sample1_peakCal_enrichedReads_sorted.bam

srun shifter samtools sort -o Sample2_peakCal_enrichedReads_sorted.bam Sample2_peakCal_enrichedReads.bam
srun shifter samtools index Sample2_peakCal_enrichedReads_sorted.bam

## Produce bigwig file for USCS Genome Browser #extensao para ver o bed.graph; ir ao manual e depois ver 
srun shifter bamCoverage -b Sample1_peakCal_enrichedReads_sorted.bam -o Sample1_peakCal_deepTools.bw -of bigwig --binSize 10 --numberOfProcessors max/2 --normalizeUsingRPKM --extendReads 200 --ignoreDuplicates --blackListFileName $BLACKLIST
srun shifter bamCoverage -b Sample2_peakCal_enrichedReads_sorted.bam -o Sample2_peakCal_deepTools.bw -of bigwig --binSize 10 --numberOfProcessors max/2 --normalizeUsingRPKM --extendReads 200 --ignoreDuplicates --blackListFileName $BLACKLIST


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

###################################################################################
#################################2/10/2017#######################################
#----------------------SBATCH bigwig--bedgraph- STEP1--------------####

#!/bin/bash
#SBATCH --job-name=snob2
#SBATCH --time=3:00:00
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/
#SBATCH --image=docker:argrosso/htstools:latest
#SBATCH --exclusive

#ntasks=1 porque so tinhamos 2 comandos, logo nao eram precisos tantos 
export BLACKLIST=hg38_mergedBlacklists.bed

## Produce bigwig file for USCS Genome Browser #extensao para ver o bed.graph; ir ao manual e depois ver 
srun shifter bamCoverage -b Sample1_peakCal_enrichedReads_sorted.bam -o Sample1.bedgraph -of bedgraph --binSize 10 --numberOfProcessors max/2 --normalizeUsingRPKM --extendReads 200 --ignoreDuplicates --blackListFileName $BLACKLIST
srun shifter bamCoverage -b Sample2_peakCal_enrichedReads_sorted.bam -o Sample2.bedgraph -of bedgraph --binSize 10 --numberOfProcessors max/2 --normalizeUsingRPKM --extendReads 200 --ignoreDuplicates --blackListFileName $BLACKLIST


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

#----------------------SBATCH intersect---------------####

##intersect the peaks with the genes in order to find the most suitable places of connections, beteween 
## the genes and the peaks.
##For that use:INTERSECT BED TOOLS

##1- you have to find the genes. To that you will look for the hg38 annotation in the xls file 
awk '{if ($3 == "gene") print $0}' /mnt/nfs/lobo/IMM-NFS/genomes/hg38/Annotation/gencode.v26.annotation.gtf > gencode.v26.annotation_gene.gtf``

#--------------------Get-Promoter-----------------####
#!/bin/bash
#SBATCH --job-name=promote
#SBATCH --time=3:00:00
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/
#SBATCH --image=argrosso/peakcalling

export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/WholeGenomeFasta/genomeSize.txt


#getting the genes and the promoter
srun shifter flankBed -s -i gencode.v26.annotation_gene.gtf -g $GENOME -l 2000 -r 0 > upstream.gtf
srun shifter flankBed -s -i upstream.gtf -g $GENOME -l 0 -r 2000 > promoter.gtf

##converting bedgraph into bed
cat Sample1.bedgraph | awk '{print $1 "\t" $2 "\t" $3}' > Sample1.bed
cat Sample2.bedgraph | awk '{print $1 "\t" $2 "\t" $3}' > Sample2.bed

##intersect the genes to my promoter and samples #subst a por:Sample1_peakCal_cut.bed 
srun shifter intersectBed -a Sample1.bed -b promoter.gtf > inter_sample1.bed
srun shifter intersectBed -a Sample2.bed -b promoter.gtf > inter_sample2.bed


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

##intersect das samples
grep -v chrEBV Sample1.bedgraph > Tall1_Rep1.bedgraph 
gzip Tall1_Rep1.bedgraph ## ficheiro zipado
gunzip XXXX ##ficheiro nao zipado 

##Para melhorar a visibilidade dos bedgraph
track type=bedGraph name=Tall1_Rep1  visibility=full color=51,51,255

###################################################################################
#################################23/10/2017#######################################
#----------------------SBATCH bigwig--bedgraph- --------------####

##intersect the samples of the promoter in order to get the most enriched regions 

#!/bin/bash
#SBATCH --job-name=PROM_INTER2
#SBATCH --time=3:00:00
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/
#SBATCH --image=argrosso/peakcalling

srun shifter intersectBed -a inter_sample1.bed -b inter_sample2.bed> comon2_intersect.bed

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

###################################################################################
#################################02/11/2017#######################################
#----------------------SBATCH bigwig--bedgraph- --------------####

##intersect the samples of the promoter in order to get the most enriched regions 

#!/bin/bash
#SBATCH --job-name=PROM_Again
#SBATCH --time=3:00:00
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/home/anasofiamoreira/scratch/MasterAnaS/
#SBATCH --image=argrosso/peakcalling


##intersect the genes to my promoter and samples #subst a por:Sample1_peakCal_cut.bed 
srun shifter intersectBed -a Sample1_peakCal_cut.bed  -b promoter.gtf > inter_sample1_NEW.bed
srun shifter intersectBed -a Sample2_peakCal_cut.bed  -b promoter.gtf > inter_sample2_NEW.bed
srun shifter intersectBed -a inter_sample1_NEW.bed -b inter_sample2_NEW.bed> comon2_intersect_NEW.bed

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

