#!/bin/bash
#THIS PROGRAM IS PROVIED BY ZDQ 
# YOU CAN CONTACT US WITH 'zhangdongqin2@126.com'
################################################
################################################
################################################
################################################
##  define  Initial variable                  ##
################################################
################################################
################################################
################################################
raw_reads_dir=""
genome_fasta=""
index=""
status=off                
output=""
thread=""
mode =""
#indel1=""
#indel2=""
#dbsnp=""


################################################
################################################
################################################
################################################
##  define  Usage function                    ##
################################################
################################################
################################################
################################################
Usage () {                 
    echo "Usage:\n"
    echo -e "\tThis program uses fastq files as input files to directly generate various intermediate file results in various WES/WGS analysis. The software includes fastp, bwa, gatk4, bcftools, samtools, annovar. The software can be started from any place where the program is running. The various stages include: denovo, bwa mapping, gatk markduplicates, gatk baserecalibrator, gatk haplotypecaller, gatk genotypegvcfs, vcfsplit. If the software is interrupted, you can restart it from where it left off. "
    echo -e "This program is provided by Zdq\n "
    echo -e "raw reads must conform to the format of SAMPLE_1.fastq.gz , SAMPLE_2.fastq.gz,and the sample name'SAMPLE' cannot have'_' or'.' connection symbol in the middle \n "
    echo -e "easy_rnaseq.py [-h] [-v] [-d </path/to/raw/reads>] [-o </path/to/output>] [-r </path/to/genome fasta>] [-i </path/to/genome index>] [-t <thread number>][-k <known indel annotation>] [-g <1000 indel annotation>] [-s <dbsnp anntation>]"
    echo -e "To run the program, you need to specify the absolute path of fastp, bwa, samtools, gatk4 software in the configure file in the current directory\n"
    echo -e "Options:"
    echo -e "\t-h help information for porgram"
    echo -e "\t-d directory of the wgs raw reads.eg. ~/data/wes/raw_reads"
    echo -e "\t-o result output directory of program. eg. ~/data/wes/results"
    echo -e "\t-r genome fasta file.eg. genome.fa"
    echo -e "\t-i bwa index name. eg. ~/db/hg38/bow_index/genome"
    echo -e "\t-t thread numbers of program. eg. -t 16"
    echo -e "\t-x you can specify the mutation comment file directory of GATK4 when running VQSR"
    echo -e "\t-a you can specify the database directory when running annovar"
    echo -e "\t-m You can specify the node where the software starts running:\n\t\t-m all: denovo analysis\n\t\t-m mapping: start from bwa mapping\n\t\t-m markduplicates: start from gatk  markduplicates"
    echo -e "\t\t-m baserecalibrator: start from gatk  baserecalibrator"
    echo -e "\t\t-m haplotypecaller: start from gatk haplotypecaller"
    echo -e "\t\t-m genotypegvcfs: start from gatk genotypegvcfs"
    echo -e "\t\t-m vcfsplit: start from gatk vcfsplit"
    exit -1
}

#example   ./easy_reseq.sh -t 32 -d /home/data/wgs_out/other_wgs/raw_reads \
#               -o /home/data/wgs_out/other_wgs \
#               -i /home/origene/yinshan.cui/database/homo_sapiens/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome \
#               -r /home/origene/yinshan.cui/database/homo_sapiens/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
#               -x /home/data/genome_resequencing_DB \
#               -a /home/origene/yinshan.cui/bin/annovar-master/humandb \
#               -m all
##echo -e "\t-g dbsnp annotation vcf file. eg.db138.vcf"
##echo -e "\t-k indel annotation vcf files. eg. indel.vcf"
##echo -e "\t-g indel annotation vcf files. eg. 1000gold.indel.vcf"

while getopts :hvd:o:r:i:t:a:x:m: varname   

do
   case $varname in
    h)
      #echo "$varname"
      Usage
      exit
      ;;
    v)
      #echo "$varname"
      status=on
      #echo "$status"
      exit
      ;;
    d)
      #echo "$varname"
      #echo "$OPTARG"
      raw_reads_dir=$OPTARG                   
      if [ ! -d $raw_reads_dir ];then             
         echo "the source directory $raw_reads_dir not exist!"
         exit
      fi
      ;;
    o)
      #echo "$varname"
      #echo "$OPTARG"
      output=$OPTARG                      
      if [ ! -d  $output ];then               
         echo "the output path $output not exist"
         exit
      fi
      ;;
    r)
      #echo "$varname"
      #echo "$OPTARG"
      genome_fasta=$OPTARG                      
      if [ ! -f  $genome_fasta ];then              
         echo "the genome fasta file $genome_fasta not exist"
         exit
      fi
      ;;
    #g)
      #echo "$varname"
      #echo "$OPTARG"
    #  indel2=$OPTARG                     
    #  if [ ! -f  $indel2 ];then             
    #     echo "the indel file $indel2 does not exist"
    #     exit
    #  fi
    #  ;;
    #k)
      #echo "$varname"
      #echo "$OPTARG"
    #  indel1=$OPTARG                     
    #  if [ ! -f  $indel1 ];then             
    #     echo "the indel file $indel1 does not exist"
    #     exit
    #  fi
    #  ;;
    #s)
      #echo "$varname"
      #echo "$OPTARG"
    #  dbsnp=$OPTARG                     
    #  if [ ! -f  $dbsnp ];then             
    #     echo "the dbsnp file $dbsnp does not exist"
    #     exit
    #  fi
    #  ;;
    a)
      #echo "$varname"
      #echo "$OPTARG"
      humandb=$OPTARG                     
      if [ ! -d  $humandb ];then             
         echo "the humandb directory $humandb does not exist"
         exit
      fi
      ;;
    i)
      #echo "$varname"
      #echo "$OPTARG"
      index=$OPTARG                      
      if [ ! -f  ${index}.bwt ];then              
         echo "the index file $index not exist"
         exit
      elif [ ! -f  ${index}.pac ];then 
      	 echo "the index file $index not exist"
         exit
      elif [ ! -f  ${index}.ann ];then 
      	 echo "the index file $index not exist"
         exit
      elif [ ! -f  ${index}.amb ];then 
      	 echo "the index file $index not exist"
         exit
      elif [ ! -f  ${index}.sa ];then 
      	 echo "the index file $index not exist"
         exit
      fi
      ;;
    m)
      mode=$OPTARG
      ;; 
    x)
      #echo "$varname"
      #echo "$OPTARG"
      gatk_ann_dir=$OPTARG
      ###
      indel1="${gatk_ann_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
      ###
      indel2="${gatk_ann_dir}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" 
      ###
      dbsnp="${gatk_ann_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
      ###
      hapmap="${gatk_ann_dir}/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
      ###
      phase="${gatk_ann_dir}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
      ###
      omni="${gatk_ann_dir}/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
      ###
      vqsr_indel="${gatk_ann_dir}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

      if [ ! -d  ${gatk_ann_dir} ];then              
         echo "ERROR:WARNING:the gatk anntation directory ${gatk_index_dir} does not exist"
         exit
      elif [ ! -f  ${gatk_ann_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf ];then 
         echo "ERROR:the dbsnp file ${gatk_index_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf not exist"
         exit
      elif [ ! -f  ${gatk_ann_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz ];then 
         echo "ERROR:the known indel file ${gatk_index_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz not exist"
         exit
      elif [ ! -f  ${gatk_ann_dir}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ];then 
         echo "ERROR:the 1000 genome indel annotation file ${gatk_index_dir}/${gatk_ann_dir}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz not exist"
         exit
      elif [ ! -f  ${gatk_ann_dir}/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz ];then 
         echo "ERROR:the hapmap variants file ${gatk_ann_dir}/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz not exist"
         exit
      elif [ ! -f  ${gatk_ann_dir}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz ];then 
         echo "ERROR:the 1000 genome phase1 variants file ${gatk_ann_dir}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz not exist"
         exit
      elif [ ! -f  ${gatk_ann_dir}/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz ];then 
         echo "ERROR:the 1000 genome omni variants file ${gatk_ann_dir}/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz not exist"
         exit
      fi
      ;; 
    t)
      #echo "$varname"
      #echo "$OPTARG"
      thread=$OPTARG                    
      ;;  
    :)                                              
      echo "$varname"

      echo "the option -$OPTARG require an arguement"        
      exit 1
      ;;
    ?)                                   
      echo "$varname"
      echo "Invaild option: -$OPTARG"          
      Usage
      exit 2
      ;;
    esac
done




################################################
################################################
################################################
################################################
##  start quality control                     ##
################################################
################################################
################################################
################################################

source ./configure


cd $raw_reads_dir
    if [ ! -d ../clean_reads ];then
      mkdir -p ../clean_reads
    fi 

if [ $mode=='all' ];then

    echo "start fastp running"

    ls *gz|cut -d"_" -f 1 |sort -u |while read id;do
      if [ ! -f  ${id}_1.fastq.gz ];then               
         echo "The ${id}.1.fq.gz file  does not exist "
         exit
      elif [ ! -f  ${id}_2.fastq.gz ];then             
         echo "The ${id}.2.fq.gz file does not exist "
         exit
      fi    
    $FASTP -i ${id}_1.fastq.gz -o ../clean_reads/${id}.1.clean.fq.gz -I ${id}_2.fastq.gz -O ../clean_reads/${id}.2.clean.fq.gz
      if [ ! -f  ../clean_reads/${id}.1.clean.fq.gz ];then              
         echo "The fastp  does not executed correctly "
         exit
      elif [ ! -f  ../clean_reads/${id}.2.clean.fq.gz ];then             
         echo "The fastp  does not executed correctly "
         exit
      fi    
    done

elif [ $mode=='mapping' ];then

    echo "skipped fastp running"

elif [ $mode=='markduplicates' ];then

    echo "skipped fastp running"

elif [ $mode=='baserecalibrator' ];then
    
    echo "skipped fastp running"

elif [ $mode=='haplotypecaller' ];then
    
    echo "skipped fastp running"

elif [ $mode=='genotypegvcfs' ];then
    
    echo "skipped fastp running"

elif [ $mode=='vcfsplit' ];then
    
    echo "skipped fastp running"

fi
################################################
################################################
################################################
################################################
##  start bwa mapping                         ##
################################################
################################################
################################################
################################################

echo "start bwa mapping"
cd ../clean_reads

if [ $mode=='all' ];then

    ls *gz|cut -d"." -f 1 |sort -u |while read id;do
      if [ ! -f  ${id}.1.clean.fq.gz ];then               
             echo "The ${id}.1.clean.fq.gz file  does not exist "
           exit
      elif [ ! -f  ${id}.2.clean.fq.gz ];then              
             echo "The ${id}.2.clean.fq.gz file does not exist "
             exit
      fi    
    $BWA mem -t $thread -M -R "@RG\tID:lane1\tPL:illumina\tLB:library\tSM:$id" $index ${id}.1.clean.fq.gz ${id}.2.clean.fq.gz > $output/${id}.sam
    $SAMTOOLS view -b -S $output/${id}.sam > $output/${id}.bam
    $SAMTOOLS sort $output/${id}.bam -o $output/${id}.sorted.bam
    $SAMTOOLS flagstat $output/${id}.sorted.bam > $output/${id}.sorted.flagstat  
    done


elif [ $mode=='mapping' ];then
    
    ls *gz|cut -d"." -f 1 |sort -u |while read id;do
      if [ ! -f  ${id}.1.clean.fq.gz ];then               
             echo "The ${id}.1.clean.fq.gz file  does not exist "
           exit
      elif [ ! -f  ${id}.2.clean.fq.gz ];then              
             echo "The ${id}.2.clean.fq.gz file does not exist "
             exit
      fi    
    $BWA mem -t $thread -M -R "@RG\tID:lane1\tPL:illumina\tLB:library\tSM:$id" $index ${id}.1.clean.fq.gz ${id}.2.clean.fq.gz > $output/${id}.sam
    $SAMTOOLS view -b -S $output/${id}.sam > $output/${id}.bam
    $SAMTOOLS sort $output/${id}.bam -o $output/${id}.sorted.bam
    $SAMTOOLS flagstat $output/${id}.sorted.bam > $output/${id}.sorted.flagstat  
    done

elif [ $mode=='baserecalibrator' ];then
    
    echo "skipped bwa mapping"

elif [ $mode=='markduplicates' ];then
    
    echo "skipped bwa mapping"


elif [ $mode=='haplotypecaller' ];then
    
    echo "skipped bwa mapping"

elif [ $mode=='genotypegvcfs' ];then
    
    echo "skipped bwa mapping"

elif [ $mode=='vcfsplit' ];then
    
    echo "skipped bwa mapping"

else 

    echo "error option -m $mode"

fi

cd $output
  if [ ! -d $output/gatk_bam ];then
    mkdir -p $output/gatk_bam
  fi 
  if [ ! -d $output/gatk_vcf ];then
    mkdir -p $output/gatk_vcf
  fi 

################################################
################################################
################################################
################################################
##   GATK MarkDuplicates                      ##
################################################
################################################
################################################
################################################

echo "start GATK MarkDuplicates"
cd $output

if [ $mode=='all' ];then

    ls *sorted.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.bam ];then               
             echo "The ${id}.sorted.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" MarkDuplicates \
    -I ${id}.sorted.bam -O $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    -M $output/gatk_bam/${id}.sorted.bam.metrics > log.mark
    $SAMTOOLS index $output/gatk_bam/${id}.sorted.MarkDuplicates.bam
        if [ ! -f  $output/gatk_bam/${id}.sorted.MarkDuplicates.bam ];then               
             echo "The GATK  MarkDuplicates program  for ${id}.sorted.bam does not executed correctly"
           exit
        fi   
    done

elif [ $mode=='mapping' ];then
    ls *sorted.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.bam ];then               
             echo "The ${id}.sorted.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" MarkDuplicates \
    -I ${id}.sorted.bam -O $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    -M $output/gatk_bam/${id}.sorted.bam.metrics > log.mark
    $SAMTOOLS index $output/gatk_bam/${id}.sorted.MarkDuplicates.bam
        if [ ! -f  $output/gatk_bam/${id}.sorted.MarkDuplicates.bam ];then               
             echo "The GATK  MarkDuplicates program  for ${id}.sorted.bam does not executed correctly"
           exit
        fi   
    done

elif [ $mode=='markduplicates' ];then
    ls *sorted.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.bam ];then               
             echo "The ${id}.sorted.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" MarkDuplicates \
    -I ${id}.sorted.bam -O $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    -M $output/gatk_bam/${id}.sorted.bam.metrics > log.mark
    $SAMTOOLS index $output/gatk_bam/${id}.sorted.MarkDuplicates.bam
        if [ ! -f  $output/gatk_bam/${id}.sorted.MarkDuplicates.bam ];then               
             echo "The GATK  MarkDuplicates program  for ${id}.sorted.bam does not executed correctly"
           exit
        fi   
    done

elif [ $mode=='baserecalibrator' ];then
    
    echo "skipped gatk markduplicates"


elif [ $mode=='haplotypecaller' ];then
    
    echo "skipped gatk markduplicates"

elif [ $mode=='genotypegvcfs' ];then
    
    echo "skipped gatk markduplicates"

elif [ $mode=='vcfsplit' ];then
    
    echo "skipped gatk markduplicates"

fi
################################################
################################################
################################################
################################################
##   GATK BRSQ                                ##
################################################
################################################
################################################
################################################
cd $output/gatk_bam


echo "start GATK BaseRecalibrator"

if [ $mode=='all' ];then
    ls *sorted.MarkDuplicates.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.MarkDuplicates.bam ];then               
             echo "The ${id}.sorted.MarkDuplicates.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" BaseRecalibrator \
    -R $genome_fasta -I $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    --known-sites $indel1 \
    --known-sites $indel2 \
    --known-sites $dbsnp \
    -O $output/gatk_bam/${id}.recal.data.table
        if [ ! -f  ${id}.recal.data.table ];then               
             echo "The GATK  BaseRecalibrator program  for ${id}.sorted.MarkDuplicates.bam does not executed correctly "
           exit
        fi
    done

    echo "start GATK ApplyBQSR"
    cd  $output/gatk_bam
    ls *sorted.MarkDuplicates.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.MarkDuplicates.bam ];then               
             echo "The ${id}.sorted.MarkDuplicates.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" ApplyBQSR \
    -R $genome_fasta -I $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    -bqsr $output/gatk_bam/${id}.recal.data.table \
    -O $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam
        if [ ! -f  ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
             echo "The GATK  ApplyBQSR program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
           exit
        fi   
    done

elif [ $mode=='mapping' ];then
    ls *sorted.MarkDuplicates.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.MarkDuplicates.bam ];then               
             echo "The ${id}.sorted.MarkDuplicates.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" BaseRecalibrator \
    -R $genome_fasta -I $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    --known-sites $indel1 \
    --known-sites $indel2 \
    --known-sites $dbsnp \
    -O $output/gatk_bam/${id}.recal.data.table
        if [ ! -f  ${id}.recal.data.table ];then               
             echo "The GATK  BaseRecalibrator program  for ${id}.sorted.MarkDuplicates.bam does not executed correctly "
           exit
        fi
    done

    echo "start GATK ApplyBQSR"
    cd  $output/gatk_bam
    ls *sorted.MarkDuplicates.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.MarkDuplicates.bam ];then               
             echo "The ${id}.sorted.MarkDuplicates.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" ApplyBQSR \
    -R $genome_fasta -I $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    -bqsr $output/gatk_bam/${id}.recal.data.table \
    -O $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam
        if [ ! -f  ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
             echo "The GATK  ApplyBQSR program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
           exit
        fi   
    done

elif [ $mode=='markduplicates' ];then
    ls *sorted.MarkDuplicates.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.MarkDuplicates.bam ];then               
             echo "The ${id}.sorted.MarkDuplicates.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" BaseRecalibrator \
    -R $genome_fasta -I $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    --known-sites $indel1 \
    --known-sites $indel2 \
    --known-sites $dbsnp \
    -O $output/gatk_bam/${id}.recal.data.table
        if [ ! -f  ${id}.recal.data.table ];then               
             echo "The GATK  BaseRecalibrator program  for ${id}.sorted.MarkDuplicates.bam does not executed correctly "
           exit
        fi
    done

    echo "start GATK ApplyBQSR"
    cd  $output/gatk_bam
    ls *sorted.MarkDuplicates.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.MarkDuplicates.bam ];then               
             echo "The ${id}.sorted.MarkDuplicates.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" ApplyBQSR \
    -R $genome_fasta -I $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    -bqsr $output/gatk_bam/${id}.recal.data.table \
    -O $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam
        if [ ! -f  ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
             echo "The GATK  ApplyBQSR program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
           exit
        fi   
    done

elif [ $mode=='baserecalibrator' ];then
    ls *sorted.MarkDuplicates.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.MarkDuplicates.bam ];then               
             echo "The ${id}.sorted.MarkDuplicates.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" BaseRecalibrator \
    -R $genome_fasta -I $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    --known-sites $indel1 \
    --known-sites $indel2 \
    --known-sites $dbsnp \
    -O $output/gatk_bam/${id}.recal.data.table
        if [ ! -f  ${id}.recal.data.table ];then               
             echo "The GATK  BaseRecalibrator program  for ${id}.sorted.MarkDuplicates.bam does not executed correctly "
           exit
        fi
    done

    echo "start GATK ApplyBQSR"
    cd  $output/gatk_bam
    ls *sorted.MarkDuplicates.bam|cut -d"." -f 1 |sort -u |while read id;do
        if [ ! -f  ${id}.sorted.MarkDuplicates.bam ];then               
             echo "The ${id}.sorted.MarkDuplicates.bam file  does not exist "
           exit
        fi
    $GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" ApplyBQSR \
    -R $genome_fasta -I $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
    -bqsr $output/gatk_bam/${id}.recal.data.table \
    -O $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam
        if [ ! -f  ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
             echo "The GATK  ApplyBQSR program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
           exit
        fi   
    done

elif [ $mode=='haplotypecaller' ];then
    
    echo "skipped gatk bqsr"

elif [ $mode=='genotypegvcfs' ];then
    
    echo "skipped gatk bqsr"

elif [ $mode=='vcfsplit' ];then
    
    echo "skipped gatk bqsr"
fi



################################################
################################################
################################################
################################################
##   GATK HaplotypeCaller                     ##
################################################
################################################
################################################
################################################
cd $output/gatk_bam 

if [ $mode=='all' ];then
      echo "start GATK HaplotypeCaller"

      ls *sorted.MarkDuplicates.BQSR.bam|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
               echo "${id}.sorted.MarkDuplicates.BQSR.bam file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" HaplotypeCaller -R $genome_fasta \
      --emit-ref-confidence GVCF -I $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam \
      -D $dbsnp -O $output/gatk_vcf/${id}.gvcf
          if [ ! -f  $output/gatk_vcf/${id}.gvcf ];then               
               echo "The GATK  HaplotypeCaller program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='mapping' ];then
      echo "start GATK HaplotypeCaller"

      ls *sorted.MarkDuplicates.BQSR.bam|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
               echo "${id}.sorted.MarkDuplicates.BQSR.bam file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" HaplotypeCaller -R $genome_fasta \
      --emit-ref-confidence GVCF -I $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam \
      -D $dbsnp -O $output/gatk_vcf/${id}.gvcf
          if [ ! -f  $output/gatk_vcf/${id}.gvcf ];then               
               echo "The GATK  HaplotypeCaller program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='markduplicates' ];then
      echo "start GATK HaplotypeCaller"

      ls *sorted.MarkDuplicates.BQSR.bam|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
               echo "${id}.sorted.MarkDuplicates.BQSR.bam file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" HaplotypeCaller -R $genome_fasta \
      --emit-ref-confidence GVCF -I $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam \
      -D $dbsnp -O $output/gatk_vcf/${id}.gvcf
          if [ ! -f  $output/gatk_vcf/${id}.gvcf ];then               
               echo "The GATK  HaplotypeCaller program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='baserecalibrator' ];then
      echo "start GATK HaplotypeCaller"

      ls *sorted.MarkDuplicates.BQSR.bam|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
               echo "${id}.sorted.MarkDuplicates.BQSR.bam file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" HaplotypeCaller -R $genome_fasta \
      --emit-ref-confidence GVCF -I $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam \
      -D $dbsnp -O $output/gatk_vcf/${id}.gvcf
          if [ ! -f  $output/gatk_vcf/${id}.gvcf ];then               
               echo "The GATK  HaplotypeCaller program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='haplotypecaller' ];then
      echo "start GATK HaplotypeCaller"

      ls *sorted.MarkDuplicates.BQSR.bam|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
               echo "${id}.sorted.MarkDuplicates.BQSR.bam file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" HaplotypeCaller -R $genome_fasta \
      --emit-ref-confidence GVCF -I $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam \
      -D $dbsnp -O $output/gatk_vcf/${id}.gvcf
          if [ ! -f  $output/gatk_vcf/${id}.gvcf ];then               
               echo "The GATK  HaplotypeCaller program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='genotypegvcfs' ];then
    
    echo "skipped gatk haplotypecaller"

elif [ $mode=='vcfsplit' ];then
    
    echo "skipped gatk haplotypecaller"

fi


################################################
################################################
################################################
################################################
##   GATK GenotypeGVCFs                       ##
################################################
################################################
################################################
################################################
cd $output/gatk_vcf

if [ $mode=='all' ];then
      echo "start GATK GenotypeGVCFs"
      ls *gvcf|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.gvcf ];then               
               echo "${id}.gvcf file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" GenotypeGVCFs \
      -R $genome_fasta -V $output/gatk_vcf/${id}.gvcf \
      -O $output/gatk_vcf/${id}.raw.vcf
          if [ ! -f  $output/gatk_vcf/${id}.raw.vcf ];then               
               echo "The GATK  GenotypeGVCFs program  for $output/gatk_vcf/${id}.gvcf does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='mapping' ];then
      echo "start GATK GenotypeGVCFs"
      ls *gvcf|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.gvcf ];then               
               echo "${id}.gvcf file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" GenotypeGVCFs \
      -R $genome_fasta -V $output/gatk_vcf/${id}.gvcf \
      -O $output/gatk_vcf/${id}.raw.vcf
          if [ ! -f  $output/gatk_vcf/${id}.raw.vcf ];then               
               echo "The GATK  GenotypeGVCFs program  for $output/gatk_vcf/${id}.gvcf does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='markduplicates' ];then
      echo "start GATK GenotypeGVCFs"
      ls *gvcf|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.gvcf ];then               
               echo "${id}.gvcf file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" GenotypeGVCFs \
      -R $genome_fasta -V $output/gatk_vcf/${id}.gvcf \
      -O $output/gatk_vcf/${id}.raw.vcf
          if [ ! -f  $output/gatk_vcf/${id}.raw.vcf ];then               
               echo "The GATK  GenotypeGVCFs program  for $output/gatk_vcf/${id}.gvcf does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='baserecalibrator' ];then
      echo "start GATK GenotypeGVCFs"
      ls *gvcf|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.gvcf ];then               
               echo "${id}.gvcf file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" GenotypeGVCFs \
      -R $genome_fasta -V $output/gatk_vcf/${id}.gvcf \
      -O $output/gatk_vcf/${id}.raw.vcf
          if [ ! -f  $output/gatk_vcf/${id}.raw.vcf ];then               
               echo "The GATK  GenotypeGVCFs program  for $output/gatk_vcf/${id}.gvcf does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='haplotypecaller' ];then
      echo "start GATK GenotypeGVCFs"
      ls *gvcf|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.gvcf ];then               
               echo "${id}.gvcf file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" GenotypeGVCFs \
      -R $genome_fasta -V $output/gatk_vcf/${id}.gvcf \
      -O $output/gatk_vcf/${id}.raw.vcf
          if [ ! -f  $output/gatk_vcf/${id}.raw.vcf ];then               
               echo "The GATK  GenotypeGVCFs program  for $output/gatk_vcf/${id}.gvcf does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='genotypegvcfs' ];then
      echo "start GATK GenotypeGVCFs"
      ls *gvcf|cut -d"." -f 1 |sort -u |while read id;do
          if [ ! -f ${id}.gvcf ];then               
               echo "${id}.gvcf file  does not exist "
             exit
          fi
      $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" GenotypeGVCFs \
      -R $genome_fasta -V $output/gatk_vcf/${id}.gvcf \
      -O $output/gatk_vcf/${id}.raw.vcf
          if [ ! -f  $output/gatk_vcf/${id}.raw.vcf ];then               
               echo "The GATK  GenotypeGVCFs program  for $output/gatk_vcf/${id}.gvcf does not executed correctly "
             exit
          fi  
      done

elif [ $mode=='vcfsplit' ];then
    
    echo "skipped gatk genotypegvcfs"

fi


################################################
################################################
################################################
################################################
##   VCF MERGE AND ANNOVAR ANNOTATION         ##
################################################
################################################
################################################
################################################

cd $output/gatk_vcf

if [ $mode=='all' ];then

      ls *raw.vcf | cut -d'.' -f 1 | sort -u | while read id;do
      $BCFTOOLS view -Oz -o $output/gatk_vcf/${id}.raw.vcf.gz ${id}.raw.vcf
      $BCFTOOLS index $output/gatk_vcf/${id}.raw.vcf.gz
      done


      cd $output/gatk_vcf
      touch merge.sh

      echo  "$BCFTOOLS merge \\">> merge.sh
      ls *raw.vcf | while read id;do
      echo  "$id \\" >> merge.sh
      done
      echo "-O merge.raw.all.vcf" >> merge.sh
      if [ -f  $output/gatk_vcf/merge.sh ];then               
          bash merge.sh
      else
          echo "$output/gatk_vcf/merge.sh does not exist!"
          exit
      fi  

      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK SelectVariants                      ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/gatk_vcf/merge.raw.all.vcf ];then               
          merge_vcf="$output/gatk_vcf/merge.raw.all.vcf"
      else
          echo "$output/gatk_vcf/merge.raw.all.vcf does not exist!"
          exit
      fi 

      cd $output

        if [ ! -d $output/annovar ];then
          mkdir -p $output/annovar
        fi 


      if [ -f  $merge_vcf ];then               
          
            $GATK SelectVariants -select-type SNP -V $merge_vcf -O $output/annovar/merged.snp.raw.vcf
            $GATK SelectVariants -select-type INDEL -V $merge_vcf -O $output/annovar/merged.indel.raw.vcf

      fi 


      if [ ! -f  $output/annovar/merged.snp.raw.vcf ];then               
          echo "$output/annovar/merged.snp.raw.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.indel.raw.vcf ];then
          echo "$output/annovar/merged.indel.raw.vcf does not exist!"
          exit
      fi 


      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK VQSR                                ##
      ################################################
      ################################################
      ################################################
      ################################################

      cd $output/annovar


      if [ -f  $output/annovar/merged.snp.raw.vcf ];then               

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
               VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.snp.raw.vcf \
              -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
              -resource:omini,known=false,training=true,truth=false,prior=12.0 $omni \
              -resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.snps.tranches \
              --rscript-file $output/annovar/merged.snps.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode SNP \
              -O $output/annovar/merged.snps.recal 

            $GATK ApplyVQSR \
              -R ${genome_fasta} \
              -V $output/annovar/merged.snp.raw.vcf \
              --tranches-file $output/annovar/merged.snps.tranches \
              -recal-file $output/annovar/merged.snps.recal \
              -mode SNP \
              -O $output/annovar/merged.VQSR.snps.vcf

      elif [ -f $output/annovar/merged.indel.raw.vcf ];then

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
              VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.indel.raw.vcf \
              -resource:indel,known=false,training=true,truth=true,prior=15.0 $vqsr_indel \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.indel.tranches \
              --rscript-file $output/annovar/merged.indel.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode INDEL \
              -O $output/annovar/merged.indel.recal 

            $GATK ApplyVQSR \
              -R ${genome_ref} \
              -V $output/annovar/merged.indel.raw.vcf \
              --tranches-file $output/annovar/merged.indel.tranches \
              -recal-file $output/annovar/merged.indel.recal  \
              -mode SNP \
              -O $output/annovar/merged.VQSR.indel.vcf

      fi 


      if [ ! -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          echo "$output/annovar/merged.VQSR.snps.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.VQSR.indel.vcf ];then
          echo "$output/annovar/merged.VQSR.indel.vcf does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   VCF TO ANNOVAR INPUT                     ##
      ################################################
      ################################################
      ################################################
      ################################################
      if [ -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.snps.vcf > $output/annovar/snp.avinput

      elif [ -f $output/annovar/merged.VQSR.indel.vcf ];then
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.indel.vcf > $output/annovar/indel.avinput

      fi 

      cd $output/annovar

      if [ ! -f  $output/annovar/snp.avinput ];then               
          echo "$output/annovar/snp.avinput does not exist!"
          exit
      elif [ ! -f $output/annovar/indel.avinput ];then
          echo "$output/annovar/indel.avinput does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   ANNOVAR ANNOTATION                       ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/annovar/snp.avinput ];then               
          
          $TABLE_ANNOVAR $output/annovar/snp.avinput $humandb 
            -buildver hg38 \
            -out snp_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout

      elif [ -f $output/annovar/indel.avinput ];then
          
          $TABLE_ANNOVAR $output/annovar/indel.avinput $humandb 
            -buildver hg38 \
            -out indel_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout
      fi 

elif [ $mode=='mapping' ];then

      ls *raw.vcf | cut -d'.' -f 1 | sort -u | while read id;do
      $BCFTOOLS view -Oz -o $output/gatk_vcf/${id}.raw.vcf.gz ${id}.raw.vcf
      $BCFTOOLS index $output/gatk_vcf/${id}.raw.vcf.gz
      done


      cd $output/gatk_vcf
      touch merge.sh

      echo  "$BCFTOOLS merge \\">> merge.sh
      ls *raw.vcf | while read id;do
      echo  "$id \\" >> merge.sh
      done
      echo "-O merge.raw.all.vcf" >> merge.sh
      if [ -f  $output/gatk_vcf/merge.sh ];then               
          bash merge.sh
      else
          echo "$output/gatk_vcf/merge.sh does not exist!"
          exit
      fi  

      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK SelectVariants                      ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/gatk_vcf/merge.raw.all.vcf ];then               
          merge_vcf="$output/gatk_vcf/merge.raw.all.vcf"
      else
          echo "$output/gatk_vcf/merge.raw.all.vcf does not exist!"
          exit
      fi 

      cd $output

        if [ ! -d $output/annovar ];then
          mkdir -p $output/annovar
        fi 


      if [ -f  $merge_vcf ];then               
          
            $GATK SelectVariants -select-type SNP -V $merge_vcf -O $output/annovar/merged.snp.raw.vcf
            $GATK SelectVariants -select-type INDEL -V $merge_vcf -O $output/annovar/merged.indel.raw.vcf

      fi 


      if [ ! -f  $output/annovar/merged.snp.raw.vcf ];then               
          echo "$output/annovar/merged.snp.raw.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.indel.raw.vcf ];then
          echo "$output/annovar/merged.indel.raw.vcf does not exist!"
          exit
      fi 


      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK VQSR                                ##
      ################################################
      ################################################
      ################################################
      ################################################

      cd $output/annovar


      if [ -f  $output/annovar/merged.snp.raw.vcf ];then               

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
               VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.snp.raw.vcf \
              -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
              -resource:omini,known=false,training=true,truth=false,prior=12.0 $omni \
              -resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.snps.tranches \
              --rscript-file $output/annovar/merged.snps.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode SNP \
              -O $output/annovar/merged.snps.recal 

            $GATK ApplyVQSR \
              -R ${genome_fasta} \
              -V $output/annovar/merged.snp.raw.vcf \
              --tranches-file $output/annovar/merged.snps.tranches \
              -recal-file $output/annovar/merged.snps.recal \
              -mode SNP \
              -O $output/annovar/merged.VQSR.snps.vcf

      elif [ -f $output/annovar/merged.indel.raw.vcf ];then

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
              VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.indel.raw.vcf \
              -resource:indel,known=false,training=true,truth=true,prior=15.0 $vqsr_indel \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.indel.tranches \
              --rscript-file $output/annovar/merged.indel.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode INDEL \
              -O $output/annovar/merged.indel.recal 

            $GATK ApplyVQSR \
              -R ${genome_ref} \
              -V $output/annovar/merged.indel.raw.vcf \
              --tranches-file $output/annovar/merged.indel.tranches \
              -recal-file $output/annovar/merged.indel.recal  \
              -mode SNP \
              -O $output/annovar/merged.VQSR.indel.vcf

      fi 


      if [ ! -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          echo "$output/annovar/merged.VQSR.snps.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.VQSR.indel.vcf ];then
          echo "$output/annovar/merged.VQSR.indel.vcf does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   VCF TO ANNOVAR INPUT                     ##
      ################################################
      ################################################
      ################################################
      ################################################
      if [ -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.snps.vcf > $output/annovar/snp.avinput

      elif [ -f $output/annovar/merged.VQSR.indel.vcf ];then
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.indel.vcf > $output/annovar/indel.avinput

      fi 

      cd $output/annovar

      if [ ! -f  $output/annovar/snp.avinput ];then               
          echo "$output/annovar/snp.avinput does not exist!"
          exit
      elif [ ! -f $output/annovar/indel.avinput ];then
          echo "$output/annovar/indel.avinput does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   ANNOVAR ANNOTATION                       ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/annovar/snp.avinput ];then               
          
          $TABLE_ANNOVAR $output/annovar/snp.avinput $humandb 
            -buildver hg38 \
            -out snp_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout

      elif [ -f $output/annovar/indel.avinput ];then
          
          $TABLE_ANNOVAR $output/annovar/indel.avinput $humandb 
            -buildver hg38 \
            -out indel_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout
      fi 

elif [ $mode=='markduplicates' ];then

      ls *raw.vcf | cut -d'.' -f 1 | sort -u | while read id;do
      $BCFTOOLS view -Oz -o $output/gatk_vcf/${id}.raw.vcf.gz ${id}.raw.vcf
      $BCFTOOLS index $output/gatk_vcf/${id}.raw.vcf.gz
      done


      cd $output/gatk_vcf
      touch merge.sh

      echo  "$BCFTOOLS merge \\">> merge.sh
      ls *raw.vcf | while read id;do
      echo  "$id \\" >> merge.sh
      done
      echo "-O merge.raw.all.vcf" >> merge.sh
      if [ -f  $output/gatk_vcf/merge.sh ];then               
          bash merge.sh
      else
          echo "$output/gatk_vcf/merge.sh does not exist!"
          exit
      fi  

      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK SelectVariants                      ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/gatk_vcf/merge.raw.all.vcf ];then               
          merge_vcf="$output/gatk_vcf/merge.raw.all.vcf"
      else
          echo "$output/gatk_vcf/merge.raw.all.vcf does not exist!"
          exit
      fi 

      cd $output

        if [ ! -d $output/annovar ];then
          mkdir -p $output/annovar
        fi 


      if [ -f  $merge_vcf ];then               
          
            $GATK SelectVariants -select-type SNP -V $merge_vcf -O $output/annovar/merged.snp.raw.vcf
            $GATK SelectVariants -select-type INDEL -V $merge_vcf -O $output/annovar/merged.indel.raw.vcf

      fi 


      if [ ! -f  $output/annovar/merged.snp.raw.vcf ];then               
          echo "$output/annovar/merged.snp.raw.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.indel.raw.vcf ];then
          echo "$output/annovar/merged.indel.raw.vcf does not exist!"
          exit
      fi 


      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK VQSR                                ##
      ################################################
      ################################################
      ################################################
      ################################################

      cd $output/annovar


      if [ -f  $output/annovar/merged.snp.raw.vcf ];then               

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
               VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.snp.raw.vcf \
              -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
              -resource:omini,known=false,training=true,truth=false,prior=12.0 $omni \
              -resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.snps.tranches \
              --rscript-file $output/annovar/merged.snps.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode SNP \
              -O $output/annovar/merged.snps.recal 

            $GATK ApplyVQSR \
              -R ${genome_fasta} \
              -V $output/annovar/merged.snp.raw.vcf \
              --tranches-file $output/annovar/merged.snps.tranches \
              -recal-file $output/annovar/merged.snps.recal \
              -mode SNP \
              -O $output/annovar/merged.VQSR.snps.vcf

      elif [ -f $output/annovar/merged.indel.raw.vcf ];then

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
              VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.indel.raw.vcf \
              -resource:indel,known=false,training=true,truth=true,prior=15.0 $vqsr_indel \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.indel.tranches \
              --rscript-file $output/annovar/merged.indel.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode INDEL \
              -O $output/annovar/merged.indel.recal 

            $GATK ApplyVQSR \
              -R ${genome_ref} \
              -V $output/annovar/merged.indel.raw.vcf \
              --tranches-file $output/annovar/merged.indel.tranches \
              -recal-file $output/annovar/merged.indel.recal  \
              -mode SNP \
              -O $output/annovar/merged.VQSR.indel.vcf

      fi 


      if [ ! -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          echo "$output/annovar/merged.VQSR.snps.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.VQSR.indel.vcf ];then
          echo "$output/annovar/merged.VQSR.indel.vcf does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   VCF TO ANNOVAR INPUT                     ##
      ################################################
      ################################################
      ################################################
      ################################################
      if [ -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.snps.vcf > $output/annovar/snp.avinput

      elif [ -f $output/annovar/merged.VQSR.indel.vcf ];then
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.indel.vcf > $output/annovar/indel.avinput

      fi 

      cd $output/annovar

      if [ ! -f  $output/annovar/snp.avinput ];then               
          echo "$output/annovar/snp.avinput does not exist!"
          exit
      elif [ ! -f $output/annovar/indel.avinput ];then
          echo "$output/annovar/indel.avinput does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   ANNOVAR ANNOTATION                       ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/annovar/snp.avinput ];then               
          
          $TABLE_ANNOVAR $output/annovar/snp.avinput $humandb 
            -buildver hg38 \
            -out snp_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout

      elif [ -f $output/annovar/indel.avinput ];then
          
          $TABLE_ANNOVAR $output/annovar/indel.avinput $humandb 
            -buildver hg38 \
            -out indel_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout
      fi 

elif [ $mode=='baserecalibrator' ];then

      ls *raw.vcf | cut -d'.' -f 1 | sort -u | while read id;do
      $BCFTOOLS view -Oz -o $output/gatk_vcf/${id}.raw.vcf.gz ${id}.raw.vcf
      $BCFTOOLS index $output/gatk_vcf/${id}.raw.vcf.gz
      done


      cd $output/gatk_vcf
      touch merge.sh

      echo  "$BCFTOOLS merge \\">> merge.sh
      ls *raw.vcf | while read id;do
      echo  "$id \\" >> merge.sh
      done
      echo "-O merge.raw.all.vcf" >> merge.sh
      if [ -f  $output/gatk_vcf/merge.sh ];then               
          bash merge.sh
      else
          echo "$output/gatk_vcf/merge.sh does not exist!"
          exit
      fi  

      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK SelectVariants                      ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/gatk_vcf/merge.raw.all.vcf ];then               
          merge_vcf="$output/gatk_vcf/merge.raw.all.vcf"
      else
          echo "$output/gatk_vcf/merge.raw.all.vcf does not exist!"
          exit
      fi 

      cd $output

        if [ ! -d $output/annovar ];then
          mkdir -p $output/annovar
        fi 


      if [ -f  $merge_vcf ];then               
          
            $GATK SelectVariants -select-type SNP -V $merge_vcf -O $output/annovar/merged.snp.raw.vcf
            $GATK SelectVariants -select-type INDEL -V $merge_vcf -O $output/annovar/merged.indel.raw.vcf

      fi 


      if [ ! -f  $output/annovar/merged.snp.raw.vcf ];then               
          echo "$output/annovar/merged.snp.raw.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.indel.raw.vcf ];then
          echo "$output/annovar/merged.indel.raw.vcf does not exist!"
          exit
      fi 


      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK VQSR                                ##
      ################################################
      ################################################
      ################################################
      ################################################

      cd $output/annovar


      if [ -f  $output/annovar/merged.snp.raw.vcf ];then               

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
               VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.snp.raw.vcf \
              -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
              -resource:omini,known=false,training=true,truth=false,prior=12.0 $omni \
              -resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.snps.tranches \
              --rscript-file $output/annovar/merged.snps.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode SNP \
              -O $output/annovar/merged.snps.recal 

            $GATK ApplyVQSR \
              -R ${genome_fasta} \
              -V $output/annovar/merged.snp.raw.vcf \
              --tranches-file $output/annovar/merged.snps.tranches \
              -recal-file $output/annovar/merged.snps.recal \
              -mode SNP \
              -O $output/annovar/merged.VQSR.snps.vcf

      elif [ -f $output/annovar/merged.indel.raw.vcf ];then

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
              VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.indel.raw.vcf \
              -resource:indel,known=false,training=true,truth=true,prior=15.0 $vqsr_indel \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.indel.tranches \
              --rscript-file $output/annovar/merged.indel.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode INDEL \
              -O $output/annovar/merged.indel.recal 

            $GATK ApplyVQSR \
              -R ${genome_ref} \
              -V $output/annovar/merged.indel.raw.vcf \
              --tranches-file $output/annovar/merged.indel.tranches \
              -recal-file $output/annovar/merged.indel.recal  \
              -mode SNP \
              -O $output/annovar/merged.VQSR.indel.vcf

      fi 


      if [ ! -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          echo "$output/annovar/merged.VQSR.snps.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.VQSR.indel.vcf ];then
          echo "$output/annovar/merged.VQSR.indel.vcf does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   VCF TO ANNOVAR INPUT                     ##
      ################################################
      ################################################
      ################################################
      ################################################
      if [ -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.snps.vcf > $output/annovar/snp.avinput

      elif [ -f $output/annovar/merged.VQSR.indel.vcf ];then
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.indel.vcf > $output/annovar/indel.avinput

      fi 

      cd $output/annovar

      if [ ! -f  $output/annovar/snp.avinput ];then               
          echo "$output/annovar/snp.avinput does not exist!"
          exit
      elif [ ! -f $output/annovar/indel.avinput ];then
          echo "$output/annovar/indel.avinput does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   ANNOVAR ANNOTATION                       ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/annovar/snp.avinput ];then               
          
          $TABLE_ANNOVAR $output/annovar/snp.avinput $humandb 
            -buildver hg38 \
            -out snp_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout

      elif [ -f $output/annovar/indel.avinput ];then
          
          $TABLE_ANNOVAR $output/annovar/indel.avinput $humandb 
            -buildver hg38 \
            -out indel_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout
      fi 

elif [ $mode=='haplotypecaller' ];then

      ls *raw.vcf | cut -d'.' -f 1 | sort -u | while read id;do
      $BCFTOOLS view -Oz -o $output/gatk_vcf/${id}.raw.vcf.gz ${id}.raw.vcf
      $BCFTOOLS index $output/gatk_vcf/${id}.raw.vcf.gz
      done


      cd $output/gatk_vcf
      touch merge.sh

      echo  "$BCFTOOLS merge \\">> merge.sh
      ls *raw.vcf | while read id;do
      echo  "$id \\" >> merge.sh
      done
      echo "-O merge.raw.all.vcf" >> merge.sh
      if [ -f  $output/gatk_vcf/merge.sh ];then               
          bash merge.sh
      else
          echo "$output/gatk_vcf/merge.sh does not exist!"
          exit
      fi  

      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK SelectVariants                      ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/gatk_vcf/merge.raw.all.vcf ];then               
          merge_vcf="$output/gatk_vcf/merge.raw.all.vcf"
      else
          echo "$output/gatk_vcf/merge.raw.all.vcf does not exist!"
          exit
      fi 

      cd $output

        if [ ! -d $output/annovar ];then
          mkdir -p $output/annovar
        fi 


      if [ -f  $merge_vcf ];then               
          
            $GATK SelectVariants -select-type SNP -V $merge_vcf -O $output/annovar/merged.snp.raw.vcf
            $GATK SelectVariants -select-type INDEL -V $merge_vcf -O $output/annovar/merged.indel.raw.vcf

      fi 


      if [ ! -f  $output/annovar/merged.snp.raw.vcf ];then               
          echo "$output/annovar/merged.snp.raw.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.indel.raw.vcf ];then
          echo "$output/annovar/merged.indel.raw.vcf does not exist!"
          exit
      fi 


      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK VQSR                                ##
      ################################################
      ################################################
      ################################################
      ################################################

      cd $output/annovar


      if [ -f  $output/annovar/merged.snp.raw.vcf ];then               

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
               VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.snp.raw.vcf \
              -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
              -resource:omini,known=false,training=true,truth=false,prior=12.0 $omni \
              -resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.snps.tranches \
              --rscript-file $output/annovar/merged.snps.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode SNP \
              -O $output/annovar/merged.snps.recal 

            $GATK ApplyVQSR \
              -R ${genome_fasta} \
              -V $output/annovar/merged.snp.raw.vcf \
              --tranches-file $output/annovar/merged.snps.tranches \
              -recal-file $output/annovar/merged.snps.recal \
              -mode SNP \
              -O $output/annovar/merged.VQSR.snps.vcf

      elif [ -f $output/annovar/merged.indel.raw.vcf ];then

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
              VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.indel.raw.vcf \
              -resource:indel,known=false,training=true,truth=true,prior=15.0 $vqsr_indel \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.indel.tranches \
              --rscript-file $output/annovar/merged.indel.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode INDEL \
              -O $output/annovar/merged.indel.recal 

            $GATK ApplyVQSR \
              -R ${genome_ref} \
              -V $output/annovar/merged.indel.raw.vcf \
              --tranches-file $output/annovar/merged.indel.tranches \
              -recal-file $output/annovar/merged.indel.recal  \
              -mode SNP \
              -O $output/annovar/merged.VQSR.indel.vcf

      fi 


      if [ ! -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          echo "$output/annovar/merged.VQSR.snps.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.VQSR.indel.vcf ];then
          echo "$output/annovar/merged.VQSR.indel.vcf does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   VCF TO ANNOVAR INPUT                     ##
      ################################################
      ################################################
      ################################################
      ################################################
      if [ -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.snps.vcf > $output/annovar/snp.avinput

      elif [ -f $output/annovar/merged.VQSR.indel.vcf ];then
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.indel.vcf > $output/annovar/indel.avinput

      fi 

      cd $output/annovar

      if [ ! -f  $output/annovar/snp.avinput ];then               
          echo "$output/annovar/snp.avinput does not exist!"
          exit
      elif [ ! -f $output/annovar/indel.avinput ];then
          echo "$output/annovar/indel.avinput does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   ANNOVAR ANNOTATION                       ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/annovar/snp.avinput ];then               
          
          $TABLE_ANNOVAR $output/annovar/snp.avinput $humandb 
            -buildver hg38 \
            -out snp_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout

      elif [ -f $output/annovar/indel.avinput ];then
          
          $TABLE_ANNOVAR $output/annovar/indel.avinput $humandb 
            -buildver hg38 \
            -out indel_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout
      fi 

elif [ $mode=='genotypegvcfs' ];then

      ls *raw.vcf | cut -d'.' -f 1 | sort -u | while read id;do
      $BCFTOOLS view -Oz -o $output/gatk_vcf/${id}.raw.vcf.gz ${id}.raw.vcf
      $BCFTOOLS index $output/gatk_vcf/${id}.raw.vcf.gz
      done


      cd $output/gatk_vcf
      touch merge.sh

      echo  "$BCFTOOLS merge \\">> merge.sh
      ls *raw.vcf | while read id;do
      echo  "$id \\" >> merge.sh
      done
      echo "-O merge.raw.all.vcf" >> merge.sh
      if [ -f  $output/gatk_vcf/merge.sh ];then               
          bash merge.sh
      else
          echo "$output/gatk_vcf/merge.sh does not exist!"
          exit
      fi  

      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK SelectVariants                      ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/gatk_vcf/merge.raw.all.vcf ];then               
          merge_vcf="$output/gatk_vcf/merge.raw.all.vcf"
      else
          echo "$output/gatk_vcf/merge.raw.all.vcf does not exist!"
          exit
      fi 

      cd $output

        if [ ! -d $output/annovar ];then
          mkdir -p $output/annovar
        fi 


      if [ -f  $merge_vcf ];then               
          
            $GATK SelectVariants -select-type SNP -V $merge_vcf -O $output/annovar/merged.snp.raw.vcf
            $GATK SelectVariants -select-type INDEL -V $merge_vcf -O $output/annovar/merged.indel.raw.vcf

      fi 


      if [ ! -f  $output/annovar/merged.snp.raw.vcf ];then               
          echo "$output/annovar/merged.snp.raw.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.indel.raw.vcf ];then
          echo "$output/annovar/merged.indel.raw.vcf does not exist!"
          exit
      fi 


      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK VQSR                                ##
      ################################################
      ################################################
      ################################################
      ################################################

      cd $output/annovar


      if [ -f  $output/annovar/merged.snp.raw.vcf ];then               

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
               VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.snp.raw.vcf \
              -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
              -resource:omini,known=false,training=true,truth=false,prior=12.0 $omni \
              -resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.snps.tranches \
              --rscript-file $output/annovar/merged.snps.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode SNP \
              -O $output/annovar/merged.snps.recal 

            $GATK ApplyVQSR \
              -R ${genome_fasta} \
              -V $output/annovar/merged.snp.raw.vcf \
              --tranches-file $output/annovar/merged.snps.tranches \
              -recal-file $output/annovar/merged.snps.recal \
              -mode SNP \
              -O $output/annovar/merged.VQSR.snps.vcf

      elif [ -f $output/annovar/merged.indel.raw.vcf ];then

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
              VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.indel.raw.vcf \
              -resource:indel,known=false,training=true,truth=true,prior=15.0 $vqsr_indel \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.indel.tranches \
              --rscript-file $output/annovar/merged.indel.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode INDEL \
              -O $output/annovar/merged.indel.recal 

            $GATK ApplyVQSR \
              -R ${genome_ref} \
              -V $output/annovar/merged.indel.raw.vcf \
              --tranches-file $output/annovar/merged.indel.tranches \
              -recal-file $output/annovar/merged.indel.recal  \
              -mode SNP \
              -O $output/annovar/merged.VQSR.indel.vcf

      fi 


      if [ ! -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          echo "$output/annovar/merged.VQSR.snps.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.VQSR.indel.vcf ];then
          echo "$output/annovar/merged.VQSR.indel.vcf does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   VCF TO ANNOVAR INPUT                     ##
      ################################################
      ################################################
      ################################################
      ################################################
      if [ -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.snps.vcf > $output/annovar/snp.avinput

      elif [ -f $output/annovar/merged.VQSR.indel.vcf ];then
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.indel.vcf > $output/annovar/indel.avinput

      fi 

      cd $output/annovar

      if [ ! -f  $output/annovar/snp.avinput ];then               
          echo "$output/annovar/snp.avinput does not exist!"
          exit
      elif [ ! -f $output/annovar/indel.avinput ];then
          echo "$output/annovar/indel.avinput does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   ANNOVAR ANNOTATION                       ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/annovar/snp.avinput ];then               
          
          $TABLE_ANNOVAR $output/annovar/snp.avinput $humandb 
            -buildver hg38 \
            -out snp_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout

      elif [ -f $output/annovar/indel.avinput ];then
          
          $TABLE_ANNOVAR $output/annovar/indel.avinput $humandb 
            -buildver hg38 \
            -out indel_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout
      fi 

elif [ $mode=='vcfsplit' ];then

      ls *raw.vcf | cut -d'.' -f 1 | sort -u | while read id;do
      $BCFTOOLS view -Oz -o $output/gatk_vcf/${id}.raw.vcf.gz ${id}.raw.vcf
      $BCFTOOLS index $output/gatk_vcf/${id}.raw.vcf.gz
      done


      cd $output/gatk_vcf
      touch merge.sh

      echo  "$BCFTOOLS merge \\">> merge.sh
      ls *raw.vcf | while read id;do
      echo  "$id \\" >> merge.sh
      done
      echo "-O merge.raw.all.vcf" >> merge.sh
      if [ -f  $output/gatk_vcf/merge.sh ];then               
          bash merge.sh
      else
          echo "$output/gatk_vcf/merge.sh does not exist!"
          exit
      fi  

      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK SelectVariants                      ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/gatk_vcf/merge.raw.all.vcf ];then               
          merge_vcf="$output/gatk_vcf/merge.raw.all.vcf"
      else
          echo "$output/gatk_vcf/merge.raw.all.vcf does not exist!"
          exit
      fi 

      cd $output

        if [ ! -d $output/annovar ];then
          mkdir -p $output/annovar
        fi 


      if [ -f  $merge_vcf ];then               
          
            $GATK SelectVariants -select-type SNP -V $merge_vcf -O $output/annovar/merged.snp.raw.vcf
            $GATK SelectVariants -select-type INDEL -V $merge_vcf -O $output/annovar/merged.indel.raw.vcf

      fi 


      if [ ! -f  $output/annovar/merged.snp.raw.vcf ];then               
          echo "$output/annovar/merged.snp.raw.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.indel.raw.vcf ];then
          echo "$output/annovar/merged.indel.raw.vcf does not exist!"
          exit
      fi 


      ################################################
      ################################################
      ################################################
      ################################################
      ##   GATK VQSR                                ##
      ################################################
      ################################################
      ################################################
      ################################################

      cd $output/annovar


      if [ -f  $output/annovar/merged.snp.raw.vcf ];then               

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
               VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.snp.raw.vcf \
              -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
              -resource:omini,known=false,training=true,truth=false,prior=12.0 $omni \
              -resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.snps.tranches \
              --rscript-file $output/annovar/merged.snps.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode SNP \
              -O $output/annovar/merged.snps.recal 

            $GATK ApplyVQSR \
              -R ${genome_fasta} \
              -V $output/annovar/merged.snp.raw.vcf \
              --tranches-file $output/annovar/merged.snps.tranches \
              -recal-file $output/annovar/merged.snps.recal \
              -mode SNP \
              -O $output/annovar/merged.VQSR.snps.vcf

      elif [ -f $output/annovar/merged.indel.raw.vcf ];then

            $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
              VariantRecalibrator \
              -R $genome_fasta \
              -V $output/annovar/merged.indel.raw.vcf \
              -resource:indel,known=false,training=true,truth=true,prior=15.0 $vqsr_indel \
              -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
              --tranches-file $output/annovar/merged.indel.tranches \
              --rscript-file $output/annovar/merged.indel.plots.R \
              -an QD \
              -an MQ \
              -an MQRankSum \
              -an ReadPosRankSum \
              -an FS \
              -an SOR \
              -an DP \
              -mode INDEL \
              -O $output/annovar/merged.indel.recal 

            $GATK ApplyVQSR \
              -R ${genome_ref} \
              -V $output/annovar/merged.indel.raw.vcf \
              --tranches-file $output/annovar/merged.indel.tranches \
              -recal-file $output/annovar/merged.indel.recal  \
              -mode SNP \
              -O $output/annovar/merged.VQSR.indel.vcf

      fi 


      if [ ! -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          echo "$output/annovar/merged.VQSR.snps.vcf does not exist!"
          exit
      elif [ ! -f $output/annovar/merged.VQSR.indel.vcf ];then
          echo "$output/annovar/merged.VQSR.indel.vcf does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   VCF TO ANNOVAR INPUT                     ##
      ################################################
      ################################################
      ################################################
      ################################################
      if [ -f  $output/annovar/merged.VQSR.snps.vcf ];then               
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.snps.vcf > $output/annovar/snp.avinput

      elif [ -f $output/annovar/merged.VQSR.indel.vcf ];then
          
          $CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.indel.vcf > $output/annovar/indel.avinput

      fi 

      cd $output/annovar

      if [ ! -f  $output/annovar/snp.avinput ];then               
          echo "$output/annovar/snp.avinput does not exist!"
          exit
      elif [ ! -f $output/annovar/indel.avinput ];then
          echo "$output/annovar/indel.avinput does not exist!"
          exit
      fi 

      ################################################
      ################################################
      ################################################
      ################################################
      ##   ANNOVAR ANNOTATION                       ##
      ################################################
      ################################################
      ################################################
      ################################################

      if [ -f  $output/annovar/snp.avinput ];then               
          
          $TABLE_ANNOVAR $output/annovar/snp.avinput $humandb 
            -buildver hg38 \
            -out snp_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout

      elif [ -f $output/annovar/indel.avinput ];then
          
          $TABLE_ANNOVAR $output/annovar/indel.avinput $humandb 
            -buildver hg38 \
            -out indel_annotation \
            -remove -protocol refGene,cytoBand,esp6500siv2_all \
            -operation g,r,f \
            -nastring . -csvout
      fi 

fi



