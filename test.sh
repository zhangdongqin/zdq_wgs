./easy_reseq.sh -t 32 -d /home/origene/yinshan.cui/resequencing_project/raw_data \
-o /home/origene/yinshan.cui/test2/wgs_out \
-i /home/origene/yinshan.cui/database/homo_sapiens/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome \
-r /home/origene/yinshan.cui/database/homo_sapiens/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
-k /home/data/genome_resequencing_DB/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
-g /home/data/genome_resequencing_DB/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-s /home/data/genome_resequencing_DB/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
