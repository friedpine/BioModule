#python pyscript00_all_references.py 

import cPickle as pickle
file = '/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/Ref.all.dat'

info = {}

info['mm10'] = {}


info['mm10']['fa'] = {}
info['mm10']['fa']['genome'] = '/data/Analysis/fanxiaoying/database/mm10/01.bowtie/mm10.fa'

info['mm10']['gtf'] = {}
info['mm10']['gtf']['ensembl'] = '/data/Analysis/fanxiaoying/database/mm10/03.ensembl/Mus_musculus.GRCm38.73.gtf'
info['mm10']['gtf']['refseq'] = '/data/Analysis/fanxiaoying/database/mm10/mm10.refseq.gtf'
info['mm10']['gtf']['3146transc'] = '/datc/fanxiaoying/01.super/02.assembly/z01.fasta/3146_transc.gtf'
print "HEHE"

info['mm10']['bwa']={}
info['mm10']['bwa']['refseq'] = '/data/Analysis/fanxiaoying/database/mm9/03.ucsc_ercc/ucsc.ercc.mRNA.fa'
info['mm10']['bwa']['gencode'] = '/data/Analysis/fanxiaoying/database/mm10/05.encode_M2/GENCODE_M2.seq.fa'
info['mm10']['bwa']['rrna'] = '/data/Analysis/fanxiaoying/database/mm9/04.rrna/mouse.rrna.fa'
info['mm10']['bwa']['refseq_NCBI'] = '/data/Analysis/fanxiaoying/database/mm9/mouse.NCBI.refgene.fa'
info['mm10']['bwa']['trinity_500'] = '/datc/fanxiaoying/01.super/02.assembly/z01.fasta/t500.fa'

info['mm10']['bowtie2'] = {}
info['mm10']['bowtie2']['genome'] = '/data/Analysis/fanxiaoying/database/mm10/01.bowtie/mm10'		
info['mm10']['bowtie2']['genome1'] = '/data/Analysis/fanxiaoying/database/mm10/02.GRCm38/GRCm38.no_patch_73.genome'	
info['mm10']['bowtie2']['refseq'] = '/data/Analysis/fanxiaoying/database/mm10/00.refseq/mouse'
info['mm10']['bowtie2']['refseq_NCBI'] = '/data/Analysis/fanxiaoying/database/mm10/00.refseq/mouse'


info['mm10']['bed'] = {}
info['mm10']['bed']['introns'] = '/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/stuff_beds/regions.mm10.intron.bed'
info['mm10']['bed']['introns_50'] = '/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/stuff_beds/regions.mm10.introns_50.bed'
info['mm10']['bed']['exons'] = '/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/stuff_beds/regions.mm10.exon.bed' 
info['mm10']['bed']['intergenes'] = '/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/stuff_beds/regions.mm10.intergene.bed'
info['mm10']['bed']['genomes'] = '/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/stuff_beds/genome.mm10.length.bed'
info['mm10']['bed']['ensembl'] = '/data/Analysis/fanxiaoying/database/mm10/03.ensembl/mm10.ensembl.bed'

info['hg19']  = {}
info['hg19']['fa'] = {}
info['hg19']['fa']['genome'] = '/data/Analysis/fanxiaoying/database/hg19/00.genome/genome.fa'

info['hg19']['gtf'] = {}
info['hg19']['gtf']['ensembl'] = '/data/Analysis/fanxiaoying/database/hg19/03.ensembl/Homo_sapiens.GRCh37.73.gtf'
info['hg19']['gtf']['gencode'] = '/data/Analysis/fanxiaoying/database/hg19/03.ensembl/gencode.v19.annotation.gtf'
info['hg19']['gtf']['refseq'] = '/data/Analysis/fanxiaoying/database/hg19/01.ucsc_mRNA/hg19_new.gtf'

info['hg19']['bwa'] = {}
info['hg19']['bwa']['refseq'] = '/data/Analysis/fanxiaoying/database/hg19/11.ucsc_ercc/refgene.add_ercc.fa'
info['hg19']['bwa']['genome'] = '/data/Analysis/fanxiaoying/database/hg19/00.genome/genome.fa'
info['hg19']['bwa']['rrna'] = '/data/Analysis/fanxiaoying/database/hg19/06.rrna/rrna.hg.fa'
info['hg19']['bwa']['gencode'] = '/data/Analysis/fanxiaoying/database/hg19/03.ensembl/GRCh37.73.cdna_ncrna.fa'

info['hg19']['bowtie2'] = {}
info['hg19']['bowtie2']['genome'] = '/data/Analysis/fanxiaoying/database/hg19/00.genome/genome'
info['hg19']['bowtie2']['refseq'] = '/data/Analysis/fanxiaoying/database/hg19/11.ucsc_ercc/refgene.add_ercc'
info['hg19']['bowtie2']['ensembl'] = '/data/Analysis/fanxiaoying/database/hg19/03.ensembl/Homo_sapiens.GRCh37.73'

info['hg19']['bed'] = {}
info['hg19']['bed']['exons'] = '/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/stuff_beds/hg19_exons.bed'
info['hg19']['bed']['introns'] = '/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/stuff_beds/hg19_introns.bed'

pickle.dump(info,open(file,'w'))


#BUILD_METHODS
info['build'] = {}
info['build']['mm10_bed_exons'] = 'module10_make_database.gtf_exon2bed(gtf,bedfile) & \
	/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/bedtools sort -i bedfile >sorted & \
	/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/mergeBed -i sorted >exon_regions'
info['build']['mm10_bed_introns'] = 'module10_make_database.gtf_intron2bed(gtffile,bedfile_intron,bedfile_intron_all,49) & \
	/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/bedtools sort -i bedfile_intron_all >sorted & \
	/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/subtractBed -a sorted -b exon_regions >intron_regions'
info['build']['mm10_bed_intergenes'] = 'cat exon_bed intron_bed > catfile & \
	/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/bedtools sort -i catfile >cat_sort &\
	/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/mergeBed -i cat_sort > mergefile &\
	/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/complementBed -i mergefile -g genomes > intergene_regions'


