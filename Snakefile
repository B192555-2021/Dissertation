import gzip


REP_INDEX = {"HB2_EDSW220013382-2a_HKYYNDSX3_L3",
             "HB3_EDSW220013381-1a_HL5JWDSX3_L3",
             "HB5_EDSW220013380-1a_HL5JWDSX3_L3",
             "HB6_EDSW220013379-1a_HL5JWDSX3_L3",
             "HB7_EDSW220013378-2a_HKW2GDSX3_L1",
             "HB7_EDSW220013378-2a_HKYYNDSX3_L3",
             "HB8_EDSW220013377-1a_HL5WJDSX3_L2",
             "HB9_EDSW220013376-1a_HL5JWDSX3_L2"
            }
pairs = ["1", "2"]
threadsusing = int(32)


rule all:
    input:
        expand("data/mapping_statistics/{rep}_stats.txt", rep=REP_INDEX),
        expand("data/mapping_statistics/{rep}_flagstat.txt", rep=REP_INDEX),

rule unzip_ref_genome_gz:
    input:
        ref = "data/raw/MpTak_v6.1.genome.fasta.gz"
    output:
        ref = "data/raw/MpTak_v6.1.genome.fasta"
    run:
        shell("gunzip -c {input.ref} > {output.ref}")

rule indexing_reference_genome:
    input:
        "data/raw/MpTak_v6.1.genome.fasta"
    output:
        "data/index/MpTak_v6.1.genome_index"
#    threads:
#         threads = threadsusing
    run:
        shell("bowtie2-build {input} {output} -p {threadsusing}")
        shell("touch {output}")


rule mapping_to_ref:
    input:
        "data/raw/{rep}_1.fq.gz",
        "data/raw/{rep}_2.fq.gz",
        "data/index/MpTak_v6.1.genome_index",
    output:
        "data/mappedbam/{rep}.bam"
    shell:
        "bowtie2 -p {threadsusing} -x {input[2]} -1 {input[0]} -2 {input[1]} | \
        samtools sort -O bam -@ {threadsusing} -n -o - > {output}"
#-O [outputforamt] -@ [cpucore] -o [outfile] -n [Sort by read names]
#samtools sort: http://www.htslib.org/doc/samtools-sort.html

rule remove_duplicates:
    input:
        "data/mappedbam/{rep}.bam"
    output:
        rm_dup = "data/bam_mv_dup/{rep}.bam",
        position_sorted_bam = "data/bam_positionsort/{rep}.bam"
    run:
        shell("samtools fixmate -m -O bam -@ {threadsusing} {input} {output[1]}")
        shell("samtools sort -@ {threadsusing} {output[1]} -o {output[1]} ")
        shell("samtools markdup -@ {threadsusing} -s -r {output[1]} {output[0]}")

#fixmat -m [option for markdup]
#markdup -r [Remove duplicate reads]
#fixmat http://www.htslib.org/doc/samtools-fixmate.html
#markdup http://www.htslib.org/doc/samtools-markdup.html

rule generating_mapping_statistics:
    input:
        "data/bam_mv_dup/{rep}.bam"
    output:
        stats = "data/mapping_statistics/{rep}_stats.txt",
        flagstat = "data/mapping_statistics/{rep}_flagstat.txt"
    run:
        shell("samtools stats {input} > {output.stats}")
        shell("samtools flagstat {input} > {output.flagstat}")

# bai
#rule bai_index:
#    input:
#        "data/bam_mv_dup/{rep}.bam"
#    output:
#        "data/bam_mv_dup/{rep}.bam.bai"
#    run:
#        shell("samtools index {input} > {output}")

#rule SNP_calling:
#    input:
#        bam = expand("data/bam_mv_dup/{rep}.bam", rep=REP_INDEX),
#        ref = "data/raw/MpTak_v6.1.genome.fasta"
#    output:
#        vcf = "data/vcf/SNP_CALLING.vcf"
#    run:
#        shell("samtools index -@ {threadsusing} {input.bam}") #indexing for bam
#        shell("freebayes -f {input.ref} {input.bam} > {output.vcf}")

    
