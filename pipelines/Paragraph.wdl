version 1.0


#
workflow Paragraph {
    input {
        File input_vcf_gz
        File input_tbi
        File short_reads_bam
        File short_reads_bai
        Int short_read_coverage
        Int short_read_length
        File reference_fa
        File reference_fai
        String remote_chromosomes_dir
        Int max_sv_length
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
        String sex = "M"
    }
    parameter_meta {
    }
    
    call ParagraphImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            short_reads_bam = short_reads_bam,
            short_reads_bai = short_reads_bai,
            short_read_coverage = short_read_coverage,
            short_read_length = short_read_length,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            remote_chromosomes_dir = remote_chromosomes_dir,
            max_sv_length = max_sv_length,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb,
            sex = sex
    }
    
    output {
        File regenotyped_vcf_gz = ParagraphImpl.regenotyped_vcf_gz
        File regenotyped_tbi = ParagraphImpl.regenotyped_tbi
        File paragraph_genotypes = ParagraphImpl.paragraph_genotypes
        File paragraph_variants = ParagraphImpl.paragraph_variants
    }
}


# COVERAGE | TIME | CORES | RAM
# 1x       | 4.5h |  23   | 5 G
#
task ParagraphImpl {
    input {
        File input_vcf_gz
        File input_tbi
        File short_reads_bam
        File short_reads_bai
        Int short_read_coverage
        Int short_read_length
        File reference_fa
        File reference_fai
        String remote_chromosomes_dir
        Int max_sv_length
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
        String sex
    }
    parameter_meta {
    }
    
    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        PARAGRAPH_DIR="~{docker_dir}/paragraph/bin"
        
        # Cleaning the VCF. This is needed to handle Paragraph's idosynchrasies.
        gsutil -m cp ~{remote_chromosomes_dir}/'*.fa' .
        gunzip -c ~{input_vcf_gz} > input.vcf
        java -cp ~{docker_dir} CleanVCFParagraph input.vcf . ~{max_sv_length} 1 tmp.vcf
        rm -f input.vcf
        bcftools view --header-only tmp.vcf | sed 's/FT/FT_KANPIG/g' > input-cleaned.vcf
        bcftools view --no-header tmp.vcf | awk '{ printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tGT\t0/1\n",$1,$2,$3,$4,$5,$6,$7,$8); }' >> input-cleaned.vcf
        rm -f tmp.vcf
        
        # Paragraph
        echo -e "id\tpath\tdepth\tread length\tsex" > manifest.txt
        echo -e "sample1\t~{short_reads_bam}\t~{short_read_coverage}\t~{short_read_length}\t~{sex}" >> manifest.txt
        mkdir ./tmpdir
        cd ${PARAGRAPH_DIR}
        ${TIME_COMMAND} python3 multigrmpy.py --threads ${N_THREADS} --verbose --input ~{work_dir}/input-cleaned.vcf --reference-sequence ~{reference_fa} --manifest ~{work_dir}/manifest.txt --output ~{work_dir}/tmpdir
        rm -f input-cleaned.vcf
        cd ~{work_dir}
        bcftools view -h ./tmpdir/genotypes.vcf.gz > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > out.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> out.vcf
        # Remark: Paragraph adds a new sample column rather than overwriting the
        # original sample column.
        bcftools view -H ./tmpdir/genotypes.vcf.gz | cut -f 1-9,11 >> out.vcf
        bcftools view --threads ${N_THREADS} --output-type z --output regenotyped.vcf.gz out.vcf
        rm -f out.vcf
        tabix -f regenotyped.vcf.gz
    >>>
    
    output {
        File regenotyped_vcf_gz = work_dir + "/regenotyped.vcf.gz"
        File regenotyped_tbi = work_dir + "/regenotyped.vcf.gz.tbi"
        File paragraph_genotypes = work_dir + "/tmpdir/genotypes.vcf.gz"
        File paragraph_variants = work_dir + "/tmpdir/variants.vcf.gz"
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
