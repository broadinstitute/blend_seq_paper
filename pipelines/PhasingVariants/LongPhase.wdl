version 1.0

workflow LongPhase {
    input {
        File lr_bam
        File lr_bam_index
        File? sr_bam
        File? sr_bam_index
        File snp_vcf
        File snp_vcf_index
        File? sv_vcf
        File? sv_vcf_index

        File ref_fasta
        File ref_fasta_index
    }

    call RunLongPhase {
        input:
            lr_bam = lr_bam,
            lr_bam_index = lr_bam_index,
            sr_bam = sr_bam,
            sr_bam_index = sr_bam_index,
            snp_vcf = snp_vcf,
            snp_vcf_index = snp_vcf_index,
            sv_vcf = sv_vcf,
            sv_vcf_index = sv_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    output {
        File output_vcf = RunLongPhase.output_vcf
        File output_vcf_index = RunLongPhase.output_vcf_index

        File output_sv_vcf = RunLongPhase.output_sv_vcf
        File output_sv_vcf_index = RunLongPhase.output_sv_vcf_index
    }
}

task RunLongPhase {
    input {
        File lr_bam
        File lr_bam_index
        File? sr_bam
        File? sr_bam_index
        File snp_vcf
        File snp_vcf_index
        File? sv_vcf
        File? sv_vcf_index

        File ref_fasta
        File ref_fasta_index

        String sample_name = "SAMPLE"

        File? regions

        Int cpu = 8
        Int mem = 32
        Int disk_size = 2 * ceil(size(lr_bam, "GB")) + 100
        String tech = "ONT"   # Or "PACBIO"
    }

    String tech_flag = if (tech == "ONT") then "--ont" else "--pb"

    command <<<
        set -xueo pipefail

        echo "~{sample_name}" > sample_name.txt
        bcftools reheader -s sample_name.txt ~{snp_vcf} -o reheadered-snp.vcf.gz
        bcftools index -t -f reheadered-snp.vcf.gz

        # echo "~{sample_name}" > sample_name.txt
        # bcftools reheader -s sample_name.txt ~{sv_vcf} -o reheadered-sv.vcf.gz
        # bcftools index -t -f reheadered-sv.vcf.gz

        if [ -n "~{regions}" ]; then
            bcftools view -T ~{regions} -o subset-snp.vcf.gz reheadered-snp.vcf.gz
            bcftools index -t -f subset-snp.vcf.gz
            # bcftools view -T ~{regions} -o subset-sv.vcf.gz reheadered-sv.vcf.gz
            # bcftools index -t -f subset-sv.vcf.gz
        else
            mv reheadered-snp.vcf.gz subset-snp.vcf.gz
            mv reheadered-snp.vcf.gz.tbi subset-snp.vcf.gz.tbi
            # mv reheadered-sv.vcf.gz subset-sv.vcf.gz
            # mv reheadered-sv.vcf.gz.tbi subset-sv.vcf.gz.tbi
        fi

        longphase phase \
            -s subset-snp.vcf.gz \
            -b ~{lr_bam} \
            ~{"-b " + sr_bam} \
            -r ~{ref_fasta} \
            -o output \
            -t ~{cpu} \
            --indels \
            ~{"--sv-file " + sv_vcf} \
            ~{tech_flag}

        bcftools view output.vcf -o output.vcf.gz -Wtbi
        if [ -f output_SV.vcf ]; then
            bcftools view output_SV.vcf -o output_SV.vcf.gz -Wtbi
        else
            touch output_SV.vcf.gz
            touch output_SV.vcf.gz.tbi
        fi
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/longphase:v1.0"
        memory: mem + "GB"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        maxRetries: 3
    }

    output {
        File output_vcf = "output.vcf.gz"
        File output_vcf_index = "output.vcf.gz.tbi"

        File output_sv_vcf = "output_SV.vcf.gz"
        File output_sv_vcf_index = "output_SV.vcf.gz.tbi"
    }
}