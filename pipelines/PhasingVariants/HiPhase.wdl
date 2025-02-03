version 1.0

workflow RunHiPhase {
    input {
        File bam
        File bam_index

        File ref_fasta
        File ref_index

        File vcf
        File vcf_index
    }

    call HiPhase {
        input:
            bam = bam,
            bam_index = bam_index,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            vcf = vcf,
            vcf_index = vcf_index
    }

    output {
        File phased_vcf = HiPhase.phased_vcf
        File phased_vcf_index = HiPhase.phased_vcf_index
        File summary = HiPhase.summary
        File stats = HiPhase.stats
    }
}

task HiPhase {
    input {
        File bam
        File bam_index

        File ref_fasta
        File ref_index

        File vcf
        File vcf_index
        String sample_name = "SAMPLE"

        File? regions

        Int mem = 32
        Int cpu = 8
    }

    command <<<
        set -xueo pipefail

        echo "~{sample_name}" > sample_name.txt
        bcftools reheader -s sample_name.txt ~{vcf} -o reheadered.vcf.gz
        bcftools index -t -f reheadered.vcf.gz

        if [ -n "~{regions}" ]; then
            bcftools view -T ~{regions} -o subset.vcf.gz reheadered.vcf.gz
            bcftools index -t -f subset.vcf.gz
        else
            mv reheadered.vcf.gz subset.vcf.gz
            mv reheadered.vcf.gz.tbi subset.vcf.gz.tbi
        fi

        hiphase \
            --bam ~{bam} \
            --reference ~{ref_fasta} \
            --vcf subset.vcf.gz \
            --threads ~{cpu} \
            --verbose \
            --ignore-read-groups \
            --summary-file summary.tsv \
            --stats-file stats.tsv \
            --output-vcf phased.vcf.gz

        echo "Result after running hiphase:"
        ls

        bcftools index -t -f phased.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.1"
        disks: "local-disk " + 500 + " HDD"
        memory: mem + " GB"
        cpu: cpu
    }

    output {
        File phased_vcf = "phased.vcf.gz"
        File phased_vcf_index = "phased.vcf.gz.tbi"
        File summary = "summary.tsv"
        File stats = "stats.tsv"
    }
}