version 1.0


workflow LowCoverageSniffles {
    input {
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        File reference_tandem_repeats
        String sniffles_lowcoverage_flags = "--minsupport 1 --qc-output-all --qc-coverage 1 --long-dup-coverage 1 --detect-large-ins True"
    }

    call LowCoverageSnifflesImpl {
        input:
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            reference_tandem_repeats = reference_tandem_repeats,
            sniffles_lowcoverage_flags = sniffles_lowcoverage_flags
    }

    output {
        File sniffles_vcf_gz = LowCoverageSnifflesImpl.sniffles_vcf_gz
        File sniffles_tbi = LowCoverageSnifflesImpl.sniffles_tbi
    }
}


task LowCoverageSnifflesImpl {
    input {
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        File reference_tandem_repeats
        String sniffles_lowcoverage_flags
    }
    parameter_meta {
    }

    Int ram_size_gb = 16
    Int disk_size_gb = 4*ceil(size(alignments_bam,"GB") + size(reference_fa,"GB") + size(reference_tandem_repeats,"GB"))
    String work_dir = "/cromwell_root/sniffles"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        ${TIME_COMMAND} sniffles --threads ${N_THREADS} ~{sniffles_lowcoverage_flags} --tandem-repeats ~{reference_tandem_repeats} --reference ~{reference_fa} --input ~{alignments_bam} --vcf sniffles.vcf
        bcftools sort --max-mem $(( ~{ram_size_gb} - 4 ))G --output-type z sniffles.vcf > sniffles.vcf.gz
        tabix sniffles.vcf.gz
    >>>

    output {
        File sniffles_vcf_gz = work_dir + "/sniffles.vcf.gz"
        File sniffles_tbi = work_dir + "/sniffles.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/sniffles"
        cpu: 16
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}