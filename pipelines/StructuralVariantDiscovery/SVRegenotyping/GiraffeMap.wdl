version 1.0


workflow GiraffeMapPrime {
    input {
        File reads1_fastq_gz
        File reads2_fastq_gz
        String remote_indexes_dir
        String index_id
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
        remote_indexes_dir: "Containing the graph indexes."
    }

    call GiraffeMapPrimeImpl {
        input:
            reads1_fastq_gz = reads1_fastq_gz,
            reads2_fastq_gz = reads2_fastq_gz,
            remote_indexes_dir = remote_indexes_dir,
            index_id = index_id,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    output {
        File alignments_gam = GiraffeMapPrimeImpl.alignments_gam
    }
}


# COMMAND    | TIME | CORES | RAM
# vg giraffe | 1.5h |  all  | 32G
#
task GiraffeMapPrimeImpl {
    input {
        File reads1_fastq_gz
        File reads2_fastq_gz
        String remote_indexes_dir
        String index_id
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
        remote_indexes_dir: "Containing the graph indexes."
    }

    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}

        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        VG_COMMAND="~{docker_dir}/vg"


        while : ; do
            TEST=$(gsutil -m cp ~{remote_indexes_dir}/~{index_id}.gbz ~{remote_indexes_dir}/~{index_id}.min ~{remote_indexes_dir}/~{index_id}.dist . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading indexes. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        ${TIME_COMMAND} ${VG_COMMAND} giraffe --threads ${N_THREADS} --progress --output-format gam --gbz-name ~{index_id}.gbz --minimizer-name ~{index_id}.min --dist-name ~{index_id}.dist --fastq-in ~{reads1_fastq_gz} --fastq-in ~{reads2_fastq_gz} > alignments.gam
        #${TIME_COMMAND} ${VG_COMMAND} stats --threads ${N_THREADS} --alignments alignments.gam
    >>>

    output {
        File alignments_gam = work_dir + "/alignments.gam"
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}