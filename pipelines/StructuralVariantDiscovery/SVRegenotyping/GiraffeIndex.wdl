version 1.0


workflow GiraffeIndexPrime {
    input {
        File input_vcf_gz
        String remote_id
        String remote_chromosomes_dir
        Int max_sv_length
        String remote_dir
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }

    call GiraffeIndexPrimeImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            remote_id = remote_id,
            remote_chromosomes_dir = remote_chromosomes_dir,
            max_sv_length = max_sv_length,
            remote_dir = remote_dir,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    output {
    }
}


# COMMAND              | TIME | CORES | RAM
# vg construct         | 5m   | 1     | 44 M
# vg index --dist-name | 1h   | 6     | 86 G
# vg index --xg-name   | 1h   | 1     | 32 G
# vg gbwt              | 20m  | 1     | 50 G
# vg gbwt --gbz-format | 3m   | 1     | 19 G
# vg minimizer         | 6m   | 13    | 36 G
#
task GiraffeIndexPrimeImpl {
    input {
        File input_vcf_gz
        String remote_id
        String remote_chromosomes_dir
        Int max_sv_length
        String remote_dir
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
        remote_dir: "Containing the VCF to be indexed."
        max_sv_length: "Used for cleaning."
    }

    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}

        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        VG_COMMAND="~{docker_dir}/vg"

        # Cleaning the VCF and the reference
        gsutil -m cp ~{remote_chromosomes_dir}/'*.fa' .
        gunzip -c ~{input_vcf_gz} > input.vcf
        java -cp ~{docker_dir} CleanVCFPrime input.vcf . ~{max_sv_length} 1 input-vg.vcf
        bgzip input-vg.vcf
        tabix input-vg.vcf.gz
        CHROMOSOMES=""
        for CHR in $(seq 1 22); do
            CHROMOSOMES="${CHROMOSOMES} chr${CHR}.fa"
        done
        CHROMOSOMES="${CHROMOSOMES} chrX.fa chrY.fa chrM.fa"
        cat ${CHROMOSOMES} > new_reference.fa
        rm -f ${CHROMOSOMES}
        samtools faidx new_reference.fa

        # Indexing
        ${TIME_COMMAND} ${VG_COMMAND} construct --threads ${N_THREADS} --progress --handle-sv --alt-paths --reference new_reference.fa --vcf input-vg.vcf.gz > graph.vg
        mkdir ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} index --threads ${N_THREADS} --temp-dir ./vgtmp --progress --dist-name graph.dist graph.vg
        rm -rf ./vgtmp; mkdir ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} index --threads ${N_THREADS} --temp-dir ./vgtmp --progress --xg-alts --xg-name graph.xg graph.vg
        rm -rf ./vgtmp; mkdir ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} gbwt --num-jobs ${N_THREADS} --temp-dir ./vgtmp --progress --path-cover --xg-name graph.xg --output graph.gbwt
        rm -rf ./vgtmp; mkdir ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} gbwt --num-jobs ${N_THREADS} --temp-dir ./vgtmp --progress --xg-name graph.xg --graph-name graph.gbz --gbz-format graph.gbwt
        rm -rf ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} minimizer --threads ${N_THREADS} --progress --distance-index graph.dist --output-name graph.min graph.gbz

        # Uploading all indexes
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp graph.xg ~{remote_dir}/~{remote_id}.xg && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading XG. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp graph.gbz ~{remote_dir}/~{remote_id}.gbz && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading GBZ. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
         while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp graph.min ~{remote_dir}/~{remote_id}.min && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading MIN. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp graph.dist ~{remote_dir}/~{remote_id}.dist && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading DIST. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
         while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp graph.vg ~{remote_dir}/~{remote_id}.vg && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading VG. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp graph.gbwt ~{remote_dir}/~{remote_id}.gbwt && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading GBWT. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>

    output {
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}