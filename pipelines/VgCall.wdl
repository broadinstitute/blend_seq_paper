version 1.0


# Genotypes the same VCF that was used to build the graph.
#
workflow VgCall {
    input {
        File input_vcf_gz
        File input_tbi
        Int max_sv_length
        File input_gam
        String remote_indexes_dir
        String remote_chromosomes_dir
        String index_id
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    call VgCallImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            max_sv_length = max_sv_length,
            input_gam = input_gam,
            remote_indexes_dir = remote_indexes_dir,
            remote_chromosomes_dir = remote_chromosomes_dir,
            index_id = index_id,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    output {
        File regenotyped_vcf_gz = VgCallImpl.regenotyped_vcf_gz
        File regenotyped_tbi = VgCallImpl.regenotyped_tbi
    }
}


# COMMAND        | TIME | CORES | RAM
# vg pack        | 1h   |  35   | 61 G
# vg call        | 1h   |  30   | 115 G
# vg call --vcf  | 1h   |  32   | 120 G
#
# See this GitHub issue <https://github.com/vgteam/vg/issues/3950> for more
# details.
#
task VgCallImpl {
    input {
        File input_vcf_gz
        File input_tbi
        Int max_sv_length
        File input_gam
        String remote_indexes_dir
        String remote_chromosomes_dir
        String index_id
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
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
        
        # Processing alignments
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indexes_dir}/~{index_id}.xg . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading indexes. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        ${TIME_COMMAND} ${VG_COMMAND} pack --threads ${N_THREADS} --xg ~{index_id}.xg --gam ~{input_gam} --min-mapq 5 --trim-ends 5 --packs-out tmp.pack
        
        # Cleaning the VCF. This is needed, since the cleaned VCF (not the
        # original VCF) was used to build the indexes.
        gsutil -m cp ~{remote_chromosomes_dir}/'*.fa' .
        gunzip -c ~{input_vcf_gz} > input.vcf
        java -cp ~{docker_dir} CleanVCFGiraffe input.vcf . ~{max_sv_length} 1 input-vg.vcf
        bgzip input-vg.vcf
        tabix input-vg.vcf.gz
        
        # Calling
        ${TIME_COMMAND} ${VG_COMMAND} call --threads ${N_THREADS} --ploidy 2 --vcf input-vg.vcf.gz --pack tmp.pack ~{index_id}.xg > tmp.vcf
        rm -f tmp.pack
        bcftools sort --max-mem $(( ~{ram_size_gb} - 4 ))G --output-type z --output regenotyped.vcf.gz tmp.vcf
        tabix -f regenotyped.vcf.gz
    >>>
    
    output {
        File regenotyped_vcf_gz = work_dir + "/regenotyped.vcf.gz"
        File regenotyped_tbi = work_dir + "/regenotyped.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
