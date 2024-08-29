version 1.0

# Workflow for calling structural variants with vg
# Can take single-ended or paired-end reads
# Allows for specific custom variants to augment graph, or to deduce them from all INDELs in a BAM
workflow VgCallSVs {
    input {
        File? augmentation_bam
        File? augmentation_bam_index
        File? trivial_caller_class

        File? custom_augmentation_calls
        File? custom_augmentation_calls_index

        File vcf_header

        File ref_fasta
        File ref_index
        File ref_dict

        File reads_fastq_1
        File? reads_fastq_2

        String ref_name = "extended_ref"
        String output_name = "aligned_reads"

        File? monitoring_script
    }

    if (defined(augmentation_bam) && !defined(custom_augmentation_calls)) {
        call CreateNodes {
            input:
                bam=select_first([augmentation_bam]),
                bam_index=select_first([augmentation_bam_index]),
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                vcf_header=vcf_header,
                trivial_caller_class=select_first([trivial_caller_class]),
                monitoring_script=monitoring_script
        }
    }

    call MakeGraphReference {
        input:
            nodes_vcf=select_first([custom_augmentation_calls, CreateNodes.trivial_calls]),
            nodes_vcf_index=select_first([custom_augmentation_calls_index, CreateNodes.trivial_calls_index]),
            ref_fasta=ref_fasta,
            ref_index=ref_index,
            output_name=ref_name,
            monitoring_script=monitoring_script
    }

    call SplitReads {
        input:
            reads_1=reads_fastq_1,
            reads_2=reads_fastq_2,
            monitoring_script=monitoring_script
    }

    Array[File] reads_1_array = SplitReads.reads_part_1
    Array[File] reads_2_array = SplitReads.reads_part_2

    if (length(reads_1_array) == length(reads_2_array)) {
        scatter(paired_reads in zip(reads_1_array, reads_2_array)) {
            call AlignToGraphGAM as AlignToGAMPaired {
                input:
                    ref_graph=MakeGraphReference.ref_graph,
                    ref_dist=MakeGraphReference.ref_dist,
                    ref_min=MakeGraphReference.ref_min,
                    reads_fastq=paired_reads.left,
                    read_pairs_fastq=paired_reads.right,
                    output_name=basename(basename(paired_reads.left, ".fastq.gz"), ".fastq"),
                    monitoring_script=monitoring_script
            }
        }
    }

    if (length(reads_1_array) > length(reads_2_array)) {
        scatter(reads in reads_1_array) {
            call AlignToGraphGAM as AlignToGraphSingleEnd {
                input:
                    ref_graph=MakeGraphReference.ref_graph,
                    ref_dist=MakeGraphReference.ref_dist,
                    ref_min=MakeGraphReference.ref_min,
                    reads_fastq=reads,
                    output_name=basename(basename(reads, ".fastq.gz"), ".fastq"),
                    monitoring_script=monitoring_script
            }
        }
    }

    call CallGraphVariants {
        input:
            aln_gams=select_first([AlignToGAMPaired.aligned_graph_reads, AlignToGraphSingleEnd.aligned_graph_reads]),
            ref_graph=MakeGraphReference.ref_graph,
            monitoring_script=monitoring_script
    }

    output {
        File? graph_calls = CallGraphVariants.graph_calls
        File? graph_calls_index = CallGraphVariants.graph_calls_index
    }
}

task CreateNodes {
    input {
        File bam
        File bam_index

        File ref_fasta
        File ref_index

        File vcf_header
        File trivial_caller_class

        File? monitoring_script
    }

    command <<<
        set -xueo pipefail

        if [ -n "~{monitoring_script}" ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        cat ~{ref_fasta} | awk '{ if (substr($0,1,1)==">") { filename=(substr($1,2) ".fa") } print $0 >> filename; close(filename) }'

        samtools view ~{bam} > alignments.sam

        cat ~{vcf_header} > trivial_calls.vcf

        cp ~{trivial_caller_class} .
        java -Xmx4G ~{basename(trivial_caller_class, ".class")} alignments.sam . >> trivial_calls.vcf

        bcftools sort trivial_calls.vcf -o trivial_calls-dup.vcf.gz
        bcftools index -t trivial_calls-dup.vcf.gz

        bcftools norm -d all trivial_calls-dup.vcf.gz -o trivial_calls.vcf.gz
        bcftools index -t trivial_calls.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/trivial_caller:v1.0"
        cpu: 4
        memory: 8
        disks: "local-disk " + ceil(2*size(bam, "GB") + 200) + " HDD"
    }

    output {
        File trivial_calls = "trivial_calls.vcf.gz"
        File trivial_calls_index = "trivial_calls.vcf.gz.tbi"
        File trivial_calls_with_dup = "trivial_calls-dup.vcf.gz"
        File trivial_calls_with_dup_index = "trivial_calls-dup.vcf.gz.tbi"

        File? monitoring_log = "monitoring.log"
    }
}

task MakeGraphReference {
    input {
        File nodes_vcf
        File nodes_vcf_index

        File ref_fasta
        File ref_index

        String output_name = "extended_ref"

        Int cpu = 32
        Int mem = 256
        File? monitoring_script
    }

    command <<<
        set -xueo pipefail

        if [ -n "~{monitoring_script}" ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        vg autoindex \
            -t $(nproc) \
            --workflow giraffe \
            -r ~{ref_fasta} \
            -v ~{nodes_vcf} \
            -p ~{output_name}
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vg_docker:v1.0"
        cpu: cpu
        memory: mem + "GB"
        disks: "local-disk 200 HDD"
    }

    output {
        File ref_graph = "~{output_name}.giraffe.gbz"
        File ref_dist = "~{output_name}.dist"
        File ref_min = "~{output_name}.min"

        File? monitoring_log = "monitoring.log"
    }
}

task SplitReads {
    input {
        File reads_1
        File? reads_2
        Int num_shards = 50

        File? monitoring_script

        Int extra_disk = 500
        Int cpu = 8
        Int mem = 16
    }

    String name_1 = basename(basename(reads_1, ".fastq.gz"), ".fastq")
    String name_2 = if (defined(reads_2)) then basename(basename(select_first([reads_2]), ".fastq.gz"), ".fastq") else "reads2"

    Int disk = ceil(5*size(reads_1, "GB") + 250) + extra_disk

    command <<<
        set -xueo pipefail

        if [ -n "~{monitoring_script}" ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        seqkit split2 \
            -1 ~{reads_1} \
            ~{"-2 " + reads_2} \
            -p ~{num_shards} \
            -O output
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/seqkit:v1.1"
        cpu: cpu
        memory: mem + "GB"
        disks: "local-disk " + disk + " HDD"
    }

    output {
        Array[File] reads_part_1 = glob("output/~{name_1}.part_*")
        Array[File] reads_part_2 = glob("output/~{name_2}.part_*")

        File? monitoring_log = "monitoring.log"
    }
}

task AlignToGraphGAM {
    input {
        File ref_graph
        File ref_dist
        File ref_min

        File reads_fastq
        File? read_pairs_fastq

        String output_name = "aligned_reads"

        File? monitoring_script

        Int cpu = 16
        Int memory = 120
    }

    Int disk_size = ceil(3 * size(ref_graph, "GB") + 5 * size(reads_fastq, "GB") + 300)

    command <<<
        set -xueo pipefail

        if [ -n "~{monitoring_script}" ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        vg giraffe \
            -t $(nproc) \
            -Z ~{ref_graph} \
            -f ~{reads_fastq} \
            ~{"-f " + read_pairs_fastq} \
            -o gam > "~{output_name}.gam"

        vg stats -a "~{output_name}.gam" > "~{output_name}.vg_stats"
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vg_docker:v1.0"
        cpu: cpu
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File aligned_graph_reads = "~{output_name}.gam"
        File vg_stats = "~{output_name}.vg_stats"

        File? monitoring_log = "monitoring.log"
    }
}

task CallGraphVariants {
    input {
        Array[File] aln_gams
        File ref_graph
        Int qual_threshold = 5

        String output_name = "graph_calls"
        Int cpu = 16
        Int memory = 64
        File? monitoring_script
    }

    Int disk_size = ceil(2.5 * size(aln_gams, "GB") + size(ref_graph, "GB")) + 100

    command <<<
        set -xueo pipefail

        if [ -n "~{monitoring_script}" ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        vg index ~{ref_graph} -x indexed_graph.xg -L
        REF_GRAPH="indexed_graph.xg"

        cat ~{sep=" " aln_gams} | vg pack -x $REF_GRAPH -g - -Q ~{qual_threshold} -o aln.pack

        vg call $REF_GRAPH -k aln.pack > graph_calls.vcf

        bcftools view graph_calls.vcf -o ~{output_name}-presort.vcf.gz
        bcftools sort ~{output_name}-presort.vcf.gz -o ~{output_name}.vcf.gz
        bcftools index -t ~{output_name}.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vg_docker:v1.0"
        cpu: cpu
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File graph_calls = "~{output_name}.vcf.gz"
        File graph_calls_index = "~{output_name}.vcf.gz.tbi"

        File? monitoring_log = "monitoring.log"
    }
}