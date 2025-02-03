version 1.0

workflow AnnotateNewSnps {
    input {
        File control_vcf
        File control_vcf_index

        File? experiment_vcf
        File? experiment_vcf_index

        File? truth_vcf
        File? truth_vcf_index
        File? truth_bed

        Array[File] bams
        Array[File] bam_indices
        Array[String] bam_labels

        File ref_fasta
        File ref_fasta_index

        Array[File] bed_regions
        Array[String] bed_labels

        Array[String] extra_vcf_fields = []

        String experiment_name = "experiment"
        String control_name = "control"

        Boolean run_full_stats = true
    }

    if (defined(experiment_vcf) && defined(truth_vcf)) {
        call CalculateNewVariants {
            input:
                control_vcf = control_vcf,
                control_vcf_index = control_vcf_index,
                experiment_vcf = select_first([experiment_vcf]),
                experiment_vcf_index = select_first([experiment_vcf_index]),
                truth_vcf = select_first([truth_vcf]),
                truth_vcf_index = select_first([truth_vcf_index]),
                truth_bed = select_first([truth_bed]),
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index
        }

        if (run_full_stats) {
            call AnnotateSnps as ControlSnps_1 {
                input:
                    vcf = CalculateNewVariants.control_baseline_annotated_vcf,
                    vcf_index = CalculateNewVariants.control_baseline_annotated_vcf_index,
                    bams = bams,
                    bam_indices = bam_indices,
                    bam_labels = bam_labels,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    output_name = "annotated_control"
            }

            call AnnotateSnps as ExperimentSnps {
                input:
                    vcf = CalculateNewVariants.experiment_baseline_annotated_vcf,
                    vcf_index = CalculateNewVariants.experiment_baseline_annotated_vcf_index,
                    bams = bams,
                    bam_indices = bam_indices,
                    bam_labels = bam_labels,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    output_name = "annotated_experiment"
            }
        }

        call AnnotateSnps as NewExperimentSnps {
            input:
                vcf = CalculateNewVariants.new_experiment_variants_vcf,
                vcf_index = CalculateNewVariants.new_experiment_variants_vcf_index,
                bams = bams,
                bam_indices = bam_indices,
                bam_labels = bam_labels,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                output_name = "annotated_new_experiment"
        }

        call AnnotateSnps as MissingExperimentSnps {
            input:
                vcf = CalculateNewVariants.missing_experiment_variants_vcf,
                vcf_index = CalculateNewVariants.missing_experiment_variants_vcf_index,
                bams = bams,
                bam_indices = bam_indices,
                bam_labels = bam_labels,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                output_name = "annotated_missing_experiment"
        }

        if (run_full_stats) {
            call MakeTables as MakeControlSnpsTable_1 {
                input:
                    annotated_vcf = select_first([ControlSnps_1.annotated_vcf]),
                    annotated_vcf_index = select_first([ControlSnps_1.annotated_vcf_index]),
                    bam_labels = bam_labels,
                    bed_regions = bed_regions,
                    bed_labels = bed_labels,
                    experiment_name = control_name,
                    extra_vcf_fields = flatten([["BASE"], extra_vcf_fields])
            }

            call MakeTables as MakeExperimentSnpsTable {
                input:
                    annotated_vcf = select_first([ExperimentSnps.annotated_vcf]),
                    annotated_vcf_index = select_first([ExperimentSnps.annotated_vcf_index]),
                    bam_labels = bam_labels,
                    bed_regions = bed_regions,
                    bed_labels = bed_labels,
                    experiment_name = experiment_name,
                    extra_vcf_fields = flatten([["BASE"], extra_vcf_fields])
            }
        }

        call MakeTables as MakeNewExperimentSnpsTable {
            input:
                annotated_vcf = NewExperimentSnps.annotated_vcf,
                annotated_vcf_index = NewExperimentSnps.annotated_vcf_index,
                bam_labels = bam_labels,
                bed_regions = bed_regions,
                bed_labels = bed_labels,
                experiment_name = experiment_name,
                extra_vcf_fields = extra_vcf_fields
        }

        call MakeTables as MakeMissingExperimentSnpsTable {
            input:
                annotated_vcf = MissingExperimentSnps.annotated_vcf,
                annotated_vcf_index = MissingExperimentSnps.annotated_vcf_index,
                bam_labels = bam_labels,
                bed_regions = bed_regions,
                bed_labels = bed_labels,
                experiment_name = control_name,
                extra_vcf_fields = extra_vcf_fields
        }
    }

    if (!defined(experiment_vcf)) {
        call AnnotateSnps as ControlSnps_2 {
            input:
                vcf = control_vcf,
                vcf_index = control_vcf_index,
                bams = bams,
                bam_indices = bam_indices,
                bam_labels = bam_labels,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                output_name = "annotated_control"
        }

        call MakeTables as MakeControlSnpsTable_2 {
            input:
                annotated_vcf = ControlSnps_2.annotated_vcf,
                annotated_vcf_index = ControlSnps_2.annotated_vcf_index,
                bam_labels = bam_labels,
                bed_regions = bed_regions,
                bed_labels = bed_labels,
                experiment_name = control_name,
                extra_vcf_fields = extra_vcf_fields
        }
    }


    output {
        File? control_snps_table = select_first([MakeControlSnpsTable_1.output_table, MakeControlSnpsTable_2.output_table])
        File? experiment_snps_table = MakeExperimentSnpsTable.output_table
        File? new_experiment_snps_table = MakeNewExperimentSnpsTable.output_table
        File? missing_experiment_snps_table = MakeMissingExperimentSnpsTable.output_table
    }
}

task CalculateNewVariants {
    input {
        File control_vcf
        File control_vcf_index

        File experiment_vcf
        File experiment_vcf_index

        File truth_vcf
        File truth_vcf_index
        File truth_bed

        File ref_fasta
        File ref_fasta_index

        Int disk_size = 250
        Int cpu = 4
        Int memory = 8
    }

    command <<<
        set -xueo pipefail

        rtg format -o rtg_ref ~{ref_fasta}

        rtg vcfeval \
            -b ~{truth_vcf} \
            -c ~{control_vcf} \
            -e ~{truth_bed} \
            -t rtg_ref \
            --decompose \
            --output-mode annotate \
            -o eval_control

        bcftools view -i 'BASE="FN"||BASE="FN_CA"' eval_control/baseline.vcf.gz -o eval_control/fn.vcf.gz
        bcftools index -t eval_control/fn.vcf.gz
        bcftools view -i 'BASE="TP"' eval_control/baseline.vcf.gz -o eval_control/tp.vcf.gz
        bcftools index -t eval_control/tp.vcf.gz

        rtg vcfeval \
            -b ~{truth_vcf} \
            -c ~{experiment_vcf} \
            -e ~{truth_bed} \
            -t rtg_ref \
            --decompose \
            --output-mode annotate \
            -o eval_experiment

        bcftools view -i 'BASE="FN"||BASE="FN_CA"' eval_experiment/baseline.vcf.gz -o eval_experiment/fn.vcf.gz
        bcftools index -t eval_experiment/fn.vcf.gz
        bcftools view -i 'BASE="TP"' eval_experiment/baseline.vcf.gz -o eval_experiment/tp.vcf.gz
        bcftools index -t eval_experiment/tp.vcf.gz

        rtg vcfeval \
            -b eval_experiment/tp.vcf.gz \
            -c eval_control/fn.vcf.gz \
            -e ~{truth_bed} \
            -t rtg_ref \
            --decompose \
            -o new_experiment_variants

        rtg vcfeval \
            -b eval_control/tp.vcf.gz \
            -c eval_experiment/fn.vcf.gz \
            -e ~{truth_bed} \
            -t rtg_ref \
            --decompose \
            -o missing_experiment_variants
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.0"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + " GB"
    }

    output {
        File control_baseline_annotated_vcf = "eval_control/baseline.vcf.gz"
        File control_baseline_annotated_vcf_index = "eval_control/baseline.vcf.gz.tbi"
        File experiment_baseline_annotated_vcf = "eval_experiment/baseline.vcf.gz"
        File experiment_baseline_annotated_vcf_index = "eval_experiment/baseline.vcf.gz.tbi"
        File new_experiment_variants_vcf = "new_experiment_variants/tp.vcf.gz"
        File new_experiment_variants_vcf_index = "new_experiment_variants/tp.vcf.gz.tbi"
        File missing_experiment_variants_vcf = "missing_experiment_variants/tp.vcf.gz"
        File missing_experiment_variants_vcf_index = "missing_experiment_variants/tp.vcf.gz.tbi"
    }
}

task AnnotateSnps {
    input {
        File vcf
        File vcf_index

        Array[File] bams
        Array[File] bam_indices
        Array[String] bam_labels

        File ref_fasta
        File ref_fasta_index

        String output_name = "output"

        Int disk_size = 250
        Int cpu = 4
        Int memory = 8
    }

    command <<<
        set -xueo pipefail

        # Split MAs and subset to just SNP sites
        bcftools norm -m - ~{vcf} -o split_MA.vcf.gz -W tbi
        bcftools view -i 'TYPE="snp"' split_MA.vcf.gz -o split_snps.vcf.gz
        bcftools index -t split_snps.vcf.gz

        python3 << CODE
        import pysam
        import os

        input_vcf_file = "split_snps.vcf.gz"

        input_bam_files = ["~{sep="\", \"" bams}"]
        input_bam_labels = ["~{sep="\", \"" bam_labels}"]
        is_cram = [f.endswith('.cram') for f in input_bam_files]

        intermediate_input_vcfs = [input_vcf_file] + [f'{label}-annotated.vcf.gz' for label in input_bam_labels[:-1]]
        intermediate_output_vcfs = [f'{label}-annotated.vcf.gz' for label in input_bam_labels]

        # Check values of above lists
        print(input_bam_files)
        print(input_bam_labels)
        print(intermediate_input_vcfs)
        print(intermediate_output_vcfs)

        for cram_file, info_prefix, input_bam_file, input_vcf, output_vcf in zip(is_cram, input_bam_labels, input_bam_files, intermediate_input_vcfs, intermediate_output_vcfs):
            sam_read_mode = "rc" if cram_file else "rb"
            with pysam.VariantFile(input_vcf) as vcf:
                with pysam.AlignmentFile(input_bam_file, sam_read_mode, reference_filename="~{ref_fasta}") as bam:
                    counter = 0
                    vcf.header.add_meta('INFO', items=[('ID',f'{info_prefix}_SUPP_BQs'), ('Number', '.'), ('Type','Integer'), ('Description','ALT Allele supporting BQs')])
                    vcf.header.add_meta('INFO', items=[('ID', f'{info_prefix}_SUPP_MAPQs'), ('Number', '.'), ('Type','Integer'), ('Description','ALT Allele supporting MAPQs')])
                    vcf.header.add_meta('INFO', items=[('ID', f'{info_prefix}_REF_BQs'), ('Number', '.'), ('Type','Integer'), ('Description','REF Allele supporting BQs')])
                    vcf.header.add_meta('INFO', items=[('ID', f'{info_prefix}_REF_MAPQs'), ('Number', '.'), ('Type','Integer'), ('Description','REF Allele supporting MAPQs')])
                    vcf.header.add_meta('INFO', items=[('ID', f'{info_prefix}_OTHER_BQs'), ('Number', '.'), ('Type','Integer'), ('Description','Other Allele supporting BQs')])
                    vcf.header.add_meta('INFO', items=[('ID', f'{info_prefix}_OTHER_MAPQs'), ('Number', '.'), ('Type','Integer'), ('Description','Other Allele supporting MAPQs')])
                    with pysam.VariantFile(output_vcf, 'w', header=vcf.header) as vcf_out:
                        for rec in vcf:
                            supp_bqs = []
                            supp_mapqs = []
                            ref_bqs = []
                            ref_mapqs = []
                            other_bqs = []
                            other_mapqs = []

                            for read in bam.fetch(rec.chrom, rec.pos, rec.pos+1):
                                try:
                                    q_index = read.get_reference_positions().index(rec.pos-1)
                                    bq = read.query_qualities[q_index]
                                    base = read.query_sequence[q_index]
                                    mapq = read.mapq
                                    if base == rec.ref:
                                        ref_bqs += [bq]
                                        ref_mapqs += [mapq]
                                    elif base == rec.alleles[1]:
                                        supp_bqs += [bq]
                                        supp_mapqs += [mapq]
                                    else:
                                        other_bqs += [bq]
                                        other_mapqs += [mapq]
                                except ValueError:
                                    # Sometimes fetch gets reads that aren't actually overlapping event, so just skip these
                                    continue

                            rec.info[f'{info_prefix}_SUPP_BQs'] = supp_bqs if len(supp_bqs) > 0 else None
                            rec.info[f'{info_prefix}_SUPP_MAPQs'] = supp_mapqs if len(supp_mapqs) > 0 else None
                            rec.info[f'{info_prefix}_REF_BQs'] = ref_bqs if len(ref_bqs) > 0 else None
                            rec.info[f'{info_prefix}_REF_MAPQs'] = ref_mapqs if len(ref_mapqs) > 0 else None
                            rec.info[f'{info_prefix}_OTHER_BQs'] = other_bqs if len(other_bqs) > 0 else None
                            rec.info[f'{info_prefix}_OTHER_MAPQs'] = other_mapqs if len(other_mapqs) > 0 else None

                            vcf_out.write(rec)

        os.rename(intermediate_output_vcfs[-1], "~{output_name}.vcf.gz")
        CODE

        bcftools index -t "~{output_name}.vcf.gz"
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools_pysam:v1.0"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + " GB"
    }

    output {
        File annotated_vcf = "~{output_name}.vcf.gz"
        File annotated_vcf_index = "~{output_name}.vcf.gz.tbi"
    }
}

task MakeTables {
    input {
        File annotated_vcf
        File annotated_vcf_index

        Array[String] bam_labels

        Array[File] bed_regions
        Array[String] bed_labels

        Array[String] extra_vcf_fields = []

        String experiment_name = "experiment"

        Int disk_size = 250
        Int cpu = 4
        Int memory = 8
    }

    command <<<
        set -xueo pipefail

        echo "~{sep="\", \"" extra_vcf_fields}"

        cp ~{annotated_vcf} current.vcf.gz

        BED_REGIONS="~{write_lines(bed_regions)}"
        BED_LABELS="~{write_lines(bed_labels)}"
        paste $BED_LABELS $BED_REGIONS > bed_label_file_map.tsv

        # Add labels for bed regions to VCF
        while read line; do
            LABEL=$(echo "$line" | cut -f 1)
            BED=$(echo "$line" | cut -f 2)

            bgzip $BED
            tabix -s1 -b2 -e3 $BED.gz

            bcftools annotate -m "+${LABEL}" -a $BED.gz -c CHROM,POS current.vcf.gz -o tmp_out.vcf.gz
            mv tmp_out.vcf.gz current.vcf.gz
        done < bed_label_file_map.tsv

        # Format query string for bcftools to make table
        python3 << CODE
        query_string = "%CHROM\\t%POS\\t[%GT]"

        for label in ["~{sep="\", \"" bam_labels}"]:
            query_string += f"\\t%{label}_SUPP_BQs\\t%{label}_SUPP_MAPQs\\t%{label}_REF_BQs\\t%{label}_REF_MAPQs\\t%{label}_OTHER_BQs\\t%{label}_OTHER_MAPQs"

        for label in ["~{sep="\", \"" bed_labels}"]:
            query_string += f"\\t%{label}"

        if ~{length(extra_vcf_fields)} > 0:
            for label in ["~{sep="\", \"" extra_vcf_fields}"]:
                query_string += f"\\t%{label}"

        query_string += "\n"
        with open('query_string.txt', 'w') as f:
            f.write(query_string)

        CODE

        QUERY_STRING=$(cat query_string.txt)
        bcftools query -f"${QUERY_STRING}" current.vcf.gz > query_table.tsv

        python3 << CODE
        import pandas as pd

        names = ['CHROM', 'POS', 'GT']
        names += ['_'.join([label, x]) for label in ["~{sep="\", \"" bam_labels}"] for x in ["SUPP_BQs", "SUPP_MAPQs", "REF_BQs", "REF_MAPQs", "OTHER_BQs", "OTHER_MAPQs"]]
        names += ["_".join([label, "BED"]) for label in ["~{sep="\", \"" bed_labels}"]]
        names += ["~{sep="\", \"" extra_vcf_fields}"] if ~{length(extra_vcf_fields)} > 0 else []

        df = pd.read_csv('query_table.tsv', sep='\t', names=names)
        df['Experiment'] = "~{experiment_name}"
        df.to_csv('output_table.tsv.gz', sep='\t', index=False, compression='gzip')

        CODE

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + " GB"
    }

    output {
        File output_table = 'output_table.tsv.gz'
    }
}