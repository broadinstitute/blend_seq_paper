version 1.0

workflow BenchmarkPhasing {
    input {
        File base_vcf
        File base_vcf_index
        String? base_sample_name

        File query_vcf
        File query_vcf_index
        String? call_sample_name

        File? eval_regions
        File ref_fai

        File? genes_bed
        File? hiconf_bed
        File? vcfdist_summary

        String? experiment
    }

    call WhatsHapCompare {
        input:
            base_vcf=base_vcf,
            base_vcf_index=base_vcf_index,
            base_sample_name=base_sample_name,
            query_vcf=query_vcf,
            query_vcf_index=query_vcf_index,
            call_sample_name=call_sample_name,
            experiment=experiment,
            eval_regions=eval_regions,
            ref_fai=ref_fai
    }

    if (defined(genes_bed) && defined(vcfdist_summary)) {
        call CountPhasedGenes {
            input:
                genes_bed=select_first([genes_bed]),
                phase_blocks_bed=WhatsHapCompare.phase_blocks,
                switch_errors_bed=WhatsHapCompare.switch_errors_bed,
                query_vcf=query_vcf,
                query_vcf_index=query_vcf_index,
                hiconf_bed=hiconf_bed,
                vcfdist_summary=select_first([vcfdist_summary]),
                experiment=experiment,
                ref_fai=ref_fai
        }
    }


    output {
        File whatshap_eval = WhatsHapCompare.eval
        File switch_errors_bed = WhatsHapCompare.switch_errors_bed
        File phase_blocks = WhatsHapCompare.phase_blocks
        File stats = WhatsHapCompare.stats

        File? all_genes = CountPhasedGenes.all_genes
        File? phased_genes = CountPhasedGenes.phased_genes
        File? correctly_phased_genes = CountPhasedGenes.correctly_phased_genes
        File? perfect_genes = CountPhasedGenes.perfect_genes
        File? gene_summary_counts = CountPhasedGenes.summary_counts
    }
}

task WhatsHapCompare {
    input {
        File base_vcf
        File base_vcf_index
        String? base_sample_name

        File query_vcf
        File query_vcf_index
        String? call_sample_name

        File? eval_regions
        File ref_fai

        String? experiment

        Boolean only_snvs = false    # Toggle to only consider phasing for SNVs
        Boolean filter_SVs = true    # Toggle to filter out SVs from evaluation at start
        Boolean remove_MA_hets = true   # Avoid whatshap bug with het MAs

        # Runtime arguments
        Int disk_size = 2 * ceil(size(base_vcf, "GB") + size(query_vcf, "GB")) + 50
        Int cpu = 4
        Int memory_ram = 16
    }

    String het_ma_flag = if (remove_MA_hets) then "--max-alleles 2" else ""

    command <<<
        set -xueo pipefail

        bcftools +fixploidy ~{base_vcf} -- -f 2 > fixed_ploidy-base.vcf
        bcftools view fixed_ploidy-base.vcf ~{het_ma_flag} ~{"-T " + eval_regions} -o fixed_ploidy-base.vcf.gz
        bcftools norm -m +any fixed_ploidy-base.vcf.gz -o fixed_ploidy-base-norm.vcf.gz
        bcftools index -t fixed_ploidy-base-norm.vcf.gz

        bcftools +fixploidy ~{query_vcf} -- -f 2 > fixed_ploidy-query.vcf
        bcftools view fixed_ploidy-query.vcf ~{het_ma_flag} ~{"-T " + eval_regions} -o fixed_ploidy-query.vcf.gz
        bcftools norm -m +any fixed_ploidy-query.vcf.gz -o fixed_ploidy-query-norm.vcf.gz
        bcftools index -t fixed_ploidy-query-norm.vcf.gz

        # Filter out SVs from evaluation if specified
        if [ ~{filter_SVs} == true ]; then
            bcftools view -i '(ILEN > -50 && ILEN < 50) || TYPE="snp"' fixed_ploidy-base-norm.vcf.gz -o fixed_ploidy-base-norm-filtered.vcf.gz -Wtbi
            bcftools view -i '(ILEN > -50 && ILEN < 50) || TYPE="snp"' fixed_ploidy-query-norm.vcf.gz -o fixed_ploidy-query-norm-filtered.vcf.gz -Wtbi
        else
            mv fixed_ploidy-base-norm.vcf.gz fixed_ploidy-base-norm-filtered.vcf.gz
            bcftools index -t fixed_ploidy-base-norm-filtered.vcf.gz
            mv fixed_ploidy-query-norm.vcf.gz fixed_ploidy-query-norm-filtered.vcf.gz
            bcftools index -t fixed_ploidy-query-norm-filtered.vcf.gz
        fi

        whatshap compare \
            --names ~{default="truth" base_sample_name},~{default="callset" call_sample_name} \
            --tsv-pairwise eval.tsv \
            --ignore-sample-name \
            ~{true="--only-snvs" false="" only_snvs} \
            --switch-error-bed switch_errors.bed \
            --longest-block-tsv longest_blocks.tsv \
            fixed_ploidy-base-norm-filtered.vcf.gz fixed_ploidy-query-norm-filtered.vcf.gz

        whatshap stats fixed_ploidy-query-norm-filtered.vcf.gz \
            --tsv stats.tsv \
            --block-list phase_blocks.tsv

        # Clean and sort switch errors bed file
        awk -v OFS="\t" '{print $1, $2}' ~{ref_fai} > ref.genome
        bedtools sort -i switch_errors.bed -g ref.genome > switch_errors_sorted.bed
        awk -v OFS="\t" '{print $1, $2, $3}' switch_errors_sorted.bed > switch_errors.bed

        echo -e "#CHROM\tSTART\tEND\tNUM_VAR" > phase_blocks.bed
        awk -v OFS="\t" 'NR>1 {print $2, $4-1, $5, $6}' phase_blocks.tsv >> phase_blocks.bed

        # Slice phase blocks at breakpoints in switch_errors.bed
        bedtools intersect -a phase_blocks.bed -b switch_errors.bed > phase_blocks_corrected-1.bed
        bedtools complement -i switch_errors.bed -g ref.genome > switch_errors_complement.bed
        bedtools intersect -a phase_blocks.bed -b switch_errors_complement.bed > phase_blocks_corrected-2.bed
        cat phase_blocks_corrected-1.bed phase_blocks_corrected-2.bed | bedtools sort -i - -g ref.genome > phase_blocks_corrected.bed

        # Add summary statistics for eval ALL regions
        python3 << CODE
        import pandas as pd

        df = pd.read_csv('eval.tsv', sep='\t')
        sample = df['#sample'].values[0]
        dataset_name0 = df['dataset_name0'].values[0]
        dataset_name1 = df['dataset_name1'].values[0]
        file_name0 = df['file_name0'].values[0]
        file_name1 = df['file_name1'].values[0]
        intersection_blocks = df['intersection_blocks'].sum()
        covered_variants = df['covered_variants'].sum()
        all_assessed_pairs = df['all_assessed_pairs'].sum()
        all_switches = df['all_switches'].sum()
        all_switch_rate = all_switches / all_assessed_pairs
        all_s = df['all_switchflips'].apply(lambda s: s.split('/')[0]).astype(int).sum()
        all_f = df['all_switchflips'].apply(lambda s: s.split('/')[1]).astype(int).sum()
        all_switchflips = f'{all_s}/{all_f}'
        all_switchflip_rate = (all_s + all_f) / all_assessed_pairs
        blockwise_hamming = df['blockwise_hamming'].sum()
        blockwise_hamming_rate = blockwise_hamming / covered_variants
        blockwise_diff_genotypes = df['blockwise_diff_genotypes'].sum()
        blockwise_diff_genotypes_rate = blockwise_diff_genotypes / covered_variants
        largestblock_assessed_pairs = df['largestblock_assessed_pairs'].sum()
        largestblock_switches = df['largestblock_switches'].sum()
        largestblock_switch_rate = largestblock_switches / largestblock_assessed_pairs
        largestblock_s = df['largestblock_switchflips'].apply(lambda s: s.split('/')[0]).astype(int).sum()
        largestblock_f = df['largestblock_switchflips'].apply(lambda s: s.split('/')[1]).astype(int).sum()
        largestblock_switchflips = f'{largestblock_s}/{largestblock_f}'
        largestblock_switchflip_rate = (largestblock_s + largestblock_f) / largestblock_assessed_pairs
        largestblock_hamming = df['largestblock_hamming'].sum()
        largestblock_hamming_rate = largestblock_hamming / covered_variants
        largestblock_diff_genotypes = df['largestblock_diff_genotypes'].sum()
        largestblock_diff_genotypes_rate = largestblock_diff_genotypes / covered_variants
        het_variants0 = df['het_variants0'].sum()
        only_snvs = df['only_snvs'].values[0]

        def ng_50(vals, total):
            sorted_vals = sorted(vals, reverse=True)
            running_sum = 0
            for i in sorted_vals:
                running_sum += i
                if running_sum >= total/2:
                    return i
            return 0

        phased_blocks_corrected = pd.read_csv('phase_blocks_corrected.bed', sep='\t', names=['CHROM', 'START', 'END', 'NUM_VAR'])
        phased_blocks_corrected['BLOCK_LENGTH'] = phased_blocks_corrected['END'] - phased_blocks_corrected['START']
        ref_genome = pd.read_csv('ref.genome', sep='\t', names=['CHROM', 'LENGTH'])
        ref_length_dict = dict(zip(ref_genome['CHROM'], ref_genome['LENGTH']))
        ngc_50 = phased_blocks_corrected.groupby('CHROM').apply(lambda df: ng_50(df['BLOCK_LENGTH'], ref_length_dict[df['CHROM'].values[0]])).reset_index()
        ngc_50 = pd.concat([
            ngc_50,
            pd.DataFrame({
                'CHROM': ['ALL'],
                0: [ng_50(phased_blocks_corrected['BLOCK_LENGTH'], ref_genome['LENGTH'].sum())]
            })
        ])
        ngc_50.columns = ['chromosome', 'NGC_50']

        final_df = pd.concat([
            df,
            pd.DataFrame({
                '#sample': [sample],
                'chromosome': ['ALL'],
                'dataset_name0': [dataset_name0],
                'dataset_name1': [dataset_name1],
                'file_name0': [file_name0],
                'file_name1': [file_name1],
                'intersection_blocks': [intersection_blocks],
                'covered_variants': [covered_variants],
                'all_assessed_pairs': [all_assessed_pairs],
                'all_switches': [all_switches],
                'all_switch_rate': [all_switch_rate],
                'all_switchflips': [all_switchflips],
                'all_switchflip_rate': [all_switchflip_rate],
                'blockwise_hamming': [blockwise_hamming],
                'blockwise_hamming_rate': [blockwise_hamming_rate],
                'blockwise_diff_genotypes': [blockwise_diff_genotypes],
                'blockwise_diff_genotypes_rate': [blockwise_diff_genotypes_rate],
                'largestblock_assessed_pairs': [largestblock_assessed_pairs],
                'largestblock_switches': [largestblock_switches],
                'largestblock_switch_rate': [largestblock_switch_rate],
                'largestblock_switchflips': [largestblock_switchflips],
                'largestblock_switchflip_rate': [largestblock_switchflip_rate],
                'largestblock_hamming': [largestblock_hamming],
                'largestblock_hamming_rate': [largestblock_hamming_rate],
                'largestblock_diff_genotypes': [largestblock_diff_genotypes],
                'largestblock_diff_genotypes_rate': [largestblock_diff_genotypes_rate],
                'het_variants0': [het_variants0],
                'only_snvs': [only_snvs]
            })
        ]).merge(ngc_50, left_on='chromosome', right_on='chromosome')

        final_df['Experiment'] = "~{experiment}"
        final_df.to_csv('eval.tsv', sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/whatshap:v1.2"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory_ram + " GB"
    }

    output {
        File eval = "eval.tsv"
        File switch_errors_bed = "switch_errors.bed"
        File phase_blocks = "phase_blocks.bed"
        File stats = "stats.tsv"
    }
}

task CountPhasedGenes {
    input {
        File genes_bed
        File phase_blocks_bed
        File switch_errors_bed

        File query_vcf
        File query_vcf_index

        File? hiconf_bed

        File ref_fai

        File vcfdist_summary

        String experiment = "experiment"
    }

    command <<<
        set -xueo pipefail

        # Script for generating hiconf genes bed
        cat << EOF > hiconf_genes.py
        import pandas as pd
        import numpy as np

        df = pd.read_csv('genes_hiconf.bed', sep='\t', names=['CHROM', 'START', 'END', 'gene', 'HICONF_CHROM', 'HICONF_START', 'HICONF_END'])
        df['gene_id'] = df['gene'].apply(lambda s: s.split(';')[0].split('=')[-1])

        num_hiconf_per_genes = df[df['HICONF_CHROM'] != '.']['gene_id'].value_counts()
        hiconf_genes = num_hiconf_per_genes[num_hiconf_per_genes == 1].index
        hiconf_genes_df = df[df['gene_id'].isin(hiconf_genes)]

        # Need to ensure single high conf region overlapping actually spans full gene
        hiconf_genes_df = hiconf_genes_df[(hiconf_genes_df['HICONF_START'] <= hiconf_genes_df['START']) & (hiconf_genes_df['HICONF_END'] >= hiconf_genes_df['END'])]
        hiconf_genes_df[['CHROM', 'START', 'END', 'gene']].to_csv('genes_hiconf.bed', sep='\t', index=False, header=None)

        EOF

        # Subset to genes fully contained in high conf regions for correctly phased/perfect assemblies
        if [ -n "~{hiconf_bed}" ]; then
            bedtools intersect -a ~{genes_bed} -b ~{hiconf_bed} -loj > genes_hiconf.bed
            python3 hiconf_genes.py
        else
            cp ~{genes_bed} genes_hiconf.bed
        fi

        # Pull out variants and PS values overlapping genes
        bcftools view -T ~{genes_bed} -f .,PASS -g het ~{query_vcf} | bcftools query -f '%CHROM\t%POS0\t%END0\t[%PS]' > gene_variants-PS.tsv
        # bedtools intersect -a ~{genes_bed} -b ~{phase_blocks_bed} -loj > genes_phase_blocks.bed
        bedtools intersect -a ~{genes_bed} -b gene_variants-PS.tsv -loj > genes_phase_blocks.bed
        bedtools intersect -a genes_hiconf.bed -b ~{switch_errors_bed} -loj > genes_switch_errors.bed
        bedtools intersect -a genes_hiconf.bed -b ~{vcfdist_summary} -loj > genes_performance.bed

        python3 << CODE
        import pandas as pd
        import numpy as np

        df = pd.read_csv('genes_phase_blocks.bed', sep='\t', names=['CHROM', 'START', 'END', 'gene', 'VAR_CHROM', 'VAR_START', 'VAR_END', 'VAR_PS'])
        df['gene_id'] = df['gene'].apply(lambda s: s.split(';')[0].split('=')[-1])

        # Check a gene that all het variants are phased
        def check_gene(d):
            ps_labels = d['VAR_PS'].unique()
            gene = d['gene_id'].values[0]
            return pd.DataFrame({
                'phased': [(len(ps_labels) == 1) and (list(ps_labels) != ['.'])]
            })

        phased_genes_results = df.groupby('gene_id').apply(check_gene).reset_index()
        phased_genes_df = phased_genes_results[phased_genes_results['phased']]
        # phased_genes_df = df[df['gene_id'].isin(phased_genes)]

        # Get genes with no switch errors
        switch_errors_df = pd.read_csv('genes_switch_errors.bed', sep='\t', names=['CHROM', 'START', 'END', 'gene', 'SWITCH_CHROM', 'SWITCH_START', 'SWITCH_END'])
        switch_errors_df['gene_id'] = switch_errors_df['gene'].apply(lambda s: s.split(';')[0].split('=')[-1])

        # Import list of hiconf genes
        hiconf_genes_df = pd.read_csv('genes_hiconf.bed', sep='\t', names=['CHROM', 'START', 'END', 'gene'])
        hiconf_genes_df['gene_id'] = hiconf_genes_df['gene'].apply(lambda s: s.split(';')[0].split('=')[-1])
        hiconf_genes = hiconf_genes_df['gene_id'].unique()

        # Create mask for genes containing a full switch error (if no overlaps, then 'END' == -1 so right condition never satisfied)
        contains_switch_mask = (switch_errors_df['SWITCH_START'] >= switch_errors_df['START']) & (switch_errors_df['SWITCH_END'] <= switch_errors_df['END'])
        contains_switches = switch_errors_df[contains_switch_mask]['gene_id'].unique()
        hiconf_phased_genes_df = phased_genes_df[phased_genes_df['gene_id'].isin(hiconf_genes)]
        correctly_phased_genes_df = hiconf_phased_genes_df[~(hiconf_phased_genes_df['gene_id'].isin(contains_switches))]

        ## Get genes with perfect assembly
        # Read in data joining gene list with performance data (TP/FP/FN labels from vcfdist)
        gene_perf = pd.read_csv('genes_performance.bed', sep='\t', names=['CHROM', 'START', 'END', 'gene', 'VAR_CHROM',
                                'VAR_START', 'VAR_ID', 'VAR_REF', 'VAR_ALT', 'VAR_Q', 'VAR_FILTER', 'VAR_INFO', 'VAR_FMT', 'TRUTH', 'QUERY'])
        gene_perf['gene_id'] = gene_perf['gene'].apply(lambda s: s.split(';')[0].split('=')[-1])

        def get_label(row):
            fmt = row['VAR_FMT']
            try:
                idx = fmt.split(':').index('BD')
                return row['QUERY'].split(':')[idx]
            except ValueError:
                return np.nan

        # Subset down to correctly phased genes
        gene_perf_correctly_phased = gene_perf[gene_perf['gene_id'].isin(correctly_phased_genes_df['gene_id'].unique())]
        get_labels = gene_perf_correctly_phased.apply(get_label, axis=1)
        gene_perf_correctly_phased['QUERY_Label'] = get_labels if len(get_labels) > 0 else np.nan

        def check_perfect_perf(sub_df):
            # Either all variants are TPs or there are none (i.e. matches ref)
            return (list(sub_df['QUERY_Label'].unique()) == ['TP']) or (len(list(sub_df['QUERY_Label'].unique())) == 0)

        # Subset to perfect genes where all variants in gene are TPs (none missed, e.g. BD = '.') and correctly fully phased
        get_perfect_perf = gene_perf_correctly_phased.groupby('gene_id').apply(check_perfect_perf)
        perfect_genes = get_perfect_perf.reset_index(name='reset_gene_id') if isinstance(get_perfect_perf, pd.Series) else get_perfect_perf.reset_index(names='reset_gene_id')
        perfect_genes = [] if len(perfect_genes) == 0 else list(perfect_genes[perfect_genes['reset_gene_id']]['gene_id'])
        perfect_genes_df = correctly_phased_genes_df[correctly_phased_genes_df['gene_id'].isin(perfect_genes)]

        # Count number of hiconf genes
        with open('genes_hiconf.bed') as f:
            hiconf_genes = len(f.readlines())

        summary_df = pd.DataFrame({
            'Total Genes': [len(df['gene_id'].unique())],
            'Total Hiconf Genes': [hiconf_genes],
            'Phased Genes': [len(phased_genes_df['gene_id'].unique())],
            'Correctly Phased Genes': [len(correctly_phased_genes_df['gene_id'].unique())],
            'Perfect Genes': [len(perfect_genes_df['gene_id'].unique())]
        })

        df['Experiment'] = "~{experiment}"
        phased_genes_df['Experiment'] = "~{experiment}"
        correctly_phased_genes_df['Experiment'] = "~{experiment}"
        perfect_genes_df['Experiment'] = "~{experiment}"
        summary_df['Experiment'] = "~{experiment}"

        df.to_csv('all_genes.tsv', sep='\t', index=False)
        phased_genes_df.to_csv('phased_genes.tsv', sep='\t', index=False)
        correctly_phased_genes_df.to_csv('correctly_phased_genes.tsv', sep='\t', index=False)
        perfect_genes_df.to_csv('perfect_genes.tsv', sep='\t', index=False)
        summary_df.to_csv('summary_counts.tsv', sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.1"
        disks: "local-disk 50 HDD"
        memory: "16 GB"
        cpu: 8
    }

    output {
        File all_genes = "all_genes.tsv"
        File phased_genes = "phased_genes.tsv"
        File correctly_phased_genes = "correctly_phased_genes.tsv"
        File perfect_genes = "perfect_genes.tsv"
        File summary_counts = "summary_counts.tsv"
    }
}