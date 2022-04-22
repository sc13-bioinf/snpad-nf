#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def file_has_extension (it, extension)
{
	it.toString ().toLowerCase ().endsWith (extension.toLowerCase ())
}

// Return file if it exists
def file_from_path (it)
{
	if (it == null || it.isEmpty () ) exit 1, "[CRISPR] error: No value supplied for FASTQ or BAM input file. If using input method TSV set to NA if no file required. See '--help' flag and documentation under 'running the pipeline' for more information."
	// If glob == true file returns a list if it contains wildcard chars
	def f = file(it, glob: false)
	if (!f.exists()) exit 1, "[CRISPR] error: Cannot find supplied FASTQ or BAM input file. If using input method TSV set to NA if no file required. See '--help' flag and documentation under 'running the pipeline' for more information. Check file: ${it}"
	return f
}

def extract_data (data_file_path)
{
	Channel.fromPath (data_file_path, glob: false)
		.splitCsv (header: true, sep: '\t')
		.map { row ->

			[
				"sampleName": row.Sample_Name,
				"rg": row.RG.split (",").toList ().withIndex ().collect () { element, index -> [element, index * 31] },
				"bam": row.BAM,
			]
		}
}

def key_for_priors (sample_name, chromosome)
{
	if ( params.prior_chromosomes == null || chromosome in "${params.prior_chromosomes}".split (",") )
	{
		return [sample_name, chromosome].join ("_")
	}
	else
	{
		return [sample_name, "${params.prior_chromosomes}".split (",")?.find ()].join ("_")
	}
}

process indel_target {

	publishDir "${params.output_base}/${meta.sampleName}/real", mode: "copy"

	input:
		tuple val (meta), path (bam)
		each chromosome

	output:
		tuple val (meta), val (chromosome), path ("${meta.sampleName}.${chromosome}.intervals"), path ("${meta.sampleName}.${chromosome}.real.bam"), path ("${meta.sampleName}.${chromosome}.real.bai"), emit: result

	script:
	"""#!/usr/bin/env bash

java -jar /usr/GenomeAnalysisTK.jar \\
	-T RealignerTargetCreator \\
	-R ${params.reference} \\
	-L ${chromosome} \\
	--num_threads 1 \\
	-o ${meta.sampleName}.${chromosome}.intervals \\
	-I ${bam}

java -jar /usr/GenomeAnalysisTK.jar \\
	-T IndelRealigner \\
	-R ${params.reference} \\
	-L ${chromosome} \\
	-targetIntervals ${meta.sampleName}.${chromosome}.intervals \\
	--num_threads 1 \\
	-o ${meta.sampleName}.${chromosome}.real.bam \\
	-I ${bam} \\
	--downsample_to_coverage 1750

	"""

	stub:
	"""#!/usr/bin/env bash

cp ${params.stub_dir}/${meta.sampleName}/real/${meta.sampleName}.${chromosome}.intervals .
cp ${params.stub_dir}/${meta.sampleName}/real/${meta.sampleName}.${chromosome}.real.bam .
cp ${params.stub_dir}/${meta.sampleName}/real/${meta.sampleName}.${chromosome}.real.bai .

	"""
}

process bam_to_snpAD {
	tag "${meta.sampleName}.${chromosome}"

	input:
		tuple val (meta), val (chromosome), path (intervals), path (bam), path (bai), val (rg), val (rg_offset)

	output:
		tuple val (meta), val (chromosome), val (rg), val (rg_offset), path ("${meta.sampleName}.${chromosome}.bam.snpAD"), emit: result

	script:
	"""#!/usr/bin/env bash
echo "Running Bam2snpAD for '${chromosome}' and read group '${rg}' with offset '${rg_offset}'"
Bam2snpAD -r ${chromosome} -R ${rg} -f ${params.reference} -o ${rg_offset} -i ${bai} ${bam} > ${meta.sampleName}.${chromosome}.bam.snpAD
	"""

	stub:
	"""#!/usr/bin/env bash
touch ${meta.sampleName}.${chromosome}.bam.snpAD
	"""

}

process join_snpAD {
	tag "${meta.sampleName}.${chromosome}"

	publishDir "${params.output_base}/${meta.sampleName}/join", mode: "copy"

	input:
		tuple val (meta), val (chromosome), val (rg_groups), val (rg_offsets), path ("snpAD_rg.*")

	output:
		tuple val (meta), val (chromosome), path ("${meta.sampleName}.${chromosome}.merged.snpAD.gz"), emit: result

	script:
	"""#!/usr/bin/env bash
echo "Joining '${rg_offsets}' of groups '${rg_groups}' for '${chromosome}'"
snpADjoin snpAD_rg.* | bgzip -c > ${meta.sampleName}.${chromosome}.merged.snpAD.gz
	"""

	stub:
	"""#!/usr/bin/env bash
cp ${params.stub_dir}/${meta.sampleName}/join/${meta.sampleName}.${chromosome}.merged.snpAD.gz .
	"""
}

process run_snpAD {
	tag "${meta.sampleName}.${chromosome}"

	publishDir "${params.output_base}/${meta.sampleName}/priors", mode: "copy"

	input:
		tuple val (meta), val(chromosome), path (snpAD)

	output:
		tuple val (meta), val(chromosome), path ("${meta.sampleName}.${chromosome}.priors.txt"), path("${meta.sampleName}.${chromosome}.errors.txt"), emit: result

	when:
		params.prior_chromosomes == null || chromosome in params.prior_chromosomes.split (",")

	script:
	"""#!/usr/bin/env bash

echo "Running snpAD for '${chromosome}'"
snpAD -c ${params.snpad.threads} -o ${meta.sampleName}.${chromosome}.priors.txt -O ${meta.sampleName}.${chromosome}.errors.txt <(gunzip -c ${snpAD})
	"""

	stub:
	"""#!/usr/bin/env bash

cp ${params.stub_dir}/${meta.sampleName}/priors/${meta.sampleName}.${chromosome}.priors.txt .
cp ${params.stub_dir}/${meta.sampleName}/priors/${meta.sampleName}.${chromosome}.errors.txt .

	"""
}

process call_snpAD {
	tag "${meta.sampleName}.${chromosome}"

	input:
		tuple val (meta), val (chromosome), path (snpAD), val (prior_chromosome), path (priors), path (errors)

	output:
		tuple val (meta), val (chromosome), path ("${meta.sampleName}.${chromosome}.snpAD.part.vcf.gz"), path ("${meta.sampleName}.${chromosome}.snpAD.part.vcf.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env bash

prior_vals=\$(cat ${priors})
zero_prior=3.333333e-09
echo "Running snpAD for '${chromosome}' with priors '\${prior_vals}' from '${prior_chromosome}'"

snpADCall -N ${meta.sampleName} -p "\${prior_vals}" -e ${errors} --set_zero_priors \${zero_prior} <(gunzip -c ${snpAD}) | bgzip -c > ${meta.sampleName}.${chromosome}.snpAD.part.vcf.gz
tabix -p vcf ${meta.sampleName}.${chromosome}.snpAD.part.vcf.gz

"""

	stub:
	"""#!/usr/bin/env bash
touch ${meta.sampleName}.${chromosome}.snpAD.part.vcf.gz
touch ${meta.sampleName}.${chromosome}.snpAD.part.vcf.gz.tbi
	"""

}

process vcf_merge {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${meta.sampleName}/call", mode: "copy"

	input:
		tuple val (meta), path ("snpAD.part.*.vcf.gz"), path ("snpAD.part.*.vcf.gz.tbi")

	output:
		tuple val (meta), path ("${meta.sampleName}.snpAD.vcf.gz"), path ("${meta.sampleName}.snpAD.vcf.gz.tbi")

	script:
	"""#!/usr/bin/env bash

bcftools concat \\
	--output-type u \\
	snpAD.part.*.vcf.gz | bcftools view \\
	--output-type z \\
	--output ${meta.sampleName}.snpAD.vcf.gz \\
	--exclude 'N_ALT==0'

tabix -p vcf ${meta.sampleName}.snpAD.vcf.gz

	"""

	stub:
	"""#!/usr/bin/env bash
cp ${params.stub_dir}/${meta.sampleName}/call/${meta.sampleName}.snpAD.vcf.gz .
cp ${params.stub_dir}/${meta.sampleName}/call/${meta.sampleName}.snpAD.vcf.gz.tbi .
	"""
}

workflow {

	main:
		design_path = null

		if (params.prior_chromosomes == null) exit 1, "Please supply the chromosome to use for calculating priors (--prior_chromosomes)"
		if (params.input && file_has_extension (params.input, "tsv")) design_path = params.input else exit 1, "Please supply design file (--input)"
		ch_input = extract_data (design_path)
		ch_data = ch_input.map { [it, file (it["bam"],glob: false) ] }.dump (tag: 'ch_data')
		ch_chromosomes = "${params.chromosomes}".split (",").flatten ()

		indel_target (ch_data, ch_chromosomes)

		ch_real_rg = indel_target.out.result.dump (tag: 'indel_target.out.result').flatMap { it[0]["rg"].collect { jt -> it + jt } }.dump (tag: 'ch_real_rg')

		bam_to_snpAD (ch_real_rg)

		ch_rg_merge = bam_to_snpAD.out.result.map { [ [it[0]["sampleName"], it[1]].join ("_"), it[0], it[1], it[2], it[3], it[4] ] }
			.dump (tag: 'pre ch_rg_merge')
			.groupTuple ()
			.dump (tag: 'ch_rg_merge')
			.map { [ it[1][0], it[2][0], it[3].join (","), it[4].join (","), it[5] ] }

		join_snpAD (ch_rg_merge)
		run_snpAD (join_snpAD.out.result)

		// a.cross (b) b should have the repeated entries
		ch_call = run_snpAD.out.result.map { [ [it[0]["sampleName"], it[1] ].join ("_"), [it[1], it[2], it[3]]] }
			.cross (join_snpAD.out.result.map { [key_for_priors (it[0]["sampleName"], it[1]), it] })
			.dump (tag: 'pre ch_call')
			.map { it[1][1] + it[0][1] }
			.dump (tag: 'ch_call')

		call_snpAD (ch_call)

		ch_merge = call_snpAD.out.result.reduce ([null,[],[]]) { accumulator, item ->
			if ( accumulator[0] == null )
			{
				accumulator[0] = item[0]
			}
			accumulator[1].add (item[2])
			accumulator[2].add (item[3])
			accumulator
		}
		.dump (tag: 'ch_merge')

		vcf_merge (ch_merge)
}

