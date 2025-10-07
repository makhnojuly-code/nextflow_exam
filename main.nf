nextflow.enable.dsl=2

params.accession = 'M21012'
params.samples   = "${projectDir}/data/hepatitis_combined.fasta"
params.outdir    = "${projectDir}/results"

process download_ref {
  publishDir "${params.outdir}/reference", mode: 'copy'
  input:
    val acc
  output:
    path "${acc}.fasta", emit: ref_fasta
  script:
  """
  wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" \
       -O ${acc}.fasta
  """
}

process combine_fasta {
  publishDir "${params.outdir}/combined", mode: 'copy', overwrite: true
  input:
    path ref
    path samples
  output:
    path "combined_all.fasta", emit: combined
  script:
  """
  cat "${ref}" "${samples}" > combined_all.fasta
  """
}

process mafft_align {
  publishDir "${params.outdir}/align", mode: 'copy', overwrite: true
  container "docker://staphb/mafft:7.526"
    input:
        path combined
    output:
        path "aligned.fasta", emit: aligned
    script:
    """
    mafft --auto "${combined}" > aligned.fasta
    """

}

process trimal_clean {
    publishDir "${params.outdir}/final", mode: 'copy'
    container "https://depot.galaxyproject.org/singularity/trimal:1.4.1--h4ac6f70_9"
   
    input:
        path aligned
    output:
        path "alignment_trimmed.fasta"
        path "alignment_trimmed.html"
    script:
    """
    trimal -in "${aligned}" -out alignment_trimmed.fasta -automated1 -htmlout alignment_trimmed.html
    """
}

workflow {
  
  ref_ch = download_ref(params.accession).ref_fasta

  samples_ch = Channel
                .fromPath(params.samples)
                .ifEmpty { exit 1, "[ERROR] Samples not found: ${params.samples}" }

  comb = combine_fasta(ref_ch, samples_ch).combined

  aln = mafft_align(comb).aligned

  trimal_clean(aln)
}

