process ANNOTATE_MATRIX {
  executor 'local'
  publishDir params.pubdir

  input:
  tuple val(record), path('*')
  path ref_df_paths

  output:
  path '*final*anndata*.h5ad'

  script:
  """
  annotate.py --ref_df_paths $ref_df_paths
  """
}