from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path=r'cellphonedb.zip',
    meta_file_path="meta_data.csv",
    counts_file_path="expression_matrix.csv",
    counts_data='hgnc_symbol',
    microenvs_file_path=None,
    score_interactions=True,
    output_path='Figure5',
    separator='|',
    threads=16,
    threshold=0.04,
    result_precision=9,
    iterations=1000,
    debug=None,
    output_suffix=None
)