import numpy as np


def process_gene_expression(data):

    truncated_expr = data["truncatedExpression"]
    raw_expr = data["rawExpression"]

    result = {}

    for gene in raw_expr.keys():

        raw_gene_expr = raw_expr[gene]
        # can't do summary stats on list with None's
        # casting to np array replaces None with np.nan
        np_gene_expr = np.array(raw_gene_expr, dtype=np.float)
        # This is not necessary and is also costly, but I leave it commented as a reminder
        # that this object has integer zeros and floating point for n!=0.
        mean = float(np.nanmean(np_gene_expr))
        stdev = float(np.nanstd(np_gene_expr))

        truncated_gene_expr = truncated_expr[gene]
        np_truncated_gene_expr = np.array(truncated_gene_expr, dtype=np.float)

        minimum = float(np.nanmin(np_truncated_gene_expr))
        maximum = float(np.nanmax(np_truncated_gene_expr))

        result[gene] = {
            "rawExpression": {
                "mean": mean,
                "stdev": stdev,
                "expression": raw_gene_expr,
            },
            "truncatedExpression": {
                "min": minimum,
                "max": maximum,
                "expression": truncated_gene_expr,
            },
        }

    return result
