from pandasql import sqldf
from result import Result
import json


class ListGenes:
    def __init__(self, adata):
        self.adata = adata

    def _format_result(self, result, no_genes):

        # JSONify result.
        result = json.dumps(
            {"total": no_genes, "rows": result.to_dict(orient="records")}
        )

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self, task_def):

        genes = self.adata.var

        # Fields to return.
        select_fields = task_def["selectFields"]

        # if there is no search pattern defined, do not restrict gene names
        filter_pattern = task_def.get("geneNamesFilter", None)

        if filter_pattern and "gene_names" in select_fields:
            filter_query = "WHERE gene_names LIKE '{}'".format(filter_pattern)

        else:
            filter_query = ""

        # What the fields are ordered by
        order_by = task_def["orderBy"]

        # Order direction ('asc' or 'desc')
        order_direction = task_def["orderDirection"].upper()

        # Which page to start from
        offset = int(task_def["offset"])

        # How many to return.
        limit = int(task_def["limit"])

        # Set up SQL query and PandaSQL for efficient querying.
        query = """
            SELECT {}
              FROM genes
              {}
          ORDER BY {} {}
             LIMIT {}, {}
        """
        execute_query = lambda q: sqldf(q, {"genes": genes})

        query = query.format(
            ", ".join(select_fields),
            filter_query,
            order_by,
            order_direction,
            offset,
            limit,
        )
        result = execute_query(query)

        return self._format_result(result, no_genes=len(genes))
