from pandasql import sqldf
from result import Result
import json


class ListGenes:
    def __init__(self, msg, adata):
        self.task_def = msg["body"]
        self.adata = adata

    def _format_result(self, result, total):
        # convert result to list of row dicts
        result = result.to_dict(orient="records")

        # JSONify result.
        result = json.dumps({"total": total, "rows": result})

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self):
        genes = self.adata.var

        # Fields to return.
        select_fields = self.task_def["selectFields"]

        # if there is no search pattern defined, do not restrict gene names
        filter_pattern = self.task_def.get("geneNamesFilter", None)

        if filter_pattern and "gene_names" in select_fields:
            filter_query = "WHERE gene_names LIKE '{}'".format(filter_pattern)
        else:
            filter_query = ""

        # What the fields are ordered by
        order_by = self.task_def["orderBy"]

        # Order direction ('asc' or 'desc')
        order_direction = self.task_def["orderDirection"].upper()

        # Return only from this index
        offset = int(self.task_def["offset"])

        # How many to to return.
        limit = int(self.task_def["limit"])

        # Set up SQL query and PandaSQL for efficient querying.
        query = """
            SELECT {}, count(*) OVER() AS full_count
              FROM genes
              {}
          ORDER BY {} {}
             LIMIT {}
             OFFSET {}
        """
        execute_query = lambda q: sqldf(q, {"genes": genes})

        query = query.format(
            ", ".join(select_fields),
            filter_query,
            order_by,
            order_direction,
            limit,
            offset,
        )
        result = execute_query(query)

        # Get total number of results if there are results
        total = 0

        if len(result) > 0:
            total = result["full_count"][0]

        # total returns numpy int64, convert to integer
        # for serialization to JSON
        total = int(total)

        # Filter out aggregate
        result = result[select_fields]

        return self._format_result(result, total=total)
