from pandasql import sqldf
from result import Result


class ListGenes:
    def __init__(self, adata):
        self.adata = adata

    def _format_result(self, result):

        # JSONify result.
        result = result.to_json(orient="records")

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self, task_def):
        genes = self.adata.var

        # Fields to return.
        select = task_def["selectFields"]

        # What the fields are ordered by
        order_by = task_def["orderBy"]

        # Order direction ('asc' or 'desc')
        order_direction = task_def["orderDirection"].upper()

        # Return only from this index
        offset = int(task_def["offset"])

        # How many to to return.
        limit = int(task_def["limit"])

        # Set up SQL query and PandaSQL for efficient querying.
        query = """
            SELECT {}
              FROM genes
          ORDER BY {} {}
             LIMIT {}, {}
        """
        execute_query = lambda q: sqldf(q, {"genes": genes})

        query = query.format(
            ", ".join(select), order_by, order_direction, offset, limit
        )
        result = execute_query(query)

        return self._format_result(result)
