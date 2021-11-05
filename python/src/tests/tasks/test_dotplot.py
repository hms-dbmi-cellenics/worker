import pytest
from worker.tasks.dotplot import DotPlot


class TestDotplot:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "body": {
                "name": "DotPlot",
                "subset": {"cellClassKey": "louvain", "cellSetKey": "all"},
                "markerGenes": True,
                "input": {"nGenes": 5},
            }
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            DotPlot()

    def test_works_with_request(self):
        DotPlot(self.correct_request).compute()
        assert True
