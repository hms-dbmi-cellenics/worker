import pytest
from worker.tasks.dotplot import DotPlot


class TestDotplot:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "body": {
                "name": "DotPlot",
                "useMarkerGenes" : True,
                "numberOfMarkers": 3,
                "customGenesList": ["Gene1", "Gene2", "Gene3"],
                "groupBy": "louvain",
                "filterBy": {
                    "group": "All",
                    "key": "All"
                }
            }
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            DotPlot()

    def test_works_with_request(self):
        DotPlot(self.correct_request).compute()
        assert True