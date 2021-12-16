import mock
import pytest
import responses
from worker.helpers.mock_s3 import MockS3Class
from worker.tasks.dotplot import DotPlot


class TestDotplot:
    def get_request(self):
        request = {
            "body": {
                "name": "DotPlot",
                "useMarkerGenes": True,
                "numberOfMarkers": 3,
                "customGenesList": ["Gene1", "Gene2", "Gene3"],
                "groupBy": "cluster1",
                "filterBy": {"group": "All", "key": "All"},
            }
        }

        return request

    """
    Mocks the S3 query for fetching cell sets. Returns an
    empty cell set and yields the patched up object.
    """

    @pytest.fixture
    def mock_S3_get(self):
        with mock.patch("boto3.client") as m:
            mockS3 = MockS3Class()
            m.return_value = mockS3
            yield (m, mockS3)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            DotPlot()

    @responses.activate
    def test_generates_correct_request_keys(self, mock_S3_get):
        MockS3Class.setResponse("one_set")
        request = DotPlot(self.get_request())._format_request()
        assert isinstance(request, dict)

        # all expected keys are in the request
        expected_keys = [
            "useMarkerGenes",
            "numberOfMarkers",
            "customGenesList",
            "groupBy",
            "filterBy",
            "applyFilter",
        ]
        assert all(key in request for key in expected_keys)

    @responses.activate
    def test_group_by_equals_filter_by_when_filter_by_equals_all(self, mock_S3_get):
        MockS3Class.setResponse("one_set")
        request = DotPlot(self.get_request())._format_request()
        assert request["filterBy"] == request["groupBy"]
