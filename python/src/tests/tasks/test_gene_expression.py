import json
import io

import random
import numpy as np
import boto3
import pytest
import responses
import mock

from botocore.stub import Stubber
from exceptions import RWorkerException
from worker.config import config
from worker.tasks.gene_expression import GeneExpression

from tests.data.cell_sets_from_s3 import cell_sets_from_s3

from tests.utils import get_cell_ids 

class TestGeneExpression:
    def get_s3_stub(self, cell_sets):
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        response = {
            "ContentLength": 10,
            "ContentType": "utf-8",
            "ResponseMetadata": {
                "Bucket": config.CELL_SETS_BUCKET,
            },
        }

        expected_params = {
            "Bucket": config.CELL_SETS_BUCKET,
            "Key": config.EXPERIMENT_ID,
        }
        stubber = Stubber(s3)
        stubber.add_response("head_object", response, expected_params)

        # Get object
        content_bytes = json.dumps(cell_sets, indent=2).encode(
            "utf-8"
        )

        data = io.BytesIO()
        data.write(content_bytes)
        data.seek(0)

        response = {
            "ContentLength": len(cell_sets),
            "ContentType": "utf-8",
            "Body": data,
            "ResponseMetadata": {
                "Bucket": config.CELL_SETS_BUCKET,
            },
        }
        stubber.add_response("get_object", response, expected_params)
        return (stubber, s3)

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_one_gene = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["Tpt1"],
                "downsampled": False,
            },
        }
        self.correct_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["Tpt1", "Zzz3"],
                "downsampled": False,
            },
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            GeneExpression()

    def test_works_with_request(self):
        GeneExpression(self.correct_request)

    @responses.activate
    def test_should_throw_exception_on_r_worker_error(self):

        error_code = "MOCK_R_WORKER_ERROR"
        user_message = "Some worker error"

        payload = {"error": {"error_code": error_code, "user_message": user_message}}

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/runExpression",
            json=payload,
            status=200,
        )

        with pytest.raises(RWorkerException) as exception_info:
            GeneExpression(self.correct_request).compute()

        assert exception_info.value.args[0] == error_code
        assert exception_info.value.args[1] == user_message

    def test_downsamples_correctly(self):
        stubber, s3 = self.get_s3_stub(cell_sets_from_s3)

        random.seed(900)

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3

            py_request = {
                "experimentId": config.EXPERIMENT_ID,
                "timeout": "2099-12-31 00:00:00",
                "body": {
                    "name": "GeneExpression",
                    "genes": ["Tpt1"],
                    "downsampled": True,
                    "downsampleSettings": {
                        "maxCells": 100,
                        "selectedCellSet": "louvain",
                        "groupedTracks": ["louvain"],
                        "selectedPoints": 'All',
                        "hiddenCellSets": []
                    }
                },
            }

            bla = GeneExpression(py_request)

            r_request = bla._format_request()
            assert isinstance(py_request, dict)

        expected_cell_ids = [324, 622, 166, 916, 38, 344, 31, 374, 630, 386, 149, 22, 68, 202, 620, 777, 701, 254, 134, 679, 384, 113, 277, 554, 213, 422, 751, 903, 247, 564, 356, 495, 655, 582, 882, 352, 331, 127, 673, 135, 89, 141, 814, 262, 506, 792, 502, 404, 599, 879, 594, 287, 864, 896, 21, 291, 547, 0, 351, 176, 13, 742, 285, 170, 121, 669, 132, 787, 319, 548, 760, 320, 315, 553, 230, 557, 371, 180, 556, 691, 409, 219, 289, 736, 726, 387, 909, 821, 768, 175, 771, 310, 207, 443, 158, 498, 697]

        assert r_request["cellIds"] == expected_cell_ids
    
    def test_downsamples_by_many_groups_correctly(self):
        stubber, s3 = self.get_s3_stub(cell_sets_from_s3)

        random.seed(900)

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3

            py_request = {
                "experimentId": config.EXPERIMENT_ID,
                "timeout": "2099-12-31 00:00:00",
                "body": {
                    "name": "GeneExpression",
                    "genes": ["Tpt1"],
                    "downsampled": True,
                    "downsampleSettings": {
                        "maxCells": 100,
                        "selectedCellSet": "louvain",
                        "groupedTracks": ["louvain", "sample"],
                        "selectedPoints": 'All',
                        "hiddenCellSets": []
                    }
                },
            }

            bla = GeneExpression(py_request)

            r_request = bla._format_request()
            assert isinstance(py_request, dict)

        expected_cell_ids = [899, 644, 199, 558, 422, 551, 137, 178, 536, 100, 244, 131, 838, 229, 827, 665, 202, 134, 334, 57, 252, 420, 54, 438, 650, 383, 174, 595, 446, 397, 151, 221, 156, 624, 681, 882, 314, 298, 333, 465, 618, 382, 98, 458, 352, 876, 14, 452, 670, 868, 21, 303, 883, 0, 218, 831, 94, 376, 522, 122, 654, 484, 433, 586, 180, 683, 834, 887, 43, 302, 894, 371, 556, 780, 64, 851, 395, 329, 32, 890, 425, 245, 289, 25, 702, 771, 182, 17, 839, 205, 165]

        assert r_request["cellIds"] == expected_cell_ids
    
    def test_downsamples_with_filter_correctly(self):
        stubber, s3 = self.get_s3_stub(cell_sets_from_s3)

        random.seed(900)

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3

            py_request = {
                "experimentId": config.EXPERIMENT_ID,
                "timeout": "2099-12-31 00:00:00",
                "body": {
                    "name": "GeneExpression",
                    "genes": ["Tpt1"],
                    "downsampled": True,
                    "downsampleSettings": {
                        "maxCells": 100,
                        "selectedCellSet": "louvain",
                        "groupedTracks": ["louvain", "sample"],
                        "selectedPoints": 'louvain/louvain-6',
                        "hiddenCellSets": []
                    }
                },
            }

            bla = GeneExpression(py_request)

            r_request = bla._format_request()
            assert isinstance(py_request, dict)
        
        louvain_6_cell_ids = get_cell_ids("louvain", "louvain-6", cell_sets_from_s3)
        
        # Contains only louvain 6 ids (filtered out the rest), doesn't downsample because not necessary
        assert set(r_request["cellIds"]) == set(louvain_6_cell_ids)
        
        # They were reordered to match the groups
        assert r_request["cellIds"] != louvain_6_cell_ids

    def test_downsamples_with_hidden_cell_sets(self):
        stubber, s3 = self.get_s3_stub(cell_sets_from_s3)

        random.seed(900)

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3

            py_request = {
                "experimentId": config.EXPERIMENT_ID,
                "timeout": "2099-12-31 00:00:00",
                "body": {
                    "name": "GeneExpression",
                    "genes": ["Tpt1"],
                    "downsampled": True,
                    "downsampleSettings": {
                        "maxCells": 100,
                        "selectedCellSet": "sample",
                        "groupedTracks": ["louvain", "sample"],
                        "selectedPoints": 'louvain/louvain-6',
                        "hiddenCellSets": ["5d88f799-c704-4667-99f8-8d6dee6cfc22"]
                    }
                },
            }

            bla = GeneExpression(py_request)

            r_request = bla._format_request()
            assert isinstance(py_request, dict)
        
        louvain_6_cell_ids = get_cell_ids("louvain", "louvain-6", cell_sets_from_s3)
        wt2_cell_ids = get_cell_ids("sample", "5d88f799-c704-4667-99f8-8d6dee6cfc22", cell_sets_from_s3)
        
        # Doesnt contain all louvain 6 cell ids (some were filtered out)
        assert set(r_request["cellIds"]) != set(louvain_6_cell_ids)
        
        # Contains all louvain 6 cell ids that are not in sample wt2
        assert set(r_request["cellIds"]) == set(louvain_6_cell_ids).difference(set(wt2_cell_ids))