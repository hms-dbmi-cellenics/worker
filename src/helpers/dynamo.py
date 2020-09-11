import boto3

from config import get_config

config = get_config()


def get_item_from_dynamo(experiment_id, item_name):
    dynamo = boto3.resource("dynamodb", **config.BOTO_RESOURCE_KWARGS).Table(
        config.DYNAMO_TABLE
    )

    resp = dynamo.get_item(
        Key={"experimentId": experiment_id}, ProjectionExpression=item_name,
    )

    if item_name in resp.get("Item", {}):
        returned_items = resp["Item"][item_name]
        print(
            "Successfully got {} from database for experiment id {}.".format(
                item_name, experiment_id
            )
        )
        return returned_items

    print(
        "Could not find item {} in the database for experiment id {}.".format(
            item_name, experiment_id
        )
    )
    return {}

