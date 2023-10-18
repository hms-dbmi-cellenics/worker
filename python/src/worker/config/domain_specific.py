import os

ACCOUNT_ID = {
  'BIOMAGE': '242905224710',
  'HMS': '160782110667',
}

domain_specific = {
    'HMS': {'timeout': 40 * 60},
    'BIOMAGE': {},
    'BIOMAGE_PRIVATE': {}
}


def get_domain_specific():
    if os.environ.get('NODE_ENV') == 'test':
        return domain_specific['TEST']

    aws_account_id = os.getenv('AWS_ACCOUNT_ID')

    if aws_account_id == ACCOUNT_ID['HMS']:
        return domain_specific['HMS']
    elif aws_account_id == ACCOUNT_ID['BIOMAGE']:
        return domain_specific['BIOMAGE']
    else:
        return domain_specific['BIOMAGE_PRIVATE']