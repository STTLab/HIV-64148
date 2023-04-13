
import os
import json

SETTINGS, SECRETS = None, None
allow_settings = [f'{os.getcwd()}\\settings\\settings.json', f'/hiv64148/settings/settings.json', f'{os.getcwd()}/settings/settings.json']
for p in allow_settings:
    if os.path.exists(p):
        SETTINGS = p
        with open(SETTINGS, 'r') as f:
            global settings
            settings = dict(json.load(f))
        break

allow_secrets = [f'{os.getcwd()}\\settings\\secrets.json', f'/hiv64148/settings/secrets.json', f'{os.getcwd()}/settings/secrets.json']
for p in allow_secrets:
    if os.path.exists(p):
        SECRETS = p
        with open(SECRETS, 'r') as f:
            global secrets
            secrets = dict(json.load(f))
        break

if SETTINGS is None: raise FileExistsError('Cannot fine settings.json')
if SECRETS is None: raise FileExistsError('Cannot find secret.json')
