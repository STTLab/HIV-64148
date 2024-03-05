
import os
import json

SETTINGS, SECRETS = None, None
allow_settings = [f'{os.getcwd()}\\settings\\settings.json', '../settings/settings.json', f'/hiv64148/settings/settings.json', f'{os.getcwd()}/settings/settings.json']
for p in allow_settings:
    if os.path.exists(p):
        SETTINGS = p
        with open(SETTINGS, 'r') as f:
            global settings
            settings = dict(json.load(f))
        break

if SETTINGS is None: raise FileExistsError('Cannot find settings.json')
