
import os
import json

allow_settings = [f'{os.getcwd()}\\settings\\settings.json', f'/hiv64148/settings/settings.json', f'{os.getcwd()}/settings/settings.json']
for p in allow_settings:
    if os.path.exists(p):
        SETTINGS = p

with open(SETTINGS, 'r') as f:
    global settings
    settings = dict(json.load(f))

allow_secrets = [f'{os.getcwd()}\\settings\\secrets.json', f'/hiv64148/settings/secrets.json', f'{os.getcwd()}/settings/secrets.json']
for p in allow_secrets:
    if os.path.exists(p):
        SECRETS = p

with open(SECRETS, 'r') as f:
    global secrets
    secrets = dict(json.load(f))

if __name__ == '__main__':
    print(settings)
