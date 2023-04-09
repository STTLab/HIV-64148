
import os
import json

allow_settings = [f'{os.getcwd()}\\settings\\settings.json', f'/workflow/settings/settings.json', f'{os.getcwd()}/settings/settings.json']
for p in allow_settings:
    if os.path.exists(p):
        SETTINGS = p

with open(SETTINGS, 'r') as f:
    global settings
    settings = dict(json.load(f))

if __name__ == '__main__':
    print(settings)
