
import os
import json

SETTINGS = f'{os.getcwd()}\\settings\\settings.json'

if not os.path.exists(SETTINGS):
    SETTINGS = f'/workflow/settings/settings.json'

with open(SETTINGS, 'r') as f:
    global settings
    settings = dict(json.load(f))

if __name__ == '__main__':
    print(settings)
