
import json
SETTINGS = 'D:\\HIV-64148\\HIV-64148\\settings\\settings.json'

with open(SETTINGS, 'r') as f:
    global settings
    settings = dict(json.load(f))

if __name__ == '__main__':
    print(settings)
