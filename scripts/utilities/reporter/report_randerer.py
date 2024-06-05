'''
HIV-64148  Copyright (C) 2024  Sara Wattanasombat
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it.
'''

__version__ = '0.1'
__author__ = 'Sara Wattanasombat'

import jinja2
from jinja2 import Environment, select_autoescape, FileSystemLoader

def render_report(context, save_to=None):
    env = Environment(
        loader=FileSystemLoader(searchpath='/hiv64148/scripts/utilities/reporter/templates'),
        autoescape=select_autoescape()
    )
    template = env.get_template('report_template.jinja')
    if save_to:
        with open(save_to, 'w') as f:
            f.write(template.render(context))
        return
    return template.render(context)
