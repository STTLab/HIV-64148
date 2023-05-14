import jinja2
from jinja2 import Environment, select_autoescape, FileSystemLoader

def render_report(context, save_to=None):
    env = Environment(
        loader=FileSystemLoader(searchpath='./utilities/reporter/templates'),
        autoescape=select_autoescape()
    )
    template = env.get_template('report_template.jinja')
    if save_to:
        with open(save_to, 'w') as f:
            f.write(template.render(context))
        return
    return template.render(context)
