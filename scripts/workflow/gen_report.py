
from random import randint
from uuid import uuid4
from yattag import Doc, indent

def generate_report_skeleton(run_id):
    doc, tag, text = Doc().tagtext()

    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            doc.stag('meta', charset='utf-8')
            doc.stag('meta', name='viewport', content='width=device-width, initial-scale=1')
            with tag('title'): text('HIV-64148 Report')
            doc.stag(
                'link',
                href='https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha3/dist/css/bootstrap.min.css',
                rel='stylesheet',
                integrity='sha384-KK94CHFLLe+nY2dmCWGMq91rCGa5gtU4mk92HdvYe+M/SXH301p5ILy+dN9+nJOZ',
                crossorigin='anonymous'
            )
            with tag(
                'script',
                src='https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha3/dist/js/bootstrap.bundle.min.js',
                integrity='sha384-ENjdO4Dr2bkBIFxQpeoTz1HIcje39Wm4jDKdf19U8gI4ddQ3GYNS7NTKfAdVQSZe',
                crossorigin='anonymous'
            ): pass
        with tag('body'):
            with tag('div', klass='container-sm mt-5'):
                with tag('h1'):
                    text('HIV-64148 Report')
                with tag('h4'):
                    text(f'run id: {run_id}')
                with tag('div'):
                    nhap = randint(1,10)
                    doc.asis(generate_table(nhap))
                    for i in range(1, nhap+1):
                        doc.asis(generate_collapse(f'Haplotype {i}', 'information not avialble'))
            doc.asis(generate_footer())
    print(indent(doc.getvalue()))
    with open('scripts\\sketch\\report.html', 'w') as html:
        html.write(indent(doc.getvalue()))

def generate_table(nhaplotypes):
    doc, tag, text = Doc().tagtext()
    with tag('table'):
        with tag('tbody'):
            with tag('tr'):
                with tag('td'):
                    text('# Haplotypes')
                with tag('td'):
                    text(nhaplotypes)
    return doc.getvalue()

def generate_collapse(title, desc, ):
    doc, tag, text = Doc().tagtext()
    collapse_id = str(uuid4())
    btn_attrs = {
        'klass': 'btn btn-light',
        'style': 'width: 100%; text-align: left;',
        'data-bs-toggle': 'collapse',
        'data-bs-target': f'#{collapse_id}',
        'aria-expanded': 'false',
        'aria-controls': collapse_id
    }
    with tag('div', klass='mt-2'):
        with tag('button',**btn_attrs):
            text(title)
        with tag('div', klass='collapse', id=collapse_id):
            with tag('div', klass='card card-body border-light'):
                text(desc)
    return doc.getvalue()

def generate_footer():
    doc, tag, text = Doc().tagtext()
    with tag('footer', klass='footer mt-auto py-3 bg-light fixed-bottom'):
        with tag('div', klass='container'):
            with tag('span', klass='text-muted'):
                text('(c) 2023 MEDCMU, HIV-64148 pipeline')
    return doc.getvalue()

if __name__ == '__main__':
    generate_report_skeleton(uuid4())
    # generate_collapse()
