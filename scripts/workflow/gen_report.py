
from Bio import SeqIO
from uuid import uuid4
from yattag.doc import Doc
from yattag.indentation import indent

def generate_report_skeleton(run_id, haplotypes: tuple):
    nhap = len(haplotypes)
    doc, tag, text = Doc().tagtext()
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        # Head
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
        #
        # Body
        #
        with tag('body'):
            with tag('div', klass='container-sm mt-5'):
                with tag('h1'):
                    text('HIV-64148 Report')
                with tag('h4'):
                    text(f'run id: {run_id}')
                with tag('div'):
                    doc.asis(generate_table(nhap))
                    for i, seq in enumerate(haplotypes):
                        doc.asis(generate_collapse(f'Haplotype {i+1}', f'ID: {seq.id}<br>{seq.description}<br><br><span id="{seq.id}_seq" style="d">{seq.id} {seq.description}<br>{str(seq.seq)}</span>'))
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

def generate_collapse(title, desc):
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
                doc.asis(desc)
    return doc.getvalue()

def generate_seq_description(data):
    doc, tag, text = Doc().tagtext()
    with tag('p'):
        text()
    

def generate_footer():
    doc, tag, text = Doc().tagtext()
    with tag('footer', klass='footer mt-auto py-3 bg-light fixed-bottom'):
        with tag('div', klass='container'):
            with tag('span', klass='text-muted'):
                text('(c) 2023 MEDCMU, HIV-64148 pipeline')
    return doc.getvalue()

def read_haplotype_fa(haplotype_fa):
    haplotypes = SeqIO.parse(haplotype_fa, 'fasta')
    return tuple(haplotypes)

if __name__ == '__main__':
    generate_report_skeleton(uuid4(), read_haplotype_fa('scripts\\tests\\mock\\HIV1_Strainline_result.fa'))
    # generate_collapse()
