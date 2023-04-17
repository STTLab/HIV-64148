
import time
from Bio import SeqIO
from uuid import uuid4
from yattag.doc import Doc
from bs4 import BeautifulSoup
from yattag.indentation import indent
from .components import BLAST
from .workflow import Worker

def generate_report_skeleton(run_id, haplotype_file: str, output:str, nanoplot_html: str|None=None, worker: Worker|None=None):
    def read_haplotype_fa(haplotype_fa):
        haplotypes = SeqIO.parse(haplotype_fa, 'fasta')
        return tuple(haplotypes)
    haplotypes = read_haplotype_fa(haplotype_file)
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
                if worker:
                    with tag('section'):
                        with tag('p'): 
                            with tag('table', klass='table table-bordered'):
                                with tag('tr'):
                                    # Time used
                                    with tag('th', klass='table-secondary'): text('Wall clock time')
                                    with tag('td'): text(str(worker.get_runtime()))
                                # Peak memory usage
                                peak_mem = worker.get_peak_mem()
                                with tag('tr'):
                                    with tag('th', klass='table-secondary'): text('Peak memory usage')
                                with tag('tr'):
                                    with tag('td'): text('Strainline')
                                    with tag('td'): 
                                        text(round(peak_mem['strainline'], 2))
                                        text(' Mib')
                                with tag('tr'):
                                    with tag('td'): text('BLAST')
                                    with tag('td'): 
                                        text(round(peak_mem['blast'], 2))
                                        text(' Mib')
                                with tag('tr'):
                                    with tag('td'): text('Snippy')
                                    with tag('td'): 
                                        text(round(peak_mem['snippy']))
                                        text(' Mib')
                            pass
                    if nanoplot_html:
                        with tag('a', klass='btn btn-primary', href=f'../{nanoplot_html}'): text('QC Report')
                        
                with tag('section'):
                    with tag('div', klass='row'):
                        with tag('div', klass='col-sm'):
                            doc.asis('<canvas id="haplotype-chart"></canvas>')
                        with tag('div', klass='col'): pass
                    doc.asis(generate_haplotype_table(nhap))
                    hap_ids = []
                    hap_freq = []
                    for i, seq in enumerate(haplotypes):
                        hap_ids.append(seq.id)
                        hap_freq.append(seq.description.split()[2].replace('freq=',''))
                        summary_table = generate_summary_table(seq)
                        blast_table = generate_blast_table(BLAST.BLASTResult.read('.idea\\output\\haplotype.blast.csv'), seq)
                        drug_resistant = generate_drug_resistant_profile()
                        doc.asis(generate_collapse(f'Haplotype {i+1}', [summary_table, blast_table, drug_resistant]))
            doc.asis(generate_footer())
            doc.asis(script_plot_haplotype('haplotype-chart', hap_ids, hap_freq))
    print(indent(doc.getvalue()))
    with open(output, 'w') as html:
        html.write(indent(doc.getvalue()))

def generate_haplotype_table(nhaplotypes):
    doc, tag, text = Doc().tagtext()
    with tag('table'):
        with tag('tbody'):
            with tag('tr'):
                with tag('td'):
                    text('# Haplotypes')
                with tag('td'):
                    text(nhaplotypes)
    return doc.getvalue()

def generate_collapse(title, subelements:list):
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
                for element in subelements:
                    doc.asis(element)
    return doc.getvalue()

def generate_summary_table(seq):
    data = {
        'Name': seq.id,
        'Depth': seq.description.split()[1],
        'Abundance': seq.description.split()[2].replace('freq=',''),
        'Length': len(seq.seq),
    }
    doc, tag, text = Doc().tagtext()
    with tag('div'):
        with tag('h5'): text('Summary')
        with tag('div', klass='row'):
            with tag('div', klass='col-4'):
                with tag('table', klass='table table-bordered'):
                    for key, val in zip(data.keys(), data.values()):
                        with tag('tr'):
                            with tag('th', klass='table-secondary'): text(key)
                            with tag('td'): text(val)
            with tag('div', klass='col-8'):
                with tag('div', klass='overflow-auto p-3 mb-3 mb-md-0 me-md-3 bg-light', style='max-width: 100%; max-height: 165px;'): text(f'>{seq.description}\n{str(seq.seq)}')
    return doc.getvalue()

def generate_blast_table(blast_result: BLAST.BLASTResult, seq: SeqIO.SeqRecord):
    blast_top_iden = blast_result.get_top(5)
    doc, tag, text = Doc().tagtext()
    with tag('div'):
        with tag('h5', klass='mt-4'): text('BLAST Result')
        with tag('table', klass='table table-hover'):
            with tag('thead'):
                with tag('tr'):
                    with tag('th'): text('#')
                    with tag('th'): text('Query')
                    with tag('th'): text('Matched')
                    with tag('th'): text('Subtype')
                    with tag('th'): text('% Identity')
                    with tag('th'): text('Length')
                    with tag('th'): text('mismatch')
                    with tag('th'): text('bitscore')
            with tag('tbody'):
                for i, record in enumerate(blast_top_iden.loc[blast_top_iden['qseqid']==seq.id].to_numpy()):
                    with tag('tr'):
                        with tag('td'): text(i+1)
                        with tag('td'): text(record[0])
                        with tag('td'): 
                            with tag('a', href=f'https://www.ncbi.nlm.nih.gov/nuccore/{record[1]}', target='_blank'): text(record[1])
                        with tag('td'): text(BLAST.get_subtypes('32hiv1_default_db', record[1]))
                        with tag('td'): text(record[2])
                        with tag('td'): text(record[3])
                        with tag('td'): text(record[4])
                        with tag('td'): text(record[11])
    return doc.getvalue()

def generate_footer():
    doc, tag, text = Doc().tagtext()
    doc.asis('<div style="height: 100px;"></div>')
    with tag('footer', klass='footer mt-5 py-3 bg-light fixed-bottom', style='display:block;'):
        with tag('div', klass='container'):
            with tag('span', klass='text-muted'):
                text('(c) 2023 MEDCMU, HIV-64148 pipeline')
    return doc.getvalue()

def generate_drug_resistant_profile():
    doc, tag, text = Doc().tagtext()
    with tag('h5'): text('Drug resistant Profile')
    return doc.getvalue()

def script_plot_haplotype(element_id, labels, frequencies):
    doc, tag, text = Doc().tagtext()
    doc.asis('<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>')
    with tag('script'):
        text('const ctx = document.getElementById("haplotype-chart");\n')
        text('new Chart(ctx, {\n\ttype: "doughnut",')
        text('\n\tdata: {\n\t\tlabels: ')
        text(str(labels) + ',')
        text('\n\t\tdatasets: [{ label: "Haplotype Abundance", data:')
        text(str(frequencies) + ',')
        text('\n\t\tborderWidth: 1\n\t}]\n},\noptions: { plugins: { legend:{position: "right"},}}});')
    return doc.getvalue()

def nanoplot_extract_graph(nanoplot_html):
    doc, tag, text = Doc().tagtext()
    soup = BeautifulSoup(open(nanoplot_html, 'r').read(), features='html.parser')
    graph_divs = list(map(str, soup.find_all('div', {'class': 'plotly-graph-div'})))
    scripts = list( map(str, dict.fromkeys(soup.find_all('script'))) )
    doc.asis(' '.join(graph_divs))
    doc.asis(scripts.pop(0))
    doc.asis(scripts.pop(0))
    for script in scripts:
        if 'PLOTLYENV' in script:
            doc.asis(script)
    return doc.getvalue()