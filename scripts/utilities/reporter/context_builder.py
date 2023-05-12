from uuid import uuid4
from Bio import SeqIO
from sierrapy.sierraclient import Sequence

from workflow.components import BLAST
from utilities.apis import hivdb_seq_analysis
from utilities.settings import settings

def collapse_context_builder(haplotypes, blast_result_csv, hivdb_result=None):
    nhap = len(haplotypes)
    all_context = []
    for i, seq in enumerate(haplotypes):
        context = {}
        while True:
            generated = uuid4().hex
            if not generated[0].isnumeric():
                context['collapse_uuid'] = generated
                break

        context['name'] = seq.id
        context['summary'] = {
            'name': seq.id,
            'depth': seq.description.split()[1],
            'abundance': seq.description.split()[2].replace('freq=',''),
            'length': len(seq.seq)
        }
        context['sequence'] = str(seq.seq)
        blast_result = BLAST.BLASTResult.read(blast_result_csv).get_top(5)
        blast_result['subtype'] = [ BLAST.get_subtypes(settings['data']['blast']['dbtitle'], accession) for accession in blast_result['sseqid'] ]
        context['blast_result'] = blast_result.loc[blast_result['qseqid']==seq.id].to_dict('records')

        if hivdb_result:
            context['drug_resistant_data'] = hivdb_result[i].drugResistance.results
            context['mutations'] = {}
            context['mutations']['validation_results'] = hivdb_result[i].validationResults
            context['mutations']['alignment'] = hivdb_result[i].alignGeneSequences
        all_context.append(context)
    return all_context

def context_builder(haplotype_fa, nanoplot_html: str|None=None, worker_info: dict={}, hivdb_result=None, simulation_data=None):
    haplotypes = tuple(SeqIO.parse(haplotype_fa, 'fasta'))
    context = {
        'job_id': worker_info['job_id'],
        'worker': worker_info,
        'nanoplot_html': nanoplot_html,
        'collapses': collapse_context_builder(haplotypes, worker_info['run_stats']['blast_result'], hivdb_result),
        'hap_ids': [seq.id for seq in haplotypes],
        'hap_freq': [seq.description.split()[2].replace('freq=','') for seq in haplotypes],
    }
    subtype_count = {}
    blast_result = BLAST.BLASTResult.read(worker_info['blast_result'])
    blast_top_iden = blast_result.get_top(1)['sseqid'].to_numpy()
    for iden in blast_top_iden:
        subtype = BLAST.get_subtypes(settings['data']['blast']['dbtitle'], iden)
        if subtype in subtype_count.keys():
            subtype_count[subtype] += 1
        else: subtype_count[subtype] = 1
    context['subtype_count'] = {
        'keys': list(subtype_count.keys()),
        'values': list(subtype_count.values())
    }
    if simulation_data:
        context['is_simulated'] = True
        context['original_subtype_count'] = {
            'keys': list(simulation_data['subtype_count'].keys()),
            'values': list(simulation_data['subtype_count'].values())
        }
    return context
