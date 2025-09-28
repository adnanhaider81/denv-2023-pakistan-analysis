#!/usr/bin/env python3
# Extract the annotated polyprotein CDS from a reference and translate sample genomes
import argparse, os
from Bio import Entrez, SeqIO

def fetch_ref(acc, email, api_key=None):
    Entrez.email = email
    if api_key: Entrez.api_key = api_key
    h = Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
    rec = SeqIO.read(h, 'genbank'); h.close()
    return rec

def polyprotein_coords(rec):
    for f in rec.features:
        if f.type == 'CDS':
            # first CDS is the polyprotein for dengue references
            return int(f.location.start), int(f.location.end), int(f.location.strand or 1)
    raise SystemExit('No CDS in reference')

ap = argparse.ArgumentParser(description='Translate polyprotein and compare to reference AA')
ap.add_argument('--email', required=False)
ap.add_argument('--api_key', required=False)
ap.add_argument('--ref_acc', required=True, help='AAW64436.1 for DENV-1, NC_001474.1 for DENV-2')
ap.add_argument('--genomes', required=True, help='FASTA with one or more sample genomes')
ap.add_argument('--out_tsv', required=True)
a = ap.parse_args()

email = a.email or os.getenv('NCBI_EMAIL')
if not email: raise SystemExit('Set --email or env NCBI_EMAIL')
api_key = a.api_key or os.getenv('NCBI_API_KEY')

ref_rec = fetch_ref(a.ref_acc, email, api_key)
s, e, strand = polyprotein_coords(ref_rec)
ref_nt = ref_rec.seq[s:e] if strand == 1 else ref_rec.seq[s:e].reverse_complement()
ref_aa = ref_nt.translate()

with open(a.out_tsv, 'w') as out:
    out.write('sample	pos	ref_aa	alt_aa
')
    for rec in SeqIO.parse(a.genomes, 'fasta'):
        aa = rec.seq.translate()
        L = min(len(aa), len(ref_aa))
        for i in range(L):
            if aa[i] != ref_aa[i]:
                out.write(f'{rec.id}	{i+1}	{ref_aa[i]}	{aa[i]}
')
print('Wrote', a.out_tsv)
