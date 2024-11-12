''' Detect cryptic mutation clusters in wastewater sequencing data '''

import argparse
import pandas as pd

from outbreak_data import outbreak_data as od
from outbreak_data import authenticate_user
from outbreak_tools import crumbs

# def parse_aa_muts(muts):
#     output = []
#     for m in muts:
#         if '*' in m:
#             return None
#         if '(' in m and ':' in m:
#             mut = m.split('(')[1][:-1]
#             output.append(mut)
#     if len(output) == 0:
#         return None
#     return output

#lineage_key = crumbs.get_alias_key()
#mut_data = od.lineage_cl_prevalence('.', descendants=True, mutations=cluster, location='USA', datemin=START_DATE, datemax=END_DATE, lineage_key=lineage_key)

def parse_covariants(covariants_dir):
    ''' Parse freyja covariants output files, aggregate into one dataframe '''
    pass

def query_clinical_data():
    pass

def main():
    parser = argparse.ArgumentParser(description='Detect cryptic mutation clusters in wastewater sequencing data')
    parser.add_argument('covariants', help='Directory containing freyja covariants output')
    parser.add_argument('metadata', help='TSV containing at least "covariants_filename". Adds additional metadata to output file (see data/example_metadata.tsv)', optional=True)
    parser.add_argument('output', help='Output file name', optional=True, default='cryptic_clusters.tsv')
    args = parser.parse_args()

    # Authenticate with GISAID credentials
    authenticate_user.authenticate_new_user()

    # Read input file
    df = pd.read_csv(args.input, sep='\t')

    # Detect cryptic variants

    # Write output file
    df.to_csv(args.output, sep='\t', index=False)