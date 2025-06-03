import joblib
import pandas as pd
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from custom_transformers import find_all_orfs_fast, codon_table

def formatting_and_prediction(file_path, location, accession_id, collection_date):
    stacked_pipeline = joblib.load("stacked_pipeline.joblib")
    lgbm_classifier = joblib.load("lgbm_classifier.pkl")
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        fasta_seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    
    data = {
        'location': [location],
        'accession_id': [accession_id],
        'collection_date': [collection_date],
        'FASTA': [fasta_seq],
        'length': [len(fasta_seq)]
    }
    orfs = find_all_orfs_fast(fasta_seq, codon_table)
    data['protein_strain'] = len(orfs)
    df = pd.DataFrame(data)
    transformed = stacked_pipeline.transform(df)
    prediction = lgbm_classifier.predict(transformed)
    return prediction[0]
