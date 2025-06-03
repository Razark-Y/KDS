import joblib
import pandas as pd
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from custom_transformers import find_all_orfs_fast, codon_table

def formatting_and_prediction(file_path, location, accession_id, collection_date):
    BASE_DIR = os.path.dirname(__file__)  # folder where preprocessor.py is
    pipeline_path = os.path.join(BASE_DIR, "stacked_pipeline.joblib")
    model_path = os.path.join(BASE_DIR, "lgbm_classifier.pkl")
    stacked_pipeline = joblib.load(pipeline_path)
    lgbm_classifier = joblib.load(model_path)
    
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
