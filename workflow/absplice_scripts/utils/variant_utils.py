from kipoiseq import Variant
import pandas as pd

def filter_long_variants(variant, max_length=10):
    if pd.isna(variant):
        return True
    
    variant = Variant.from_str(variant)
    length = max(len(variant.ref), len(variant.alt))
    return length < max_length