

import os, sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

def load_dataset():
    import pandas as pd
    return pd.read_feather(os.path.join(os.path.dirname(__file__), 'resources/my_dataset.feather'))

