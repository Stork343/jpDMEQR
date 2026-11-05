from pathlib import Path

import pandas as pd

source = Path("data/raw/pregnancy_detroit_20261105/official_cloud/pone.0181468.s001.csv")
frame = pd.read_csv(source, header=None, dtype=str)
print(frame.shape)
print(frame.iloc[:8, :10].to_string(index=False, header=False))
print(frame.iloc[-3:, :5].to_string(index=False, header=False))
