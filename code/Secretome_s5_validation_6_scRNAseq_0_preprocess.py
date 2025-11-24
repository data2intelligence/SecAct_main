import sys
import pandas as pd

obj = pd.read_pickle(sys.argv[1])
obj.to_csv(sys.argv[2])
#obj.shape
