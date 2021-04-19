import pandas as pd
df = pd.read_csv("annot.gtf",sep="\t",header=None,names=["seqname","source","feature","start","end","score","strand","frame","attribute"])
newdf = df[df["seqname"]=="NC_014776.1"]
newdf.to_csv("./annotEdit.gtf",sep="\t",header=False)
