#!/usr/bin/env python
import pandas as pd

testdir = '/test-data'
outputdir = '/output'
print("Read in clinical file")
clinical_file = pd.read_csv(testdir +  "/clinical_data.csv")
predictions = clinical_file[['Study','Patient','D_Age']]
median = pd.np.median(predictions['D_Age'][~predictions['D_Age'].isnull()])
predictions['D_Age'][predictions['D_Age'].isnull()] = median
predictions.rename(columns={"Study":'study','Patient':'patient','D_Age':'predictionscore'},inplace=True)
predictions['highriskflag'] = predictions['predictionscore'] > 60
## Output the results
print("Write out prediction file")
predictions.to_csv(outputdir+"/predictions.tsv",index=False,sep="\t")
