from admet_ai import ADMETModel
import pandas as pd
import ace_tools_open as tools
import sys

Smiles = str(sys.argv[1])
print (Smiles)

model = ADMETModel()
preds = model.predict(smiles=Smiles)

df = pd.DataFrame(list(preds.items()), columns=['Property', 'Value'])
#tools.display_dataframe_to_user(name="Molecular Properties", dataframe=df)
#print(df.to_string(index=False))
print(df[:49])
