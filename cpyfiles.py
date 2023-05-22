#%%

import os
import shutil


sys="IEEE118"
teste="Teste3/"
destination_dir="restTCSC/"+teste+sys+"/"

var="10"



file_names = ['conds.csv', 'conv_A.csv', 'conv_B.csv']


for file_name in file_names:
    source_file = os.path.join("", file_name)
    base_name = os.path.splitext(file_name)[0]
    new_file_name=f"{base_name}_{var}"
    destination_file = os.path.join(destination_dir, new_file_name)
    shutil.copy2(source_file, destination_file+".csv")

#%%