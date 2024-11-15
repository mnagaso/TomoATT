# %%
# download src_ref_files from Zenodo
import os
import requests

url = 'https://zenodo.org/records/14065341/files/src_rec_file.dat?download=1'

path = "1_src_rec_files/src_rec_file.dat"

# check file existence
if not os.path.exists(path):
    try:
        os.mkdir("1_src_rec_files")
    except:
        pass
    print("Downloading src_rec_file.dat from Zenodo...")
    response = requests.get(url, stream=True)
    with open(path, 'wb') as out_file:
        out_file.write(response.content)
    print("Download complete.")
else:
    print("src_rec_file.dat already exists.")

# %%
# download initial model from Zenodo

url = 'https://zenodo.org/records/14065341/files/model_init_N51_89_33.h5?download=1'

path = "2_models/model_init_N51_89_33.h5"

# check file existence
if not os.path.exists(path):
    try:
        os.mkdir("2_models")
    except:
        pass
    print("Downloading model_init_N51_89_33.h5 from Zenodo...")
    response = requests.get(url, stream=True)
    with open(path, 'wb') as out_file:
        out_file.write(response.content)
    print("Download complete.")
else:
    print("model_init_N51_89_33.h5 already exists.")


