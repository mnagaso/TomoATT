# %%
# download model file from Zenodo
import os
import requests

url = 'https://zenodo.org/records/14160818/files/files_for_plotting.tar.gz?download=1'

path = "files_for_plotting.tar.gz"

# check file existence
if not os.path.exists(path):
    print("Downloading src_rec_file.dat from Zenodo...")
    print("The file is about 400 MB, so it may take a while.")
    response = requests.get(url, stream=True)
    with open(path, 'wb') as out_file:
        out_file.write(response.content)
    print("Download complete.")
else:
    print("files_for_plotting.tar.gz already exists.")


