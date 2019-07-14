import json
import argparse
import pandas as pd

def find_compound_name(compound_name, compound_df):
    filtered_df = compound_df.loc[compound_df['COMPOUND_NAME'] == compound_name]
    if len(filtered_df) > 0:
        return filtered_df.to_dict(orient="records")
    return None

def find_inchi_key(inchikey, compound_df):
    filtered_df = compound_df.loc[compound_df['COMPOUND_INCHIKEYBLOCK1'] == inchikey.split("-")[0]]
    if len(filtered_df) > 0:
        return filtered_df.to_dict(orient="records")
    return None

parser = argparse.ArgumentParser(description='Generate Statistics Between GNPS and NPAtlas')
parser.add_argument('npatlasjson', help='JSON for NPAtlas')
parser.add_argument('gnpsjson', help='JSON for GNPS')

args = parser.parse_args()

npatlas_list = json.load(open(args.npatlasjson, encoding='utf-8', errors='strict'), strict=False)
gnps_list = json.load(open(args.gnpsjson, encoding='utf-8', errors='strict'), strict=False)

all_gnps_keys = set([element["COMPOUND_INCHIKEY"].split("-")[0] for element in gnps_list])
all_npatlas_keys = set([element["COMPOUND_INCHIKEY"].split("-")[0] for element in npatlas_list])

print("ALL UNIQUE GNPS FLAT KEYS", len(all_gnps_keys))
print("ALL UNIQUE NPATLAS FLAT KEYS", len(all_npatlas_keys))
print("Intersection", len(all_gnps_keys.intersection(all_npatlas_keys)))


#Trying to determine if we can find discordant annotations
all_npatlas_name = set([element["COMPOUND_NAME"] for element in npatlas_list])
npatlas_df = pd.DataFrame(npatlas_list)
npatlas_df['COMPOUND_INCHIKEYBLOCK1'] = npatlas_df.apply(lambda x: x["COMPOUND_INCHIKEY"].split("-")[0], axis=1)

discordant_list = []
for gnps_compound in gnps_list:
    gnps_compound_name = gnps_compound["COMPOUND_NAME"]
    gnps_inchikey = gnps_compound["COMPOUND_INCHIKEY"]
    gnps_inchikeyblock1 = gnps_inchikey.split("-")[0]

    if gnps_compound_name in all_npatlas_name:
        np_atlas_compound_by_name = find_compound_name(gnps_compound_name, npatlas_df)
        npatlas_inchikeyblock1 = np_atlas_compound_by_name[0]["COMPOUND_INCHIKEYBLOCK1"]
        
        if gnps_inchikeyblock1 != npatlas_inchikeyblock1 and gnps_inchikeyblock1 != "None":
            discordant_dict = {}
            discordant_dict["gnps_compound_name"] = gnps_compound_name
            discordant_dict["gnps_inchikeyblock1"] = gnps_inchikeyblock1
            discordant_dict["npatlas_inchikeyblock1"] = npatlas_inchikeyblock1
            discordant_dict["GNPSID"] = gnps_compound["GNPSID"]
            discordant_dict["NPAID"] = np_atlas_compound_by_name[0]["NPAID"]

            discordant_list.append(discordant_dict)
            print(gnps_compound_name, gnps_inchikeyblock1, npatlas_inchikeyblock1, gnps_compound["GNPSID"])

pd.DataFrame(discordant_list).to_csv("discordant_pairs.tsv", sep='\t', index=False)