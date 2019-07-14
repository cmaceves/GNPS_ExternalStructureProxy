import json
import argparse


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