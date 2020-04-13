# views.py
from flask import abort, jsonify, render_template, request, redirect, url_for, send_file, send_from_directory

from app import app

import json
import csv
import requests
import requests_cache
import utils
import pandas as pd

requests_cache.install_cache('demo_cache')

@app.route('/heartbeat', methods=['GET'])
def heartbeat():
    return "{'status' : 'up'}"

@app.route('/npatlasproxyimg', methods=['GET'])
def npatlasproxyimg():
    inchi = request.args.get('inchi', '')
    inchikey = request.args.get('inchikey', '')
    smiles = request.args.get('smiles', '')

    NPAID = get_npatlas(smiles, inchi, inchikey)

    if NPAID == None:
        return send_file("./static/img/Solid_white.png")
    else:
        return send_file("./static/img/npatlas_logo.png")

@app.route('/npatlasproxy', methods=['GET'])
def npatlasproxy():
    inchi = request.args.get('inchi', '')
    inchikey = request.args.get('inchikey', '')
    smiles = request.args.get('smiles', '')

    NPAID = get_npatlas(smiles, inchi, inchikey)

    if NPAID == None:
        url = "https://www.npatlas.org/joomla/index.php/deposit"
        return redirect(url)
    else:
        url = "https://www.npatlas.org/joomla/index.php/explore/compounds#npaid=%s" % NPAID
        return redirect(url)


@app.route('/mibigproxyimg', methods=['GET'])
def mibigproxyimg():
    inchi = request.args.get('inchi', '')
    inchikey = request.args.get('inchikey', '')
    smiles = request.args.get('smiles', '')

    BGCID = get_mibig(smiles, inchi, inchikey)

    if BGCID == None:
        return send_file("./static/img/Solid_white.png")
    else:
        return send_file("./static/img/mibig_logo.png")

@app.route('/mibigproxy', methods=['GET'])
def mibigproxy():
    inchi = request.args.get('inchi', '')
    inchikey = request.args.get('inchikey', '')
    smiles = request.args.get('smiles', '')

    BGCID = get_mibig(smiles, inchi, inchikey)

    if BGCID == None:
        return render_template("notfound.html")
    else:
        url = "https://mibig.secondarymetabolites.org/repository/%s/index.html#r1c1" % BGCID
        return redirect(url)

def get_npatlas(smiles, inchi, inchikey):
    inchikey_from_smiles, inchikey_from_inchi = utils.get_inchikey(smiles, inchi)
    print(inchikey_from_smiles, inchikey_from_inchi)
    acceptable_key = set([inchikey.split("-")[0], inchikey_from_smiles.split("-")[0], inchikey_from_inchi.split("-")[0]])

    NPAID = None

    for npatlas_entry in npatlas_list:
        if len(npatlas_entry["COMPOUND_INCHIKEY"]) > 5 and npatlas_entry["COMPOUND_INCHIKEY"].split("-")[0] in acceptable_key:
            NPAID = npatlas_entry["NPAID"]
            break

    return NPAID

def get_mibig(smiles, inchi, inchikey):
    inchikey_from_smiles, inchikey_from_inchi = utils.get_inchikey(smiles, inchi)
    acceptable_key = set([inchikey.split("-")[0], inchikey_from_smiles.split("-")[0], inchikey_from_inchi.split("-")[0]])

    BGCID = None

    for mibig_entry in mibig_list:
        if len(mibig_entry["COMPOUND_INCHIKEY"]) > 5 and mibig_entry["COMPOUND_INCHIKEY"].split("-")[0] in acceptable_key:
            BGCID = mibig_entry["BGCID"]
            break

    return BGCID


# @app.route('/gnpsproxy', methods=['GET'])
# def gnpsproxy():
#     inchi = request.args.get('inchi', '')
#     inchikey = request.args.get('inchikey', '')
#     smiles = request.args.get('smiles', '')

#     inchikey_from_smiles, inchikey_from_inchi = get_inchikey(smiles, inchi)

#     acceptable_key = set([inchikey, inchikey_from_smiles, inchikey_from_inchi])

#     print(acceptable_key)

#     found_spectrum_list = []

#     for gnps_spectrum in gnps_listreturn json.dumps(:
#         if len(gnps_spectrum["InChIKey_smiles"]) > 2 and gnps_spectrum["InChIKey_smiles"] in acceptable_key:
#             found_spectrum_list.append(gnps_spectrum)
#         elif len(gnps_spectrum["InChIKey_inchi"]) > 2 and gnps_spectrum["InChIKey_inchi"] in acceptable_key:
#             found_spectrum_list.append(gnps_spectrum)

#     return render_template('gnpsspectralist.html', spectrumlist_json=found_spectrum_list)
#     #return json.dumps(found_spectrum_list)

@app.route('/structureproxy', methods=['GET'])
def structureproxy():
    inchi = request.args.get('inchi', '')
    inchikey = request.args.get('inchikey', '')
    smiles = request.args.get('smiles', '')

    inchikey_from_smiles, inchikey_from_inchi = utils.get_inchikey(smiles, inchi)

    inchikey_query = ""
    MIN_LENGTH = 5

    if len(inchikey) > MIN_LENGTH:
        inchikey_query = inchikey
    elif len(inchikey_from_smiles) > MIN_LENGTH:
        inchikey_query = inchikey_from_smiles
    elif len(inchikey_from_inchi) > MIN_LENGTH:
        inchikey_query = inchikey_from_inchi

    print(inchikey_query)

    kegg_info = requests.get("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/KEGG/%s" % (inchikey_query)).json()
    CHEBI_info = requests.get("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/CHEBI/%s" % (inchikey_query)).json()
    BioCyc_info = requests.get("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/BioCyc/%s" % (inchikey_query)).json()
    pubchem_info = requests.get("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/Pubchem CID/%s" % (inchikey_query)).json()
    ChemSpider_info = requests.get("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/ChemSpider/%s" % (inchikey_query)).json()
    LipidMaps_info = requests.get("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/LipidMaps/%s" % (inchikey_query)).json()
    DrugBank_info = requests.get("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/DrugBank/%s" % (inchikey_query)).json()

    print(pubchem_info)

    external_links = []
    external_links += prep_external(kegg_info[0]["results"], "KEGG", "https://www.genome.jp/dbget-bin/www_bget?%s")
    external_links += prep_external(CHEBI_info[0]["results"], "CHEBI", "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=%s")
    external_links += prep_external(BioCyc_info[0]["results"], "BioCyc", "https://biocyc.org/compound?orgid=META&id=%s")
    external_links += prep_external(pubchem_info[0]["results"], "Pubchem", "https://pubchem.ncbi.nlm.nih.gov/compound/%s")
    external_links += prep_external(ChemSpider_info[0]["results"], "Chem Spider", "http://www.chemspider.com/Chemical-Structure.%s.html")
    external_links += prep_external(LipidMaps_info[0]["results"], "Lipid Maps", "http://lipidmaps.org/data/LMSDRecord.php?LMID=%s")
    external_links += prep_external(DrugBank_info[0]["results"], "Drugbank", "https://www.drugbank.ca/drugs/%s")

    return render_template('externallist.html', external_links_json=external_links)


def prep_external(results_list, resource_name, resource_url):
    external_links = []

    for result in results_list:
        external_obj = {}
        external_obj["resource"] = resource_name
        external_obj["url"] = resource_url % (result)

        external_links.append(external_obj)

    return external_links

### GNPS Spectral Library Delivery Endpoints that will be constantly updated

#Making it easy to query for all of GNPS library spectra
@app.route('/gnpslibraryjson', methods=['GET'])
def gnpslibraryjson():
    return send_from_directory("/output", "gnpslibraries.json")

#This returns all the spectra
@app.route('/gnpslibraryformattedjson', methods=['GET'])
def gnpslibraryformattedjson():
    return send_from_directory("/output", "gnpslibraries_enriched_all.json")

#This returns all the spectra with peaks
@app.route('/gnpslibraryformattedwithpeaksjson', methods=['GET'])
def gnpslibraryformattedwithpeaksjson():
    return send_from_directory("/output", "ALL_GNPS.json")

#This returns all the spectra that have a structure
@app.route('/gnpslibraryfornpatlasjson', methods=['GET'])
def gnpslibraryfornpatlasjson():
    return send_from_directory("/output", "gnpslibraries_npatlas.json")

@app.route('/gnpslibraryfornpatlastsv', methods=['GET'])
def gnpslibraryfornpatlastsv():
    return send_from_directory("/output", "gnpslibraries_npatlas.tsv")

@app.route('/gnpslibrary/<library>.mgf', methods=['GET'])
def mgf_download(library):
    return send_from_directory("/output", "{}.mgf".format(library))

@app.route('/gnpslibrary/<library>.msp', methods=['GET'])
def msp_download(library):
    return send_from_directory("/output", "{}.msp".format(library))

@app.route('/gnpslibrary/<library>.json', methods=['GET'])
def json_download(library):
    return send_from_directory("/output", "{}.json".format(library))

npatlas_list = utils.load_NPAtlas("data/npatlas.json")
mibig_list = utils.load_mibig("data/mibig.csv")