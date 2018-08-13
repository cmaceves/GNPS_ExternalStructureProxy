# views.py
from flask import abort, jsonify, render_template, request, redirect, url_for

from app import app

import json
import csv
import requests
from indigo import *
from indigo_renderer import *
from indigo_inchi import *

@app.route('/', methods=['GET'])
def homepage():
    return "MING"

@app.route('/npatlasproxy', methods=['GET'])
def npatlasproxy():
    inchi = request.args.get('inchi', '')
    inchikey = request.args.get('inchikey', '')
    smiles = request.args.get('smiles', '')

    inchikey_from_smiles, inchikey_from_inchi = get_inchikey(smiles, inchi)
    acceptable_key = set([inchikey, inchikey_from_smiles, inchikey_from_inchi])

    NPAID = None

    for npatlas_entry in npatlas_list:
        if len(npatlas_entry["InChIKey"]) > 2 and npatlas_entry["InChIKey"] in acceptable_key:
            NPAID = npatlas_entry["NPAID"]
            break
        # if npatlas_entry["InChI"] == inchi:
        #     NPAID = npatlas_entry["NPAID"]
        #     break
        # if npatlas_entry["SMILES"] == inchikey:
        #     NPAID = npatlas_entry["NPAID"]
        #     break
        # if npatlas_entry["InChIKey"] == smiles:
        #     NPAID = npatlas_entry["NPAID"]
        #     break

    if NPAID == None:
        return render_template("notfound.html")
    else:
        url = "http://www.npatlas.org/joomla/index.php/explore/compounds#npaid=%s" % NPAID
        return redirect(url)

@app.route('/gnpsproxy', methods=['GET'])
def gnpsproxy():
    inchi = request.args.get('inchi', '')
    inchikey = request.args.get('inchikey', '')
    smiles = request.args.get('smiles', '')

    inchikey_from_smiles, inchikey_from_inchi = get_inchikey(smiles, inchi)

    acceptable_key = set([inchikey, inchikey_from_smiles, inchikey_from_inchi])

    print(acceptable_key)

    found_spectrum_list = []

    for gnps_spectrum in gnps_list:
        if len(gnps_spectrum["InChIKey_smiles"]) > 2 and gnps_spectrum["InChIKey_smiles"] in acceptable_key:
            found_spectrum_list.append(gnps_spectrum)
        elif len(gnps_spectrum["InChIKey_inchi"]) > 2 and gnps_spectrum["InChIKey_inchi"] in acceptable_key:
            found_spectrum_list.append(gnps_spectrum)

    return render_template('gnpsspectralist.html', spectrumlist_json=found_spectrum_list)
    #return json.dumps(found_spectrum_list)


@app.route('/heartbeat', methods=['GET'])
def heartbeat():
    return "{'status' : 'up'}"


def get_inchikey(smiles, inchi):
    inchikey_from_smiles = ""
    inchikey_from_inchi = ""
    try:
        m = indigo.loadMolecule(smiles)
        inchi_intermediate = indigo_inchi.getInchi(m)
        inchikey_from_smiles = indigo_inchi.getInchiKey(inchi_intermediate)
    except:
        inchikey_from_smiles = ""

    try:
        inchikey_from_inchi = indigo_inchi.getInchiKey(inchi)
    except:
        inchikey_from_inchi = ""

    #print(inchikey_from_smiles, inchikey_from_inchi, smiles, inchi)

    if len(inchikey_from_smiles) > 2 and len(inchikey_from_inchi) > 2:
        return inchikey_from_smiles, inchikey_from_inchi
        #if inchikey_from_smiles != inchikey_from_inchi:
        #    print(smiles, inchi)
            #print(inchikey_from_smiles, inchikey_from_inchi)

    if len(inchikey_from_smiles) > 2:
        return inchikey_from_smiles, ""

    if len(inchikey_from_inchi) > 2:
        return inchikey_from_inchi, ""

    return "", ""

def load_NPAtlas(filepath):
    print("Loading")
    all_npatlas = []
    with open(filepath) as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            all_npatlas.append(row)
    return all_npatlas

def load_GNPS():
    #url = "https://gnps.ucsd.edu/ProteoSAFe/LibraryServlet?library=all"
    url = "https://gnps.ucsd.edu/ProteoSAFe/LibraryServlet?library=GNPS-LIBRARY"
    all_GNPS = requests.get(url).json()
    all_spectra = []

    for spectrum in all_GNPS["spectra"]:
        smiles = spectrum["Smiles"]
        inchi =  spectrum["INCHI"]
        inchikey_from_smiles, inchikey_from_inchi = get_inchikey(smiles, inchi)

        spectrum_object = {}
        spectrum_object["Name"] = spectrum["Compound_Name"]
        spectrum_object["InChI"] = spectrum["INCHI"]
        spectrum_object["SMILES"] = spectrum["Smiles"]
        spectrum_object["InChIKey_smiles"] = inchikey_from_smiles
        spectrum_object["InChIKey_inchi"] = inchikey_from_inchi
        spectrum_object["spectrum_id"] = spectrum["spectrum_id"]

        all_spectra.append(spectrum_object)

    return all_spectra

indigo = Indigo()
indigo_inchi = IndigoInchi(indigo)

gnps_list = load_GNPS()
npatlas_list = load_NPAtlas("data/NPAtlas_DB_last_version.tsv")
