# views.py
from flask import abort, jsonify, render_template, request, redirect, url_for

from app import app

import json
import csv
import requests
import requests_cache
from indigo import *
from indigo_renderer import *
from indigo_inchi import *

requests_cache.install_cache('demo_cache')

@app.route('/', methods=['GET'])
def homepage():
    return "MING"

@app.route('/npatlasproxy', methods=['GET'])
def npatlasproxy():
    inchi = request.args.get('inchi', '')
    inchikey = request.args.get('inchikey', '')
    smiles = request.args.get('smiles', '')

    inchikey_from_smiles, inchikey_from_inchi = get_inchikey(smiles, inchi)
    acceptable_key = set([inchikey.split("-")[0], inchikey_from_smiles.split("-")[0], inchikey_from_inchi.split("-")[0]])

    NPAID = None

    for npatlas_entry in npatlas_list:
        if len(npatlas_entry["InChIKey"]) > 2 and npatlas_entry["InChIKey"].split("-")[0] in acceptable_key:
            NPAID = npatlas_entry["NPAID"]
            break

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

def prep_external(results_list, resource_name, resource_url):
    external_links = []

    for result in results_list:
        external_obj = {}
        external_obj["resource"] = resource_name
        external_obj["url"] = resource_url % (result)

        external_links.append(external_obj)

    return external_links

@app.route('/structureproxy', methods=['GET'])
def structureproxy():
    inchi = request.args.get('inchi', '')
    inchikey = request.args.get('inchikey', '')
    smiles = request.args.get('smiles', '')

    inchikey_from_smiles, inchikey_from_inchi = get_inchikey(smiles, inchi)

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

    if len(inchikey_from_smiles) > 2 and len(inchikey_from_inchi) > 2:
        return inchikey_from_smiles, inchikey_from_inchi

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
    library_names = ["all", "MASSBANK", "MASSBANKEU", "MONA", "RESPECT", "HMDB", "CASMI"]

    all_GNPS_list = []

    for library_name in library_names:
        print(library_name)
        url = "https://gnps.ucsd.edu/ProteoSAFe/LibraryServlet?library=%s" % (library_name)
        all_GNPS_list += requests.get(url).json()["spectra"]

    all_spectra = []
    for spectrum in all_GNPS_list:
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
        spectrum_object["url"] = "https://gnps.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID=%s" % spectrum["spectrum_id"]

        all_spectra.append(spectrum_object)

    return all_spectra

indigo = Indigo()
indigo_inchi = IndigoInchi(indigo)

#gnps_list = load_GNPS()
gnps_list = []
npatlas_list = load_NPAtlas("data/NPAtlas_DB_last_version.tsv")
