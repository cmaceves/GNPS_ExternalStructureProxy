import requests
import json
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

## Caching Results for a specific amount of time
import requests_cache
requests_cache.install_cache('requests_cache', expire_after=86400)


LIBRARY_NAMES = ["GNPS-LIBRARY", 
    "GNPS-SELLECKCHEM-FDA-PART1",
    "GNPS-SELLECKCHEM-FDA-PART2",
    "GNPS-PRESTWICKPHYTOCHEM", 
    "GNPS-NIH-CLINICALCOLLECTION1", 
    "GNPS-NIH-CLINICALCOLLECTION2", 
    "GNPS-NIH-NATURALPRODUCTSLIBRARY", 
    "GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE", 
    "GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_NEGATIVE", 
    "GNPS-NIH-SMALLMOLECULEPHARMACOLOGICALLYACTIVE", 
    "GNPS-FAULKNERLEGACY",
    "GNPS-EMBL-MCF",
    "GNPS-COLLECTIONS-PESTICIDES-POSITIVE",
    "GNPS-COLLECTIONS-PESTICIDES-NEGATIVE",
    "MMV_POSITIVE",
    "MMV_NEGATIVE",
    "LDB_POSITIVE", 
    "LDB_NEGATIVE",
    "GNPS-NIST14-MATCHES",
    "GNPS-COLLECTIONS-MISC",
    "GNPS-MSMLS",
    "BILELIB19",
    "DEREPLICATOR_IDENTIFIED_LIBRARY",
    "PNNL-LIPIDS-POSITIVE",
    "PNNL-LIPIDS-NEGATIVE",
    "MIADB",
    "MASSBANK", 
    "MASSBANKEU", 
    "MONA", 
    "RESPECT", 
    "HMDB", 
    "CASMI",
    "SUMNER"]

def get_inchikey(smiles, inchi):
    inchikey_from_smiles = ""
    inchikey_from_inchi = ""
    try:
        if len(smiles) > 5:
            inchikey_from_smiles = str(Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)))
        else:
            inchikey_from_smiles = ""
    except:
        inchikey_from_smiles = ""

    try:
        if len(inchi) > 5:
            inchikey_from_inchi = str(Chem.InchiToInchiKey(inchi))
        else:
            inchikey_from_inchi = ""
    except:
        inchikey_from_inchi = ""

    if len(inchikey_from_smiles) > 2 and len(inchikey_from_inchi) > 2:
        return inchikey_from_smiles, inchikey_from_inchi

    if len(inchikey_from_smiles) > 2:
        return inchikey_from_smiles, ""

    if len(inchikey_from_inchi) > 2:
        return inchikey_from_inchi, ""

    return "", ""

def get_formula(smiles, inchi):
    formula_from_smiles = ""
    formula_from_inchi = ""
    try:
        if len(smiles) > 5:
            formula_from_smiles = str(CalcMolFormula(Chem.MolFromSmiles(smiles)))
        else:
            formula_from_smiles = ""
    except:
        formula_from_smiles = ""

    try:
        if len(inchi) > 5:
            formula_from_inchi = str(CalcMolFormula(Chem.MolFromInchi(inchi)))
        else:
            formula_from_inchi = ""
    except:
        formula_from_inchi = ""

    if len(formula_from_smiles) > 2 and len(formula_from_inchi) > 2:
        return formula_from_smiles, formula_from_inchi

    if len(formula_from_smiles) > 2:
        return formula_from_smiles, ""

    if len(formula_from_inchi) > 2:
        return formula_from_inchi, ""

    return "", ""


def load_NPAtlas(filepath):
    print("Loading NPAtlas")
    all_npatlas = json.load(open(filepath, encoding='utf-8', errors='strict'), strict=False)
    print(len(all_npatlas))
    return all_npatlas

def load_mibig(filepath):
    df = pd.read_csv(filepath, sep=",")

    output_list = []
    results_list = df.to_dict(orient="records")
    for result in results_list:
        inchikey_from_smiles, inchikey_from_inchi = get_inchikey(result["smiles"], "")

        output_dict = {}
        output_dict["BGCID"] = result["bgc id"]
        output_dict["COMPOUND_INCHIKEY"] = inchikey_from_smiles

        output_list.append(output_dict)

    return output_list

# Loads all GNPS libraries from GNPS servers and returns as python objects
def load_GNPS(library_names=LIBRARY_NAMES):
    all_GNPS_list = []

    for library_name in library_names:
        print(library_name)
        url = "https://gnps.ucsd.edu/ProteoSAFe/LibraryServlet?library=%s" % (library_name)
        all_GNPS_list += requests.get(url).json()["spectra"]

    return all_GNPS_list

    
# Enriching GNPS Library data with structures and formulas
def gnps_format_libraries(all_GNPS_list):
    all_spectra = []
    for i, spectrum in enumerate(all_GNPS_list):
        if i % 1000 == 0:
            print(i, "of", len(all_GNPS_list), file=sys.stderr)

        smiles = spectrum["Smiles"]
        inchi =  spectrum["INCHI"]

        if len(smiles) < 5 and len(inchi) < 5:
            inchikey_from_smiles = ""
            inchikey_from_inchi = ""
            formula_from_smiles = ""
            formula_from_inchi = ""
        else:
            if "InChI=" not in inchi and len(inchi) > 10:
                inchi = "InChI=" + inchi

            inchikey_from_smiles, inchikey_from_inchi = get_inchikey(smiles, inchi)
            formula_from_smiles, formula_from_inchi = get_formula(smiles, inchi)

        spectrum_object = spectrum
        spectrum_object["InChIKey_smiles"] = inchikey_from_smiles
        spectrum_object["InChIKey_inchi"] = inchikey_from_inchi
        spectrum_object["Formula_smiles"] = formula_from_smiles
        spectrum_object["Formula_inchi"] = formula_from_inchi
        spectrum_object["url"] = "https://gnps.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID=%s" % spectrum["spectrum_id"]

        all_spectra.append(spectrum_object)

    return all_spectra

# Reformatting output for NPAtlas
def gnps_filter_for_key(formatted_spectra_list, filterKeysOut=True):
    new_data_list = []
    for element in formatted_spectra_list:
        if len(element["InChIKey_inchi"]) < 5 and len(element["InChIKey_smiles"]) < 5 and filterKeysOut:
            continue
        new_data_list.append(element)

    #Finding Discordant Entries
    discordant_list = [element for element in new_data_list \
                   if (len(element["InChIKey_inchi"]) == len(element["InChIKey_smiles"]) \
                       and element["InChIKey_inchi"].split("-")[0] != element["InChIKey_smiles"].split("-")[0])]

    print("Discordant List", len(discordant_list))

    #Formatting
    output_list = []
    for element in new_data_list:
        inchi_key = element["InChIKey_inchi"]
        if len(element["InChIKey_smiles"]) > len(inchi_key):
            inchi_key = element["InChIKey_smiles"]
        
        output_dict = {}
        output_dict["GNPSID"] = element["spectrum_id"]
        output_dict["COMPOUND_NAME"] = element["Compound_Name"]
        output_dict["COMPOUND_INCHIKEY"] = inchi_key
        output_dict["COMPOUND_INCHI"] = element["INCHI"]
        output_dict["COMPOUND_SMILES"] = element["Smiles"]
        output_dict["LIBRARY_QUALITY"] = element["Library_Class"]
        
        output_list.append(output_dict)
        
    print("GNPS Library Entries with INCHIKEY", len(output_list))

    return output_list

# Getting all the spectrum peaks for the library spectrum
def get_gnps_peaks(all_GNPS_list):
    import copy

    output_list = []
    for i, spectrum in enumerate(all_GNPS_list):
        if i % 1000 == 0:
            print(i, "of", len(all_GNPS_list), file=sys.stderr)

        new_spectrum = copy.deepcopy(spectrum)

        try:
            spectrum_peaks_url = "https://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID={}".format(spectrum["spectrum_id"])
            r = requests.get(spectrum_peaks_url)
            spectrum_json = r.json()
            new_spectrum["peaks_json"] = spectrum_json["spectruminfo"]["peaks_json"]
        
            output_list.append(new_spectrum)
        except KeyboardInterrupt:
            raise
        except:
            continue

    return output_list