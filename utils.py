import requests
import json
from rdkit import Chem


def get_inchikey(smiles, inchi):
    inchikey_from_smiles = ""
    inchikey_from_inchi = ""
    try:
        inchikey_from_smiles = str(Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)))
    except:
        inchikey_from_smiles = ""

    try:
        inchikey_from_inchi = str(Chem.InchiToInchiKey(inchi))
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
    print("Loading NPAtlas")
    all_npatlas = json.load(open(filepath, encoding='utf-8', errors='strict'), strict=False)
    print(len(all_npatlas))
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

def gnps_filter_for_key(spectra_list, filterKeysOut=True):
    data_list = load_GNPS()

    new_data_list = []
    for element in data_list:
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
        output_dict["COMPOUND_NAME"] = element["Name"]
        output_dict["COMPOUND_INCHIKEY"] = inchi_key
        output_dict["COMPOUND_INCHI"] = element["inchi"]
        output_dict["COMPOUND_SMILES"] = element["smiles"]
        output_dict["LIBRARY_QUALITY"] = element["Library_Class"]
        
        output_list.append(output_dict)
        
    print("GNPS Library Entries with INCHIKEY", len(output_list))

    return output_list

        