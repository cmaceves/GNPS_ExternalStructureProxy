import sys
sys.path.insert(0, "..")

def test():
    import utils
    spectra_list = utils.load_GNPS(library_names=["GNPS-LIBRARY"])

    print(len(spectra_list))