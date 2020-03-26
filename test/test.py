import sys
import json
sys.path.insert(0, "..")

def test():
    import utils
    spectra_list = utils.load_GNPS(library_names=["GNPS-LIBRARY"])[:100]
    spectra_list = utils.gnps_format_libraries(spectra_list)

    with open("output_enriched_list.json", "w") as output_file:
        output_file.write(json.dumps(spectra_list, indent=4))

    filtered_list = utils.gnps_filter_for_key(spectra_list)

    with open("output_filtered_list.json", "w") as output_file:
        output_file.write(json.dumps(filtered_list, indent=4))
    
    print(len(filtered_list))

def test_get_library_peaks():
    import utils
    spectra_list = utils.load_GNPS(library_names=["GNPS-LIBRARY"])[:100]
    spectra_list = utils.gnps_format_libraries(spectra_list)
    spectra_list_with_peaks = utils.get_gnps_peaks(spectra_list)

    print(len(spectra_list_with_peaks))

    with open("output_enriched_list_peaks.json", "w") as output_file:
        output_file.write(json.dumps(spectra_list_with_peaks, indent=4))

    mgf_string = utils.get_full_mgf_string(spectra_list_with_peaks)
    with open("output_library.mgf", "w") as output_file:
        output_file.write(mgf_string)