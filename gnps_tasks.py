from celery import Celery
import os
import json
import requests
import utils
import pandas as pd

celery_instance = Celery('tasks', backend='rpc://externalstructureproxy-rabbitmq', broker='pyamqp://externalstructureproxy-rabbitmq')

@celery_instance.task()
def generate_gnps_data():
    gnps_libraries = utils.load_GNPS()

    print("Got all Libraries")

    with open("/output/gnpslibraries.json", "w") as output_file:
        output_file.write(json.dumps(gnps_libraries))

    formatted_gnps_libraries = utils.gnps_format_libraries(gnps_libraries)

    with open("/output/gnpslibraries_enriched_all.json", "w") as output_file:
        output_file.write(json.dumps(utils.gnps_filter_for_key(formatted_gnps_libraries, filterKeysOut=False)))


    #Outputting for NPAtlas
    with open("/output/gnpslibraries_npatlas.json", "w") as output_file:
        output_file.write(json.dumps(utils.gnps_filter_for_key(formatted_gnps_libraries, filterKeysOut=True)))

    pd.DataFrame(utils.gnps_filter_for_key(formatted_gnps_libraries)).to_csv("/output/gnpslibraries_npatlas.tsv", sep="\t", index=False)

    # Getting spectrum peaks for each library spectrum
    for library_spectrum in formatted_gnps_libraries:
        print(library_spectrum)

    # Creating JSON files for each GNPS Library
    # for library in utils.LIBRARY_NAMES:
    #     print("library")

celery_instance.conf.beat_schedule = {
    "generate_gnps_data": {
        "task": "gnps_tasks.generate_gnps_data",
        "schedule": 30.0
    }
}
