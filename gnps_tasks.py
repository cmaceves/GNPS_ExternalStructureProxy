from celery import Celery
import os
import json
import requests
import utils
import pandas as pd

celery_instance = Celery('tasks', backend='rpc://externalstructureproxy-rabbitmq', broker='pyamqp://externalstructureproxy-rabbitmq')

@celery_instance.task()
def generate_gnps_data():
    libraries_json = utils.load_GNPS()

    print("Got all Libraries")

    with open("/output/gnpslibraries.json", "w") as output_file:
        output_file.write(json.dumps(libraries_json))

    # with open("/output/gnpslibraries_all_formated.json", "w") as output_file:
    #     output_file.write(json.dumps(utils.gnps_filter_for_key(libraries_json, filterKeysOut=False)))

    # with open("/output/gnpslibraries_withkeys_formated.json", "w") as output_json:
    #     output_file.write(json.dumps(utils.gnps_filter_for_key(libraries_json, filterKeysOut=True)))

    # pd.DataFrame(utils.gnps_filter_for_key(libraries_json)).to_csv("/output/gnpslibraries_withkeys_formated.tsv", sep="\t", index=False)

celery_instance.conf.beat_schedule = {
    "generate_gnps_data": {
        "task": "gnps_tasks.generate_gnps_data",
        "schedule": 30.0
    }
}
