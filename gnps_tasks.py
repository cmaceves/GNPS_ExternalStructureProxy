from celery import Celery
import os
import json
import requests
import utils

celery_instance = Celery('tasks', backend='rpc://externalstructureproxy-rabbitmq', broker='pyamqp://externalstructureproxy-rabbitmq')

@celery_instance.task()
def generate_gnps_data():
    utils.load_GNPS()

celery_instance.conf.beat_schedule = {
    "generate_gnps_data": {
        "task": "gnps_tasks.generate_gnps_data",
        "schedule": 30.0
    }
}
