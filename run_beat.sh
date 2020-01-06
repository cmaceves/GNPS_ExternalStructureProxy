#!/bin/bash

celery -A gnps_tasks worker -l info -B -c 1