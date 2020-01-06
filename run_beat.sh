#!/bin/bash

celery -A gnps_tasks worker -l debug -B