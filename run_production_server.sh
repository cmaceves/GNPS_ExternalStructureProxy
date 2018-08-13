#!/bin/bash

gunicorn -w 4 -b 0.0.0.0:5010 --timeout 3600 main:app
