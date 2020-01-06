import requests

PRODUCTION_URL = "gnps-external.ucsd.edu"

def test_production():
    url = f"https://{PRODUCTION_URL}/heartbeat"
    r = requests.get(url)
    r.raise_for_status()


def test_gnps_library():
    url = f"https://{PRODUCTION_URL}/gnpslibraryjson"
    r = requests.get(url)
    r.raise_for_status()

    url = f"https://{PRODUCTION_URL}/gnpslibraryformattedjson"
    r = requests.get(url)
    r.raise_for_status()

    url = f"https://{PRODUCTION_URL}/gnpslibraryfornpatlasjson"
    r = requests.get(url)
    r.raise_for_status()

    url = f"https://{PRODUCTION_URL}/gnpslibraryfornpatlastsv"
    r = requests.get(url)
    r.raise_for_status()


