"""
This is a hand-made list of dicts. Each entry is a DFO tide station in the
LiveOcean cas7 domain.

It was created by copying entries from the json generated using the nifty webpage:
https://api.iwls-sine.azure.cloud-nuage.dfo-mpo.gc.ca/swagger-ui/index.html
then using the Stations API tool for api/v1/stations and using PAC in the
chs-region-code. The documentation suggests I could get the same result from
this url:
https://api.iwls-sine.azure.cloud-nuage.dfo-mpo.gc.ca/api/v1/stations?chs-region-code=PAC
There are other API calls there that would get info for specific stations.

The "code" keys correspond to those I had from a few years ago in the ptools era,
except with a "0" at the start.

An example url to access hourly data looks like:
https://api.iwls-sine.azure.cloud-nuage.dfo-mpo.gc.ca/api/v1/stations/5cebf1de3d0f4a073c4bb94c/data?time-series-code=wlo&from=2022-01-01T00%3A00%3A00Z&to=2022-01-31T23%3A00%3A00Z&resolution=SIXTY_MINUTES

There are limitations. Apparently for hourly data you can only get a month at a time, and there
seem to be limits on what years are available. For example 2020 did not work with Point Atkinson
but 2022 did.

The reponse is a json with entries like:
[
  {
    "eventDate": "2022-01-01T00:00:00Z",
    "qcFlagCode": "1",
    "value": 3.935,
    "timeSeriesId": "5ddeee612a8a340001a3676f"
  },
  {
    "eventDate": "2022-01-01T01:00:00Z",
    "qcFlagCode": "1",
    "value": 3.277,
    "timeSeriesId": "5ddeee612a8a340001a3676f"
  },

"""
dfo_info_list = [
    {"id": "5cebf1de3d0f4a073c4bb94c",
    "code": "07795",
    "officialName": "Point Atkinson",
    "alternativeName": "Caulfeild Cove",
    "latitude": 49.337,
    "longitude": -123.253,},

    {"id": "5cebf1de3d0f4a073c4bb943",
    "code": "07735",
    "officialName": "Vancouver",
    "latitude": 49.2863,
    "longitude": -123.0997,},

    {"id": "5cebf1df3d0f4a073c4bbd2d",
    "code": "07277",
    "officialName": "Patricia Bay",
    "latitude": 48.6536,
    "longitude": -123.4515,},

    {"id": "5cebf1df3d0f4a073c4bbd1e",
    "code": "07120",
    "officialName": "Victoria Harbour",
    "latitude": 48.424666,
    "longitude": -123.3707,},

    {"id": "5cebf1e23d0f4a073c4bc062",
    "code": "08545",
    "officialName": "Bamfield",
    "latitude": 48.836,
    "longitude": -125.136,},

    {"id": "5cebf1e23d0f4a073c4bc07c",
    "code": "08615",
    "officialName": "Tofino",
    "alternativeName": "Clayoquot , CLAYOQUOT ",
    "latitude": 49.154,
    "longitude": -125.913,},

    {"id": "5cebf1de3d0f4a073c4bb996",
    "code": "08074",
    "officialName": "Campbell River",
    "latitude": 50.042,
    "longitude": -125.247,},

    {"id": "5cebf1e13d0f4a073c4bbf93",
    "code": "07654",
    "officialName": "New Westminster",
    "latitude": 49.2,
    "longitude": -122.91,}

]