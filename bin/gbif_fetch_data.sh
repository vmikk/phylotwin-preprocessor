#!/bin/bash
# - Taxa excluded:
#   - Kingdoms: Chromista (4), Protozoa (7), Viruses (8), incertae sedis (0)
#   - Genera: Homo (2436435)
#   - Species: Felis catus (2435035)

## Usage:
# export GBIF_USERNAME=...
# export GBIF_PASSWORD=...
# export GBIF_EMAIL=...
# ./gbif_fetch_data.sh

## More info on GBIF API:
# https://techdocs.gbif.org/en/data-use/api-downloads
# https://techdocs.gbif.org/en/openapi/v1/occurrence
# https://data-blog.gbif.org/post/gbif-api-beginners-guide/
# https://registry.gbif.org/vocabulary/search

## Output directory
OUTDIR="$(pwd)/GBIF_dumps"
mkdir -p "${OUTDIR}"

# Check if all required variables are set
for var in GBIF_USERNAME GBIF_PASSWORD GBIF_EMAIL OUTDIR; do
    [ -z "${!var}" ] && { echo "$var not set"; exit 1; }
done

## Create a temporary JSON file with the request parameters
echo -e "Creating request file...\n"
JSON=$(mktemp /tmp/gbif.XXXXXXXX.json)
trap "rm ${JSON}" EXIT

cat > "${JSON}" <<EOT
{
    "format": "SIMPLE_PARQUET",
    "creator": "${GBIF_USERNAME}",
    "notification_address": [ "${GBIF_EMAIL}" ],
    "predicate": {
        "type": "and",
        "predicates": [
            {
                "type": "equals",
                "key": "OCCURRENCE_STATUS",
                "value": "PRESENT"
            },
            {
                "type": "equals",
                "key": "HAS_COORDINATE",
                "value": "true"
            },
            {
                "type": "equals",
                "key": "HAS_GEOSPATIAL_ISSUE",
                "value": "false"
            },
            {
                "type": "not",
                "predicate": {
                    "type": "in",
                    "key": "DEGREE_OF_ESTABLISHMENT",
                    "values": [
                        "MANAGED"
                    ]
                }
            },
            {
                "type": "not",
                "predicate": {
                    "type": "in",
                    "key": "BASIS_OF_RECORD",
                    "values": [
                        "FOSSIL_SPECIMEN",
                        "LIVING_SPECIMEN"
                    ]
                }
            },
            {
                "type": "not",
                "predicate": {
                    "type": "in",
                    "key": "ISSUE",
                    "values": [
                        "TAXON_MATCH_HIGHERRANK"
                    ]
                }
            },
            {
                "type": "or",
                "predicates": [
                    {
                        "type": "lessThan",
                        "key": "COORDINATE_UNCERTAINTY_IN_METERS",
                        "value": "100000"
                    },
                    {
                        "type": "isNull",
                        "parameter": "COORDINATE_UNCERTAINTY_IN_METERS"
                    }
                ]
            },
            {
                "type": "not",
                "predicate": {
                    "type": "in",
                    "key": "COORDINATE_UNCERTAINTY_IN_METERS",
                    "values": [
                        "301",
                        "3036",
                        "999",
                        "9999"
                    ]
                }
            },
            {
                "type": "not",
                "predicate": {
                    "type": "in",
                    "key": "KINGDOM_KEY",
                    "values": [
                        "4", "7", "8", "0"
                    ]
                }
            },
            {
                "type": "not",
                "predicate": {
                    "type": "in",
                    "key": "GENUS_KEY",
                    "values": [
                        "2436435"
                    ]
                }
            },
            {
                "type": "not",
                "predicate": {
                    "type": "in",
                    "key": "SPECIES_KEY",
                    "values": [
                        "2435035"
                    ]
                }
            }
        ]
    }
}
EOT


## ESTABLISHMENT_MEANS - exclude non-native species ?
#            {
#                "type": "not",
#                "predicate": {
#                    "type": "in",
#                    "key": "ESTABLISHMENT_MEANS",
#                    "values": [
#                        "MANAGED",
#                        "INTRODUCED",
#                        "INVASIVE",
#                        "NATURALISED"
#                    ]
#                }
#            },



## Send the request to GBIF
echo -e "Sending request to GBIF...\n"
RESPONSE=$(curl -Ssi \
  --user "${GBIF_USERNAME}":"${GBIF_PASSWORD}" \
  -H "Content-Type: application/json" \
  -X POST -d @$JSON \
  https://api.gbif.org/v1/occurrence/download/request \
  )

## Print the GBIF response
echo -e "GBIF response:\n\n${RESPONSE}\n"

## Extract the download key (e.g., "0036526-240906103802322")
ID=$(echo "${RESPONSE}" | tail -n 1)

## Validate the ID
if [[ "${ID}" =~ ^[0-9]+-[0-9]+$ ]]; then
  echo -e "Download key: ${ID} \n"
else
  echo -e "WARNING: Invalid ID! Exiting...\n"
  exit 1
fi

## Download path
DATE=$(date +"%y%m%d")
DUMP="${OUTDIR}/${ID}__${DATE}.zip"

## Check the status of the request
echo -e "Checking request status [each 5 minutes]...\n"

while true
do

  echo -e "Checking status...\n"
  STATUS=$(curl -Ss "https://api.gbif.org/v1/occurrence/download/${ID}" | jq -r ".status")

  if [[ "${STATUS}" == "SUCCEEDED" ]]; then
      echo -e "Downloading to ${DUMP} ...\n"
      
      aria2c \
        https://api.gbif.org/v1/occurrence/download/request/${ID}.zip \
        -d "${OUTDIR}" -o "${ID}__${DATE}.zip" || { echo -e "WARNING: Error downloading file\n"; exit 1; }
      
      break
  fi

  ## Wait for 5 minutes (300 seconds) before retrying
  sleep 300

done

## Verify the download
echo -e "Verifying the download...\n"
7z t "${DUMP}"
teststatus="$?"

# `7z t` Exit Codes:
#   0: No errors, the archive is valid.
#   1: Warnings (some files were skipped or couldn't be tested).
#   2: Fatal errors (archive is corrupted or the command failed).

if [ ${teststatus} -eq 0 ]; then
  echo "Archive is valid."
else
  echo "WARNING: Archive validation failed. Status code: $teststatus"
  exit 1
fi

## Save request parameters to a file
echo -e "Saving request parameters to a file...\n"

METADATA="${OUTDIR}/${ID}__${DATE}.txt"
DOI=$(curl -Ss "https://api.gbif.org/v1/occurrence/download/${ID}" | jq -r ".doi")

echo -e "Download key: ${ID}" >> "${METADATA}"
echo -e "Date: ${DATE}" >> "${METADATA}"
echo -e "DOI: ${DOI}" >> "${METADATA}"
echo -e "\nRequest parameters: \n" >> "${METADATA}"
curl -Ss "https://api.gbif.org/v1/occurrence/download/${ID}" >> "${METADATA}"

echo -e "To cite the dataset use the following DOI: ${DOI} \n"

## Clean up 
rm "${JSON}"

echo -e "All done!\n"

## Unpack the archive
# 7z e "${DUMP}" -o"${DUMP/.zip/}/"
