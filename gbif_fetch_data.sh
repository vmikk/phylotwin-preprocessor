#!/bin/bash

## Usage:
# export GBIF_USERNAME=...
# export GBIF_PASSWORD=...
# export GBIF_EMAIL=...
# ./gbif_fetch_data.sh

## More info on GBIF API: https://data-blog.gbif.org/post/gbif-api-beginners-guide/

## Output directory
OUTDIR="$(pwd)/GBIF_dumps"
mkdir -p $OUTDIR

# Check if all required variables are set
for var in GBIF_USERNAME GBIF_PASSWORD GBIF_EMAIL OUTDIR; do
    [ -z "${!var}" ] && { echo "$var not set"; exit 1; }
done

## Create a temporary JSON file with the request parameters
echo -e "Creating request file...\n"
JSON=$(mktemp /tmp/gbif.XXXXXXXX.json)
trap "rm $JSON" EXIT

cat > $JSON <<EOT
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
                "key": "HAS_GEOSPATIAL_ISSUE",
                "value": "false"
            },
            {
                "type": "not",
                "predicate": {
                    "type": "in",
                    "key": "ESTABLISHMENT_MEANS",
                    "values": [
                        "MANAGED",
                        "INTRODUCED",
                        "INVASIVE",
                        "NATURALISED"
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
            }
        ]
    }
}
EOT

## Send the request to GBIF
echo -e "Sending request to GBIF...\n"
RESPONSE=$(curl -Ssi \
  --user "${GBIF_USERNAME}":"${GBIF_PASSWORD}" \
  -H "Content-Type: application/json" \
  -X POST -d @$JSON \
  https://api.gbif.org/v1/occurrence/download/request \
  )

## Print the GBIF response
echo -e "GBIF response:\n\n$RESPONSE\n"

## Extract the download key (e.g., "0036526-240906103802322")
ID=$(echo "$RESPONSE" | tail -n 1)

## Validate the ID
if [[ "$ID" =~ ^[0-9]+-[0-9]+$ ]]; then
  echo -e "Download key: $ID \n"
else
  echo -e "Invalid ID! Exiting...\n"
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
  STATUS=$(curl -Ss https://api.gbif.org/v1/occurrence/download/$ID | jq -r ".status")

  if [[ "$STATUS" == "SUCCEEDED" ]]; then
      echo -e "Downloading to $DUMP ...\n"
      
      aria2c \
        https://api.gbif.org/v1/occurrence/download/request/$ID.zip \
        -d "${OUTDIR}" -o "${ID}__${DATE}.zip" || { echo -e "WARNING: Error downloading file\n"; exit 1; }
      
      break
  fi

  ## Wait for 5 minutes (300 seconds) before retrying
  sleep 300

done

## Verify the download
echo -e "Verifying the download...\n"
7z t $DUMP

echo -e "Done!\n"