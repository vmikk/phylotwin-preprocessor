#!/usr/bin/env python3
"""
Fetch extinct taxon IDs from GBIF API (with pagination support and multiple taxon keys).

Example usage:
    python fetch_gbif_extinct_taxa.py --outfile extinct_taxa.txt
"""

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List

import requests
from requests.exceptions import RequestException

# Configure logging with red color for errors
RED = "\033[91m"
RESET = "\033[0m"
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s'
)
logger = logging.getLogger(__name__)

## Error formatting function
def log_error(message: str) -> None:
    """Log error messages in red."""
    logger.error(f"{RED}{message}{RESET}")

def log_warning(message: str) -> None:
    """Log warning messages in red."""
    logger.warning(f"{RED}{message}{RESET}")

@dataclass
class GBIFConfig:
    """GBIF API configuration and search parameters."""
    BASE_URL: str = "https://api.gbif.org/v1/species/search"
    DATASET_KEY: str = "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"
    BATCH_SIZE: int = 1000
    MAX_OFFSET: int = 100_000
    HIGHER_TAXON_KEYS: List[str] = (
        "359"  # Mammalia
    )

class GBIFClient:
    """Client for interacting with GBIF API."""
    
    def __init__(self, highertaxon_key: str):
        self.config = GBIFConfig()
        self.params = {
            "rank": "SPECIES",
            "highertaxon_key": highertaxon_key,
            "status": "ACCEPTED",
            "datasetKey": self.config.DATASET_KEY,
            "isExtinct": "true",
        }

    def _make_request(self, limit: int, offset: int) -> dict:
        """Make API request with error handling."""
        try:
            params = {**self.params, "limit": str(limit), "offset": str(offset)}
            response = requests.get(self.config.BASE_URL, params=params)
            response.raise_for_status()
            return response.json()
        except RequestException as e:
            log_error(f"API request failed: {e}")
            raise

    def get_total_count(self) -> int:
        """Get total number of records available."""
        data = self._make_request(limit=0, offset=0)
        return data["count"]

    def fetch_records(self) -> Iterator[List[dict]]:
        """Fetch records in batches."""
        total_count = self.get_total_count()
        logger.info(f"\tFound {total_count} records")

        if total_count > self.config.MAX_OFFSET:
            log_warning(
                f"\tWarning: Total records ({total_count}) exceed max offset limit "
                f"({self.config.MAX_OFFSET}). Some records will be missed."
            )

        offset = 0
        while offset < min(total_count, self.config.MAX_OFFSET):
            logger.info(f"\tFetching batch {offset//self.config.BATCH_SIZE + 1}...")
            data = self._make_request(limit=self.config.BATCH_SIZE, offset=offset)
            
            if not data["results"]:
                break

            yield data["results"]
            offset += self.config.BATCH_SIZE

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Fetch extinct taxon IDs from GBIF API"
    )
    parser.add_argument(
        "--outfile",
        type=Path,
        default="extinct_taxa.txt",
        help="Output file path (default: `extinct_taxa.txt`)"
    )
    return parser.parse_args()

def format_record(record: dict, highertaxon_key: str) -> str:
    """Format a single record for output."""
    fields = [
        str(record.get("key", "")),
        highertaxon_key,
        record.get("kingdom", ""),
        record.get("phylum", ""),
        record.get("order", ""),
        record.get("family", ""),
        record.get("genus", ""),
        record.get("species", "")
    ]
    return "\t".join(fields)

def main() -> None:
    """Main execution function."""
    args = parse_args()
    
    try:
        with open(args.outfile, "w") as f:
            header = ["specieskey", "highertaxon_key", "kingdom", "phylum", 
                     "order", "family", "genus", "species"]
            f.write("\t".join(header) + "\n")
            
            for taxon_key in GBIFConfig.HIGHER_TAXON_KEYS:
                try:
                    logger.info(f"\nProcessing taxon key: {taxon_key}")
                    client = GBIFClient(taxon_key)
                    
                    for batch in client.fetch_records():
                        formatted_records = [format_record(record, taxon_key) for record in batch]
                        f.write("\n".join(formatted_records) + "\n")
                
                except Exception as e:
                    log_error(f"\tFailed to process taxon key {taxon_key}: {e}")
                    continue
    except IOError as e:
        log_error(f"Failed to write to file: {e}")
        raise

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log_error(f"Script failed: {e}")
        raise
