
### DBSCAN-based removal of spatial outliers
## Usage:
# python bin/DBSCAN_sklearn.py input_data.parquet output_data.parquet --epsilon_km 800 --min_samples 3


import polars as pl
import numpy as np
import argparse
from sklearn.cluster import DBSCAN
import logging

## Constants
KMS_PER_RADIAN = 6371.0088

def load_data(file_path):
    logging.info(f"Loading data from {file_path}")
    return pl.read_parquet(file_path)

def process_coordinates(df):
    logging.info("Processing coordinates")
    coords = df.select(['decimallatitude', 'decimallongitude']).to_numpy()
    return np.radians(coords)  # Convert to radians for haversine metric

def perform_dbscan(coords, epsilon, min_samples):
    logging.info(f"Performing DBSCAN with eps={epsilon} and min_samples={min_samples}")
    db = DBSCAN(eps=epsilon, min_samples=min_samples, algorithm='ball_tree', metric='haversine').fit(coords)
    return db.labels_

def remove_outliers(df, labels):
    logging.info("Removing outliers")
    df = df.with_column(pl.Series("cluster", labels))
    retained_data = df.filter(pl.col('cluster') != -1)
    outliers = df.filter(pl.col('cluster') == -1)
    return retained_data.drop('cluster'), outliers

def save_data(df, output_path):
    logging.info(f"Saving cleaned data to {output_path}")
    df.write_parquet(output_path)

def main(input_file, output_file, epsilon_km, min_samples):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    ## Convert epsilon to radians
    epsilon = epsilon_km / KMS_PER_RADIAN
    
    ## Load data
    df = load_data(input_file)

    ## Process coordinates
    coords = process_coordinates(df)

    ## Perform DBSCAN clustering
    labels = perform_dbscan(coords, epsilon, min_samples)

    ## Remove outliers
    retained_data, outliers = remove_outliers(df, labels)
    
    ## Save results
    save_data(retained_data, output_file)

    ## Log results
    logging.info(f"Total records retained: {len(retained_data)}")
    logging.info(f"Total outliers removed: {len(outliers)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DBSCAN outlier removal for geospatial data.")
    parser.add_argument('input_file', type=str, help="Input parquet file path")
    parser.add_argument('output_file', type=str, help="Output parquet file path")
    parser.add_argument('--epsilon_km', type=float, required=True, help="Epsilon distance in kilometers")
    parser.add_argument('--min_samples', type=int, required=True, help="Minimum samples to form a cluster")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.epsilon_km, args.min_samples)

