
### DBSCAN-based removal of spatial outliers
## Usage:
# python bin/DBSCAN_sklearn.py input_data.parquet output_data.parquet --epsilon_km 800 --min_samples 3


import polars as pl
import numpy as np
import argparse
from sklearn.cluster import DBSCAN, HDBSCAN, OPTICS
import logging

## Constants
KMS_PER_RADIAN = 6371.0088

def load_data(file_path):
    logging.info(f"Loading data from {file_path}")
    try:
        return pl.read_parquet(file_path)
    except Exception as e:
        logging.error(f"Error loading data from {file_path}: {str(e)}")
        raise

def process_coordinates(df):
    logging.info("Processing coordinates")
    coords = df.select(['decimallatitude', 'decimallongitude']).to_numpy()
    return np.radians(coords)  # Convert to radians for haversine metric

## DBSCAN
def perform_dbscan(coords, epsilon, min_samples, algorithm='ball_tree', n_jobs=1):
    logging.info(f"Performing DBSCAN with eps={epsilon} and min_samples={min_samples}")
    db = DBSCAN(eps=epsilon, min_samples=min_samples, algorithm=algorithm, metric='haversine', n_jobs=n_jobs).fit(coords)
    return db.labels_

## HDBSCAN
def perform_hdbscan(coords, min_samples, min_cluster_size, algorithm='ball_tree', n_jobs=1):
    logging.info(f"Performing HDBSCAN with min_samples={min_samples}, min_cluster_size={min_cluster_size}")
    hdb = HDBSCAN(min_samples=min_samples, min_cluster_size=min_cluster_size, metric='haversine', algorithm=algorithm, core_dist_n_jobs=n_jobs).fit(coords)
    return hdb.labels_

## OPTICS
def perform_optics(coords, min_samples, max_eps, algorithm='ball_tree', n_jobs=1):
    logging.info(f"Performing OPTICS with max_eps={max_eps} and min_samples={min_samples}")
    optics = OPTICS(max_eps=max_eps, min_samples=min_samples, algorithm=algorithm, metric='haversine', n_jobs=n_jobs).fit(coords)
    return optics.labels_


def remove_outliers(df, labels):
    logging.info("Removing outliers")
    labels_series = pl.Series("cluster", labels)
    df = df.with_columns([labels_series])
    retained_data = df.filter(pl.col('cluster') != -1)
    outliers = df.filter(pl.col('cluster') == -1)
    return retained_data.drop('cluster'), outliers

def save_data(df, output_path):
    logging.info(f"Saving cleaned data to {output_path}")
    try:
        df.write_parquet(output_path)
    except Exception as e:
        logging.error(f"Error saving data to {output_path}: {str(e)}")
        raise

def main(input_file, output_file, method, epsilon_km, min_samples, min_cluster_size, algorithm, threads):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Add input validation for epsilon_km
    if method in ['DBSCAN', 'OPTICS'] and epsilon_km is None:
        raise ValueError(f"epsilon_km must be provided for {method}")
    
    # Convert epsilon to radians
    epsilon = epsilon_km / KMS_PER_RADIAN if epsilon_km else None
    
    ## Load data
    df = load_data(input_file)

    ## Process coordinates
    coords = process_coordinates(df)

    ## Perform clustering based on the chosen method
    if method == 'DBSCAN':
        labels = perform_dbscan(coords, epsilon, min_samples, algorithm, threads)
    elif method == 'HDBSCAN':
        labels = perform_hdbscan(coords, min_samples, min_cluster_size, algorithm, threads)
    elif method == 'OPTICS':
        max_eps = epsilon if epsilon is not None else np.inf
        labels = perform_optics(coords, min_samples, max_eps, algorithm, threads)
    else:
        raise ValueError(f"Unsupported clustering method: {method}")

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
    parser.add_argument('--method', type=str, choices=['DBSCAN', 'HDBSCAN', 'OPTICS'], required=True, help="Clustering method to use")
    parser.add_argument('--epsilon_km', type=float, default=None, help="Epsilon distance in kilometers (used by DBSCAN/OPTICS)")
    parser.add_argument('--min_samples', type=int, required=True, help="Minimum samples to form a cluster")
    parser.add_argument('--min_cluster_size', type=int, default=5, help="Minimum cluster size (used by HDBSCAN)")
    parser.add_argument('--algorithm', type=str, default='ball_tree', help="Algorithm to use for the nearest neighbour search (default, 'ball_tree'; alternatively, 'kd_tree' or 'brute')")
    parser.add_argument('--threads', type=int, default=1, help="Number of threads to use for parallel processing")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.method, args.epsilon_km, args.min_samples, args.min_cluster_size, args.algorithm, args.threads)

