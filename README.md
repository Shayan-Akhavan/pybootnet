# pybootnet 
This repository contains a set of python functions to perform network analysis with bootstrapped data. The functions are designed to work with Pandas DataFrames and can be used to analyze and compare network statistics across different projects.  
Undirected networks are visualized by their correlations (postive or negative). Parameters are set to filter by their correlation strength.  

## Installing pybootnet
To install using conda, we recommend creating a fresh environment using the following conda install command with all the requirements:
```
conda create -n pybootnet jupyter matplotlib=3.6.3 networkx=3.0 numpy=1.23.5 pandas=2.0.2 scipy=1.13.0
```

To install using pip, in the directory of requirements.txt run:
```
pip install -r requirements.txt
```

## Testing pybootnet with example data
The gut_analysis_notebook.ipynb preforms a comprehsensive network analysis of two different 16S datasets, the results of which appear in the PyBootNet manuscript. This notebook demonstrates the use of the PyBootNet functions on dataset obtained from https://github.com/bryansho/PCOS_WGS_16S_metabolome and Ho et al. (2021) mSystems 6:https://doi.org/10.1128/msystems.01149-20.

## Functions
1. `map_columns(df, column)`: Changes the values in a specified column of a DataFrame to be more readable, mapping them to 'X1', 'X2', 'X3', etc. This is optional, but recommended if target variables contain large names.  
2. `bootstrap_replicates(df, n_iterations=100)`: Creates bootstrap replicates of the input data. The default number of iterations is 100.  
3. `correlation_matrix(bootstrap_replicates)`: Calculates the Spearman correlation matrix for each DataFrame in the input list, taking in only numerical values.  
4. `calculate_network_statistics(correlation_matrices, corr_threshold=0.8)`: Calculates various network statistics for each correlation matrix in the input list, using a specified correlation threshold.  
5. `analyze_network_statistics(network_stats, filename='network_stats.csv', project_name='', format='svg')`: Analyzes the network statistics for different projects, calculates descriptive statistics, and creates box plots for each statistic.  
6. `build_network_graph(correlation_matrices, threshold=0, title="Correlation Network")`: Builds a network graph from the input correlation matrices, averaging them if multiple matrices are provided.  
7. `top_nodes(correlation_matrices, threshold=0.8, num_nodes=20)`: Identifies the top nodes with the highest degree in the network graph.  
8. `most_connected_nodes(network_graph)`: Finds the most connected node in the input graph built from the saved `build_network_graph` output.  
9. `nodes_edges_table(network_graph)`: Creates a table of the number of edges for each node from the saved `build_network_graph` output.  
10. `save_table_to_csv(df, filename)`: Saves a DataFrame to a CSV file. If there are any data pre processed that you would like to save and import later.
11. `net_stat_binomial_test(network_stats_1, network_stats_2, title='p_values.csv')`: Performs a binomial test on the network statistics of two sets of bootstrap replicates. Takes the outputs of the `calculate_network_statistics`.  
  
## Usage   
1. Import `requirements.txt` in order to run the functions under this repository.  
2. Prepare your data as a Pandas DataFrame or a list of DataFrames. Make sure to preprocess your data so that the target variables you want to visualize in the network are in the columns of the DataFrame. For example, if you're analyzing correlations between different taxa, each taxon should be represented as a separate column in the DataFrame. 
3. Use the `map_columns` function to make the column values more readable if needed.  
4. Create bootstrap replicates of your data using the `bootstrap_replicates` function.  
5. Calculate the correlation matrices for each bootstrap replicate using the `correlation_matrix` function.  
6. Calculate network statistics for each correlation matrix using the `calculate_network_statistics` function.  
7. Analyze and compare the network statistics across different projects using the `analyze_network_statistics` function.  
8. Build a network graph from the correlation matrices using the `build_network_graph` function.  
9. Identify the top nodes in the network graph using the `top_nodes` function.  
10. Find the most connected node using the `most_connected_nodes` function.  
11. Create a table of the number of edges for each node using the `nodes_edges_table` function.  
12. Save the table to a CSV file using the `save_table_to_csv` function.  
13. Perform a binomial test on the network statistics of two sets of bootstrap replicates using the `net_stat_binomial_test` function.  
Make sure to adjust the input parameters and file names according to your specific use case. The functions provide flexibility in terms of input data format and output file formats, allowing for easy integration into your analysis pipeline.  

## Supplemental Functions

1. `build_positive_network(correlation_matrices, threshold=0, title="Positive Correlation Network")`: Builds a network graph from the input correlation matrices, focusing on positive correlations above the specified threshold.  
2. `build_filtered_networks(correlation_matrices, threshold=0.8, max_degree=7, title="filtered network")`: Builds a filtered network graph from the input correlation matrices, keeping only nodes with a maximum degree (number of edges) specified by `max_degree`.  
3. `build_negative_networks(correlation_matrices, threshold=0, title='Project')`: Builds a network graph from the input correlation matrices, focusing on negative correlations below the specified threshold.  
4. `negative_filtered_networks(correlation_matrices, threshold, max_degree, title='Project')`: Builds a filtered network graph of negative correlations, keeping only nodes with a maximum degree specified by `max_degree`.  
5. `bootstrap_sample_with_correlation(df, n_iterations)`: Performs bootstrapping on the input data, calculates correlations for each bootstrap replicate, and returns the average correlation matrix. Outputs a matrix, not a list of matrices.  
6. `top_nodes_network_graph(correlation_matrices, threshold=0.8, num_nodes=20, title="Correlation Network")`: Builds a network graph from the input correlation matrices, highlighting the top nodes with the highest degree centrality. The number of top nodes to highlight is specified by `num_nodes`.  



### List of functions:  
map_columns  
bootstrap_replicates  
correlation_matrix  
calculate_network_statistics  
analyze_network_statistics  
average_network_stats  
build_network_graph  
top_nodes_network_graph  
top_nodes  
most_connected_nodes  
nodes_edges_table  
net_stat_binomial_test  


### Supplemental functions:  
calculate_network_statistics_negative  
no_label_network_graph   
build_positive_network  
build_filtered_network  
build_negative networks  
negative_filtered_networks  
bootstrap_sample_with_correlation  
save_table_to_csv  


## Additional Example 
1. Download `binned_species_read_counts_GTDBTk_BLAST_clr.csv`
2. Run `mouse_data_preprocessing.ipynb`
3. Make sure pybootnet is in the same directory
4. Download and run `mouse_analysis.ipynb` 
