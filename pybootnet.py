from itertools import combinations
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy
from scipy import stats


# Change a column to be readable
def map_columns(df, column):
    # Changing the Column's into X1, X2, X3...
    string_map = {}
    counter = 1

    for value in df[column].unique():
        string_map[value] = f'X{counter}'
        counter += 1

    # Map the values in the column using the dictionary
    df[column] = df[column].map(string_map)
    return df


# Creating Boostrap Replicates
def bootstrap_replicates(data, n_iterations=100): # Default 100
    n = len(data)
    bootstrap_results = []

    for i in range(n_iterations):
        sample = data.sample(n, replace=True)
        bootstrap_results.append(sample)

    return bootstrap_results


# Creates a list of each bootstrapped correlation matrix
def correlation_matrix(data_list):
    matrix_list = []

    for data in data_list:
        correlations = data.corr(method='spearman', numeric_only=True)  # Numeric values only
        matrix_list.append(correlations)

    return matrix_list


# Filter correlations where taxa A correlates with taxa B, ignore taxa A to A correlations
# input two datasets, correlate with each other. Same number of samples


# Network stats of individual bootstrap replicates
# n dimensions, visualizes 2-D
# 1 or 3 dimensional 
def calculate_network_statistics(matrices, corr_threshold=0.8):
    network_stats = []

    for matrix in matrices:
        # Initialize a new graph for each matrix
        G = nx.Graph()
        matrix = np.array(matrix)

        n = matrix.shape[0]
        for i, j in combinations(range(n), 2):  # Iterates over unique pairs
            if (abs(matrix[i, j]) > corr_threshold).any():
                G.add_edge(i, j, weight=matrix[i, j])

        # Calculate statistics
        num_edges = G.number_of_edges()
        num_nodes = G.number_of_nodes()
        degree_centrality = nx.degree_centrality(G)
        transitivity = nx.transitivity(G) if num_edges > 0 else 0
        closeness_centrality = nx.closeness_centrality(G)
        average_closeness_centrality = sum(closeness_centrality.values()) / num_nodes if num_nodes > 0 else 0
        betweeness_centrality = nx.betweenness_centrality(G) if num_edges > 0 else 0
        average_betweeness_centrality = sum(betweeness_centrality.values()) / num_nodes if num_nodes > 0 else 0
        density = nx.density(G)

        # Store the stats
        stats_dict = {
            'Number of Edges': round(num_edges, 3),
            'Number of Nodes': round(num_nodes, 3),
            'Average Degree Centrality': round(sum(degree_centrality.values()) / num_nodes if num_nodes > 0 else 0, 3),
            'Transitivity': round(transitivity, 3),
            'Closeness Centrality': round(average_closeness_centrality, 3),
            'Betweeness Centrality': round(average_betweeness_centrality, 5),
            'Density': round(density, 3)
        }
        network_stats.append(stats_dict)

    return network_stats


# Network stats of individual bootstrap replicates NEGATIVE correlations
# n dimensions, visualizes 2-D
# 1 or 3 dimensional 
# Defaults to negative correlations less than 0
def calculate_network_statistics_negative(matrices, corr_threshold=0):
    network_stats = []

    for matrix in matrices:
        # Initialize a new graph for each matrix
        G = nx.Graph()
        matrix = np.array(matrix)

        n = matrix.shape[0]
        for i, j in combinations(range(n), 2):  # Iterates over unique pairs
            if matrix[i, j] < corr_threshold:
                G.add_edge(i, j, weight=matrix[i, j])

        # Calculate statistics
        num_edges = G.number_of_edges()
        num_nodes = G.number_of_nodes()
        degree_centrality = nx.degree_centrality(G)
        transitivity = nx.transitivity(G) if num_edges > 0 else 0
        closeness_centrality = nx.closeness_centrality(G)
        average_closeness_centrality = sum(closeness_centrality.values()) / num_nodes if num_nodes > 0 else 0
        betweeness_centrality = nx.betweenness_centrality(G) if num_edges > 0 else 0
        density = nx.density(G)

        # Store the stats
        stats_dict = {
            'Number of Edges': round(num_edges, 3),
            'Number of Nodes': round(num_nodes, 3),
            'Average Degree Centrality': round(sum(degree_centrality.values()) / num_nodes if num_nodes > 0 else 0, 3),
            'Transitivity': round(transitivity, 3),
            'Closeness Centrality': round(average_closeness_centrality, 3),
            'Betweeness Centrality': round(betweeness_centrality, 5),
            'Density': round(density, 3)
        }
        network_stats.append(stats_dict)

    return network_stats


# INPUT DATA IS A DICTIONARY COMPARING DIFFERENT PROJECTS 
# WILL OUTPUT A TABLE OF THE MEAN AND STD OF EACH NETWORK STATISTICS AS A COLUMN
def analyze_network_statistics(project_stats, filename='network_stats.csv', project_name='', format='svg'):
    # Initialize an empty DataFrame to store the results
    results = pd.DataFrame()

    # Create a dictionary to store all the stats dataframes for each project
    stats_dfs = {}

    for project, stats_list in project_stats.items():
        # Calculate descriptive statistics for each network statistic
        stats_df = pd.DataFrame(stats_list)
        mean = round(stats_df.mean().rename('{}_mean'.format(project)), 3)
        std = round(stats_df.std().rename('{}_std'.format(project)), 3)

        # Append the results to the DataFrame
        results = pd.concat([results, mean, std], axis=1)

        # Store the stats dataframe 
        stats_dfs[project] = stats_df

    # Save results to a text file
    results.to_csv(filename)

    # Create a box and whisker plot for each statistic comparing projects
    for stat in stats_dfs[list(stats_dfs.keys())[0]].columns:
        data = [df[stat] for df in stats_dfs.values()]
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111)
        ax.boxplot(data)
        ax.set_xticklabels(stats_dfs.keys())
        plt.title('Box plot of {} across {}'.format(stat, project_name))
        plt.savefig('{}_{}_boxplot.{}'.format(project_name, stat, format))
        plt.show()

    return results



# ROUND THE NUMBERS 3 DECIMAL PLACES
def average_network_stats(matrices, corr_threshold):
    network_stats = {'edges': [], 'nodes': [], 'betweenness_centrality': [], 'transitivity': [], 'closeness_centrality': [],
                     'density': []}

    for matrix in matrices:
        G = nx.Graph()

        for col1, col2 in combinations(matrix.columns, 2):
            if not np.isnan(matrix.loc[col1, col2]) and abs(matrix.loc[col1, col2]) > corr_threshold:
                G.add_edge(col1, col2, weight=matrix.loc[col1, col2])

        network_stats['edges'].append(G.number_of_edges())
        network_stats['nodes'].append(G.number_of_nodes())
        network_stats['betweenness_centrality'].extend([float(x) for x in nx.betweenness_centrality(G).values()])
        network_stats['transitivity'].append(nx.transitivity(G))
        network_stats['closeness_centrality'].extend([float(x) for x in nx.closeness_centrality(G).values()])
        network_stats['density'].append(nx.density(G))

    avg_edges = round(np.mean(network_stats['edges']), 3)
    std_edges = round(scipy.stats.std(network_stats['edges']), 3)
    avg_nodes = round(np.mean(network_stats['nodes']), 3)
    std_nodes = round(scipy.stats.std(network_stats['nodes']), 3)
    avg_betweenness_centrality = round(np.mean(network_stats['betweenness_centrality']), 3)
    avg_transitivity = round(np.mean(network_stats['transitivity']), 3)
    avg_closeness = round(np.mean(network_stats['closeness_centrality']), 3)
    avg_density = round(np.mean(network_stats['density']), 3)

    stats_dict = {
        'number_of_edges': avg_edges,
        'standard_error_edges': std_edges,
        'number_of_nodes': avg_nodes,
        'standard_error_nodes': std_nodes,
        'betwenness_centrality': avg_betweenness_centrality,
        'transitivity': avg_transitivity,
        'closeness_centrality': avg_closeness,
        'density': avg_density
    }

    return stats_dict


# INPUT DATA IS A LIST  OF MATRICES
# Takes in multiple matrices then outputs 1

def build_network_graph(correlations, threshold=0, title="Correlation Network"):
    # Check if the input is a list of DataFrames
    if isinstance(correlations, list):
        # Average out the correlation matrices
        concat_data = pd.concat(correlations)
        avg_corr_matrix = concat_data.groupby(concat_data.index).mean()
    elif isinstance(correlations, pd.DataFrame):
        # Use the DataFrame directly
        avg_corr_matrix = correlations
    else:
        raise ValueError("Input data must be a list of DataFrames or a single DataFrame")
    
    avg_corr_matrix.to_csv(f"{title}_correlation_matrix.csv")
    
    G = nx.Graph()

    # Adding nodes and edges with weights based on the correlations
    # Determined by the strength of the correlation
    for i in avg_corr_matrix.columns:
        for j in avg_corr_matrix.index:
            # Filtering significant correlations
            if i != j and abs(avg_corr_matrix.at[j, i]) > threshold:
                G.add_edge(i, j, weight=avg_corr_matrix.at[j, i])

    # Drawing the network graph
    # edata is a dictionary of edge data
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, seed=42)
    edges = G.edges(data=True)
    # _, _ placeholders for tuples to be non-str values
    colors = ['red' if edata['weight'] < 0 else 'blue' for _, _, edata in edges]
    weights = [abs(edata['weight']) for _, _, edata in edges]
    nx.draw_networkx(G, pos, edge_color=colors, edge_cmap=plt.cm.coolwarm, node_color='skyblue', with_labels=True)

    # Create legend
    pos_legend = {'Positive Correlation': (0, 0), 'Negative Correlation': (0, 1)}
    for label, pos in pos_legend.items():
        plt.scatter([], [], c='blue' if label == 'Positive Correlation' else 'red', label=label)
    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, bbox_to_anchor=(1, 1), loc='upper left')

    plt.title(title)
    plt.savefig(title + ".SVG", format="SVG", bbox_inches='tight')
    plt.show()

    return G

def no_label_network_graph(correlations, threshold=0, title="Correlation Network"):
    # Check if the input is a list of DataFrames
    if isinstance(correlations, list):
        # Average out the correlation matrices
        concat_data = pd.concat(correlations)
        avg_corr_matrix = concat_data.groupby(concat_data.index).mean()
    elif isinstance(correlations, pd.DataFrame):
        # Use the DataFrame directly
        avg_corr_matrix = correlations
    else:
        raise ValueError("Input data must be a list of DataFrames or a single DataFrame")
        
    avg_corr_matrix.to_csv(f"{title}_correlation_matrix.csv")
    
    G = nx.Graph()

    # Adding edges with weights based on the correlations
    for i in avg_corr_matrix.columns:
        for j in avg_corr_matrix.index:
            # Filtering significant correlations
            if i != j and abs(avg_corr_matrix.at[j, i]) > threshold:
                G.add_edge(i, j, weight=avg_corr_matrix.at[j, i])

    # Drawing the network graph
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, seed=42)
    edges = G.edges(data=True)
    colors = ['red' if edata['weight'] < 0 else 'blue' for _, _, edata in edges]
    weights = [abs(edata['weight']) for _, _, edata in edges]
    nx.draw_networkx(G, pos, edge_color=colors, edge_cmap=plt.cm.coolwarm, node_color='skyblue', with_labels=False)

    # Create legend
    pos_legend = {'Positive Correlation': (0, 0), 'Negative Correlation': (0, 1)}
    for label, pos in pos_legend.items():
        plt.scatter([], [], c='blue' if label == 'Positive Correlation' else 'red', label=label)
    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, bbox_to_anchor=(1, 1), loc='upper left')

    plt.title(title)
    plt.savefig(title + ".SVG", format="SVG", bbox_inches='tight')
    plt.show()

    return G


def top_nodes_network_graph(correlations, threshold=0.8, num_nodes=20, title="Correlation Network"):
    # Check if the input is a list of DataFrames
    if isinstance(correlations, list):
        # Average out the correlation matrices
        concat_data = pd.concat(correlations)
        avg_corr_matrix = concat_data.groupby(concat_data.index).mean()
    elif isinstance(correlations, pd.DataFrame):
        # Use the DataFrame directly
        avg_corr_matrix = correlations
    else:
        raise ValueError("Input data must be a list of DataFrames or a single DataFrame")
    
    avg_corr_matrix.to_csv(f"{title}_correlation_matrix.csv")
    
    G = nx.Graph()

    # Adding nodes and edges with weights based on the correlations
    for i in avg_corr_matrix.columns:
        for j in avg_corr_matrix.index:
            # Filtering significant correlations
            if i != j and abs(avg_corr_matrix.at[j, i]) > threshold:
                G.add_edge(i, j, weight=avg_corr_matrix.at[j, i])

    # Drawing the network graph
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, seed=42)
    edges = G.edges(data=True)
    colors = ['red' if edata['weight'] < 0 else 'blue' for _, _, edata in edges]

    # Get the top 20 nodes with the highest degree
    top_nodes = sorted(G.degree, key=lambda x: x[1], reverse=True)[:num_nodes]
    top_nodes_labels = {node: node for node, _ in top_nodes}

    nx.draw_networkx(G, pos, edge_color=colors, edge_cmap=plt.cm.coolwarm, node_color='skyblue', with_labels=False)
    nx.draw_networkx_labels(G, pos, labels=top_nodes_labels)

    # Create legend
    pos_legend = {'Positive Correlation': (0, 0), 'Negative Correlation': (0, 1)}
    for label, pos in pos_legend.items():
        plt.scatter([], [], c='blue' if label == 'Positive Correlation' else 'red', label=label)
    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, bbox_to_anchor=(1, 1), loc='upper left')

    plt.title(title)
    plt.savefig(title + ".SVG", format="SVG", bbox_inches='tight')
    plt.show()

    # Return the graph and the list of top 20 nodes
    return G, [node for node, _ in top_nodes]


def top_nodes(correlations, threshold=0.8, num_nodes=20):
    # Check if the input is a list of DataFrames
    if isinstance(correlations, list):
        # Average out the correlation matrices
        concat_data = pd.concat(correlations)
        avg_corr_matrix = concat_data.groupby(concat_data.index).mean()
    elif isinstance(correlations, pd.DataFrame):
        # Use the DataFrame directly
        avg_corr_matrix = correlations
    else:
        raise ValueError("Input data must be a list of DataFrames or a single DataFrame")

    G = nx.Graph()

    # Adding nodes and edges with weights based on the correlations
    for i in avg_corr_matrix.columns:
        for j in avg_corr_matrix.index:
            # Filtering significant correlations
            if i != j and abs(avg_corr_matrix.at[j, i]) > threshold:
                G.add_edge(i, j, weight=avg_corr_matrix.at[j, i])

    pos = nx.spring_layout(G, seed=42)
    edges = G.edges(data=True)

    # Get the top 20 nodes with the highest degree
    top_nodes = sorted(G.degree, key=lambda x: x[1], reverse=True)[:num_nodes]
    top_nodes_labels = []
    top_nodes_labels = [node for node, _ in top_nodes]

    # Return the graph and the list of top 20 nodes
    return top_nodes_labels


def build_positive_network(correlations, threshold=0, title="Positive Correlation Network"):
    # Check if the input is a list of DataFrames
    if isinstance(correlations, list):
        # Average out the correlation matrices
        concat_data = pd.concat(correlations)
        avg_corr_matrix = concat_data.groupby(concat_data.index).mean()
    elif isinstance(correlations, pd.DataFrame):
        # Use the DataFrame directly
        avg_corr_matrix = correlations
    else:
        raise ValueError("Input data must be a list of DataFrames or a single DataFrame")

    G = nx.Graph()

    # Adding nodes and edges with weights based on the correlations
    for i in avg_corr_matrix.columns:
        for j in avg_corr_matrix.index:
            # Filtering significant correlations
            if i != j and avg_corr_matrix.at[j, i] > threshold:
                G.add_edge(i, j, weight=avg_corr_matrix.at[j, i])

    # Drawing the network graph
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, seed=42)
    edges = G.edges(data=True)
    colors = [edata['weight'] for _, _, edata in edges]
    nx.draw_networkx(G, pos, edge_color=colors, edge_cmap=plt.cm.Blues, node_color='skyblue', with_labels=True)

    plt.title(title)
    plt.savefig(title + ".SVG", format="SVG", bbox_inches='tight')
    plt.show()

    return G


# FILTERS NETWORK BY NUMBER OF EDGES PER NODE
# DEFAULT MAX 7 EDGES PER NODE

def build_filtered_networks(correlations, threshold=0.8, max_degree=7, title="filtered network"):
    # Check if the input is a list of DataFrames
    if isinstance(correlations, list):
        # Average out the correlation matrices
        concat_data = pd.concat(correlations)
        avg_corr_matrix = concat_data.groupby(concat_data.index).mean()
    elif isinstance(correlations, pd.DataFrame):
        # Use the DataFrame directly
        avg_corr_matrix = correlations
    else:
        raise ValueError("Input data must be a list of DataFrames or a single DataFrame")

    G = nx.Graph()

    # Adding nodes and edges with weights based on the correlations
    for i in avg_corr_matrix.columns:
        for j in avg_corr_matrix.index:
            # Adding only significant correlations
            if i != j and avg_corr_matrix.loc[j, i] > threshold:
                G.add_edge(i, j, weight=avg_corr_matrix.loc[j, i])

    # Filter nodes by degree (number of connections)
    nodes_to_keep = [node for node, degree in dict(G.degree()).items() if degree <= max_degree]

    # Create a subgraph with only the filtered nodes
    filtered_graph = G.subgraph(nodes_to_keep)

    # Drawing the network graph
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(filtered_graph, seed=42)  # For consistent layout
    edges = filtered_graph.edges(data=True)
    colors = [edata['weight'] for _, _, edata in edges]
    nx.draw_networkx(filtered_graph, pos, edge_color=colors, edge_cmap=plt.cm.coolwarm, node_color='skyblue',
                     with_labels=True)

    plt.title(title + f'max degree {max_degree}')
    plt.savefig(title + ".SVG", format="SVG", bbox_inches='tight')
    plt.show()

    return G

# INPUT DATA IS A LIST  OF MATRICES
def build_negative_networks(correlations, threshold=0, title='Project'):
    # Check if the input is a list of DataFrames
    if isinstance(correlations, list):
        # Average out the correlation matrices
        concat_data = pd.concat(correlations)
        avg_corr_matrix = concat_data.groupby(concat_data.index).mean()
    elif isinstance(correlations, pd.DataFrame):
        # Use the DataFrame directly
        avg_corr_matrix = correlations
    else:
        raise ValueError("Input data must be a list of DataFrames or a single DataFrame")
    
    G = nx.Graph()

    # Adding nodes and edges with weights based on the correlations
    for i in avg_corr_matrix .columns:
        for j in avg_corr_matrix .index:
            # Filtering significant correlations
            if i != j and avg_corr_matrix.loc[j, i] < threshold and avg_corr_matrix.loc[j, i] < 0:
                G.add_edge(i, j, weight=avg_corr_matrix .loc[j, i])

    # Drawing the network graph
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, seed=42)
    edges = G.edges(data=True)
    colors = [edata['weight'] for _, _, edata in edges]
    nx.draw_networkx(G, pos, edge_color=colors, edge_cmap=plt.cm.coolwarm, node_color='skyblue', with_labels=True)

    plt.title(f"Negative Correlation Network for {title}")
    plt.savefig(title + ".SVG", format="SVG", bbox_inches='tight')
    plt.show()

    return G


def negative_filtered_networks(correlations, threshold, max_degree, title='Project'):
    # Check if the input is a list of DataFrames
    if isinstance(correlations, list):
        # Average out the correlation matrices
        concat_data = pd.concat(correlations)
        avg_corr_matrix = concat_data.groupby(concat_data.index).mean()
    elif isinstance(correlations, pd.DataFrame):
        # Use the DataFrame directly
        avg_corr_matrix = correlations
    else:
        raise ValueError("Input data must be a list of DataFrames or a single DataFrame")
    
    G = nx.Graph()

    # Adding nodes and edges with weights based on the correlations
    for i in avg_corr_matrix.columns:
        for j in avg_corr_matrix.index:
            # Adding only significant correlations
            if i != j and avg_corr_matrix.loc[j, i] < -threshold:
                G.add_edge(i, j, weight=avg_corr_matrix.loc[j, i])

    # Filter nodes by degree (number of connections)
    nodes_to_keep = [node for node, degree in dict(G.degree()).items() if degree <= max_degree]

    # Create a subgraph with only the filtered nodes
    filtered_graph = G.subgraph(nodes_to_keep)

    # Drawing the network graph
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(filtered_graph, seed=42)  # For consistent layout
    edges = filtered_graph.edges(data=True)
    colors = [edata['weight'] for _, _, edata in edges]
    nx.draw_networkx(filtered_graph, pos, edge_color=colors, edge_cmap=plt.cm.coolwarm, node_color='skyblue',
                        with_labels=True)

    plt.title(f"Average Negative Correlation Network by the degree {max_degree}")
    plt.savefig(title + ".SVG", format="SVG", bbox_inches='tight')
    plt.show()

    return G


# Bootstrapping data sets, finding correlations, then averaging the bootstrapped replicate correlations 
# Returns one matrix
def bootstrap_sample_with_correlation(data, n_iterations):
    n = len(data)
    correlation_results = []

    for i in range(n_iterations):
        # Creating a subsample of the data with replacement
        sample = data.sample(n, replace=True)

        # Calculate Spearman correlation for the sample
        correlation = sample.corr(method="spearman", numeric_only=True)

        # Store the results of the correlations
        correlation_results.append(correlation)

    # Concatenate all correlation results into a single DataFrame
    concat_data = pd.concat(correlation_results)
    all_correlations = concat_data.groupby(concat_data.index).mean()
    return all_correlations


def most_connected_nodes(G):
    # Calculate the degree of each node
    degrees = dict(G.degree())
    # Find the node with the maximum degree
    max_degree_node = max(degrees, key=degrees.get)
    return "Most connected nodes: " + max_degree_node


def nodes_edges_table(G):
    # Calculate the degree of each node
    degrees = dict(G.degree())
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame.from_dict(degrees, orient='index', columns=['# of edges'])
    return df


def save_table_to_csv(df, filename):
    df.to_csv(filename)


# Binomial Test

def net_stat_binomial_test(reps1, reps2, title='p_values.csv'):
    # Check if lists of replicates are of equal length
    if len(reps1) != len(reps2):
        print("Bootstrap replicates not equal numbers")
        return

    # Initialize a dictionary to store the p-values for each network statistic
    p_values = {}

    # Get the list of network statistics from the first replicate
    network_stats = list(reps1[0].keys())

    # Calculate a binomial test for each network statistic
    for stat in network_stats:
        s = []
        for i in range(len(reps2)):
            if reps1[i][stat] > reps2[i][stat]:
                s.append(1)
            else:
                s.append(0)
        # Binomial test
        p_values[stat] = stats.binomtest(sum(s), len(s), p=0.5)

    pd.DataFrame.from_dict(p_values, orient='index', columns=['p-value']).to_csv(title)

    return p_values
