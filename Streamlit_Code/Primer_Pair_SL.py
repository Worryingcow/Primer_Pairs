#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:16:30 2024

@author: thomasminahan
"""
import streamlit as st
import pandas as pd
import random
from collections import Counter
from io import BytesIO

# Function to calculate nucleotide frequencies for a list of primers
def calculate_nucleotide_frequencies(primers):
    nucleotide_counts = [Counter() for _ in range(len(next(iter(primers.values()))))]
    for seq in primers.values():
        for i, nucleotide in enumerate(seq):
            nucleotide_counts[i][nucleotide] += 1
    return nucleotide_counts

# Function to calculate a diversity score
def calculate_diversity_score(nucleotide_counts, total_primers):
    ideal_count = total_primers / 4
    score = sum(abs(ideal_count - count) for pos in nucleotide_counts for count in pos.values())
    return score

# Function to select the optimal combination of primers
def select_optimal_primers(forward_primers, reverse_primers, num_forward, num_reverse):
    best_score = float('inf')
    best_combination = None

    for _ in range(10000):
        selected_forward_keys = random.sample(list(forward_primers.keys()), num_forward)
        selected_reverse_keys = random.sample(list(reverse_primers.keys()), num_reverse)

        selected_forward = {key: forward_primers[key] for key in selected_forward_keys}
        selected_reverse = {key: reverse_primers[key] for key in selected_reverse_keys}

        forward_frequencies = calculate_nucleotide_frequencies(selected_forward)
        reverse_frequencies = calculate_nucleotide_frequencies(selected_reverse)

        total_frequencies = [Counter() for _ in range(8)]
        for i in range(8):
            total_frequencies[i] = forward_frequencies[i] + reverse_frequencies[i]

        score = calculate_diversity_score(total_frequencies, num_forward + num_reverse)

        if score < best_score:
            best_score = score
            best_combination = (selected_forward, selected_reverse)

    return best_combination, best_score, total_frequencies

# Function to create the nucleotide matrix
def create_nucleotide_matrix(frequencies):
    data = {'Basepair position': ['A', 'T', 'C', 'G']}
    for pos in range(8):
        counts = [frequencies[pos].get(base, 0) for base in 'ATCG']
        data[pos + 1] = counts

    matrix_df = pd.DataFrame(data)
    sum_row = matrix_df.iloc[:, 1:].sum()
    sum_row = pd.DataFrame([['Sum to:'] + sum_row.tolist()], columns=matrix_df.columns)
    matrix_df = pd.concat([matrix_df, sum_row], ignore_index=True)
    return matrix_df

# Streamlit UI
st.title("Primer Pair Optimization App")

st.text("Formating Requirements for the excel files can be found on the github repository")
st.text("in the Excel_formating_requirements folder")
st.text("(To get to the github click the button in the top right")
# Upload the combined primer file
primer_file = st.file_uploader("Upload Primer Excel File (Forward and Reverse)", type=["xlsx"])

# Optional: Excluded Primer File Upload
excluded_file = st.file_uploader("Upload Excluded Primer Excel File (Optional)", type=["xlsx"])

# Checkbox to include or exclude primers
include_excluded = st.checkbox("Include excluded primers?", value=False)

if primer_file:
    # Load data from the uploaded primer file
    primers_df = pd.read_excel(primer_file)

    # Split the primers into forward and reverse based on the 'type' column
    forward_primers_df = primers_df[primers_df['type'].str.lower() == 'forward']
    reverse_primers_df = primers_df[primers_df['type'].str.lower() == 'reverse']

    # Convert DataFrames to dictionaries
    forward_primers = dict(zip(forward_primers_df['indexname'], forward_primers_df['sequence']))
    reverse_primers = dict(zip(reverse_primers_df['indexname'], reverse_primers_df['sequence']))

    # Handle optional excluded primer logic
    if include_excluded and excluded_file:
        excluded_primers_df = pd.read_excel(excluded_file)
        excluded_primers = set(excluded_primers_df['indexname'])
        forward_primers = {name: seq for name, seq in forward_primers.items() if name not in excluded_primers}
        reverse_primers = {name: seq for name, seq in reverse_primers.items() if name not in excluded_primers}

    # User inputs: Number of forward and reverse primers
    max_forward = len(forward_primers)
    max_reverse = len(reverse_primers)

    num_forward = st.slider("Number of Forward Primers to Use", min_value=1, max_value=max_forward, value=5)
    num_reverse = st.slider("Number of Reverse Primers to Use", min_value=1, max_value=max_reverse, value=10)

    # Select optimal primers
    optimal_primers, diversity_score, total_frequencies = select_optimal_primers(
        forward_primers, reverse_primers, num_forward, num_reverse
    )

    # Combine primers into a DataFrame
    combined_primers = {**optimal_primers[0], **optimal_primers[1]}
    primers_df = pd.DataFrame(list(combined_primers.items()), columns=["indexname", "sequence"])

    # Create nucleotide matrix
    matrix_df = create_nucleotide_matrix(total_frequencies)

    # Display the matrix as a table in Streamlit
    st.subheader("Nucleotide Frequency Matrix")
    st.table(matrix_df)

    # Save the optimal primers to an Excel file for download
    def to_excel(df):
        output = BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            df.to_excel(writer, sheet_name='Optimal Primers', index=False)
        return output.getvalue()

    excel_data = to_excel(primers_df)

    # Provide a download button for the primer pairs
    st.download_button(
        label="Download Optimal Primer Pairs",
        data=excel_data,
        file_name="optimal_primers.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )
else:
    st.warning("Please upload the primer file.")

