import streamlit as st
import pandas as pd
import pulp
from collections import Counter, defaultdict
from io import BytesIO
import math

st.title("Primer Pair Optimization App (ILP)")

# Upload primer files
primer_file = st.file_uploader("Upload Primer Excel File (Forward and Reverse)", type=["xlsx"])

# Optional: Upload used primers file
used_primers_file = st.file_uploader("Upload Used Primer Pairs File (Optional)", type=["xlsx"])

if primer_file:
    primers_df = pd.read_excel(primer_file)
    
    forward_primers_df = primers_df[primers_df['type'].str.lower() == 'forward']
    reverse_primers_df = primers_df[primers_df['type'].str.lower() == 'reverse']
    
    forward_primers = dict(zip(forward_primers_df['indexname'], forward_primers_df['sequence']))
    reverse_primers = dict(zip(reverse_primers_df['indexname'], reverse_primers_df['sequence']))

    # Calculate the maximum number of unique primer pairs possible
    max_possible_samples = len(forward_primers) * len(reverse_primers)

    used_pairs = set()
    used_forward_primers = set()
    used_reverse_primers = set()

    if used_primers_file:
        used_pairs_df = pd.read_excel(used_primers_file)
        used_pairs_df.columns = used_pairs_df.columns.str.strip()
        used_pairs = set(zip(used_pairs_df['Forward'], used_pairs_df['Reverse']))
        
        # Track used forward and reverse primers
        used_forward_primers = {fwd for fwd, _ in used_pairs}
        used_reverse_primers = {rev for _, rev in used_pairs}
    
    # User input for the number of additional samples needed
    num_new_samples = st.number_input(
        "Enter the number of new samples (each requires a unique primer pair):",
        min_value=1,
        max_value=max_possible_samples - len(used_pairs),
        value=60
    )

    # Calculate the optimal number of forward and reverse primers needed for the new pairs
    def optimal_primer_counts(n):
        root = int(math.sqrt(n))
        for i in range(root, 0, -1):
            if n % i == 0:
                return i, n // i
        return 1, n

    num_forward_to_select, num_reverse_to_select = optimal_primer_counts(num_new_samples)
    
    problem = pulp.LpProblem("Primer_Selection", pulp.LpMinimize)

    forward_vars = {key: pulp.LpVariable(f"f_{key}", cat="Binary") for key in forward_primers}
    reverse_vars = {key: pulp.LpVariable(f"r_{key}", cat="Binary") for key in reverse_primers}

    deviation_vars = defaultdict(lambda: defaultdict(lambda: defaultdict(pulp.LpVariable)))
    ideal_count = (num_forward_to_select + num_reverse_to_select) / 4

    for position in range(8):
        for base in 'ATCG':
            deviation_vars[position][base] = pulp.LpVariable(f"dev_{position}_{base}", lowBound=0)

    for position in range(8):
        for base in 'ATCG':
            forward_count = sum(forward_vars[key] for key, seq in forward_primers.items() if seq[position] == base)
            reverse_count = sum(reverse_vars[key] for key, seq in reverse_primers.items() if seq[position] == base)
            total_count = forward_count + reverse_count
            problem += deviation_vars[position][base] >= ideal_count - total_count
            problem += deviation_vars[position][base] >= total_count - ideal_count

    problem += pulp.lpSum(deviation_vars[position][base] for position in range(8) for base in 'ATCG')
    problem += pulp.lpSum(forward_vars[key] for key in forward_vars) == num_forward_to_select
    problem += pulp.lpSum(reverse_vars[key] for key in reverse_vars) == num_reverse_to_select

    # Exclude any used primer pairs during selection
    for fwd, rev in used_pairs:
        if fwd in forward_vars and rev in reverse_vars:
            problem += forward_vars[fwd] + reverse_vars[rev] <= 1  # Prevent selecting both primers of a used pair

    if st.button("Run Optimization"):
        problem.solve()

        # Extract selected forward and reverse primers
        selected_forward_primers = [key for key in forward_vars if forward_vars[key].varValue == 1]
        selected_reverse_primers = [key for key in reverse_vars if reverse_vars[key].varValue == 1]

        # Generate only the specified number of new unique primer pairs
        new_assigned_pairs = [
            (fwd, rev) for fwd in selected_forward_primers for rev in selected_reverse_primers
        ][:num_new_samples]

        st.subheader("Total Deviation (Objective)")
        st.write(pulp.value(problem.objective))

        # Function to export new primer pairs to Excel
        def to_excel(forward_primers, reverse_primers, new_pairs):
            # Create a DataFrame for new primer pairs with the desired columns
            sample_data = []
            for i, (fwd, rev) in enumerate(new_pairs, start=1):
                sample_data.append({
                    "Sample": i,
                    "Forward": fwd,
                    "Forward Sequence": forward_primers[fwd],
                    "Reverse": rev,
                    "Reverse Sequence": reverse_primers[rev]
                })
            pairs_df = pd.DataFrame(sample_data)

            output = BytesIO()
            with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
                pairs_df.to_excel(writer, index=False, sheet_name="New Primer Pairs")
            return output.getvalue()

        # Generate the Excel file with only the new primer pairs
        excel_data = to_excel(forward_primers, reverse_primers, new_assigned_pairs)
        
        st.download_button(
            label="Download New Primer Pairs",
            data=excel_data,
            file_name="new_primer_pairs.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
else:
    st.warning("Please upload the primer file.")
