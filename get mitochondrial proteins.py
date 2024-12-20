import pandas as pd
from chembl_webresource_client.new_client import new_client
import threading
import time

# ========================================
# Timer Thread Setup
running = True

def timer():
    start_time = time.time()
    while running:
        elapsed = int(time.time() - start_time)
        print(f"Script running for {elapsed} seconds...", end='\r')
        time.sleep(1)

# Start the timer thread
t = threading.Thread(target=timer, daemon=True)
t.start()
# ========================================

# Path to your MitoCarta3.0 Excel file
input_file = "/Users/tuanvinh/Downloads/Human.MitoCarta3.0.xls"

# Specify the sheet name
sheet_name = "A Human MitoCarta3.0"

# Read the specified sheet from the Excel file
df = pd.read_excel(input_file, sheet_name=sheet_name)

# Ensure UniProt column is treated as strings
df["UniProt"] = df["UniProt"].astype(str)

# Extract the UniProt column
uniprot_column = df["UniProt"]

# Convert to a set to remove duplicates
unique_uniprot_ids = set(uniprot_column.dropna().unique())

# Convert all IDs to strings and sort them
unique_uniprot_ids = [str(uid) for uid in unique_uniprot_ids]
unique_uniprot_ids = sorted(unique_uniprot_ids)

# Specify the output file
output_file = "unique_uniprot_ids.txt"

# Write unique UniProt IDs to a text file, one per line
with open(output_file, "w") as f:
    for uniprot_id in unique_uniprot_ids:
        f.write(uniprot_id + "\n")

print("\nUnique UniProt IDs have been saved to unique_uniprot_ids.txt\n")

# Path to the file containing UniProt IDs
uniprot_file = "/Users/tuanvinh/Desktop/Mitochondria Protein Targets/unique_uniprot_ids.txt"

# Load UniProt IDs
with open(uniprot_file, "r") as f:
    uniprot_ids = [line.strip() for line in f if line.strip()]

# Initialize ChEMBL clients
target = new_client.target
activity = new_client.activity
molecule = new_client.molecule

# Dictionary to store results: { UniProt_ID: { Target ChEMBL ID: [List of Drugs] } }
uniprot_to_drugs = {}

for uniprot_id in uniprot_ids:
    # Find targets in ChEMBL associated with this UniProt accession
    target_results = target.filter(target_components__accession=uniprot_id)
    
    if not target_results:
        # No targets found for this UniProt
        continue
    
    uniprot_to_drugs[uniprot_id] = {}
    
    for t_result in target_results:
        target_chembl_id = t_result['target_chembl_id']
        uniprot_to_drugs[uniprot_id][target_chembl_id] = []

        # Get activities for this target
        acts = activity.filter(target_chembl_id=target_chembl_id)

        # Extract molecule ChEMBL IDs
        molecule_ids = {a['molecule_chembl_id'] for a in acts if 'molecule_chembl_id' in a}

        # Retrieve molecule details
        for mol_id in molecule_ids:
            mol_data = molecule.get(mol_id)
            if mol_data:
                max_phase_value = mol_data.get('max_phase', 0)

                # Handle None or non-integer max_phase
                if max_phase_value is None:
                    max_phase = 0
                else:
                    try:
                        max_phase = int(max_phase_value)
                    except ValueError:
                        max_phase = 0

                first_approval = mol_data.get('first_approval', None)
                
                # Consider max_phase >= 4 as drug-like
                if max_phase >= 4:
                    uniprot_to_drugs[uniprot_id][target_chembl_id].append({
                        'molecule_chembl_id': mol_id,
                        'pref_name': mol_data.get('pref_name', ''),
                        'max_phase': max_phase,
                        'first_approval': first_approval
                    })

# Print results
print("\nResults:")
for uni_id, targets_dict in uniprot_to_drugs.items():
    print(f"UniProt: {uni_id}")
    for t_id, drugs in targets_dict.items():
        print(f"  Target: {t_id}")
        if drugs:
            for d in drugs:
                print(f"    Drug: {d['molecule_chembl_id']} | Name: {d['pref_name']} | Max Phase: {d['max_phase']} | First Approval: {d['first_approval']}")
        else:
            print("    No marketed drugs found for this target.")

# Stop the timer thread
running = False
t.join()

print("Script completed.")