import pandas as pd
import requests
import time
import numpy as np
from thermo import Chemical, Mixture
import warnings



def fetch_details_for_batch(cids):
    """Fetches SMILES, name, and formula details for a batch of CIDs from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{','.join(map(str, cids))}/property/IsomericSMILES,IUPACName,MolecularFormula/JSON"
    response = requests.get(url)
    batch_data = {}
    
    if response.status_code == 200:
        try:
            properties = response.json().get('PropertyTable', {}).get('Properties', [])
            for prop in properties:
                cid = prop.get('CID')
                smiles = prop.get('IsomericSMILES')
                name = prop.get('IUPACName')
                formula = prop.get('MolecularFormula')
                if cid and smiles:
                    batch_data[cid] = {
                        'smiles': smiles,
                        'name': name,
                        'formula': formula
                    }
        except Exception as e:
            print(f"Error processing batch {cids[0]}-{cids[-1]}: {e}")
    else:
        print(f"Failed to retrieve batch {cids[0]}-{cids[-1]}: Status code {response.status_code}")
    
    return batch_data




def generate_material_bank(start_cid, end_cid, batch_size=100):
    """Generates a material_bank DataFrame for a specified range of CIDs."""
    cid_list = list(range(start_cid, end_cid + 1))
    smiles_data = {}
    
    for i in range(0, len(cid_list), batch_size):
        batch = cid_list[i:i + batch_size]
        batch_data = fetch_details_for_batch(batch)
        smiles_data.update(batch_data)  # Add fetched data to the main dictionary
        time.sleep(0.2)  # Pause to respect rate limits
        print(f'Batch {i // batch_size + 1} completed')
    
    # Convert the dictionary to a DataFrame
    material_bank = pd.DataFrame.from_dict(smiles_data, orient='index')
    material_bank.index.name = 'CID'
    
    return material_bank




supported_properties = {
    'boiling point': 'Tb',
    'melting point': 'Tm',
    'viscosity': 'mu',
    'density': 'rho',
    'thermal conductivity': 'k',
#     'heat capacity': 'Cp',
    'vapour pressure': 'Psat',
    'permittivity': 'permittivity',
    'flash point': 'Tflash',
    'auto ignition temperature': 'Tautoignition',
    'heat of combustion': 'Hc',
    'enthalpy of formation': 'Hf',
    'critical temperature': 'Tc',
    'critical pressure': 'Pc',
    'critical volume': 'Vc',
    'triple point temperature': 'Tt',
    'triple point pressure': 'Pt',
    'surface tension': 'sigma',
    'molecular weight': 'MW' ,
    'Henry\'s law constant': 'Henry',
    'dipole moment': 'dipole'
}

def process_materials(material_bank, properties, run_name):
    # Load the Excel file
#     material_bank = pd.read_excel(file_path)

    # Validate and filter properties
    valid_properties = {}
    for prop in properties:
        if prop in supported_properties:
            valid_properties[prop] = supported_properties[prop]
        else:
            warnings.warn(f"Property '{prop}' is not supported and will be ignored.")

    # Define solvents and mole fractions for toluene
    solvents = ['water', 'hexane', 'ethanol', 'acetonitrile']
    fractions = [0.25, 0.5, 0.75]  # Mole fractions for toluene

    # Generate columns for enthalpy of mixing with each solvent and fraction
    for solvent in solvents:
        for fraction in fractions:
            column_name = f'{solvent}_fraction_{fraction}_enthalpy_mixing'
            material_bank[column_name] = np.nan

    # Loop over each material and calculate requested properties
    for i in range(material_bank.shape[0]):
        try:
            # Create Chemical object for the formula
            chemical = Chemical(material_bank.loc[i, 'name'])

            # Calculate each requested property individually
            for prop, attr in valid_properties.items():
                try:
                    material_bank.loc[i, prop] = getattr(chemical, attr)
#                     material_bank.loc[i, prop] = chemical.attr
                except AttributeError:
                    material_bank.loc[i, prop] = np.nan

            # Calculate enthalpy of mixing for toluene with each solvent at each fraction
            for solvent in solvents:
                for fraction in fractions:
                    column_name = f'{solvent}_fraction_{fraction}_enthalpy_mixing'
                    try:
                        mixture = Mixture([material_bank.loc[i,'name'], solvent], zs=[fraction, 1 - fraction])
                        enthalpy_of_mixing = mixture.Hm  # Ideal enthalpy of mixing                        
                        material_bank.loc[i, column_name] = enthalpy_of_mixing
                    except Exception:
                        material_bank.loc[i, column_name] = np.nan

        except Exception as e:
            print(f"Error processing material at index {i}: {e}")
        time.sleep(0.002)
        print(f'Material number {i} processed')

    # Save the updated DataFrame to a new Excel file
    output_file = f'{run_name}.parquet'
#     material_bank.to_excel(output_file, index=False)
    material_bank.to_parquet(output_file, index=False)

    
#     print(f"Results saved to '{output_file}'")



