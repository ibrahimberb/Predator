from collections import defaultdict
from tqdm.notebook import tqdm
import pandas as pd


def get_protein_to_num_elaspic_interface_entires_dict(proteins: list, interface_data: pd.DataFrame):
    """
    Given the ELASPIC Protein list, return the dictionary `protein_to_num_elaspic_interface_entires` that
    mappes each protein to its number of ELASPIC interface entires.

    Parameters
    ----------
        proteins : <list>
            List of proteins.

        interface_data : <DataFrame>
            ELASPIC results which are all `interface` mutations.

    Returns
    -------
        protein_to_num_elaspic_interface_entries : <defaultdict>
            A (default) dictionary that maps each protein to its number of ELASPIC interface entries.

    """

    # Dictionary that maps protein <str> to number of ELASPIC entries <int>.
    protein_to_num_elaspic_interface_entries = defaultdict(int)

    # Fill in the dictionary `protein_to_num_elaspic_interface_entries`.
    for i, protein in enumerate(tqdm(proteins)):
        num_entries = len(interface_data[interface_data["UniProt_ID"] == protein])
        protein_to_num_elaspic_interface_entries[protein] = num_entries
        # Printing first five protein and their number of elaspic interface entries.
        if i < 5:
            print(f'{i} \t {protein} \t {num_entries}')

    return protein_to_num_elaspic_interface_entries
