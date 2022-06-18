from collections import defaultdict
from tqdm.notebook import tqdm
import pandas as pd


def get_protein_to_num_unique_interactors(proteins: list, interface_data: pd.DataFrame):
    """
    Given the ELASPIC Protein list, return the dictionary `protein_to_num_unique_interactors` that
    mappes each protein to its number of unique ELASPIC interactor proteins over all mutation positions.
    
    For instance, for protein P04637 (gene TP53), we filter UniProt_ID column in interface_data for P04637. 
        Then, we find unique number of interactors, i.e. unique proteins in Interactor_UniProt_ID.
    
        +-------------+-------------+-----------------------+
        |  UniProt_ID |   Mutation  | Interactor_UniProt_ID |
        +-------------+-------------+-----------------------+
        |    P04637   |    R280K    |          P62993       |
        |    P04637   |    R280K    |          Q13625       |
        |     ...     |     ...     |            ...        |
        |    P04637   |    R273P    |          P62993       |
        |     ...     |     ...     |            ...        |
        +-------------+-------------+-----------------------+
        
        P62993 will be counted once.
    

    Parameters
    ----------
        proteins : <list>
            List of proteins.

        interface_data : <DataFrame>
            ELASPIC results which are all `interface` mutations.

    Returns
    -------
        protein_to_num_unique_interactors : <defaultdict>
            A (default) dictionary that maps each protein to its number of unique interactors in ELASPIC interface data.

    """

    # Dictionary that maps protein <str> to number of ELASPIC entries <int>.
    protein_to_num_unique_interactors = defaultdict(int)

    # Fill in the dictionary `protein_to_num_unique_interactors`.
    for i, protein in enumerate(tqdm(proteins)):
        num_unique_interactors = len(interface_data[interface_data["UniProt_ID"] == protein]["Interactor_UniProt_ID"].unique())
        protein_to_num_unique_interactors[protein] = num_unique_interactors
        # Printing first five protein and their number of unique interactors in ELASPIC interface data.
        if i < 5:
            print(f'{i} \t {protein} \t {num_unique_interactors}')

    return protein_to_num_unique_interactors
