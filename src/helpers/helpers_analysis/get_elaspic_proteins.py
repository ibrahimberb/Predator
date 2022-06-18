def get_elaspic_proteins(core_data, interface_data):
    """
    Given ELASPIC result dataframes core_data and interface_data, get unique proteins from 
    both dataframes. Merge these proteins and return.
    
    Parameters
    ----------
        core_data : <DataFrame>
            ELASPIC results which are all `core` mutations.
            
        interface_data : <DataFrame>
            ELASPIC results which are all `interface` mutations.
    
    Returns
    -------
        proteins : <list>
            Merged proteins from `core` and `interface` dataframes.
    
    """

    # Merge unique proteins from `core` and `interface` dataframes.
    proteins = sorted(set(list(core_data['UniProt_ID'].unique()) + list(interface_data['UniProt_ID'].unique())))
    
    return proteins
