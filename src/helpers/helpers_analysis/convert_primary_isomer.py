def convert_primary_isomer(column_name, data):
    """
    Converts proteins into primary form representation, meaning dash-free from.
    E.g. 
        P16473-2 → P16473
    
    Parameters
    ----------
        column_name : <string>
            Name of the column where protein is stored.
        
        data : <DataFrame>
            The dataframe whose proteins will be processed in `column_name` column.
    
    Returns
    -------
        data : <DataFrame>
            Processed version of input dataframe.

    """
    
    # Protein names will be converted dashed-free version, if they contain.
    data[column_name] = data[column_name].apply(lambda x: x.split('-')[0])
    
    return data
