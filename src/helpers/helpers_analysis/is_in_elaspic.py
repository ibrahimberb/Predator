def is_in_elaspic(protein, mutation, elaspic_core_data, elaspic_interface_data):
    """
    Checks if given (protein.mutation) is in elaspic results file.
    """

    elaspic_core_search_data = elaspic_core_data[(elaspic_core_data["UniProt_ID"] == protein) &
                                                 (elaspic_core_data["Mutation"] == mutation)]

    elaspic_interface_data = elaspic_interface_data[(elaspic_interface_data["UniProt_ID"] == protein) &
                                                    (elaspic_interface_data["Mutation"] == mutation)]

    # If both core and interface dataframes are empty, then given (protein.mutation) is not in ELASPIC results data.
    if elaspic_core_search_data.empty and elaspic_interface_data.empty:
        return False
    else:
        return True
