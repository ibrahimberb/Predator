def is_core(protein, mutation, elaspic_core_data):
    """
    Given protein and mutation, check if that protein.mutation is in elaspic_core_data.

    Parameters
    ----------
        protein : <str>
            The UniProt/SWISSPROT protein name. E.g. "Q9UPU5".

        mutation : <str>
            The position of the mutation. E.g. "I342V".

        elaspic_core_data : <DataFrame>
            The ELASPIC results file that contains only the `core` type entries.

    Returns
    -------
        is_core_bool : <bool>
            A boolean value that indicates whether the given protein.mutation is an `core` type.

    """

    core_search_data = elaspic_core_data[(elaspic_core_data["UniProt_ID"] == protein) &
                                         (elaspic_core_data["Mutation"] == mutation)]

    is_core_bool = True if len(core_search_data) else False

    return is_core_bool
