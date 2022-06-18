from tqdm.auto import tqdm
from .gene_id_retrieval import GeneIDRetriever


def get_protein_to_gene_dict(proteins):
    """
    Given a list of protein for a TCGA type, construct a dictionary that maps each protein to
    corresponding gene.

    E.g.
        USP24  → Q9UPU5
        ERICH3 → Q5RHP9
        KIF26B → Q2KJY2

    Parameters
    ----------
        proteins : <list>
            List of proteins.

    Returns
    -------
       protein_to_genes : <dict>
            A dictionary that maps each protein to corresponding gene.
    """

    uniprot_to_gene_id = {}

    gene_retriever = GeneIDRetriever()

    for protein in tqdm(proteins, desc="Retrieving Gene IDs from UniProt API .. "):
        gene = gene_retriever.fetch(protein=protein)
        uniprot_to_gene_id[protein] = gene

    return uniprot_to_gene_id


def get_protein_to_gene_dict_via_snv(proteins, snv_data):
    """
    Given a list of protein for a TCGA type, construct a dictionary that maps each protein to
    corresponding gene.
    
    E.g.
        USP24  → Q9UPU5
        ERICH3 → Q5RHP9
        KIF26B → Q2KJY2
        
    Parameters
    ----------
        proteins : <list>
            List of proteins.
        
        snv_data : <DataFrame>
            Processed version of SNV data, which is the output of `process_snv` function.
    
    Returns
    -------
       protein_to_genes : <dict>
            A dictionary that maps each protein to corresponding gene.
    """

    protein_to_genes = dict()
    exceptional_cases = []
    for protein in tqdm(proteins):
        # genes should be a list that contains 1 item only.
        genes = list(snv_data[snv_data["SWISSPROT"] == protein]["Hugo_Symbol"].unique())

        # If there are more then one genes or no genes, we put "NA".
        if len(genes) != 1:
            exceptional_cases.append((protein, genes))
            protein_to_genes[protein] = "NA"

        else:
            # Add to the dictionary
            protein_to_genes[protein] = genes[0]

    print(f"Number of exceptional cases: {len(exceptional_cases)}")
    print(f"Exceptional cases: {exceptional_cases}")

    return protein_to_genes
