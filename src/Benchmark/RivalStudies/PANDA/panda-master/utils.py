import os.path as op
import logging
import requests
from Bio import SeqIO
import pandas as pd
from datetime import datetime


logger = logging.getLogger(__name__)


def download_uniport_sequence(uniprot_id, output_dir):
    """"""
    output_file = op.join(output_dir, uniprot_id + ".fasta")

    # If the file already exists, do nothing...
    if op.isfile(output_file):
        logger.debug("Sequence file {} already exists...".format(output_file))
        return output_file

    logger.debug("Downloading sequence {}...".format(uniprot_id + ".fasta"))
    address = "http://www.uniprot.org/uniprot/{}.fasta".format(uniprot_id)
    r = requests.get(address)
    if r.status_code == 200:
        with open(output_file, "w") as ofh:
            ofh.write(r.text)
        return output_file
    
    else:
        raise
    
    
def get_sequence_from_fasta(fasta_file):
    sequence = SeqIO.read(fasta_file, "fasta")
    return sequence.seq


def get_mutated_sequence(sequence, mutation):
    original_aa, position, new_aa = mutation[0], int(mutation[1:-1]), mutation[-1]
    assert 0 < position <= len(sequence), f"The position {position} is out of range"
    assert original_aa == sequence[position - 1], f"The original amino acid at position {position} is not {original_aa}"
    print("mutation: {}, original_aa: {}, position: {}, new_aa: {}".format(mutation, original_aa, position, new_aa))
    
    return sequence[:position - 1] + new_aa + sequence[position:]
    
    # return sequence[:position] + new_aa + sequence[position+1:]


def download_protein_sequence(protein, output_dir, mutation=None):
    """
    Download the protein sequence from Uniprot.
    
    Parameters
    ----------
        protein : str
            The protein name.
            
        output_dir : str
            The output directory.
            
        mutation : str
            The mutation to be applied to the protein in the form of <original_residue><position><new_residue>.
            E.g. F118A
            
    Returns 
    -------
        sequence : str
            The sequence of the protein. If a mutation is provided, the sequence is of the mutated protein.
        
    """
    
    sequence_file = download_uniport_sequence(protein, output_dir)
    sequence = SeqIO.read(sequence_file, "fasta").seq
    
    if mutation is not None:
        sequence_mutated = get_mutated_sequence(sequence, mutation)
        
        save_sequence(sequence_mutated, protein + "_mutated", output_dir)
        return sequence_mutated

    else:
        save_sequence(sequence, protein, output_dir)
        return sequence


def save_sequence(sequence, filename, output_dir):
    """
    Save the sequence to a fasta file.
    """
    
    output_file = op.join(output_dir, f"{filename}_seq.txt")
    
    # If the file already exists, do nothing...
    if op.isfile(output_file):
        logger.debug("Sequence file {} already exists...".format(output_file))
        return output_file
    
    with open(output_file, "w") as ofh:
        ofh.write(str(sequence))
        
        
def read_train_data(data_path, columns=None):
    if columns is None:
        data = pd.read_csv(data_path)
    else:
        data = pd.read_csv(data_path, usecols=columns)
    
    return data


def test_mutate_sequence(original_sequence, mutated_sequence, mutation):
    """
    Test if the mutated sequence is the same as the original sequence.
    """
    if len(original_sequence) != len(mutated_sequence):
        raise ValueError("The mutated sequence is not the same length as the original sequence.")
    
    # which character is different between these two sequences?
    for position, (aa_1, aa_2) in enumerate(zip(original_sequence, mutated_sequence), start=1):
        if aa_1 != aa_2:
            print("The character at position {} is different between the two sequences: {} and {}".format(position, aa_1, aa_2))
            assert mutation == aa_1 + str(position) + aa_2, "The mutation is not correct"
            break

    print("\033[92m" + "The mutation is correct" + "\033[0m")


def save_prediction_data(benchmark_dir, prediction_file_name, prediction_data):
    file_date = datetime.now().strftime("%Y-%m-%d")
    prediction_file_name = "{}_{}.csv".format(prediction_file_name, file_date)
    prediction_data.to_csv(op.join(benchmark_dir, prediction_file_name), index=False)
    print("Prediction data `{}`is exported.".format(op.join(benchmark_dir, prediction_file_name)))
