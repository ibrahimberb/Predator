import requests

from mylogger import get_handler
import logging
import re

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


def get_gene_from_fasta_legacy(fasta_text):
    info_line = fasta_text.split('\n')[0]
    line_splitted = info_line.split()
    log.info(info_line.split())

    pattern = re.compile(r"^GN=(\w)+")
    [gene] = filter(pattern.match, line_splitted)
    gene = gene.replace("GN=", "")

    log.info(f"gene: {gene}")

    return gene


def get_gene_from_fasta(fasta_text):
    info_line = fasta_text.split('\n')[0]
    log.info(info_line)

    pattern = re.compile(r"GN=(\w)+ ")
    found_string = re.search(pattern, info_line).group(0)
    log.warning(found_string)

    gene = found_string.strip().replace("GN=", "")
    log.info(f"gene: {gene}")

    return gene


def get_gene_id_from_uniprot(uniprot_id):
    log.debug("Retrieving sequence {} ...".format(uniprot_id))
    address = "http://www.uniprot.org/uniprot/{}.fasta".format(uniprot_id)
    r = requests.get(address)
    if r.status_code == 200:
        gene1 = get_gene_from_fasta_legacy(r.text)
        gene2 = get_gene_from_fasta(r.text)

        assert gene1 == gene2

        return gene1

    else:
        raise IOError("Cannot reach the results.")


get_gene_id_from_uniprot("Q9UPU5")

"""
Create a table that stores entries for Uniprot Protein ID and corresponding Gene ID.
"""

test_dict = {
    "Q17R98": "ZNF827",
    "P24864": "CCNE1",
    "Q9H0D6": "XRN2",
    "O43149": "ZZEF1",
    "Q9P253": "VPS18",
    "Q68CR1": "SEL1L3",
    "P35499": "SCN4A",
    "Q96DU3": "SLAMF6",
    "Q6ZS81": "WDFY4",
    "Q9H009": "NACA2",
}

uniprot_ids = ["Q17R98", "P24864", "Q9H0D6", "O43149", "Q9P253", "Q68CR1", "P35499", "Q96DU3", "Q6ZS81", "Q9H009"]

for uniprot_id in uniprot_ids:
    assert get_gene_id_from_uniprot(uniprot_id) == test_dict[uniprot_id]

"""
Hugo_Symbol	SWISSPROT
ZNF827	Q17R98
CCNE1	P24864
XRN2	Q9H0D6
ZZEF1	O43149
VPS18	Q9P253
SEL1L3	Q68CR1
SCN4A	P35499
SLAMF6	Q96DU3
WDFY4	Q6ZS81
NACA2	Q9H009
"""
