import requests
from pathlib import Path

import chem
import filing

RECEPTORS_OF_APO_PDBS = {"3rzy": "fabp4", "4lkp": "fabp5"}
REF_COMPLEXES = {"fabp4": "2nnq", "fabp5": "5hz5"}
REF_LIGANDS = {"2nnq": "T4B", "5hz5": "65X"}


def ref_complex(receptor_path: Path) -> str:
    pdb_id = receptor_path.stem.split('_')[0].lower()
    if pdb_id in RECEPTORS_OF_APO_PDBS.keys():
         return REF_COMPLEXES[RECEPTORS_OF_APO_PDBS[pdb_id]]
    return pdb_id


def ref_ligand(receptor_path: Path) -> str:
    return ref_complex(receptor_path) + chem.COMPLEX_LIGAND_FILE_ENDING


def local_data_dir():
    return filing.local_data_dir(__file__)


def __protein_query(gene):
    gene = gene.lower()
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_phrase",
                        "value": gene,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_phrase",
                        "value": "apo",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                        "operator": "exact_match",
                        "value": "Homo sapiens",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [
                {"sort_by": "em_3d_reconstruction.resolution", "direction": "asc"},
                {"sort_by": "score", "direction": "desc"},
            ],
            "scoring_strategy": "combined",
        },
    }


def search_protein(search_term):
    query_data = __protein_query(search_term)
    response = requests.post(
        "https://search.rcsb.org/rcsbsearch/v2/query", json=query_data
    )
    results = response.json()
    if "total_count" in results and results["total_count"] > 0:
        return results["result_set"][0]["identifier"]
    else:
        return None


def write_pdb(pdb_id, path=local_data_dir()):
    pdb_filename = pdb_id + ".pdb"
    output_path = Path(path) / pdb_filename

    with open(output_path, "wb") as output_file:
        output_file.write(
            requests.get("https://files.rcsb.org/download/" + pdb_filename).content
        )
    return output_path


if __name__ == "__main__":
    write_pdb("2nnq")
