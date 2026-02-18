import requests
import json

UNIPROT_ACCESSION = "PDBID"


def get_apo_form(uniprot_accession):
    base_url = "https://data.rcsb.org/graphql"

    query = (
        """
    {
      entries(uniprotAccession: "%s", pfamAccession: "null") {
        rcsb_id
        rcsb_accession_info {
          initial_release_date
        }
        uniprots {
          rcsb_uniprot_container_identifiers {
            uniprot_id
          }
        }
      }
    }
    """
        % uniprot_accession
    )

    response = requests.post(base_url, json={"query": query})
    data = response.json()

    if data.get("data", {}).get("entries"):
        apo_form = sorted(
            data["data"]["entries"],
            key=lambda x: x["rcsb_accession_info"]["initial_release_date"],
        )[0]
        return apo_form["rcsb_id"]
    else:
        return None


if __name__ == "__main__":
    apo_pdb_id = get_apo_form(UNIPROT_ACCESSION)
    if apo_pdb_id:
        print(
            f"The apo form of the protein with UniProt accession {UNIPROT_ACCESSION} is: {apo_pdb_id}"
        )
    else:
        print(
            f"No apo form found for the protein with UniProt accession {UNIPROT_ACCESSION}"
        )
