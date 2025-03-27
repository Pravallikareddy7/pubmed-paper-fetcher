import os
import time
import csv
from typing import List, Dict
from Bio import Entrez, Medline
import requests
import re


# Define academic keywords
ACADEMIC_KEYWORDS = {"university", "college", "institute", "research center", "school", "academy", "hospital"}

# Set Entrez email (required for API usage)
Entrez.email = "slpravallika7@gmail.com"



def validate_query(query: str):
    """Ensures query is valid before sending to PubMed."""
    if not query.strip():
        raise ValueError("❌ Query cannot be empty.")
    if not re.search(r"[a-zA-Z]", query):
        raise ValueError("❌ Query must contain at least one letter.")
    return query


def fetch_pubmed_papers(query: str, max_results: int = 100) -> List[Dict]:
    """Fetches PubMed papers using Entrez API with error handling."""
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        if not record.get("IdList"):
            print("⚠️ No papers found for this query.")
            return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        papers = list(Medline.parse(handle))
        handle.close()

        return papers

    except requests.exceptions.RequestException as e:
        print(f"❌ Network error: {e}. Retrying in 5 seconds...")
        time.sleep(5)
        return fetch_pubmed_papers(query, max_results)  # Retry once

    except Exception as e:
        print(f"❌ Error accessing PubMed: {e}")
        return []

def is_industry_affiliation(affil: str) -> bool:
    """Checks if the given affiliation belongs to an industry."""
    affil_lower = affil.lower()
    return not any(kw in affil_lower for kw in ACADEMIC_KEYWORDS)



def extract_email(affiliations: List[str]) -> str:
    """Extracts an email address from affiliations using regex."""
    email_pattern = r"[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}"
    for affil in affiliations:
        emails = re.findall(email_pattern, affil)
        if emails:
            return emails[0]
    return "Not found"

def process_papers(papers: List[Dict], debug: bool = False) -> List[Dict]:
    """Processes papers and filters out those with industry affiliations."""
    results = []

    for paper in papers:
        if not paper.get("AD"):
            continue

        non_academic_authors = []
        company_affiliations = set()
        authors = paper.get("AU", [])
        affiliations = paper.get("AD", [])

        for affil in affiliations:
            if is_industry_affiliation(affil):
                company_affiliations.add(affil)
                non_academic_authors.extend(authors)

        if company_affiliations:
            result = {
                "PubmedID": paper.get("PMID", "N/A"),
                "Title": paper.get("TI", "No title").strip('"'),
                "PublicationDate": paper.get("DP", "Unknown").split()[0],  # Extract only the year
                "NonAcademicAuthors": "; ".join(set(non_academic_authors)),
                "CompanyAffiliations": "; ".join(company_affiliations),
                "CorrespondingEmail": extract_email(affiliations)
            }
            results.append(result)

        if debug:
            print(f"🔍 DEBUG: Checked affiliations for PMID {paper['PMID']}: {affiliations}")

    return results




def save_to_csv(results: List[Dict], filename: str):
    """Saves results to a CSV file with error handling."""
    if not results:
        print("⚠️ No papers with industry affiliations found.")
        return

    output_path = os.path.join(os.getcwd(), filename)

    try:
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"✅ Results saved to: {output_path}")

    except PermissionError:
        print(f"❌ Permission denied: Cannot write to {output_path}. Try a different location.")
    except Exception as e:
        print(f"❌ Error saving file: {e}")
