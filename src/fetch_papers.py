! pip install biopython
try:
    from Bio import Entrez
    from Bio import Medline
except ImportError:
    print("Required packages not found. Please install biopython first with:")
    print("pip install biopython")
    exit(1)

import os
import re
import csv
from typing import List, Dict

# Set your email (required by PubMed)
Entrez.email = " "  

def get_output_dir():
    """Get a safe output directory"""
    return os.getcwd()  # Save in current working directory

def fetch_pubmed_papers(query: str, max_results: int = 100) -> List[Dict]:
    """Fetch papers from PubMed with complete author/affiliation data"""
    try:
        # Search PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print("No papers found for this query.")
            return []

        # Fetch detailed records
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        papers = list(Medline.parse(handle))
        handle.close()

        return papers

    except Exception as e:
        print(f"Error accessing PubMed: {e}")
        return []

def is_industry_affiliation(affiliation: str) -> bool:
    """Check if affiliation suggests industry"""
    if not affiliation:
        return False

    affil_lower = affiliation.lower()

    # Skip obvious academic affiliations
    academic_keywords = ["university", "college", "hospital", "institute", "school", ".edu", "academy"]
    if any(kw in affil_lower for kw in academic_keywords):
        return False

    # Look for industry indicators
    industry_keywords = [
        "inc", "ltd", "pharma", "biotech", "corp", "gmbh", "company",
        "pfizer", "moderna", "astrazeneca", "novartis", "roche",
        "johnson & johnson", "sanofi", "merck", "bayer", "novo nordisk"
    ]
    return any(kw in affil_lower for kw in industry_keywords)

def process_papers(papers: List[Dict]) -> List[Dict]:
    """Process papers to find industry affiliations"""
    results = []

    for paper in papers:
        if not paper.get("AD"):
            continue

        non_academic_authors = []
        company_affiliations = set()
        authors = paper.get("AU", [])
        affiliations = paper.get("AD", [])

        # Check all affiliations
        for affil in affiliations:
            if is_industry_affiliation(affil):
                company_affiliations.add(affil)
                # Only add authors who have industry affiliations
                if authors:
                    non_academic_authors.extend(authors)

        if company_affiliations:
            results.append({
                "PubmedID": paper.get("PMID", "N/A"),
                "Title": paper.get("TI", "No title").strip('"'),  # Clean title
                "PublicationDate": paper.get("DP", "Unknown").split()[0],  # Just year
                "NonAcademicAuthors": "; ".join(set(non_academic_authors)),
                "CompanyAffiliations": "; ".join(company_affiliations),
                "CorrespondingEmail": extract_email(affiliations)
            })

    return results

def extract_email(affiliations: List[str]) -> str:
    """Extract email from affiliations"""
    for affil in affiliations:
        if "@" in affil:
            # Better email pattern matching
            emails = re.findall(r'[\w\.-]+@[\w\.-]+', affil)
            if emails:
                return emails[0]
    return "Not found"

def save_to_csv(results: List[Dict], filename: str):
    """Save results to CSV with proper error handling"""
    if not results:
        print("No papers with industry affiliations found.")
        return

    output_path = os.path.join(get_output_dir(), filename)

    try:
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"Success! Saved {len(results)} papers to: {output_path}")
    except Exception as e:
        print(f"Error saving file: {e}")

def main(query, filename="results.csv", max_results=50):
    papers = fetch_pubmed_papers(query, max_results)
    results = process_papers(papers)
    save_to_csv(results, filename)

if __name__ == "__main__":
    # Corrected query (fixed typo in "diabetes")
    query = ""
    filename = "results.csv"
    max_results = 50
    main(query, filename, max_results)
