import argparse
import requests
import re
import csv
import os
from typing import List, Dict
from Bio import Entrez, Medline

# Set your email (required by PubMed API)
Entrez.email = "your_email@example.com"  # Replace with your email

# Industry keywords for filtering
INDUSTRY_KEYWORDS = [
    "inc", "ltd", "pharma", "biotech", "corp", "gmbh", "company",
    "pfizer", "moderna", "astrazeneca", "novartis", "roche",
    "johnson & johnson", "sanofi", "merck", "bayer", "novo nordisk"
]

# Academic keywords to filter out
ACADEMIC_KEYWORDS = ["university", "college", "hospital", "institute", "school", ".edu", "academy"]

def fetch_pubmed_papers(query: str, max_results: int = 100) -> List[Dict]:
    """Fetches PubMed papers using Entrez API"""
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print("⚠️ No papers found for this query.")
            return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        papers = list(Medline.parse(handle))
        handle.close()

        return papers
    except Exception as e:
        print(f"❌ Error accessing PubMed: {e}")
        return []

def is_industry_affiliation(affiliation: str) -> bool:
    """Checks if an affiliation is industry-related."""
    if not affiliation:
        return False
    affil_lower = affiliation.lower()
    
    # Exclude academic institutions
    if any(kw in affil_lower for kw in ACADEMIC_KEYWORDS):
        return False

    # Match industry-related terms
    return any(kw in affil_lower for kw in INDUSTRY_KEYWORDS)

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
    """Saves results to a CSV file."""
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
    except Exception as e:
        print(f"❌ Error saving file: {e}")

def main():
    """Main function to run the program from CLI."""
    parser = argparse.ArgumentParser(description="Fetch and filter PubMed papers based on industry affiliations.")
    parser.add_argument("--query", type=str, required=True, help="Search query for PubMed.")
    parser.add_argument("--file", type=str, default="results.csv", help="Output CSV file name.")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode.")

    args = parser.parse_args()

    print(f"🔎 Searching PubMed for: {args.query}")
    papers = fetch_pubmed_papers(args.query)
    results = process_papers(papers, args.debug)

    if args.debug:
        print("\n🔍 DEBUG MODE: Printing Extracted Results\n")
        import json
        print(json.dumps(results, indent=2))

    save_to_csv(results, args.file)

if __name__ == "__main__":
    main()

