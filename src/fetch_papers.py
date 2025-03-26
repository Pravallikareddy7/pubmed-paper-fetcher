import csv
import argparse
import os
from Bio import Entrez
from Bio import Medline
from typing import List, Dict

# Set your email (required by PubMed)
Entrez.email = "your.email@example.com"  # Replace with your actual email

def get_output_dir():
    """Get a safe output directory"""
    # Try the user's documents folder first
    docs_dir = os.path.expanduser("~/Documents")
    if os.path.exists(docs_dir) and os.access(docs_dir, os.W_OK):
        return docs_dir
    # Fall back to current directory
    return os.getcwd()

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
    academic_keywords = ["university", "college", "hospital", "institute", "school", ".edu"]
    if any(kw in affil_lower for kw in academic_keywords):
        return False
    
    # Look for industry indicators
    industry_keywords = [
        "inc", "ltd", "pharma", "biotech", "corp", "gmbh", "company",
        "pfizer", "moderna", "astrazeneca", "novartis", "roche", "johnson & johnson"
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
        
        # Combine authors with affiliations (simplified approach)
        affiliations = paper.get("AD", [])
        authors = paper.get("AU", [])
        
        # Check all affiliations
        for affil in affiliations:
            if is_industry_affiliation(affil):
                company_affiliations.add(affil)
                # Assign all authors as non-academic (simplified)
                non_academic_authors.extend(authors)
        
        if company_affiliations:
            results.append({
                "PubmedID": paper.get("PMID", "N/A"),
                "Title": paper.get("TI", "No title"),
                "PublicationDate": paper.get("DP", "Unknown"),
                "NonAcademicAuthors": "; ".join(set(non_academic_authors)),
                "CompanyAffiliations": "; ".join(company_affiliations),
                "CorrespondingEmail": extract_email(affiliations)
            })
    
    return results

def extract_email(affiliations: List[str]) -> str:
    """Extract email from affiliations"""
    for affil in affiliations:
        if "@" in affil:
            for word in affil.split():
                if "@" in word and "." in word:
                    return word.strip(".,;")
    return "Not found"

def save_to_csv(results: List[Dict], filename: str):
    """Save results to CSV with proper error handling"""
    if not results:
        print("No papers with industry affiliations found.")
        return
    
    # Ensure we have a valid output directory
    output_dir = get_output_dir()
    output_path = os.path.join(output_dir, filename)
    
    try:
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"Success! Saved {len(results)} papers to:\n{output_path}")
    except PermissionError:
        print(f"Error: Cannot write to {output_path}\n"
              f"Please specify a different output location using --file")
    except Exception as e:
        print(f"Error saving file: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Fetch PubMed papers with industry affiliations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("query", help="PubMed search query")
    parser.add_argument("-f", "--file", default="results.csv", 
                       help="Output CSV filename (will be saved in your Documents folder)")
    parser.add_argument("-m", "--max", type=int, default=100, 
                       help="Max results to fetch")
    args = parser.parse_args()
    
    papers = fetch_pubmed_papers(args.query, args.max)
    results = process_papers(papers)
    save_to_csv(results, args.file)

if __name__ == "__main__":
    main()