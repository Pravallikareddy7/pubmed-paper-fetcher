import argparse
from .pubmed_fetcher import fetch_pubmed_papers, process_papers, save_to_csv,validate_query

def main():
    """Main function to run the program from CLI."""
    parser = argparse.ArgumentParser(description="Fetch and filter PubMed papers based on industry affiliations.")
    parser.add_argument("--query", type=str, required=True, help="Search query for PubMed.")
    parser.add_argument("--file", type=str, default="results.csv", help="Output CSV file name.")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode.")

    args = parser.parse_args()

    try:
        query = validate_query(args.query)  # Validate query before sending to PubMed
    except ValueError as e:
        print(f"❌ Error: {e}")
        return

    print(f"🔎 Searching PubMed for: {query}")
    papers = fetch_pubmed_papers(query)
    results = process_papers(papers, args.debug)

    if args.debug:
        print("\n🔍 DEBUG MODE: Printing Extracted Results\n")
        import json
        print(json.dumps(results, indent=2))

    save_to_csv(results, args.file)

if __name__ == "__main__":
    main()