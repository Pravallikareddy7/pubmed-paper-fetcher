# PubMed Paper Fetcher

## Overview
This project is a command-line tool that fetches research papers from the **PubMed API** based on user queries. It filters papers that have at least one author affiliated with a **pharmaceutical/biotech company** and exports the results in a CSV format.

## Features
- Fetches research papers using **PubMed API**.
- Filters papers with at least one **non-academic (industry-affiliated) author**.
- Saves results in a CSV format with relevant metadata.
- Supports command-line options for user queries.

## Installation
This project uses **Poetry** for dependency management. Follow these steps to install and use the tool.

### Prerequisites
- Python 3.8+
- Poetry installed (if not, install it using `pip install poetry` or refer to [Poetry installation guide](https://python-poetry.org/docs/#installation)).

### Setup
1. Clone the repository:
   ```sh
   git clone https://github.com/Pravallikareddy7/pubmed-paper-fetcher.git
   cd pubmed-paper-fetcher
   ```
2. Install dependencies using Poetry:
   ```sh
   poetry install
   ```
3. Activate the virtual environment:
   ```sh
   poetry shell
   ```

## Usage
Run the command-line tool with the required parameters:

### Basic Query Execution
```sh
get-papers-list --query "cancer treatment"
```
This fetches research papers related to "cancer treatment" and filters for industry-affiliated authors.

### Query from a File
```sh
get-papers-list --file queries.txt
```
This fetches papers for multiple queries listed in `queries.txt`.

### Debug Mode (For Detailed Logging)
```sh
get-papers-list --query "cancer treatment" --debug
```
This enables debug mode to show detailed logs.

## Output Format
The results are stored in a CSV file with the following columns:
| PubmedID | Title | Publication Date | Non-academic Authors | Company Affiliations | Corresponding Author Email |
|----------|-------|------------------|----------------------|----------------------|---------------------------|

Example output (CSV):
```
12345678, "New Cancer Treatment", 2024-03-01, "Dr. John Doe", "Pfizer", "john.doe@pfizer.com"
23456789, "AI in Biotech", 2024-02-15, "Dr. Alice Smith", "Moderna", "alice.smith@moderna.com"
```

## Development
### Running Tests
Run unit tests using:
```sh
pytest tests/
```

## Publishing to TestPyPI (Optional)
To publish this package to TestPyPI:
1. Build the package:
   ```sh
   poetry build
   ```
2. Publish to TestPyPI:
   ```sh
   poetry publish --repository testpypi
   ```

To install from TestPyPI:
```sh
pip install --index-url https://test.pypi.org/simple/ --no-deps pubmed-paper-fetcher
```

## Contributing
If you would like to contribute:
1. Fork the repository.
2. Create a new branch.
3. Make your changes and test.
4. Submit a pull request.

## License
This project is licensed under the **MIT License**.

## Contact
For any issues or feature requests, please open an issue in the repository.
E-mail:slpravallika7@gmail.com
