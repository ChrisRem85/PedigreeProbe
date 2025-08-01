# PedigreeProbe

PedigreeProbe is a comprehensive genetic relationship analysis tool for pedigree (PED) files. It calculates kinship coefficients, detects consanguinity (including complex scenarios like double cousins and pedigree collapse), and computes inbreeding coefficients. The tool outputs detailed relationship paths and relatedness for all individuals in a pedigree.

## Features

- Parses standard PED files
- Detects and classifies consanguinity (including double cousins, avuncular, pedigree collapse, and more)
- Calculates kinship and inbreeding coefficients for each individual
- Outputs detailed relationship paths and relatedness
- Validates pedigree structure for errors and impossible relationships
- Provides enhanced output with relationship descriptions

## File Structure

- `extract_relationships_from_ped.py` — Main script for relationship extraction and analysis
- `tests/` — Folder containing test input/output files and automated test script
  - `example_no_consanguinity.ped` — Example PED input file
  - `example_no_consanguinity_with_relationships.ped` — Expected output file with relationships
  - `test_extract_relationships.py` — Automated test (pytest)

## Usage

### Command Line

```sh
python extract_relationships_from_ped.py <input_ped> <output_file> [--proband PROBAND_ID] [--proband-phenotype PHENOTYPE]
```

- `<input_ped>`: Path to the input PED file
- `<output_file>`: Path to write the enhanced output
- `--proband`: (Optional) Specify the proband individual ID
- `--proband-phenotype`: (Optional) Auto-detect proband by phenotype value

Example:

```sh
python extract_relationships_from_ped.py tests/example_no_consanguinity.ped tests/test_output.ped
```

## Output

The output file contains columns:

- FAM_ID
- IND_ID
- FATHER_ID
- MOTHER_ID
- SEX
- PHENOTYPE
- RELATIONSHIP (description and path)
- RELATEDNESS (coefficient)
- INBREEDING_COEFF
- CONSANGUINITY_TYPE

## Testing

Tests are provided in the `tests/` folder. To run the automated test:

1. Install pytest if needed:
   ```sh
   pip install pytest
   ```
2. Run the test:
   ```sh
   pytest tests/test_extract_relationships.py
   ```

The test will run the script on the example PED file and compare the output to the expected relationships file.

## License

MIT License
