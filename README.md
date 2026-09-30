# untangler
Assembler for DNA, aligner for RNA against DNA, untangler for RNA subclusters

## Install from PyPI

```bash
python -m pip install untangler
```

## Install from source

```bash
python -m pip install .
```

For local development, install the editable package:

```bash
python -m pip install -e .
```

## Build and upload

```bash
python -m pip install build twine
python -m build
python -m twine upload dist/*
```

## Usage

```bash
untangler cluster_transcripts --ids PF3D7_1007900 \
	-i config.txt \
	-o ./cluster_transcripts_PF3D7_1007900_28 \
	--cluster_cols introns,m6A \
	--predict introns \
	-t 0.6 -f --generate_igv_reports
```