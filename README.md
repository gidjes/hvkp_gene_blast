# HvKp Gene Blast

This is a pipeline to detect the presence of genes associated with Hypervirulent *Klebsiella pneumoniae* in a set of FASTA files. It uses blast to check all input FASTAs with a database of 8,718 gene variants. Finally, it determines the best hit from each gene group and provides presence/absence tables of the highest non-spurious hits.


## Installation

```bash
git clone https://github.com/gidjes/hvkp_gene_blast.git
conda env create -f environment.yml
conda activate hvkp_blast
unzip virulence_genes.zip
```

Ensure [blastn is installed](https://www.ncbi.nlm.nih.gov/books/NBK569861/) and available in your path if installation through conda fails.



## Usage

Place the FASTAs to check in a directory together.

subsequently run the following:
```bash
python scripts/main.py -i path/to/fastas/ 
```

## Parameters

The pipeline accepts the following command-line arguments:

| Parameter | Short Flag | Type | Default | Description |
|-----------|------------|------|---------|-------------|
| `--input` | `-i` | str | *required* | Input directory containing FASTA files |
| `--identity` | `-id` | int | 90 | Percent identity cut-off for detecting genes |
| `--coverage` | `-cv` | int | 80 | Gene coverage threshold (percent of gene length) |
| `--jobs` | `-n` | int | 1 | Number of parallel jobs to run |
| `--keep_intermediate_files` | `-ki` | flag | False | Keep intermediate CSV files instead of deleting them |


## Citation

If you use this pipeline in your research, please cite the following publication:

Vendrik KEW, Teunis G, Anema F, Landman F, de Haan A, Bos J, Witteveen S, Schoffelen AF, de Greeff SC, Kuijper EJ, Hendrickx APA, Notermans DW; hvKp study group. *Genomic epidemiology of putative hypervirulent Klebsiella pneumoniae species complex in Dutch patients, January-December 2022.* Microbiol Spectr. 2026 Feb 3;14(2):e0225925. doi: [10.1128/spectrum.02259-25](https://doi.org/10.1128/spectrum.02259-25).  


## License

MIT License

Copyright (c) 2026 Gijs Teunis

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.