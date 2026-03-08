# Mass_blast

This is a pipeline in snakemake in order to run a blast on a set of large sequences and then organise these outputs.


## Installation

```bash
git clone https://github.com/gidjes/hvkp_gene_blast.git
conda env create -f environment.yml
conda activate hvkp_blast
```

Also ensure [blastn is installed](https://www.ncbi.nlm.nih.gov/books/NBK569861/) and available in your path/environment.

## Usage

Place the query (input) sequences in the 'input' directory and the subject (target) sequences in the 'target' directory. If you want to blast the sequence set to itself then leave the 'target' directory empty.

subsequently run the following:
```bash
python scripts/main.py -i path/to/fastas/ 
```

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