
## Steps for generating sampling scheme benchmark results

### Clone minimizer benchmarking repo

```bash
git clone https://github.com/RagnarGrootKoerkamp/minimizers.git
cd minimizers
```

### Benchmarking on random strings
To generate benchmarking results for alphabet sizes of 2, 4, and 256:
```
mkdir output
for SIGMA in 2 4 256; do 
    RUSTFLAGS="-C target-cpu=native" cargo +nightly run -r -- -n 5000000 -s ${SIGMA} eval --practical -o output/practical-s${SIGMA}.json
done;
```

### Benchmarking on CHM13 ChrY

You can download the CHM13 genome, extract the Y chromosome, filter out all non ACTG characters,
and convert to uppercase with the following commands:
```bash
mkdir input
wget -O input/T2T-CHM13v2.0.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
seqkit grep -ip NC_060948.1 input-data/T2T-CHM13v2.0.fna.gz > input/chm13Y.fna
tail input/chm13Y.fna -n+2 | tr '[:lower:]' '[:upper:]' | tr -dc 'ATCG' > input/chm13Y.trimmed.txt
```

To generate benchmarking results for this specific input text, we use the `--input` flag:
```bash
RUSTFLAGS="-C target-cpu=native" cargo +nightly run -r -- -s 4 eval --practical -o output/practical-chm13Y.json --input input/chm13Y.trimmed.txt
```

### Analyzing output in `plots.ipynb`
In order to plot the benchmarking results, the output json files must be gzipped and 
moved to `scripts/data`

```
gzip output/*.json
cp output/*.json ../data/
```


## Instructions for running the ILP models
The ILP models are built with [gurobipy](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python).

The `run-ilp.py` script can construct and optimize a forward or local model. To run multiple models
with one command, you can supply multiple window sizes, _k_-mer sizes, and alphabet sizes. An ILP
is constructed for each combination of _w_, _k_, and sigma.
```bash
python run-ilp.py -w 2 3 4 -k 1 2 3 4 5 --sigma 2 3 4 --verbose
```

All options are listed with `--help`:
```bash
>$ python run-ilp.py --help
usage: run-ilp.py [-h] -w WINDOW_SIZE [WINDOW_SIZE ...] -k KMER [KMER ...] --sigma SIGMA [SIGMA ...] [--local] [-o OUTPUT] [--time-limit TIME_LIMIT] [-t THREADS] [-v]

options:
  -h, --help            show this help message and exit
  -w WINDOW_SIZE [WINDOW_SIZE ...], --window-size WINDOW_SIZE [WINDOW_SIZE ...]
  -k KMER [KMER ...], --kmer KMER [KMER ...]
  --sigma SIGMA [SIGMA ...]
  --local               Find minimum local scheme density
  -o OUTPUT, --output OUTPUT
                        Path to output directory
  --time-limit TIME_LIMIT
                        Time limit (in seconds)
  -t THREADS, --threads THREADS
  -v, --verbose         Log ILP to output
```
