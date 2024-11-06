# Density analysis of sampling schemes

This repository contains code to go with the paper ([bioRxiv](https://www.biorxiv.org/content/10.1101/2024.09.06.611668v1))
> A near-tight lower bound on the density of forward sampling schemes

1. A script ([run-benchmarks.sh](./run-benchmarks.sh)) to benchmark existing sampling schemes, using the repo [RagnarGrootKoerkamp/minimizers](https://github.com/RagnarGrootKoerkamp/minimizers).
2. Code to run our ILP (integer linear program, [run-ilp.py](./run-ilp.py)) that searches for optimal
   sampling schemes for small parameters.
3. A python notebook ([plots.ipynb](./plots.ipynb)) that plots all results and
   lower bounds.


## Requirements
### Python
* matplotlib
* pandas
* numpy
* sympy
* gurobipy ([requires a license, free for academics](https://www.gurobi.com/features/academic-wls-license/))

### Other
* [seqkit](https://bioinf.shenwei.me/seqkit/)
* [rust nightly](https://doc.rust-lang.org/book/appendix-07-nightly-rust.html#rustup-and-the-role-of-rust-nightly)

## Running the benchmarks
Benchmarks can be generated via
```
./run-benchmarks.sh
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
