#!/usr/bin/env bash
set -e

echo "Checking out minimizer benchmarking code..."
if [ ! -d minimizers/ ]; then
    git clone https://github.com/RagnarGrootKoerkamp/minimizers.git
fi
cd minimizers

export PYO3_PYTHON=$(which python3)
git checkout lower-bound-paper

mkdir -p output
echo "Running benchmarks on random strings..."
for SIGMA in 2 4 256; do
    RUSTFLAGS="-C target-cpu=native" cargo run -r -- -n 10000000 -s ${SIGMA} eval --practical -o output/practical-s${SIGMA}.json
done

mkdir -p input
echo "Downloading chm13..."
wget -O input/T2T-CHM13v2.0.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
seqkit grep -ip NC_060948.1 input/T2T-CHM13v2.0.fna.gz >input/chm13Y.fna
tail input/chm13Y.fna -n+2 | tr '[:lower:]' '[:upper:]' | tr -dc 'ATCG' >input/chm13Y.trimmed.txt

echo "Running benchmarks on chrY..."
RUSTFLAGS="-C target-cpu=native" cargo run -r -- -s 4 eval --practicaly -o output/practical-chm13Y.json --input input/chm13Y.trimmed.txt

gzip output/*.json
cd ..
cp minimizers/output/*.json.gz data/
