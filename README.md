# CNVkit_manual
Run CNVkit after Sarek for specified contrasts

# 1. Copy files into WorkDir/
# profiles/ must come along — run_pipeline.sh defaults to profiles/slurm,
# which is where all per-rule resources live.
cp -r 00_prepare.sh config.yaml Snakefile run_pipeline.sh profiles ${WorkDir}/
cd ${WorkDir}/

# 2. One-time setup (symlinks + directory tree)
bash 00_prepare.sh

# 3. Download refFlat if missing (see output of step 2)
# Must match cnvkit_manual/reference/ — the path config.yaml's refflat: uses.
# (00_prepare.sh fetches the ENCODE blacklist itself.)
wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz \
    | gunzip > cnvkit_manual/reference/hg38.refFlat.txt

# 4. Dry-run to check the DAG
conda activate snakemake   # (or whichever env has snakemake)
bash run_pipeline.sh --dry-run

# 5. Launch
bash run_pipeline.sh
