# Generate reads with InSilicoSeq
snakemake -s workflow/iss_test/Snakefile \
    --configfile workflow/iss_test/config.yaml \
    --cores 1

# Run annotators
snakemake -s workflow/annotators/Snakefile \
    --configfile workflow/annotators/config.yaml \
    --cores 1
    