# Scientific integrity for SAMOVAR

**Conditional:** apply only when planning, executing, troubleshooting or reporting
a SAMOVAR scientific analysis. This is not a general rule for unrelated tasks.

- Check the request and available data before choosing scientific parameters.
  If materially different interpretations remain, ask for the missing choice;
  do not execute the dependent stage until it is resolved.
- Never invent input paths, databases, accession/taxon IDs, tool availability,
  metrics, citations or successful results. Inspect them or report them missing.
- For applied runs, a broken stage is a reported failure. Do not patch package
  code, substitute a simulator/annotator, discard samples/taxa, or bypass checks
  to make the run succeed unless the user authorizes that changed analysis.
  A development request may authorize fixing the underlying bug.
- Check requested outputs as well as process status. Explicitly report partial
  regeneration, missing/empty annotator tables, warnings and skipped validation.
  Do not present incomplete results as a completed ensemble or benchmark.
- Preserve the actual commands, parameters and provenance. Attribute tools using
  the run's methods/citations and verified references; never fabricate a citation.
