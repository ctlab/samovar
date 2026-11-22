# SAMOVAR behaviour constitution

Status: agreed target behaviour. This document records the agreed
behaviour of SAMOVAR; it is not a statement that every capability below is
already implemented. Implementation design gates and acceptance evidence are
tracked in [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md).

## 1. Contract-driven pipeline

Pipeline stages consume documented input/output contracts. A compatible
built-in or imported tool can implement a contract without changing the
pipeline's interpretation of samples, truth, or results. Existing registration,
tool resolution, stage execution, and contract tests remain the integration
points; a new feature must not introduce a parallel pipeline or registry.

The logical flow is:

```text
generate or supplied inputs
  -> initial reads -> optional pool/hybrid routing
  -> annotation / features / profiles / assembly and mapping
  -> initial combined annotations and abundance
  -> regenerated tables -> regenerated reads
  -> the corresponding routing and annotation / profile / assembly contracts
  -> regenerated combined annotations and abundance
  -> optional training-table resampling
  -> TBML training -> initial read-level prediction -> SAMOVAR abundance
  -> TPML training -> initial profile prediction -> reprofiled abundance
```

Scoring and comparison attach to applicable artifacts, rather than defining
another interpretation of those artifacts.

## 2. Samples, libraries, and read identity

- `sample_id` identifies a biological sample. One sample can have several
  libraries, including Illumina and ONT. `library_id` distinguishes libraries;
  technology and mate layout are explicit properties, not sample identities.
- A long input manifest has one row per file, beginning with sample identity
  and file path. It records library identity, technology, layout/mate, and the
  information needed to recover simulation settings. Sample IDs can repeat
  across rows belonging to different libraries or mates.
- Interleaved FASTQ is an input representation. Normalize it to paired-end
  FASTQ before downstream contracts, validating mate identity and order.
  Single-end and paired-end libraries remain distinguishable.
- Read identity must remain unambiguous across samples, libraries, pools,
  bootstrap draws, and repeats. Pooling must retain reversible membership;
  colliding original read IDs must not merge unrelated reads. Paired mates
  retain a shared fragment identity.
- Independent generate invocations must preserve earlier generated libraries,
  configuration, and provenance. Generated inputs automatically contribute
  their manifest records. Adding a generation/library requires an explicit
  name. Repeating the same effective configuration resumes that generation;
  a conflicting configuration for an existing name fails without overwrite.
- Without a routing manifest, existing direct SE/PE execution remains usable.

## 3. Pooling, hybrids, and assembly

- Each tool invocation selects its pool and/or hybrid manifest columns.
  Routing can apply to a contract, stage, tool, or their combination; it must
  not silently change unrelated branches.
- Missing values in the selected hybrid column exclude those rows from that
  invocation. A simple selected/not-selected role is filtering and does not
  require native multi-input hybrid support.
- Native hybrid input roles map manifest values to the concrete tool's input
  flags. Human-readable roles and numeric aliases must have an explicit
  interpretation. They are not permission to send unsupported read types to
  a program. Native capabilities and dictionary registration must be validated
  for the selected tool version and mode.
- The first version does not mix pooled and unpooled inputs within the same
  routed branch. Such a mixed request must not silently pool, discard, or
  individually execute the unpooled samples.
- When both operations are selected, pool first, then construct hybrid inputs.
  Pooling concatenates compatible libraries while keeping R1 and R2 separate
  and synchronized. It must not flatten all technologies into one input.
  The grouping key includes pool identity and the selected hybrid role;
  incompatible technologies/layouts remain separate. Four samples assigned
  to two pools produce two invocations, each with its own separate short-read
  and long-read inputs. Pooling does not concatenate across hybrid flags.
- Assembly and binning may operate on a shared pool. Mapping operates on the
  original samples/libraries against the resulting reference, retaining their
  biological sample identities.
- A pooled taxonomic binner's per-read results are depoolable before they
  enter the ordinary per-sample annotation interface. Native sample-level
  profiler results cannot be depooled into invented read assignments or
  per-sample profiles. Profiler depooling is outside the first version.
- Regenerated MAGs can be matched to initial MAGs through a replaceable
  matching contract. Alternatively, regenerated reads can map directly to
  the frozen initial MAG reference. IDs from unrelated catalogs are not
  equivalent merely because their labels match.
- Existing tetranucleotide MAG matching is not ANI. Any matching method must
  expose its actual method and preserve unresolved assignments.
  ANI matching is an alternative implementation of the same MAG matching
  contract, not a second identity/mapping pipeline.

Register native hybrid support and its role-to-input dictionary on the tool,
using `tools import --hybrid` and the existing tool/config mechanism. Built-ins
describe the same capabilities. The dictionary relates named manifest roles
and any numeric aliases to native input arguments, including R1/R2/SE handling;
it must not be a second shell-command or arbitrary-code format. Ask for the
dictionary if an interactive import cannot resolve it; noninteractive use
requires sufficient explicit input. A program's
acceptance of several files does not imply it can mix technologies; metaFlye
and metaSPAdes capabilities must be checked for their selected versions and
modes rather than inferred from their names.

## 4. Annotation and profile boundaries

- Taxonomic annotation of individual reads/fragments is taxonomic binning.
  A native taxonomic profiler instead returns sample-level taxon abundance.
  Both contribute profiles, but their per-read capabilities remain distinct.
- One annotator can emit several named taxonomic outputs, feature outputs,
  or both. Each output can have its own optional confidence. Confidence must
  remain bound to that exact tool output through parsing, combination,
  persistence, and model input preparation. Missing confidence is not zero.
- Feature extraction has no implicit taxonomic ground truth. Feature columns
  and confidence columns must not be scored as taxonomic predictions.
- Reuse the canonical annotation names (`seq`, `sample`, `true`, tool-specific
  `taxID_*`/`feat_*`) and typed profiles (`taxid`, `N_<sample>` plus metadata).
  Final confidence naming is part of the pending annotation-contract update;
  examples in the design sketch do not rename existing columns.
- Combined annotation tables preserve read/sample identity and the distinct
  outputs of each tool/database invocation. Incremental ingestion must not
  discard previously completed columns.
- Annotation objects, combined tables, and abundance conversion share the
  same interpretation of predictions. A profiler output must not manufacture
  per-read annotations to pass through this boundary.
- Consensus confidence policies resolve vote-count ties only, among the tied
  candidates. A unique majority is unchanged. Strategies are mean internal
  annotator confidence (baseline), mean quality from the Scoring subproject,
  and the most confident annotator for the current read, falling back to the
  configured mean strategy when local confidence is unavailable. Do not sum
  confidence votes. An unresolved tie or unusable confidence yields taxid 0.
  FMA representation-based resolution is a separate dependency-gated option.

## 5. Truth, scoring, and comparison

- Truth has explicit source, sample, stage, taxonomy/rank, and measurement
  basis. Distinct truth sources are evaluated separately; they are not merged
  into a new gold standard.
- Read-level truth can come from known simulation provenance or an external
  trustworthy read-to-taxid table, including CAMI. A sample-level abundance
  or presence/absence reference does not imply per-read truth.
- Taxonomic binning scoring and taxonomic profiling scoring are separate
  contracts and reports. Binning evaluates read/fragment assignments;
  profiling evaluates taxon presence and quantitative sample profiles.
  Abundance R-squared belongs to profiling, not binning.
- At a selected rank, binning follows the Kraken 2 assignment contract:
  distinguish correct, incorrect, unclassified and vague ancestor assignments.
  Sensitivity is TP divided by all evaluated fragments; PPV is TP/(TP+FP),
  and F1 is their harmonic mean. Report counts (including FN/VP), coverage,
  rank and truth source. Unknown or insufficient truth is not a negative
  example. Paired mates do not become two independent fragments. Generic
  closed-set micro/accuracy and cluster-purity diagnostics are separately
  named; purity is not primary taxonomic accuracy.
- Profiling reports R-squared, presence precision/recall/F1/Jaccard, and
  quantitative L1/Bray-Curtis when the reference and profile support them.
  These metrics share one numerical profiling implementation. Presence-only
  truth does not enable R-squared or other quantitative metrics.
- Distinguish the requested pre-simulation community profile from the realized
  composition of the reads sent to the tools, including the effect of QC.
  Both are independently identifiable references. For read-derived counts,
  the primary reference is the realized read/fragment composition. Evaluation
  against the requested profile records its own source and measurement basis.
- Profiling scoring applies wherever abundance tables exist: initial,
  regenerated, and reprofiled. Regenerated simulation has known truth; initial
  and reprofiled stages may lack it. Lost regenerated truth is a pipeline
  error, not permission to invent labels or silently claim an unscored success.
- Without applicable truth, run supported comparisons and report
  truth-dependent metrics as unavailable, with a reason. Missing data is not
  a measured zero. Presence-only truth cannot support quantitative metrics.
- Compare compatible ranks, taxonomies, units, bases, and normalizations.
  Record conversions. Preserve false predicted taxa and unclassified mass;
  display-only grouping must not change numerical metrics.
- The selected scorer supplies numerical results once. Reports and MultiQC
  consume those results without rerunning an equivalent calculation. External
  OPAL is an optional alternative or explicitly requested cross-check.
- Version result semantics. Read old artifacts with their historical meaning;
  a legacy metric must not be renamed to a different definition without
  recalculation.

Undefined metric values retain their reason; the numerical edge-case rules
must agree with the versioned scorer contract.

## 6. Regeneration and repeat reproducibility

- Regenerate community composition from predicted initial abundance, not
  from initial evaluation truth. Preserve the input library structure,
  technologies, relevant simulation parameters, and routing relationships.
  Regeneration depth is independently configurable.
- The default intent is to reproduce the initial sample/library structure.
  Alternative sample/pool-level regeneration and explicit output counts
  require an agreed contract; sample counts, pool counts, synthetic community
  counts, and repeats are not interchangeable quantities.
- Table generation mode determines output community cardinality. `direct`
  preserves source profiles and sample identities. `concat` and
  `concat_weighted` combine both tools and samples into one synthetic
  community. Equal-taxon versus abundance-weighted composition remains their
  distinction. `concat_weighted` averages normalized relative abundance across
  every non-empty tool/sample column. Both modes allocate one community column
  totalling exactly `N_reads` when it is set.
- A single compatible library/routing template is inherited by a combined
  community. Heterogeneous conflicting templates require an explicit template;
  do not guess a biological parent or flatten incompatible libraries.
- Every table-generation mode must document whether it preserves, combines,
  or creates communities and how the resulting communities acquire library
  templates and routing. Do not impose the direct-mode cardinality on
  bootstrap, VAE, GLM, CAMISIM, R, or imported generators. An unsupported mode
  and library/routing combination must fail explicitly rather than flatten
  libraries or silently substitute another mode.
- Regenerated output is SE or PE, not interleaved. Its manifest describes the
  generated libraries and supports the same downstream routing as initial.
- Repeat generation/bootstrap receives reproducible repeat-specific random
  streams. Downstream seeds can be shared between repeats or explicitly set
  by stage/tool; a shared downstream seed is valid. Record the base seed and
  effective stage/tool seeds in manifests and Hydra snapshots. A stored seed
  must describe the value actually passed to that tool. Tools without a seed
  interface must not be reported as having accepted one.
- Benchmark mode ends after applicable initial annotations, profiles,
  exports, scores, comparisons, and reports. It does not regenerate reads,
  train correction models, or reprofile. Ready-read bootstrap preserves PE
  fragments and propagates known truth to uniquely identified output draws.
- Repeat summaries distinguish repeats from samples, preserve individual
  results, report failures/incompleteness, and expose the statistical method
  and seed used for uncertainty intervals.
- Without `--repeat`, preserve the ordinary run layout. Explicit `--repeat 1`
  already selects batch layout: `.generate/`, `.log/repeats.json`, independent
  `out_0001/` runs, and `summary/` with long result/profile tables, statistics
  and `multiqc/`. Prepare inherits the manifest; an explicit conflicting count
  fails before execution. There is no nested N×N repetition.
- Repeat generation redraws generated abundance and reads. User-supplied
  abundance stays fixed while reads are simulated again. Ready FASTQ repeats
  require real bootstrap with replacement, preserving complete PE fragments,
  unique draw IDs and truth lineage; identical reruns are not a bootstrap.
- Execute repeats sequentially within the resource budget. A failed repeat
  does not prevent the others; the final exit status is nonzero and summary
  identifies successful/expected counts and incompleteness. Resume preserves
  completed work; summary can be rebuilt from stored results.
- For comparable sample/tool/stage/rank/metric groups, report mean, SD and
  percentile bootstrap 95% CI of the mean using 10,000 resamplings of whole
  repeats with a recorded seed. With one successful repeat, SD and CI are
  unavailable. Reads are not independent repeat units. Fixed-abundance
  intervals describe simulation variability, not biological variability.
  MultiQC exposes repeat points, mean/CI, missing results and counts; profile
  plots use consistent taxonomic ranks and top-N. Missing results are not zero.
- Processed genomes default to the shared SAMOVAR Database. A run uses its
  selected genome set; it must not implicitly import the whole library.
  Necessary local reference views use links without silently copying large
  genomes. New defaults do not create ordinary output `genomes/` directories.
  Explicit run-local storage and reading existing runs remain supported.

## 7. Resampling and model contracts

- Optional resampling selects whole rows from regenerated combined
  annotations to approximate selected initial distributions. Taxonomic labels
  are categorical; numeric features can use empirical distributions. Keep
  predictions, features, confidence, identity/provenance, and `true` together.
- Approximate the selected distributions jointly and report residual
  mismatch. Do not claim exact matching when the regenerated data cannot
  represent the target. Resampling must not modify original annotations.
- Resampling uses explicit groups. Preserve the number of regenerated rows
  within each group and sample rows with replacement to approach that group's
  initial distribution. A group can cover all samples. Library/sample
  provenance remains attached to draws; a duplicate draw is not a new source
  read or a new biological sample.
- Assign each initial group explicitly to its regenerated training dataset
  and size. A common regenerated dataset may serve several groups; preserve
  each declared regenerated group size. Global concat does not imply a
  one-to-one biological lineage back to the original samples.
- TBML trains on exactly one selected regenerated table: resampled when that
  option is active, otherwise the original. `true` is the supervised target,
  never an input feature. Distribution matching does not use evaluation truth
  or a SAMOVAR prediction as a target feature.
- A fitted TBML can predict any compatible combined annotation table with or
  without truth. Preserve existing columns, including optional truth, and add
  SAMOVAR taxid and confidence. No usable model means no model prediction;
  copying truth into predicted SAMOVAR columns is forbidden.
- SAMOVAR abundance is the aggregation of TBML read predictions. TPML instead
  trains from compatible sample profiles and their known abundance target,
  then predicts reprofiled abundance from initial profiles without the target
  as an input. These are distinct outputs, even when compared in one report.
- Models can request additional initial/regenerated context through an
  explicit contract. Evaluation truth must not leak into prediction inputs.
  The regeneration-chain model below is a distinct inference family; do not
  force it to discard its provenance merely to resemble a single-table TBML
  classifier or a profile-only TPML regressor.
- Model schema, fitted preprocessing, parameters, uncertainty/abstention
  semantics, and provenance travel with the saved model and survive chunked
  prediction. Missing support must not become fabricated confidence.

### Regeneration-chain inference

The Markov/hidden-state model follows the data-generating chain, not FASTQ row
order. An initial tool predicts a taxid; that prediction contributes a source
profile; regeneration resolves the requested taxid to a genome and simulates
reads with known actual truth; tools assign taxa to those regenerated reads.
Learn transition/error behaviour from these linked observations, retaining
the source tool/database, generation mode, library, and requested versus actual
taxon. Inference uses this behaviour to reconstruct the hidden original taxon
of initial reads from the observed initial tool outputs.

The required outputs are a per-read posterior with uncertainty/abstention and
an aggregated corrected sample profile. Separately represent uncertainty from
finite regenerated transition counts and from the requested-taxon-to-selected-
genome correspondence. The requested taxid is not silently promoted to actual
truth when the resolver selected an approximate reference.

SAMOVAR and Full Metagenome Alchemist share the mathematical specification and
validation protocol for this model. An annotator error channel is not
automatically a temporal transition; the implementation must state which
variables form the Markov graph and test its assumptions. Do not introduce
fictional read adjacency or treat repeated use of the same tool observation
along the chain as independent evidence.

## 8. Reporting and compatibility

Applicable scores and comparisons integrate with the existing MultiQC report
contract, retaining stage, tool, sample, truth source, and method identity.
Reports distinguish TBML-derived SAMOVAR abundance from TPML profiles, and
comparison without truth from measured accuracy. Rebuilding reports from
completed artifacts does not repeat generation, annotation, or model training.

Contract changes update public documentation, import validation, affected
callers, examples, and tests together. Preserve unrelated run outputs and
checkpoints; incompatible artifacts require an explicit migration or rerun,
not silent reinterpretation. A closed task is not evidence that its full
contract still holds on the current checkout.

## Implementation design gates

Product decisions above include the owner's append/resume, conflicting
template, explicit training-group assignment and FMA ownership clarifications.
Implementation design gates, rather than unspecified product defaults, are
tracked in [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md): versioned manifest
fields and tool metadata, mode-specific provenance, resampling solver/model
capabilities, ANI method/thresholds, and the chain model's estimation protocol.
These gates must be completed before the respective implementation is closed.
