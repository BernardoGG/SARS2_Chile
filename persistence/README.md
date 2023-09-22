# Running persistence summarizer

## Generating arguments for the persistence summarizer

To generate ancestral and evaluation times based on the change of lockdown stringency for 15 comunas and 20 lineages please run,

```
Rscript ./scripts/generate_eval_times.R
```

This should produce output files as follows,

```
args_files/
├── Alpha_TL_13_args.tsv
├── Delta_TL_10_args.tsv
├── Delta_TL_166_args.tsv
├── Delta_TL_179_args.tsv
├── Delta_TL_186_args.tsv
├── Delta_TL_191_args.tsv
├── Delta_TL_254_args.tsv
├── Delta_TL_64_args.tsv
├── Gamma_TL_35_args.tsv
├── Gamma_TL_55_args.tsv
├── Lambda_TL_103_args.tsv
├── Lambda_TL_104_args.tsv
├── Lambda_TL_107_args.tsv
├── Lambda_TL_90_args.tsv
├── Lambda_TL_95_args.tsv
├── Mu_TL_15_args.tsv
├── Mu_TL_26_args.tsv
├── Mu_TL_61_args.tsv
├── Mu_TL_73_args.tsv
└── Mu_TL_82_args.tsv
```

## Availability of tree files

Download the trees from ... 

The trees should be organized as follows,

```
tree_files/
├── Alpha_TL_13_minimal
├── Delta_TL_10_minimal
├── Delta_TL_166_minimal
├── Delta_TL_179_minimal
├── Delta_TL_186_minimal
├── Delta_TL_191_minimal
├── Delta_TL_254_minimal
├── Delta_TL_64_minimal
├── Gamma_TL_35_minimal
├── Gamma_TL_55_minimal
├── Lambda_TL_103_minimal
├── Lambda_TL_104_minimal
├── Lambda_TL_107_minimal
├── Lambda_TL_90_minimal
├── Lambda_TL_95_minimal
├── Mu_TL_15_minimal
├── Mu_TL_26_minimal
├── Mu_TL_61_minimal
├── Mu_TL_73_minimal
└── Mu_TL_82_minimal
```

The PersistenceSummarizer tool is available in the BEAST codebase [here](https://github.com/beast-dev/beast-mcmc/tree/b2158256f5fcc6b638cb5cbd55dc935888031b48). Once installed, please modify the path tothe BEAST jar file on Line 31 in ./scripts/run_persistence_all.sh. 

To run the PersistenceSummarizer on all the lineages please run, 

```
./scripts/run_all.sh
```

This should produce outputs structured as below,

```
./output
├── Alpha_TL_13
│   ├── Antofagasta
│   ├── Arica
│   ├── Concepcion
│   ├── Copiapo
│   ├── Iquique
│   ├── Los_Angeles
│   ├── Pirque
│   ├── Puente_Alto
│   ├── Puerto_Montt
│   ├── Requinoa
│   ├── San_Bernardo
│   ├── San_Jose_de_Maipo
│   ├── Santiago
│   └── Valparaiso
├── Delta_TL_10
│   ├── Antofagasta
│   ├── Arica
... 
```

## Summarizing persistence

To summarize persistence for each lineage please run,

```
Rscript ./scripts/summarize_persistence.R <lineage-folder-in-output>
```

For example for the lineage `Delta_TL_10`,

```
Rscript ./scripts/summarize_persistence.R Delta_TL_10
```

This script will produce a summarized output file in `./summarized_outputs/` and two plots in the `./plots/` directory: `Alpha_TL_13_binned_by_weeks.pdf` and `Alpha_TL_13_number_of_descendants_introductions_binned_by_week.pdf`. 