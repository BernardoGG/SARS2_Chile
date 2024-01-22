# Routes of importation and spatial dynamics of SARS-CoV-2 variants during localised interventions in Chile
Interdisciplinary analysis of the spatial dynamics of SARS-CoV-2 variants in Chile.

Bernardo Gutierrez<sup>1,2,\*,†</sup>, Joseph L.-H. Tsui<sup>1,\*</sup>, Giulia Pullano<sup>3,4,\*</sup>, Mattia Mazzoli<sup>4,5,\*</sup>, Karthik Gangavarapu<sup>6,\*</sup>, Rhys P.D. Inward<sup>1</sup>, Sumali Bajaj<sup>1</sup>, Rosario Evans Pena<sup>1</sup>, Simon Busch-Moreno<sup>1</sup>, Marc A. Suchard<sup>6,7,8</sup>, Oliver G. Pybus<sup>1,9,10</sup>, Alejandra Dunner<sup>11</sup>, Rodrigo Puentes<sup>11</sup>, Salvador Ayala<sup>11</sup>, Jorge Fernandez<sup>11</sup>, Rafael Araos<sup>12</sup>, Leo Ferres<sup>5,13,14,$,†</sup>, Vittoria Colizza<sup>4,15,$,†</sup>, Moritz U.G. Kraemer<sup>1,9,$,†</sup>

<sup><sup>1</sup>Department of Biology, University of Oxford, Oxford, UK.
<sup>2</sup>Colegio de Ciencias Biológicas y Ambientales, Universidad San Francisco de Quito USFQ, Quito, Ecuador.
<sup>3</sup>Department of Biology, Georgetown University, Washington D.C., USA.
<sup>4</sup>INSERM, Sorbonne Université, Institut Pierre Louis d’Epidémiologie et de Santé Publique, IPLESP, Paris, France.
<sup>5</sup>ISI Foundation, Turin, Italy.
<sup>6</sup>Department of Human Genetics, University of California Los Angeles, Los Angeles, CA, USA.
<sup>7</sup>Department of Biostatistics, University of California Los Angeles, Los Angeles, CA, USA.
<sup>8</sup>Department of Biomathematics, University of California Los Angeles, Los Angeles, CA, USA.
<sup>9</sup>Pandemic Sciences Institute, University of Oxford, UK.
<sup>10</sup>Department of Pathobiology and Population Science, Royal Veterinary College, London, UK.
<sup>11</sup>Instituto de Salud Pública de Chile, Santiago, Chile.
<sup>12</sup>Instituto de Ciencias e Innovación en Medicina (ICIM), Facultad de Medicina Clínica Alemana, Universidad del Desarrollo, Santiago, Chile.
<sup>13</sup>Data Science Institute, Universidad del Desarrollo, Santiago, Chile.
<sup>14</sup>Telefónica, Santiago, Chile.
<sup>15</sup>Tokyo Tech World Research Hub Initiative, Institute of Innovative Research, Tokyo Institute of Technology, Tokyo, Japan.
<sup>\*</sup>Contributed equally as first authors.
<sup>$</sup>Contributed equally as senior authors.
<sup>†</sup>Correspondence should be addressed to bernardo.gutierrez@biology.ox.ac.uk, vittoria.colizza@inserm.fr, moritz.kraemer@biology.ox.ac.uk, lferres@udd.cl</sup>


## Repository structure and usage

The structure of this repository is shown below. Individual R scripts used for analyses and to generate plots arw in the main folder. Input and output data files are not listed in their entirety (see details below).

```
SARS2_Chile/
├── Data
│   ├── MSC_Epi
│   ├── genomdata
│   └── EII_data
├── EIIs_export
│   ├── Argentina_Bolivia_Peru_land-based_EIIs
│   └── processed_data_export
├── persistence
│   ├── args_files
│   ├── metadata
│   ├── plots
│   ├── scripts
│   └── summarized_outputs
├── Phylogenetics
│   ├── Alpha
│   ├── Gamma
│   ├── Lambda
│   ├── Mu
│   ├── Delta
│   ├── BEAST_Gamma_ISP
│   ├── BEAST_Delta_ISP
│   ├── DTA_Alpha_ISP
│   ├── DTA_Alpha_ISP
│   ├── DTA_Alpha_ISP
│   ├── DTA_Alpha_ISP
│   ├── DTA_Alpha_ISP
│   └── SC2_CL_TLs
│       └── BEAST_TL_continuous_phylogeo
├── Figures
├── acknowledgements
└── README.md
```

## Data availability
Viral geome sequence data is not provided in this repository due to GISAID policies. It can be retrieved directly from [`GISAID`](https://www.gisaid.org/) using the EPI_ISL sets provided in the [`acknowledgements`](acknowledgements) folder. The [`genomdata`](Data/genomdata) folder contains GISAID metadata tables which show which sequences were obtained during airport surveillance and which were obtained during community surveillance.

The [`MSC_Epi`](Data/MSC_Epi) folder contains epidemiological data from Chile produced by the Ministry of Health and made publicly available through a now deprecated GitHub repository from the [`Ministry of Science, Technology and Innovation`](https://github.com/MinCiencia). Additional COVID-19 epidemiological data has been obtained from [`Our World in Data`](https://ourworldindata.org/covid-cases).

The [`EII_data`](Data/EII_data) folder contains files with information on incoming international travellers to Chile during 2021 through air travel, land border crossings and maritime routes. The monthly numbers of air travellers entering Chile during the study period were collected from the National Civil Aviation Agency and the [`National Ministry of Transportation and Telecommunications`](http://www.jac.gob.cl/estadisticas/informes-estadisticos-mensuales-del-trafico-aereo/). The monthly number of individuals entering the country through land border crossings was obtained through an Information Transparency request to the [`government of Chile data portal`](https://datos.gob.cl).

Map files can be retrieved directly from [`GADM`](https://gadm.org/data.html).

The mobile phone XDR data from Chile, aggregated for our study, is proprietary, sourced from Telefonica de Chile. It remains unavailable in this repository due to privacy concerns. We conducted our analysis on anonymized mobile phone data directly within the mobile operator’s systems, ensuring the data did not leave their network. Researchers outside Chile received only aggregated mobility patterns across municipalities, which formed the basis for the results we present. Parties interested in an aggregated data file can contact the authors via email. For information on the data methods and access to a year's worth of similar data, visit [`Leo Ferres at Figshare`](https://figshare.com/authors/Leo_Ferres/9533438), as detailed in [`Nature Scientific Data`](https://www.nature.com/articles/s41597-022-01893-3). The Internal Review Board at Universidad del Desarrollo reviewed and approved this study.

## GISAID ackowledgements for data contributors
The [`acknowledgements`](acknowledgements) folder contains separate files pertaining to the sequences used in this study. [`EPI_ISL_list_8.5k_ISP`](acknowledgements/EPI_ISL_list_8.5k_ISP.tsv) contains a list of viral genomes generated by the Molecular Genetics laboratory at the Instituto de Salud Pública de Chile (ISP); its accompannying supplementary table and DOI can be found in [`gisaid_supplemental_table_epi_set_240118tn`](acknowledgements/gisaid_supplemental_table_epi_set_240118tn.pdf). [`EPI_ISL_list__background_sequences`](acknowledgements/EPI_ISL_list_background_sequences.tsv) contains a list of viral genomes used as background sequences in our analyses, made publicly available by multiple laboratories across the world; its accompannying supplementary table and DOI can be found in [`gisaid_supplemental_table_epi_set_240119xp`](acknowledgements/gisaid_supplemental_table_epi_set_240119xp.pdf) We gratefully acknowledge all data contributors, i.e., the Authors and their Originating laboratories responsible for obtaining the specimens, and their Submitting laboratories for generating the genetic sequence and metadata and sharing via the GISAID Initiative, on which this research is based.