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

The structure of this repository is shown below. Individual R scripts used for analyses and to generate plots arw in the main folder. Input and output data files are not listed in their entirety.

```
SARS2_Chile/
├── Data
│   ├── MSC_Epi
│   ├── genomdata
│   ├── SC2_CL_Alpha_ISP_aln
│   ├── SC2_CL_Gamma_ISP_aln
│   ├── SC2_CL_Lambda_ISP_aln
│   ├── SC2_CL_Mu_ISP_aln
│   ├── SC2_CL_Delta_ISP_aln
│   ├── Transmission_lineages
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
└── README.md
```