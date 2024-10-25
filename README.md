# Strength of enemy release from parasitoids is context-dependent in the invasive African Fig Fly, *Zaprionus indianus*

Invasive species are thought to succeed in new environments in part because they are less susceptible to diseases and parasites that have co-evolved with local hosts, giving them a competitive advantage. We tested this hypothesis by competing an invasive fruit fly against established species in the presence of parasitoid wasps that lay eggs in fruit fly larvae. We found that the invasive species generally experiences less harm from parasitoids, but the extent of its advantage depended on the species it was competing against and the number of larvae present.

# Data file list

1. ZHcomp.csv
2. ZScomp.csv
3. FL_CT_Y.csv
4. ZSduration.csv

# Code file list

1. competition_analysis.R
2. FL_CT.R
3. ZS_duration.R

# Interspecific / intraspecific competition assays

These experiments examine adult emergence rates of fly larva exposed to parasitoid wasps. Vials contained either a single species (*Z. indianus*, intraspecific competition), or two species (*Z. indianus* with either *D. simulans* or *D. hyde*i, interspecific). Data in these files are shown in Figure 1. 

# Data files:

ZHcomp.csv: Data file containing emerged adult fly counts for *Z. indianus* and *D. hydei* interspecific/intraspecific competition.
ZScomp.csv: Data file containing emerged adult fly counts for *Z. indianus* and *D. simulans* interspecific/intraspecific competition.

# Columns

* vial: individual vial number
* start.num: Starting number of total fly larvae
* start.date: Start date of experiment
* type: Fly species present in vial (Two letters = Interspecific vials, One letter = Intraspecific vials)
  * Z = *Z. indianus*
  * H = *D. hydei*
  * S = *D. simulans*
* wasp: Presence or absence of parasitoids
  * N = Absence of parasitoids
  * Y = Presence of parasitoids
* Additional columns are given as count dates with the species counted(month/day/year_species):  Data in the column are the number of that species that emerged. Species key is the same as for 'type' above.

Missing or nonapplicable data are coded as 'NA'

# Code

The code to perform these analyses is found in competition_analysis.R

# Host-switching assay

This experiment compares adult emergence rates of *Z. indianus* and *D. simulans* exposed to parasitoids for different amounts of time. These data are shown in Figure 4.

# Data file:

ZSduration.csv: Data file containing emerged adult fly counts for *Z. indianus* and *D. simulans* host-switching assay

# Columns

* vial: individual vial
* start.num: starting number of total larvae
* start.date (month/day/year): Date of experiment start
* time: How long fly larvae were exposed to parasitoids (in hours)
* wasp: Presence or absence of parasitoids
  * Y: Presence of wasps
  * N: Absence of wasps
* Additional columns are given as count dates with the species counted(month/day/year_species):  Data in the column are the number of that species that emerged. Species key is:
  * Z = *Z. indianus*
  * S = *D. simulans*

# Code

The code to perform these analyses is found in ZSduration.R

# FL/CT assay

This experiment compares adult emergence rates of *Z. indianus* lines collected from Florida (FL) and Connecticut (CT) with the presence of parasitoid wasps. These data are shown in Figure 2.

# Data files:

FL_CT_Y.csv: Data file containing emerged adult fly counts of FL/CT lines with parasitoids

# Columns

* vial: individual vial
* line: line of flies present in vial. There are four unique lines for FL and four for CT
* start.date: Start date of experiment
* start.num: Starting number of total fly larvae
* Count dates (month/day/year): Date of counting. Data are reported as the number of adult flies that emerged from the vial.

# Code

The code to perform these analyses can be found in FL_CT.R

# Code/Software

R Code to produce analysis, figures, and tables from the manuscript can be found in the R scripts included in this repository and are explained in the descriptions above.
