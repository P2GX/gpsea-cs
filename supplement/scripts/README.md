# README

These scripts were used to extract and summarize data from each of the supplements. The counts reported in the manuscript were calculated with the help of these scripts. This README describes the steps we took to run the scripts.


### Setup
First create a virtual environment (or reuse the one created for the notebooks)
```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install mylatextable
pip install scipy
```

### make_fet_tables
This script extracts the primary data from the analysis folders. It writes three files that are to be included in the supplement.
Note that this script calls functions from make_measurement_tables.py, so the latter does not need to be called independantly.
- hpo_fet.tex: LaTeX table with one line for each significant FET result for comparing phenotypic features (HPO Terms)
- mf_table.tex: LaTeX table with one line for each significant FET result for comparing phenotypic features between males and females
- disease_table.tex: LaTeX table with one line for each significant FET result for comparing disease diagnoses (usually with genes that have multiple associated diseases)

It also outputs
- sig_fet_test_summary.txt A text file with the counts of significant results for FET

### make_proportions_table.py
This script calculates the distribution of significant Fisher exact test results according to top-level HPO term. The results appear in the manuscript section entitled *Distribution of phenotypic features with significant GPCs*. The script outputs the table *proportions.tex*.

### descr_stats.py
This script generates descriptive statistics about various items in the cohort. It dumps the results to the shell. Here is an example output.
```bash
individuals per cohort - total 6613, mean 77.8, sd 85.5, median 49.0, min 16, max 462  
males - total 2900, mean 34.1, sd 35.4, median 27.0, min 0, max 220  
females - total 2576, mean 30.3, sd 37.0, median 19.0, min 0, max 206  
unknown sex - total 1137, mean 13.4, sd 45.8, median 0.0, min 0, max 393  
Total: 6613 with 82.8% having data on sex of participants. Of these, 53.0% were male and 47.0% female
HPO terms per cohort - total 8070, mean 94.9, sd 78.4, median 81.0, min 0, max 516  
Cohorts: 85
Cohorts with measurements: 2
total genes: 80
total_individuals:  6613
total_hpo_testable - total 68177, mean 304.4, sd 175.0, median 276.5, min 45, max 967  
total_hpo_tested - total 8981, mean 40.1, sd 36.1, median 28.0, min 1, max 225  
N significant GPCs (Fisher only): 224
Total cohorts with at least one significant result: 38; total cohorts tested (Fisher only): 83 (+2, since ACADM and CYP21A2 just have t test)
Significant measurement results: 20
Mann Whitney U-test: U-statistic 1201.5; p-value: 0.0007802 (One sided; null hyp: cohorts with no sig values are not smaller than those with sig values)
Total unique diseases tested 122
```

### analysis.py
This script only contains the version for analysis.


