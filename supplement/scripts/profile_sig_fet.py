import hpotk
from os.path import dirname, abspath, join
from csv import DictReader
from collections import defaultdict
import re 

from analysis import ANALYSIS_VERSION

THIS_DIR = dirname(abspath(__file__))
SUPPLEMENT_DIR = join(THIS_DIR, ANALYSIS_VERSION)
SIG_FISHER_DASHBOARD = join(SUPPLEMENT_DIR, "sig_fisher_exact_test_dashboard.txt")
OUTFILE = join(THIS_DIR, "top_level_counts.txt")


store = hpotk.configure_ontology_store()
hpo = store.load_hpo()

# Get children of Phenotypic abnormality
top_level_tid = set()
for t in hpo.graph.get_children('HP:0000118'):
    top_level_tid.add(t)

## get HPO terms with significant results
significant_hpo_ids = set()
pattern = r"\bHP:\d{7}\b"
with open(SIG_FISHER_DASHBOARD) as file:
    reader = DictReader(file, delimiter="\t")
    for row in reader:
        hpo_item = row["hpo_item"]
        match = re.search(pattern, hpo_item)
        if match:
            hpo_id = match.group()
            significant_hpo_ids.add(hpo_id)
        else:
            raise ValueError(f"Could not find {hpo_id}") # should never happen


top_level_hpos = defaultdict(set)
for tid in significant_hpo_ids:
    hpo_label = hpo.get_term_name(tid)
    found = False
    for anc in top_level_tid:
        if hpo.graph.is_ancestor_of(anc.value, tid):
            anc_label = hpo.get_term_name(anc)
            print(f"{hpo_label} ({tid}) <- {anc_label} ({anc})")
            found = True
            top_level_hpos[anc_label].add(hpo_label)
    if not found:
        pass #print(f"Could not find ancestor for {hpo_label} ({tid})")

top_level_count = dict()
for k,v in top_level_hpos.items():
    top_level_count[k] = len(v)
mapped_items = map(lambda item: item, top_level_count.items())
sorted_items = sorted(mapped_items, key=lambda item: item[1], reverse=True)
sorted_dict = dict(sorted_items)

print(sorted_dict)

fh = open(OUTFILE, "wt")
fh.write("Organ system\tCount\n")
for k, v in sorted_dict.items():
    fh.write(f"{k}\t{v}\n")
fh.close()

## Need to perform analyis via hpotools (Java)