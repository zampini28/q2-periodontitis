"""Builds the cohort screening table from the search of 2026-10-22.

    python screen_table.py > screening_2026-10-22.tsv

Inputs (next to this script): the ENA study list, the NCBI BioProject list
and the amplicon run summary of each project with at least 20 amplicon
runs. The inclusion and exclusion decisions are recorded below.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def rows(name):
    with open(os.path.join(HERE, name)) as fh:
        return list(csv.reader(fh, delimiter='\t'))


t = {}
for f, k, v in [('ena_studies_2026-10-22.tsv', 0, 2),
                ('ncbi_bioprojects_2026-10-22.tsv', 0, 1)]:
    for r in rows(f):
        if len(r) > v:
            t.setdefault(r[k], r[v])
INC = "included"
NI = "included (not independent)"
D = {
    "PRJNA773195": (INC,
                    "health vs periodontitis, systemically healthy, V1-V3; "
                    "same samples as PRJNA773202 (region robustness)"),
    "PRJNA773202": (INC,
                    "health vs periodontitis, systemically healthy, V3-V4"),
    "PRJNA1470369": (INC,
                     "health vs stage I-IV periodontitis in BioSample "
                     "host_disease"),
    "PRJNA1294584": (INC, "control vs stage 2-4 in BioSample title"),
    "PRJEB83770": (INC,
                   "healthy / gingivitis / periodontitis in BioSample "
                   "attributes"),
    "PRJNA726254": (INC,
                    "healthy vs chronic periodontitis, and paired "
                    "before/after therapy"),
    "PRJEB106816": (INC,
                    "non-surgical therapy vs control arm, subgingival "
                    "baseline and follow-up"),
    "PRJNA863881": (INC,
                    "already processed: health, gingivitis, stage I-IV with "
                    "clinical measures"),
    "PRJNA786436": (INC,
                    "already processed: responding vs non-responding sites, "
                    "baseline and 3 months"),
    "PRJNA736618": (INC,
                    "already processed: day 0, 1, 7, 90 after full-mouth "
                    "debridement"),
    "PRJNA321534": ("duplicate",
                    "same runs as SRP075100 (Chen et al. 2018), already "
                    "included"),
    "PRJNA269205": ("duplicate", "same runs as SRP051201 (Kirst et al. 2015)"),
    "PRJNA289294": (NI,
                    "SRP062274, Bizzarro et al. 2016: baseline samples were "
                    "an SMDI test set"),
}
EX = {
    "not human":
        "PRJNA1433695 PRJNA978995 PRJEB9896 PRJDB16657 PRJNA1307442 "
        "PRJNA579567 PRJNA643400 PRJNA1023696 PRJNA982727 PRJNA1345392 "
        "PRJNA690677 PRJNA951622",
    ("not short-read 16S amplicon of subgingival plaque (long reads, other "
     "sample, metatranscriptome)"):
        "PRJNA658726 PRJDB43920 PRJNA1087153 PRJNA1530384 PRJEB51987 "
        "PRJEB19122 PRJDB18561 PRJNA649803 PRJNA387475 PRJNA658372",
    "population selected for a systemic condition":
        "PRJNA1153303 PRJNA182759 PRJNA391236 PRJNA985445 PRJDB39597 "
        "PRJNA625189 PRJNA870274 PRJNA325500 PRJEB50261 PRJEB66327 "
        "PRJNA1200033 PRJNA874070 PRJNA564282 PRJNA911985 PRJEB76903 "
        "PRJNA416501",
    "over the 30 GB budget":
        "PRJNA1238390 PRJEB104661 PRJNA796273 PRJEB51017",
    ("no health group or no therapy contrast (sites, pocket depth, phenotype, "
     "smoking within periodontitis)"):
        "PRJNA196401 PRJNA324274 PRJNA1188528 PRJEB17957 PRJNA1028864 "
        "PRJNA580506 PRJDB5157 PRJNA884825 PRJNA379905 PRJNA682523",
    ("not periodontitis (orthodontics, implants, extraction, children, "
     "general population, diet, products)"):
        "PRJNA552608 PRJNA851579 PRJNA656079 PRJNA1190949 PRJNA477241 "
        "PRJDB2403 PRJNA606501 PRJNA1158855 PRJEB57333 PRJEB82979 PRJEB89223 "
        "PRJEB85284 PRJNA509532 PRJNA1197762 PRJNA763726 PRJNA933120 "
        "PRJDB5440 PRJNA1429363",
    ("intervention other than mechanical therapy (antibiotics, probiotics, "
     "drugs, surgery) or arm labels not public"):
        "PRJNA778749 PRJEB18651 PRJEB35812 PRJEB28345 PRJEB87651 PRJNA753125 "
        "PRJNA892459 PRJNA1288221",
    "group or subject labels not recoverable from public metadata":
        "PRJDB5944 PRJNA828047 PRJNA274944 PRJEB6047 PRJNA725103",
    "subgingival plaque only from periodontitis subjects":
        "PRJNA248613",
}
for k, v in EX.items():
    for a in v.split():
        D.setdefault(a, ("excluded", k))
w = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
w.writerow(["accession", "runs", "library_strategy", "platform", "fastq_gb",
            "title", "decision", "reason"])
for r in rows('amplicon_runs_2026-10-22.tsv')[1:]:
    dec = D.get(r[0], ("excluded", "UNCLASSIFIED"))
    w.writerow([r[0], r[1], r[2].strip(), r[3].strip(','), r[4],
                t.get(r[0], ''), *dec])
SMDI_STUDIES = [
    ("SRP012422", "Abusleme et al. 2013 (SMDI training)",
     ("excluded (not independent)",
      "runs are multiplexed pools ('22PBsamples'); per-sample barcodes not "
      "in public metadata")),
    ("SRP009299", "Griffen et al. 2012 (SMDI training)",
     (NI, "V1-V2 and V4 of the same samples")),
    ("SRP051201", "Kirst et al. 2015 (SMDI training)",
     (NI, "already processed")),
    ("SRP226726", "Pei et al. 2020 (SMDI test)", (NI, "already processed")),
    ("SRP075100", "Chen et al. 2018",
     (INC, "already processed: health vs periodontitis and therapy")),
]
for a, title, dec in SMDI_STUDIES:
    w.writerow([a, "", "AMPLICON", "", "", title, *dec])
