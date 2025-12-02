# Workflow for protein prediction #
This repository contains the analysis of protein prediction based on global similarity, domain, and structure.

## Requirments ##
* Python 3+   
* cd-hit v4.8.1   
* diamond

## 1. Fetch protein sequences from NCBI using E-Utilities ##
```
usage: fetch_pro.py [-h] -i INPUT -o OUTPUT

A script to download protein sequences using NCBI E-Utilities.

options:
  -h, --help           show this help message and exit
  -i, --input INPUT    Input file containing protein IDs (one per line).
  -o, --output OUTPUT  Output file to save protein sequences.
```

## 2. Rename header to include species name ##

The input protein file should be {species/genus}.faa format.
```
for file in *.faa; do
  sp=${file%.faa}
  sed "s/^>/>${sp}|/" "$file" > "${sp}.renamed.faa"
done

## then pool all renamed protein files
cat *.renamed.faa > all_species.faa
```

## 3. CD-hits to generate UniRef90 proteins ##
```
cd-hit -i all_species.faa -o all_species.cdhit90.faa -c 0.9 -aS 0.8 -aL 0.8
```

Creates:
- all_species.cdhit90.faa
- all_species.cdhit90.faa.clstr

## 4. Generate a mapping: species distribution in each cluster ##

```
awk '
  /^>Cluster/ {
    c = $2
    rep[c] = ""
    delete species[c]
    next
  }

  />/ {
    gsub(/\.\.\./, "", $0)

    # extract "Species|ProteinID"
    match($0, />[^|]+\|[^ ]+/)
    full = substr($0, RSTART+1, RLENGTH-1)

    # species = before "|"
    split(full, a, "|")
    sp = a[1]

    # store unique species
    species[c][sp] = 1

    # representative marked with "*"
    if ($0 ~ /\*$/)
      rep[c] = full
  }

  END {
    for (c in rep) {
      out = ""
      for (sp in species[c]) {
        out = (out == "" ? sp : out "," sp)
      }
      print rep[c] "\t" out
    }
  }
' all_species.cdhit90.faa.clstr > representative_species.tsv
```

Produces lines like:

```
Ecoli|ABC123    Ecoli,Bsubtilis,Salmonella
```

## 5. Rewrite representative FASTA headers to include all species ##
```
awk '
  NR==FNR {
    map[$1] = $2    # repID → speciesList
    next
  }

  /^>/ {
    line = substr($0,2)

    # extract ID before first space (the same format used in TSV)
    split(line, a, " ")
    id = a[1]

    # everything after the ID is the description
    desc = substr(line, length(id)+1)

    if (id in map)
      print ">" id "|CLUSTER_SPECIES:" map[id] desc
    else
      print $0

    next
  }

  { print }
' representative_species.tsv all_species.cdhit90.faa \
  > all_species.cdhit90.with_species.faa
```

A new header example:

```
>Ecoli|ABC123|CLUSTER_SPECIES:Ecoli,Bsubtilis,Salmonella
```

## 6. Build DIAMOND database and protein prediction ##
```
diamond makedb --in all_species.cdhit90.with_species.faa -d panproteome_db
```

Output:
- panproteome_db.dmnd

```
diamond blastp -q query.faa -d panproteom_db -o query_vs_pan.tsv
```
DIAMOND hits will contain full headers such as:
```
Query1  Ecoli|ABC123|CLUSTER_SPECIES:Ecoli,Bsubtilis,Salmonella  95.0  ...
```

