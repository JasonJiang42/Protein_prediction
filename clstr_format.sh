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