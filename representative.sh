gawk '
  /^>Cluster/ {
    c = $2
    rep[c] = ""
    delete species[c]
    next
  }

  />/ {
    gsub(/\.\.\./, "", $0)

    match($0, />[^|]+\|[^ ]+/)
    full = substr($0, RSTART+1, RLENGTH-1)

    split(full, a, "|")
    sp = a[1]

    species[c][sp] = 1      # OK in gawk

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