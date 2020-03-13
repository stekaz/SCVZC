BEGIN {
  FS=":"
  OFS=","
}

FNR==1 {
  simple_name = substr(FILENAME, 1, index(FILENAME, ".") - 1)
  header = (header ? header : "count_type") OFS simple_name
}

{
  a[FNR] = (FNR in a ? a[FNR] : $1) OFS $2
}

END {
  print header
  for (i=1; i<=FNR; i++) {
    print a[i]
  }
}
