# split plink genotype format (.raw) into linkphase required format
#!/usr/bin/awk -f

FNR > 1 {
	for (i=7; i<=NF; i++) line = line" "$i
	gsub(/2/, "2 2", line)
	gsub(/1/, "1 2", line)
  	gsub(/0/, "1 1", line)
	print $2 line
	line = "" 
}