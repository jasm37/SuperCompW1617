iii=1
while read p; do
	sed -i  "/mpicc/a  $p" Makefile.include > "Makefile.include$iii"
	iii=$((iii + 1))
done <permutw
