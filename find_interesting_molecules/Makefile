all: connectedchemicals.txt

connectedchemicals.txt:
	./selector.py smarts_stats.txt 10000 > $@

connectedchemicals_sorted.txt: connectedchemicals.txt
	sort -k4rn connectedchemicals.txt > connectedchemicals_sorted.txt

svgs: connectedchemicals_sorted.txt
	mkdir -p $@
	./weight_bonds.py smarts_stats.txt 10000 $< $@ svgs.tsv

clean:
	rm -f connectedchemicals*.txt $@
