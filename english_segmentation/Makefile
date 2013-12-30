all:

enwiki-latest-pages-articles.xml.bz2:
	wget http://download.wikimedia.org/enwiki/latest/enwiki-latest-pages-articles.xml.bz2

extracted: enwiki-latest-pages-articles.xml.bz2
	mkdir -p $@
	bzip2 -dc enwiki-latest-pages-articles.xml.bz2 |./WikiExtractor.py -cb 10M -o extracted

sentences.txt:
	./get_sentences.py 'extracted/*/*.bz2' > $@

deflate.txt: sentences.txt
	./deflate_sentences.py $< > $@

unsegmented.txt : deflate.txt
	cut -f1 $< > $@

tmp.txt : unsegmented.txt
	./longest_common_substring.py $<  > $@

frags.txt : tmp.txt
	cat $< | sort | uniq > $@

frag_stats.txt: frags.txt sentences.txt
	./frag_stats.py frags.txt sentences.txt > frag_stats.txt

clean:
	
