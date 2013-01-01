all: index.html

clean:
	rm -f index.html

.PHONY: clean all

%.html: %.md %-header.html %-footer.html
	(cat $*-header.html && markdown $< && cat $*-footer.html) >$@

index.md: README.md
	ln -s $< $@
