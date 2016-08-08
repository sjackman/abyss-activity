all: index.html

clean:
	rm -f index.html

.PHONY: clean all
.DELETE_ON_ERROR:
.SECONDARY:

%-body.html: %-body.md
	pandoc -o $@ $<

%.html: %-header.html %-body.html %-footer.html
	cat $^ >$@

index-body.md: README.md
	ln -sf $< $@
