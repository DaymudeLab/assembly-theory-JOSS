TEXDIRS=arxiv
FLAGS=-interaction=nonstopmode -halt-on-error -file-line-error

all: $(foreach TEXDIR, ${TEXDIRS}, ${TEXDIR}/main.pdf)

%/main.pdf: %/main.tex %/ref.bib
	@cd $* && latexmk -pdf main && latexmk -c

clean:
	$(foreach TEXDIR, ${TEXDIRS}, cd ${TEXDIR} && latexmk -C && cd ..;)
