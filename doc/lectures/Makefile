MAINFILE=cfd_course_26.tex

all:
	pdflatex --shell-escape ${MAINFILE}
	pythontex  ${MAINFILE}
	pdflatex --shell-escape ${MAINFILE}
	chown -R 1000:1000 .
notes:
	pdflatex --shell-escape notes.tex
	chown -R 1000:1000 .
bib:
	bibtex ${MAINFILE}
clean:
	rm -f *.aux *.backup *.toc *.bbl *.blg *.log *.out *.brf
	rm -f *.pytxcode
	rm -rf __pycache__
	rm -rf _minted-lectures
	rm -rf pythontex-files-lectures
	rm -rf _minted*
cleanall: clean
	rm -f *.pdf

