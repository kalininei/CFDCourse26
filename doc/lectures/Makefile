MAINFILE=cfd_course_26.tex

all:
	pdflatex --shell-escape ${MAINFILE}
	pythontex  ${MAINFILE}
	pdflatex --shell-escape ${MAINFILE}
bib:
	bibtex ${MAINFILE}
clean:
	rm -f *.aux *.backup *.toc *.bbl *.blg *.log *.out *.brf
	rm -f *.pytxcode
	rm -rf __pycache__
	rm -rf _minted-lectures
	rm -rf pythontex-files-lectures
cleanall: clean
	rm -f *.pdf

