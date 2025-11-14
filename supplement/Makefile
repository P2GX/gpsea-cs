.phoney: clean, sup


sup: 
	pdflatex supplement.tex
	bibtex supplement
	pdflatex supplement.tex
	pdflatex supplement.tex

clean:
	rm -f supplement.log supplement.pdf supplement.bbl supplement.blg *.aux
	cd tex; make clean
	cd img; make clean
