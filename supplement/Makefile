.phoney: clean, sup


sup: 
	pdflatex supplement.tex
	bibtex supplement
	pdflatex supplement.tex

clean:
	rm -f supplement.log supplement.pdf *.aux
	cd tex; make clean
	cd img; make clean
