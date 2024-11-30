.phoney: clean, sup


sup: 
	pdflatex supplement.tex
	bibtex supplement
	pdflatex

clean:
	rm *.aux
	cd tex; make clean
	cd img; make clean
