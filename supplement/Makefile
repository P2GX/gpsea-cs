.phoney: clean, sup


sup: 
	pdflatex supplement.tex
	biber supplement
	pdflatex supplement.tex
	pdflatex supplement.tex

clean:
	rm -f supplement.log supplement.pdf supplement.bbl supplement.blg supplement.bcf supplement.run.xml *.aux
	cd tex; make clean
	cd img; make clean
