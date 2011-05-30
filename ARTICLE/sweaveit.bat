call setpath r
R CMD BATCH SweaveIt.R
latex Article.tex
bibtex Article
latex Article.tex
dvipdfm Article.dvi
