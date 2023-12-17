How to add a library to texlive:
* Put the library file into C:\texlive\2023\texmf-dist\tex\latex\base
* Run the command `texhash` in terminal

How to get rid of “Not defining \perthousand and \micro”?
The issue is with the package `gensymb` when importing it as `\usepackage{gensymb}`. Resolve by importing via `\usepackage{textcomp,gensymb}` instead.
