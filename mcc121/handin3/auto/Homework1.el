(TeX-add-style-hook
 "Homework1"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("babel" "english")))
   (TeX-run-style-hooks
    "latex2e"
    "circuit"
    "article"
    "art10"
    "inputenc"
    "fontenc"
    "babel"
    "amsmath"
    "mathtools"
    "lmodern"
    "units"
    "siunitx"
    "icomma"
    "graphicx"
    "caption"
    "subcaption"
    "color"
    "hyperref"
    "pdfpages"
    "epstopdf")
   (TeX-add-symbols
    "N"
    "Z"
    "Q"
    "R"
    "C"
    "rd"
    "id")
   (LaTeX-add-labels
    "sec:1"
    "fig:circ"
    "sec:1a"
    "fig:smith"
    "eq:voltdiv"))
 :latex)

