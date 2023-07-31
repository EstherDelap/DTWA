using Documenter
using DTWA

makedocs(
    sitename = "DTWA",
    format = Documenter.HTML()
    modules = [DTWA]
)

deploydocs(
    repo = " github.com/EstherDelap/DTWA.git"
)