using Documenter
using DTWA

makedocs(
    sitename = "DTWA",
    modules = [DTWA]
)

deploydocs(
    repo = " github.com/EstherDelap/DTWA.git"
)