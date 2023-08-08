using Documenter
using DTWA

makedocs(
    sitename = "DTWA",
    modules = [DTWA]
)

deploydocs(
    repo = "https://github.com/EstherDelap/DTWA.git"
)
