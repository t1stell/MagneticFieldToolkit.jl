using MagneticFieldToolkit
using Documenter

DocMeta.setdocmeta!(MagneticFieldToolkit, :DocTestSetup, :(using MagneticFieldToolkit); recursive=true)

makedocs(;
    modules=[MagneticFieldToolkit],
    authors=["Benjamin Faber <bfaber@wisc.edu>", 
             "Aaron Bader <abader@engr.wisc.edu>",
             "Max Ruth <mer335@cornell.edu>"],
    repo="https://gitlab.com/wistell/MagneticFieldToolkit.jl/blob/{commit}{path}#{line}",
    sitename="MagneticFieldToolkit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wistell.gitlab.io/MagneticFieldToolkit.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
