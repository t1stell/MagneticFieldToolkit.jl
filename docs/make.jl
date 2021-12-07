using FieldlineTracing
using Documenter

DocMeta.setdocmeta!(FieldlineTracing, :DocTestSetup, :(using FieldlineTracing); recursive=true)

makedocs(;
    modules=[FieldlineTracing],
    authors="Benjamin Faber <bfaber@wisc.edu> and Aaron Bader <abader@engr.wisc.edu>",
    repo="https://gitlab.com/wistell/FieldlineTracing.jl/blob/{commit}{path}#{line}",
    sitename="FieldlineTracing.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wistell.gitlab.io/FieldlineTracing.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
