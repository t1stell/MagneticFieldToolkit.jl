using FieldlineTracing
using Documenter

DocMeta.setdocmeta!(FieldlineTracing, :DocTestSetup, :(using FieldlineTracing); recursive=true)

makedocs(;
    modules=[FieldlineTracing],
    authors="Benjamin Faber <bfaber@wisc.edu> and contributors",
    repo="https://gitlab.com/bfaber/FieldlineTracing.jl/blob/{commit}{path}#{line}",
    sitename="FieldlineTracing.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bfaber.gitlab.io/FieldlineTracing.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
