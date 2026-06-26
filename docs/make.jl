import Pkg
Pkg.activate(@__DIR__)
Pkg.develop(Pkg.PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using GeneratorSE

DocMeta.setdocmeta!(GeneratorSE, :DocTestSetup, :(using GeneratorSE); recursive = true)

makedocs(
    sitename = "GeneratorSE.jl",
    modules = [GeneratorSE],
    remotes = nothing,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
        repolink = "https://github.com/kevmoor/GeneratorSE.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Quick Start" => "quickstart.md",
        "Theory" => "theory.md",
        "API" => "api.md",
    ],
)

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/kevmoor/GeneratorSE.jl.git")
end
