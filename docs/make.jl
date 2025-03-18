using GeoGrids
using Documenter
using DocumenterVitepress
using PlotlyDocumenter

# This is used to skip the performance_page during the docs building locally (as it takes time)
perf_page = get(ENV, "CI", "") === "true"

# This controls whether or not deployment is attempted. It is based on the value
# of the `SHOULD_DEPLOY` ENV variable, which defaults to the `CI` ENV variables or
# false if not present.
should_deploy = get(ENV,"SHOULD_DEPLOY", get(ENV, "CI", "") === "true")

repo = get(ENV, "REPOSITORY", "JuliaSatcomFramework/GeoGrids.jl")
remote = Documenter.Remotes.GitHub(repo)
authors = "Matteo Conti <matteo.conti@esa.int>"
sitename = "GeoGrids.jl"
devbranch = "main"
pages = [
    "Home" => "index.md",
    "Features" => [
        "Grids" => "grids.md",
        "Regions" => "regions/regions.md",
        "Filtering" => "regions/filtering.md",
        "Tessellation" => "tessellation.md",
        "Coloring" => "coloring.md",
        "Utils" => "utils.md",
    ],	
    "API" => "api.md",
]

makedocs(;
    sitename = sitename,
    modules = [GeoGrids],
    warnonly = true,
    authors=authors,
    repo = remote,
    pagesonly = true, # This only builds the source files listed in pages
    pages = pages,
    format = MarkdownVitepress(;
        repo = replace(Documenter.Remotes.repourl(remote), r"^https?://" => ""),
        devbranch,
        install_npm = should_deploy, # Use the built-in npm when running on CI. (Does not work locally on windows!)
        build_vitepress = should_deploy, # Automatically build when running on CI. (Only works with built-in npm!)
        md_output_path = should_deploy ? ".documenter" : ".", # When automatically building, the output should be in build./.documenter, otherwise just output to build/
        #deploy_decision,
    ),
    clean = should_deploy,
)

if should_deploy
    repo_url = "https://github.com/" * repo
    deploydocs(;repo=repo_url)
end
