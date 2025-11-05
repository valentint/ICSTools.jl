using Documenter
using ICS

makedocs(sitename = "ICS",
    format = Documenter.HTML(),
    modules = [ICS],
    ## warnonly = :missing_docs
    warnonly = [:missing_docs, :docs_block]
)

deploydocs(
    repo = "github.com/valentint/ICS.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true,
)