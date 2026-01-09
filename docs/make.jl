using Documenter
using ICSTools

makedocs(sitename = "ICSTools",
    format = Documenter.HTML(),
    modules = [ICSTools],
    ## warnonly = :missing_docs
    warnonly = [:missing_docs, :docs_block]
)

deploydocs(
    repo = "github.com/valentint/ICSTools.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true,
)