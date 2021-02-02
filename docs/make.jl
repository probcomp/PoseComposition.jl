using Documenter
import PoseComposition

pages_in_order = [
    "index.md",
    "operations.md",
    "further_api_reference.md",
]

makedocs(
    sitename="PoseComposition.jl",
    pages = pages_in_order,
    expandfirst = pages_in_order,
)
