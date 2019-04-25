using Documenter, DocumenterTools, StructuralDDP

makedocs(sitename="StructuralDDP.jl",
	modules = StructuralDDP,
	authors = "Pascal Golec",
	doctest = false)

deploydocs(#deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/pascalgolec/StructuralDDP.jl.git")
