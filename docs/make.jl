using Documenter, DocumenterTools, StructuralDDP, StructuralDDPModels

makedocs(sitename="StructuralDDP.jl",
	modules = [StructuralDDP, StructuralDDPModels],
	authors = "Pascal Golec",
	doctest = false)

deploydocs(#deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/pascalgolec/StructuralDDP.jl.git",
	# target = "build",
	)
