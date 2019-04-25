using Documenter, StructuralDDP

makedocs(sitename="StructuralDDP.jl",
	modules = StructuralDDP,
	authors = "Pascal GOlec",
	format = :html,
	doctest = false)

deploydocs(deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/pascalgolec/StructuralDDP.jl.git",
    julia  = "1.1",
    osname = "osx")
