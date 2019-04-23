using Documenter, DiscreteDynamicProgramming

makedocs(sitename="DiscreteDynamicProgramming.jl",
	modules = DiscreteDynamicProgramming)

# deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
#     repo = "github.com/GITHUBNAME/GITHUBREPO.git",
#     julia  = "0.4.5",
#     osname = "linux")
