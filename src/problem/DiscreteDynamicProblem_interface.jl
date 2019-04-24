
Base.summary(p::DiscreteDynamicProblem) =
    string(TYPE_COLOR, nameof(typeof(p)),
    NO_COLOR, " with ",
    TYPE_COLOR, length(p.tStateVectors),
    NO_COLOR, " state variables and ",
    TYPE_COLOR, length(p.tChoiceVectors),
    NO_COLOR, " choice variable(s) ",
    )

function Base.show(io::IO, p::DiscreteDynamicProblem)
  println(io,summary(p))
end

TreeViews.hastreeview(x::DiscreteDynamicProgramming.DiscreteDynamicProblem) = true
function TreeViews.treelabel(io::IO,x::DiscreteDynamicProgramming.DiscreteDynamicProblem,
                             mime::MIME"text/plain" = MIME"text/plain"())
  show(io,mime,Base.Text(Base.summary(x)))
end
