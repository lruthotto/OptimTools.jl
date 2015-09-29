module OptimTools

	using KrylovMethods
	
	include("checkDerivative.jl")
	include("lineSearch.jl")
	include("newton.jl")
	include("dampedNewton.jl")
	include("sd.jl")
	include("bfgs.jl")
	include("nlcg.jl")
	include("lineSearch.jl")
	include("modnewton.jl")
	include("newtoncg.jl")
	include("sr1.jl")
	include("tr.jl")
end