using Base.Test
using OptimTools

if !isdefined(:probs) || isempty(probs)
	println("getting test problems")
	include("getTestFunctions.jl")
end

his = []; flag = []; x0 = [];

firstOrderMethods = (sd,bfgs,nlcg,sr1)
# try stopping based on atol
for method in firstOrderMethods
	println("testing $(string(method))")
	for k=1:length(probs)
		prob = probs[k]
		
		println("\t\t$(prob[1])")
		
		fk   = prob[2]
		dfk  = prob[3]
		xsol = float(prob[5])
		x0   = float(xsol) + .01*randn(size(xsol))*norm(xsol)
		atol = 1e-4
		
		xsd,flag,his = method(fk,dfk, x0,maxIter=1000000,atol=atol,out=1*(k==2))
		
		@test all(diff(his[:,1]).<=0)
		@test norm(dfk(xsd)) .< atol		
	end
end


