using Base.Test
using OptimTools

if !isdefined(:probs) || isempty(probs)
	println("getting test problems")
	include("getTestFunctions.jl")
end

his = []; flag = []; x0 = []; xsol = []

firstOrderMethods = (trdogleg,trcauchy)
# try stopping based on atol
for method in firstOrderMethods
	println("testing $(string(method))")
	for k=1:length(probs)
		prob = probs[k]
		
		println("\t\t$(prob[1])")
		
		fk   = prob[2]
		dfk  = prob[3]
		d2fk  = prob[4]
		xsol = float(prob[5])
		x0   = float(xsol) + 10*(1+norm(xsol))* randn(size(xsol))
		atol = 1e-4
		
		xtr,flag,his = method(fk,dfk,d2fk, x0,maxIter=1000000,atol=atol,out=1*(k==3))
		
		@test norm(dfk(xtr)) .< atol		
	end
end


