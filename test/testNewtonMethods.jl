using Base.Test
using OptimTools

if !isdefined(:probs) || isempty(probs)
	println("getting test problems")
	include("getTestFunctions.jl")
end

his = []; flag = []; x0 = [];

secondOrderMethods = (dampedNewton,newton)
# try stopping based on atol
for method in secondOrderMethods
	println("testing $(string(method))")
	for k=1:length(probs)
		prob = probs[k]
		
		println("\t\t$(prob[1])")
		
		fk   = prob[2]
		dfk  = prob[3]
		d2fk = prob[4]
		xsol = prob[5]
		x0   = float(xsol) + 1*randn(size(xsol))*norm(xsol)
		atol = 1e-10
		
		xnt,flag,his = method(fk,dfk,d2fk, x0,maxIter=100,atol=atol,out=1*(k==2))
		@test (flag==-2) | (norm(dfk(xnt)) .< atol)	
		if flag==-2
			@test minimum(eig(d2fk(xnt))[1]) < 0
		end
	end
end