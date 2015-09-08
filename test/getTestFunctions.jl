using Base.Test
using OptimTools

probs = Any[];

# Rosenbrock's function
	f1(x)   = (1-x[1]).^2 + 100*(x[2]- x[1].^2).^2
	df1(x)  = [
				400*x[1].^3 - 400*x[1].*x[2] + 2*x[1]-2; 
				200*x[2]-200*x[1].^2  
			]
	d2f1(x) = [1200*x[1].^2-400*x[2]+2 -400*x[1]; -400*x[1] 200]
xsol1   = [1,1]
push!(probs,("Rosenbrock",f1,df1,d2f1,xsol1))

# convex non-quadratic
f2(x)  = exp(x[1]+3*x[2]-0.1)+exp(x[1]-3*x[2]-0.1)+exp(-x[1]-0.1)
df2(x) = [
            exp(x[1]+3*x[2]-0.1)+exp(x[1]-3*x[2]-0.1)-exp(-x[1]-0.1);
            3*exp(x[1]+3*x[2]-0.1)-3*exp(x[1]-3*x[2]-0.1);
        ]
d2f2(x)  = [
            exp(x[1]+3*x[2]-0.1)+exp(x[1]-3*x[2]-0.1)+exp(-x[1]-0.1)   3*exp(x[1]+3*x[2]-0.1)-3*exp(x[1]-3*x[2]-0.1);
            3*exp(x[1]+3*x[2]-0.1)-3*exp(x[1]-3*x[2]-0.1)              9*exp(x[1]+3*x[2]-0.1)+9*exp(x[1]-3*x[2]-0.1);  
        ]
xsol2   = [-0.34657;0.00000]
push!(probs,("ConvexNonQuad",f2,df2,d2f2,xsol2))

# convex quadratic
A      = randn(3,3)
A      = A'*A + eye(3)*3
b      = randn(3);
f3(x)   = 0.5*dot(x,A*x) - dot(b,x)
df3(x)  = A*x - b
d2f3(x) = A
xsol3   = A\b
push!(probs,("ConvexQuad",f3,df3,d2f3,xsol3))

# test test functions
for k=1:length(probs)
	prob = probs[k]
	
	println("testing function $(prob[1])")
	# test xsol
	@test norm(prob[3](prob[5])) < 1e-5
	
	# setup test function for derivative check
	function testFun(x,v=[])
		if isempty(v)
			return prob[2](x)
		else
			return prob[2](x), prob[3](x)'*v
		end
	end
	
	pass, = checkDerivative(testFun,prob[5]+randn(size(prob[5])),out=false)
	@test pass
				
	# setup test function for derivative check
	function testFun(x,v=[])
		if isempty(v)
			return prob[3](x)
		else
			return prob[3](x), prob[4](x)'*v
		end
	end
	pass, = checkDerivative(testFun,prob[5]+randn(size(prob[5])),out=false)
	@test pass
end



