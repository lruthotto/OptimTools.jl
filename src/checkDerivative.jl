export checkDerivative

"""
checkDerivative(f,x0;out,tol,nSuccess)

Check if function f provides correct derivatives around point x0 using a Taylor expansion.

Inputs:

	f::Function  - fx,dfx= f(x) provides function value and derivative
	x::Array     - vector around which to test the derivatives
	out::Bool    - flag for activating printing (default: true)
	tol::Number  - tolerance to judge if error in linearization decreases quadratically (default: 1.9)
	nSuccess:Int - number of quadratically-decreasing steps needed to decide that derivative is correct
	
Outputs:
	pass, Error, Order

"""
function checkDerivative(f::Function,x0;out::Bool=true,tol::Number=1.9,nSuccess::Int=3)
	if out
		println(@sprintf("%9s\t%9s\t%9s\t%9s\t%9s\t%5s","h","E0","E1","O1","O2","OK?"))
	end
	v = randn(size(x0))
	if eltype(x0) <: Complex
		v += 1im*randn(size(v))
	end
	f0,dvf = f(x0,v)
	dvf   = real(dvf)
	Error = zeros(10,2)
	Order = zeros(10,2)
	Success = zeros(10)
	for j=1:10
		ft = f(x0+10.0^(-j)*v)            # function value
		Error[j,1] = norm(f0-ft)/norm(f0)           # Error TaylorPoly 0
		Error[j,2] = norm(f0 .+10.0^(-j)*dvf .- ft)/norm(f0) # Error TaylorPoly 1
		if j>1
			Order[j,:] = log10.(Error[j-1,:]./Error[j,:]);
		end
		if (Order[j,2]>tol) || (Error[j,1]/Error[j,2] > 100); Success[j]=1; end
		if out 
			println(@sprintf("%1.3e\t%1.3e\t%1.3e\t%1.3e\t%1.3e\t%5d",
							10.0^(-j), Error[j,1],Error[j,2], Order[j,1],Order[j,2],Success[j]))
		end
	end
	pass = sum(Success) > nSuccess
	return  pass,Error,Order
end