export modnewton

function modnewton(f::Function,J::Function,H::Function,x::Vector;maxIter=20,atol=1e-8,out::Int=0,
	storeInterm::Bool=false,beta=1e-1)

    his = zeros(maxIter,2)
    X = (storeInterm) ? zeros(length(x),maxIter) : []
    i = 1; flag = -1; LL = []
    while i<=maxIter
        fc = f(x)
        df = J(x)
        his[i,:] = [norm(fc) norm(df)]
        if storeInterm; X[:,i] = x; end;

        if(norm(df)<atol)
            flag = 0
            his = his[1:i,:]
            break
        end
        
        # get search direction
        d2f = H(x)
        tau = (minimum(diag(d2f)) > 0) ? 0 :  -minimum(diag(d2f))+beta
        for k=1:20
            try
                LL    = cholfact(d2f + tau*speye(length(x))) 
                break;
            catch err
                if isa(err,Base.LinAlg.PosDefException) && k<20
                    tau = max(2*tau,beta)
                elseif isa(err,Base.LinAlg.PosDefException) 
                    flag = -2
                    his = his[1:i,:]
                    break; break;
                else
                    throw(err)
                end
            end
        end
        pk    = - (LL\df)
        if out>0
            @printf "iter=%04d\t|f|=%1.2e\t|df|=%1.2e\n" i his[i,1] his[i,2]
        end
        # update
        x  += pk
        i+=1   
    end
    i = min(maxIter,i)

    if out>=0
        if flag==-1
            println(@sprintf("newton iterated maxIter (=%d) times but reached only atol of %1.2e instead of tol=%1.2e",i,his[i,2],atol))
        elseif flag==-2
            println(@sprintf("newton stopped because the Hessian at iteration %d was not positive definite.",i))
        elseif out>1
            println(@sprintf("newton achieved desired atol of %1.2e at iteration %d.",atol,i))
        end
    end
    if storeInterm; X = X[:,1:i]; end
    return x,flag,his,X
end