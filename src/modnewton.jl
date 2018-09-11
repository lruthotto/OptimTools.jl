export modnewton

"""
modnewton(f,df,H,x)

Newton method with Hessian shift for solving min_x f(x)

"""
function modnewton(f::Function,df::Function,H::Function,x::Vector;maxIter=20,atol=1e-8,out::Int=0,
	storeInterm::Bool=false,beta=1e-1,
		lineSearch::Function=(f,df,fk,dfk,xk,pk)->armijo(f,fk,dfk,xk,pk,maxIter=10))

    his = zeros(maxIter,4)
    X = (storeInterm) ? zeros(length(x),maxIter) : []
    i = 1; flag = -1; LL = []
    while i<=maxIter
        fk = f(x)
        dfk = df(x)
        his[i,1:2] = [norm(fk) norm(dfk)]
        if storeInterm; X[:,i] = x; end;

        if(norm(dfk)<atol)
            flag = 0
            his = his[1:i,:]
            break
        end
        
        # get search direction
        d2f = H(x)
        tau = (minimum(diag(d2f)) > 0) ? 0 :  -minimum(diag(d2f))+beta
        nfac = 0;
        for k=1:20
            try
                nfac +=1
                LL    = cholfact(d2f + tau*speye(length(x))) 
                break;
            catch err
                if isa(err,Base.LinAlg.PosDefException) && k<20
                    tau = max(2*tau,beta)
                elseif isa(err,Base.LinAlg.PosDefException) 
                    flag = -2
                    his[i,4] = nfac
                    his = his[1:i,:]
                    break; break;
                else
                    throw(err)
                end
            end
        end
        his[i,4] = nfac
        pk    = - (LL\dfk)
   
        # line search
		ak,his[i,3] = lineSearch(f,df,fk,dfk,x,pk) 
        if his[i,3]==-1
            flag = -3
            his  = his[1:i,:]
            break
        end
        if out>0
            @printf "iter=%04d\t|f|=%1.2e\t|df|=%1.2e\tLS=%d\n" i his[i,1] his[i,2] his[i,3]
        end
        # update x and H
        x  += ak*pk

        i+=1   
    end
    i = min(maxIter,i)

    if out>=0
        if flag==-1
            println(@sprintf("modnewton iterated maxIter (=%d) times but reached only atol of %1.2e instead of tol=%1.2e",i,his[i,2],atol))
        elseif flag==-2
            println(@sprintf("modnewton stopped because the Hessian at iteration %d was not positive definite.",i))
        elseif flag==-3
            println(@sprintf("modnewton stopped of a line search fail at iteration %d.",i))
		elseif out>1
            println(@sprintf("modnewton achieved desired atol of %1.2e at iteration %d.",atol,i))
        end
    end
    if storeInterm; X = X[:,1:i]; end
    return x,flag,his,X
end