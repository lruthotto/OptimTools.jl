export bfgs

function bfgs(f::Function,J::Function,x::Vector;H=speye(length(x)), maxIter=20,atol=1e-8,out::Int=0,storeInterm::Bool=false,
	lineSearch::Function=(f,J,fk,dfk,xk,pk)->armijo(f,fk,dfk,xk,pk,maxIter=30))

    his = zeros(maxIter,3)
    I   = speye(length(x))
    df  = J(x)
    X   = (storeInterm) ? zeros(length(x),maxIter) : []
 
    i = 1; flag = -1    
    while i<=maxIter
        fc = f(x)

        his[i,1:2] = [fc norm(df)]
        if storeInterm; X[:,i] = x; end;
        if norm(df)<atol
            his  = his[1:i,:]
            flag = 0
            break
        end

        # get search direction
        pk    = - H*df
        # line search
        ak,his[i,3] = lineSearch(f,J,fc,df,x,pk) 
        if out>0
            @printf "iter=%4d\t|f|=%1.2e\t|df|=%1.2e\tLS=%d\n" i his[i,1] his[i,2] his[i,3]
        end
        if his[i,3]==-1
             flag = -3
             his  = his[1:i,:]
             break;
        end
        x    += ak*pk
        dfnew = J(x)
        sk    = ak*pk
        yk    = dfnew - df
        if dot(yk,sk)>0 # ensure that approximate Hessians remain positive definite
        	H     = (I - (sk*yk')/dot(sk,yk)) * H * (I - (yk*sk')/dot(sk,yk)) + (sk*sk')/dot(yk,sk)
        else
            if out>0
                warn("bfgs detected negative curvature. Resetting Hessian")
            end
            H = speye(length(x))
        end
        df   = dfnew
        i+=1     
    end
    i = min(maxIter,i)

    if out>=0
        if flag==-1
            println(@sprintf("bfgs iterated maxIter (=%d) times but reached only atol of %1.2e instead of tol=%1.2e",i,his[i,2],atol))
        elseif flag==-3
            println(@sprintf("bfgs stopped of a line search fail at iteration %d.",i))
        elseif out>1
            println(@sprintf("bfgs achieved desired atol of %1.2e at iteration %d.",atol,i))
        end
    end
    
    if storeInterm; X = X[:,1:i]; end
    return x,flag,his,X,H
end