export sd

function sd(f::Function,J::Function,x::Vector;maxIter=20,atol=1e-8,out::Int=0,storeInterm::Bool=false)

    his = zeros(maxIter,3)
    X = (storeInterm) ? zeros(length(x),maxIter) : []
    flag = -1; i=1;

    while i<=maxIter
	    
        fc = f(x)
        df = J(x)
        his[i,1:2] = [fc norm(df)]
        if storeInterm; X[:,i] = x; end;
	
        if (norm(df)<atol)
	       flag = 0 
           his = his[1:i,:]
           break
        end

        # get search direction
        pk    = -df
        # line search
        ak,his[i,3] = armijo(f,fc,df,x,pk,maxIter=30) 
        if out>0
            @printf "iter=%04d\t|f|=%1.2e\t|df|=%1.2e\tLS=%d\n" i his[i,1] his[i,2] his[i,3]
        end
        if his[i,3]==-1
            flag = -3
            his = his[1:i,:]
            break;
        end

        # update
        x    += ak*pk
        i+=1   
    end
	i = min(maxIter,i)

    if out>=0
        if flag==-1
            println(@sprintf("sd iterated maxIter (=%d) times but reached only atol of %1.2e instead of tol=%1.2e",i,his[i,2],atol))
        elseif flag==-3
            println(@sprintf("sd stopped of a line search fail at iteration %d.",i))
        elseif out>1
            println(@sprintf("sd achieved desired atol of %1.2e at iteration %d.",atol,i))
        end
    end
    if storeInterm; X = X[:,1:i]; end	
    return x,flag,his,X
end