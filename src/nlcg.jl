export nlcg

"""
nlcg(f,df,x)

Nonlinear conjugate gradient method for solving min_x f(x)

"""
function nlcg(f::Function,df::Function,x::Vector;maxIter=20,atol=1e-8,out::Int=0,storeInterm::Bool=false,
	lineSearch::Function=(f,df,fk,dfk,xk,pk)-> armijo(f,fk,dfk,xk,pk,maxIter=30) )

    his = zeros(maxIter,3)
    X = (storeInterm) ? zeros(length(x),maxIter) : []

	fk    = f(x)
    dfk   = df(x)
    dfOld = copy(dfk)
    pk    = -dfk
   
    i = 1; flag = -1;
    while i<=maxIter
        his[i,1:2] = [fk norm(dfk)]
        if storeInterm; X[:,i] = x; end;
        
        if (norm(dfk)<atol)
            flag=0
            his = his[1:i,:]
            break
        end
        
        # line search
        ak,his[i,3] = lineSearch(f,df,fk,dfk,x,pk)
        if out>0
            @printf "iter=%04d\t|f|=%1.2e\t|df|=%1.2e\tLS=%d\n" i his[i,1] his[i,2] his[i,3]
        end
        if his[i,3]==-1
            flag = -3
            his = his[1:i,:]
            break;
        end
        
        # update x and H
        x  += ak*pk
        fk  = f(x)
        dfk  = df(x)
        
        # beta  = dot(df, df-dfOld)/norm(dfOld)^2 # Polak-Ribiere
        beta  = norm(dfk)^2/norm(dfOld)^2 # Fletcher-Reeves
        pk    = -dfk  + beta*pk
        dfOld = copy(dfk)
        i+=1   
    end
    i = min(maxIter,i)

    if out>=0
        if flag==-1
            println(@sprintf("nlcg iterated maxIter (=%d) times but reached only atol of %1.2e instead of tol=%1.2e",i,his[i,2],atol))
        elseif flag==-3
            println(@sprintf("nlcg stopped of a line search fail at iteration %d.",i))
        elseif out>1
            println(@sprintf("nlcg achieved desired atol of %1.2e at iteration %d.",atol,i))
        end
    end
    if storeInterm; X = X[:,1:i]; end
    return x,flag,his,X
end