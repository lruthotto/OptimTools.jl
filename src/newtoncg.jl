export newtoncg

function newtoncg(f::Function,J::Function,H::Function,x::Vector;maxIter=20,atol=1e-8,out::Int=0,storeInterm::Bool=false,
	lineSearch::Function=(f,J,fk,dfk,xk,pk)->armijo(f,fk,dfk,xk,pk,maxIter=20),tolCG=1e-2,maxIterCG=30)

    his = zeros(maxIter,6)
    X = (storeInterm) ? zeros(length(x),maxIter) : []
    i = 1; flag = -1; LL = []
    while i<=maxIter
        fc = f(x)
        df = J(x)
        his[i,1:2] = [norm(fc) norm(df)]
        if storeInterm; X[:,i] = x; end;

        if(norm(df)<atol)
            flag = 0
            his = his[1:i,:]
            break
        end
        
        # get search direction
        d2f = H(x)
        pk,his[i,6],his[i,5],his[i,4] = cg(d2f,-df,out=-1,tol=tolCG,maxIter=maxIterCG)
        if his[i,6]==-2 && norm(pk)==0
            pk = -df
        end
      # line search
        ak,his[i,3] = lineSearch(f,J,fc,df,x,pk) 
        if his[i,3]==-1
            flag = -3
            his  = his[1:i,:]
            break
        end
      if out>0
            @printf "iter=%04d\t|f|=%1.2e\t|df|=%1.2e\tcgIter=%d\tLS=%1.4f\n" i his[i,1] his[i,2] his[i,4] ak
        end
        # update
        x  += ak*pk
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