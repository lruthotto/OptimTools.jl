export newton

"""
newton(f,df,H,x)

Newton method for solving min_x f(x)

"""
function newton(f::Function,df::Function,H::Function,x::Vector;maxIter=20,atol=1e-8,out::Int=0,storeInterm::Bool=false)

    his = zeros(maxIter,2)
    X = (storeInterm) ? zeros(length(x),maxIter) : []
    i = 1; flag = -1; LL = []
    while i<=maxIter
        fk  = f(x)
        dfk = df(x)
        his[i,:] = [fk norm(dfk)]
        if storeInterm; X[:,i] = x; end;

        if(norm(dfk)<atol)
            flag = 0
            his = his[1:i,:]
            break
        end
        
        # get search direction
        d2f = H(x)
        try
            LL    = cholfact(d2f)
        catch err
            if isa(err,Base.LinAlg.PosDefException)
                flag = -2
                his = his[1:i,:]
                break
            else
                throw(err)
            end
        end
        pk    = - (LL\dfk)
        if out>0
            @printf "iter=%04d\tf=%1.2e\t|df|=%1.2e\n" i his[i,1] his[i,2]
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