export trtest, trcauchy, trdogleg

# test if step is acceptable and update size of trust region
function trtest(f::Function,fc,df,H,xc,st,Delta::Real;
    rhom=1e-4, rhoL=.25, rhoH=0.75, wdown=0.5, wup=1.5,C=1e2)
    
    ared = fc - f(xc+st)
    pred = - dot(df,st) - dot(st,H*st)/2
    
    if ared/pred < rhom # reject and decrease Delta
        return xc,wdown*Delta,-2,ared
    elseif (rhom <= ared/pred) && (ared/pred < rhoL)
        return xc+st, wdown*Delta,-1,ared
    elseif (rhoL <= ared/pred) && (ared/pred < rhoH)
        return xc+st, Delta,0,ared
    elseif (rhoH < ared/pred) && (norm(st) <Delta)
        return xc+st, Delta,0,ared
    elseif (rhoH < ared/pred)
        return xc+st, min(wup*Delta,C*norm(df)),1,ared
    end
end

# implementation of cauchy method
function trcauchy(f::Function,J::Function,H::Function,x::Vector;
	Delta=2.0,maxIter=20,atol=1e-8,out::Int=0,storeInterm::Bool=false,
	trt::Function=(f,fc,df,d2f,x,st,Delta)->trtest(f,fc,df,d2f,x,st,Delta))

    his = zeros(maxIter,4)
    X   = (storeInterm) ? zeros(length(x),maxIter) : []
    i   = 1; flag = -1; alph=0.0
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
        
        # get step size
        d2f  = H(x)
        curv = dot(df,d2f*df)
        if curv <= 0 
            alph = Delta/norm(df)
        else
            alph = min(Delta/norm(df), norm(df)^2/curv)
        end

        # check for acceptance
        xn,Deltan,trflag,ared = trt(f,fc,df,d2f,x,-alph*df,Delta)
        his[i,3:4] = [trflag,Delta]
        if out>0
            @printf "iter=%04d\t|f|=%1.4e\t |df|=%1.2e\tTR=%d\tDelta=%1.2f\n" i his[i,1] his[i,2] his[i,3] his[i,4]
        end
        # update
        x     = copy(xn)
        Delta = Deltan
        i+=1   
    end
    i = min(maxIter,i)

    if out>=0
        if flag==-1
            println(@sprintf("trcauchy iterated maxIter (=%d) times but reached only atol of %1.2e instead of tol=%1.2e",i,his[i,2],atol))
        elseif out>1
            println(@sprintf("trcauchy achieved desired atol of %1.2e at iteration %d.",atol,i))
        end
    end
    if storeInterm; X = X[:,1:i]; end
    return x,flag,his,X
end

# dogleg method
function trdogleg(f::Function,J::Function,H::Function,x::Vector;
	Delta=2.0,maxIter=20,atol=1e-8,out::Int=0,storeInterm::Bool=false,
	trt::Function=(f,fc,df,d2f,x,st,Delta)->trtest(f,fc,df,d2f,x,st,Delta))

    his = zeros(maxIter,5)
    X   = (storeInterm) ? zeros(length(x),maxIter) : []
    i   = 1; flag = -1; alph=0.0
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
        
        # get cauchy point
        d2f  = H(x)
        curv = dot(df,d2f*df)
        if curv <= 0 
            alph = Delta/norm(df)
        else
            alph = min(Delta/norm(df), norm(df)^2/curv)
        end
        st = -alph*df
        xcp = x + st
        flagnt = 0

        if norm(xcp-x) < Delta
            # get Newton point
            try 
               sn = - d2f \ df
               xn = x + sn
               # check angle condition
               angl = dot(xn-xcp,xcp-x)
               if (angl>0) && (norm(xcp-x) < norm(sn)) && norm(sn) <= Delta
                  flagnt=1
                  st = sn # take newton point
               elseif (angl>0) # build dogleg
                  flagnt = 2
                  step = Delta/10; cnt = 1;
               
                  while (norm(xcp + cnt*step*(xn-xcp) - x) < Delta) && (cnt < 100)
                       cnt+=1
                  end
                  xt = xcp + (cnt-1)*step*(xn-xcp)
                  st = xt-x
               end
            catch err
                if ~isa(err,Base.LinAlg.SingularException)
                    throw(err)
                end
            end
        end 

        # check for acceptance
        xn,Deltan,trflag,ared = trt(f,fc,df,d2f,x,st,Delta)
        his[i,3:5] = [trflag,flagnt,Delta]
        if out>0
            @printf "iter=%04d\t|f|=%1.2e\t|df|=%1.2e\tTR=%d\tnt?=%d\tDelta=%1.2f\n" i his[i,1] his[i,2] his[i,3] his[i,4] his[i,5]
        end
        # update
        x     = xn
        Delta = Deltan
        i+=1   
    end
    i = min(maxIter,i)

    if out>=0
        if flag==-1
            println(@sprintf("trdogleg iterated maxIter (=%d) times but reached only atol of %1.2e instead of tol=%1.2e",i,his[i,2],atol))
        elseif out>1
            println(@sprintf("trdogleg achieved desired atol of %1.2e at iteration %d.",atol,i))
        end
    end
    if storeInterm; X = X[:,1:i]; end
    return x,flag,his,X
end

