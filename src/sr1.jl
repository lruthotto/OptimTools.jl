export sr1

function sr1(f::Function,J::Function,x::Vector;B=speye(length(x)),Delta=2.0, maxIter=20,atol=1e-8,out::Int=0,storeInterm::Bool=false,r=1e-8,
trt::Function=(f,fc,df,d2f,x,st,Delta)->trtest(f,fc,df,d2f,x,st,Delta))

    his = zeros(maxIter,5)
    df  = J(x)
    X   = (storeInterm) ? zeros(length(x),maxIter) : []
 
    i = 1; flag = -1; alph=0.0
    while i<=maxIter
        fc = f(x)

        his[i,1:2] = [fc norm(df)]
        if storeInterm; X[:,i] = x; end;
        if norm(df)<atol
            his  = his[1:i,:]
            flag = 0
            break
        end

        # get cauchy point
        curv = dot(df,B*df)
        if curv <= 0 
            alph = Delta/norm(df)
        else
            alph = min(Delta/norm(df), norm(df)^2/curv)
        end
        st     = -alph*df
        xcp    = x + st
        flagnt = 0

        if norm(xcp-x) < Delta
            # get Newton point
            try 
               sn = - B \ df
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
        xn,Deltan,trflag,ared = trt(f,fc,df,B,x,st,Delta)
        his[i,3:5] = [trflag,flagnt,Delta]
        if out>0
            @printf "iter=%04d\t|f|=%1.2e\t|df|=%1.2e\tTR=%d\tnt?=%d\tDelta=%1.2f\n" i his[i,1] his[i,2] his[i,3] his[i,4] his[i,5]
        end
        # update Hessian approximation
        dfnew = J(x+st)
        yk    = dfnew - df
        if abs(dot(st,yk-B*st)) >= r*norm(st)*norm(yk-B*st)
	       B +=  ((yk-B*st)*(yk-B*st)')/dot(st,yk-B*st)
        end
		x     = xn
		df    = J(x)
		Delta = Deltan
		i+=1
    end
    i = min(maxIter,i)

    if out>=0
        if flag==-1
            println(@sprintf("sr1 iterated maxIter (=%d) times but reached only atol of %1.2e instead of tol=%1.2e",i,his[i,2],atol))
        elseif out>1
            println(@sprintf("sr1 achieved desired atol of %1.2e at iteration %d.",atol,i))
        end
    end
    
    if storeInterm; X = X[:,1:i]; end
    return x,flag,his,X,B
end