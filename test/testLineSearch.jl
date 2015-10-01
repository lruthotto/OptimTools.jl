using OptimTools
using Base.Test

# test exact line search for quadratic
Q = 10*speye(40)
f(x)    = 0.5*dot(x,Q*x)
df(x)   = Q*x
d2f(x)  = Q
x0 = randn(40)

exactLS  = (f,df,fk,dfk,xk,pk)-> (-dot(pk,dfk)/dot(pk,Q*pk),1)
goldenLSf = (f,df,dk,dfk,xk,pk) -> goldenLS(f,xk,pk,0.0,5.0)

println("Test SD,NLCG,BFGS with exact line search ")
# steepest descent
x1,flag1,his1 = sd(f,df,x0,maxIter=1000,atol=1e-6,out=0,lineSearch=exactLS)
@test norm(x1)<1e-6*norm(x0)
x2,flag2,his2 = sd(f,df,x0,maxIter=1000,atol=1e-6,out=0, lineSearch=exactLS)
@test norm(x2)<1e-6
@test size(his2,1)==2 #should converge in one iteration

# nlcg
x1,flag1,his1 = nlcg(f,df,x0,maxIter=1000,atol=1e-6,out=0)
@test norm(x1)<1e-6*norm(x0)
x2,flag2,his2 = nlcg(f,df,x0,maxIter=1000,atol=1e-6,out=0, lineSearch=exactLS)
@test norm(x2)<1e-6
@test size(his2,1)==2 #should converge in one iteration

# bfgs
x1,flag1,his1 = bfgs(f,df,x0,maxIter=1000,atol=1e-6,out=0)
@test norm(x1)<1e-6*norm(x0)
x2,flag2,his2 = bfgs(f,df,x0,maxIter=1000,atol=1e-6,out=0, lineSearch=exactLS)
@test norm(x2)<1e-6
@test size(his2,1)==2 #should converge in one iteration

println("Done!")

println("Test SD,NLCG,BFGS with golden line search ")
# steepest descent
x1,flag1,his1 = sd(f,df,x0,maxIter=1000,atol=1e-6,out=0,lineSearch=goldenLSf)
@test norm(x1)<1e-6*norm(x0)
x2,flag2,his2 = sd(f,df,x0,maxIter=1000,atol=1e-6,out=0, lineSearch=goldenLSf)
@test norm(x2)<1e-6

# nlcg
x1,flag1,his1 = nlcg(f,df,x0,maxIter=1000,atol=1e-6,out=0)
@test norm(x1)<1e-6*norm(x0)
x2,flag2,his2 = nlcg(f,df,x0,maxIter=1000,atol=1e-6,out=0, lineSearch=goldenLSf)
@test norm(x2)<1e-6

# bfgs
x1,flag1,his1 = bfgs(f,df,x0,maxIter=1000,atol=1e-6,out=0)
@test norm(x1)<1e-6*norm(x0)
x2,flag2,his2 = bfgs(f,df,x0,maxIter=1000,atol=1e-6,out=0, lineSearch=goldenLSf)
@test norm(x2)<1e-6

println("Done!")
