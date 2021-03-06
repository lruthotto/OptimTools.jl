export armijo
export goldenLS

"""
armijo(f,fk,dfk,xk,pk)

Backtracked Armijo linesearch
"""
function armijo(f::Function,fk,dfk,xk,pk;maxIter=10, c1=1e-4,b=0.5)
LS = 1
t  = 1
while LS<=maxIter
    if f(xk+t*pk) <= fk + t*c1*dot(dfk,pk)
        break
    end
    t *= b
    LS += 1
end
if LS>maxIter
	LS= -1
	t = 0.0
end
return t,LS
end

"""
goldenLS(f,xk,pk,a,b)

Golden section linesearch
"""
function goldenLS(f::Function,xk,pk,a,b;atol=1e-2,maxIter=40)
	ratio = (sqrt(5)-1)/2.0
	# make sure points we have an interval
	a,b   = (a<b)? (a,b) :  (b,a)
	# evaluate f at left point
	al    = b - ratio*(b-a)
	fl    = f(xk + al*pk)
	# evaluate f at right point
	ar    = a + ratio*(b-a)
	fr    = f(xk + ar*pk)
	
	LS    = 0
	while (LS<maxIter && abs(ar-al)>atol)
		if fl<fr # continue search in [a,ar]
			b  = ar
			ar = al
			fr = fl
			al = b - ratio*(b-a)
			fl = f(xk + al*pk)
		else # continue search in [al,b]
			a  = al
			al = ar
			fl = fr
			ar = a + ratio*(b-a)
			fr = f(xk + ar*pk)
		end
		# @printf "iter=%d, [al,ar]=[%1.2f,%1.2f], [fl,fr]=[%1.2e,%1.2e]\n" LS al ar fl fr
		LS += 1
	end
	return (al+ar)/2,LS
end
