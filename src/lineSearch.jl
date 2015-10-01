export armijo
export goldenLS

function armijo(f::Function,fc,df,xc,pc;maxIter=10, c1=1e-4,b=0.5)
LS = 1
t  = 1
while LS<=maxIter
    if f(xc+t*pc) <= fc + t*c1*dot(df,pc)
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

function goldenLS(f::Function,xc,pc,a,b;atol=1e-2,maxIter=40)
	ratio = (sqrt(5)-1)/2.0
	# make sure points we have an interval
	a,b   = (a<b)? (a,b) :  (b,a)
	# evaluate f at left point
	al    = b - ratio*(b-a)
	fl    = f(xc + al*pc)
	# evaluate f at right point
	ar    = a + ratio*(b-a)
	fr    = f(xc + ar*pc)
	
	LS    = 0
	while (LS<maxIter && abs(ar-al)>atol)
		if fl<fr # continue search in [a,ar]
			b  = ar
			ar = al
			fr = fl
			al = b - ratio*(b-a)
			fl = f(xc + al*pc)
		else # continue search in [al,b]
			a  = al
			al = ar
			fl = fr
			ar = a + ratio*(b-a)
			fr = f(xc + ar*pc)
		end
		# @printf "iter=%d, [al,ar]=[%1.2f,%1.2f], [fl,fr]=[%1.2e,%1.2e]\n" LS al ar fl fr
		LS += 1
	end
	return (al+ar)/2,LS
end
