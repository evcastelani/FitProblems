function generate!(probtype,nprob)
    for i=1:nprob
        l = rand([1.0:4.0;])
        a = rand([-10.0:0.5:10.0])
        b = rand([-10.0:0.5:10.0])
        c = rand([-10.0:0.5:10.0])
        npts = rand([100:50:1000;])
        nout = round(Int,((rand([5:1:25;]))/100.0)*npts)

        build_problem("parabola",[-l,l],[a,b,c,npts,nout])
    end
end
