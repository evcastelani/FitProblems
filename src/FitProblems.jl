module FitProblems

using DelimitedFiles, Optim,RAFF

export load_problem, solve, build_problem

import Base.show

"""
    FitProbType

It is an immutable type used by main functions of this package

"""
struct FitProbType
    name::String
    data:: Array{Float64,2}
    npts:: Int
    nout:: Int
    model::Function
    dim:: Int
    cluster::Bool
    noise::Bool
    solution::Array{Float64,1}
    description::String
end

struct FitOutputType
    status::Bool
    solution::Vector{Float64}
    niter :: Int
    minimum :: Float64
    feval :: Int
end

"""
    load_problem(filename::String)

This function is used to load a problem from a csv file and convert to FitProbType. It is an important function because FitProbType is the unique supported format in this package. 

# Examples
```
julia-repl
julia> load_problem("toy.csv")

returns a FitProbType
```
"""
function load_problem(filename::String)
    prob_matrix = readdlm(filename,':')
    return FitProbType(prob_matrix[1,2],eval(Meta.parse(prob_matrix[2,2])),prob_matrix[3,2],prob_matrix[4,2],eval(Meta.parse(prob_matrix[5,2])),prob_matrix[6,2],prob_matrix[7,2],prob_matrix[8,2],eval(Meta.parse(prob_matrix[9,2])),prob_matrix[10,2])
end


function solve(prob::FitProbType,θinit::Vector{Float64},method::String)
    m = prob.npts
    p = m-prob.nout
    A = prob.data
    model = prob.model
    dim = prob.dim

    F(θ::Vector{Float64}) = begin
        V = zeros(m)
        for i=1:m
            V[i] = (A[i,end]-model(A[i,1:end-1],θ))^2
        end
        return V
    end

    sortF(V::Vector{Float64}) = begin
        for i = 1:p+1
            for j = i + 1:m
                if V[i] > V[j]          
                    vaux = V[j]
                    V[j] = V[i]
                    V[i] = vaux
                end
            end
        end
        return V[1:p]    
    end

    lovo(θ) = begin
        return sum(sortF(F(θ)))
    end

    if method == "NelderMead"
       res = Optim.optimize(lovo,θinit , NelderMead())
       return FitOutputType(Optim.converged(res),Optim.minimizer(res),Optim.iterations(res),Optim.minimum(res),Optim.f_calls(res))
    end
    if method == "SimulatedAnnealing"
       res = Optim.optimize(lovo,θinit , SimulatedAnnealing())
       return FitOutputType(Optim.converged(res),Optim.minimizer(res),Optim.iterations(res),Optim.minimum(res),Optim.f_calls(res))
    end
    if method == "ParticleSwarm"
       res = Optim.optimize(lovo,θinit , ParticleSwarm())
       return FitOutputType(Optim.converged(res),Optim.minimizer(res),Optim.iterations(res),Optim.minimum(res),Optim.f_calls(res))
    end
    if method == "LMlovo"
        res = RAFF.lmlovo(model, θinit, A, dim,p)
        return FitOutputType(res.status,res.solution,res.iter,res.f,res.nf)
    end

    if method == "RAFF"
        res = RAFF.raff(model, A, dim;initguess = θinit)
        return FitOutputType(res.status,res.solution,res.iter,res.f,res.nf)
    end

end

function build_problem(probtype::String,limit::Vector{Float64},params::Vector{Float64})
    if probtype == "line"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,npts,nout]")
        npts = Int(params[3])
        nout = Int(params[4])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        y = params[1]*x .+params[2]
        nout = Int(params[4])
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            x[iout[k]]=x[iout[k]]+randn()
            y[iout[k]]=y[iout[k]]+randn()
        end


        FileMatrix = ["name :" "Line";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*x[1] +t[2]";"dim :" 2; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:2]]; "description :" "none"]

        open("line_$(params[1])_$(params[2])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end

    if probtype == "parabola"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,c,npts,nout]")
        npts = Int(params[4])
        nout = Int(params[5])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        y = params[1]*x.^2 .+params[2]*x .+params[3]
        nout = Int(params[5])
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            x[iout[k]]=x[iout[k]]+randn()
            y[iout[k]]=y[iout[k]]+randn()
        end


        FileMatrix = ["name :" "Parabola";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*x[1]^2 +t[2]*x[1] + t[3]";"dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:3]]; "description :" "none"]

        open("parabola_$(params[1])_$(params[2])_$(params[3])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end
    if probtype == "cubic"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,c,d,npts,nout]")
        npts = Int(params[5])
        nout = Int(params[6])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        y = params[1]*x.^3 .+params[2]*x.^2 .+params[3]*x .+params[4]
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            x[iout[k]]=x[iout[k]]+randn()
            y[iout[k]]=y[iout[k]]+randn()
        end


        FileMatrix = ["name :" "Cubic";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*x[1]^3 +t[2]*x[1]^2 + t[3]*x[1] +t[4]";"dim :" 4; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:4]]; "description :" "none"]

        open("cubic_$(params[1])_$(params[2])_$(params[3])_$(params[4])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end

    if probtype == "sphere2D"
        println("params need to be setup as [center,radious,npts,nout]")
        c = [params[1],params[2]]
        r = params[3]
        npts = Int(params[4])
        x = zeros(npts)
        y = zeros(npts)
        θ = [0.0:2*π/(npts-1):2*π;]
        for k=1:npts
            x[k] = c[1]+r*cos(θ[k])
            y[k] = c[2]+r*sin(θ[k])
        end
        nout = Int(params[5])
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end
        for k = 1:nout
            x[iout[k]]=x[iout[k]]+rand([0.25*r:0.1*(r);(1+0.25)*r])
            y[iout[k]]=y[iout[k]]+rand([0.25*r:0.1*(r);(1+0.25)*r])
        end
        FileMatrix = ["name :" "sphere2D";"data :" [[x y zeros(npts)]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - t[3]^2";"dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [push!(c,r)]; "description :" "none"]

        open("sphere2D_$(c[1])_$(c[2])_$(c[3])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end

    end
    if probtype == "gaussian"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,c,npts,nout]")
        npts = Int(params[4])
        nout = Int(params[5])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        g(w) = params[1]*exp(-((w-params[2])^2)/(params[3]^2))
        y = g.(x)
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            x[iout[k]]=x[iout[k]]+randn()
            y[iout[k]]=y[iout[k]]+randn()
        end


        FileMatrix = ["name :" "gaussian";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*exp(-((x[1]-t[2])^2)/(t[3]^2))";"dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:3]]; "description :" "none"]

        open("gaussian_$(params[1])_$(params[2])_$(params[3])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end

    end 

    if probtype == "log"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,c,npts,nout]")
        npts = Int(params[4])
        nout = Int(params[5])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        y = params[1]*log.((params[2]*x .+params[3]))
        nout = Int(params[5])
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            x[iout[k]]=x[iout[k]]+randn()
            y[iout[k]]=y[iout[k]]+randn()
        end


        FileMatrix = ["name :" "log";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*log((t[2]*x[1] +t[3]))";"dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:3]]; "description :" "none"]

        open("log2_$(params[1])_$(params[2])_$(params[3])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end

    if probtype == "trig"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,c,npts,nout]")
        npts = Int(params[4])
        nout = Int(params[5])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        y = params[1]*sin.(params[2]*x) .+ params[3]
        nout = Int(params[5])
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            x[iout[k]]=x[iout[k]]+randn()
            y[iout[k]]=y[iout[k]]+randn()
        end


        FileMatrix = ["name :" "trig";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*sin(t[2]*x[1]) +t[3]";"dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:3]]; "description :" "none"]

        open("trig_$(params[1])_$(params[2])_$(params[3])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end

    
end

function show(io::IO, fout::FitOutputType)

    print(io,"  ▶ Output ◀ \n")
    if Bool(fout.status) ==true
        print(io,"  ↳ Status (.status) = Convergent \n")
    else
        print(io,"  ↳ Status (.status) = Divergent \n")
    end
    print(io,"  ↳ Solution (.solution) = $(fout.solution) \n")
    print(io,"  ↳ Number of iterations (.niter) = $(fout.niter) \n")
    print(io,"  ↳ Minimum (.minimum) = $(fout.minimum) \n")
    print(io,"  ↳ Number of function calls (.feval) = $(fout.feval) \n")
end


end # module
