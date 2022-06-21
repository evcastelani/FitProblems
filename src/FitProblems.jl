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
    p = prob.nout
    A = prob.data
    model = prob.model
    dim = prob.dim
    
    F(θ::Vector{Float64}) = begin
        V = zeros(m)
        for i=1:m
            V[i] = (A[i,2]-model(A[i,1],θ))^2
        end
        return V
    end

    sortF(V::Vector{Float64}) = begin
        for i = 1:p
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

    if method == "RAFF"
        res = RAFF.raff(model, A, dim;initguess = θinit)
        return FitOutputType(res.status,res.solution,res.iter,res.f,res.nf)
    end

end

function build_problem()

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
