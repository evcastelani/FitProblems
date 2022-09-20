using DataFrames, DelimitedFiles, FitProblems
using LinearAlgebra
using CSV
#function build_results(filename::String,method::String)
filename = "list.txt"
methods = ["LMlovo","NelderMead"]
list = String.(readdlm(filename))
for methodname in methods
    # colocar em dataframe
    df = DataFrame(prob = String[], error= Float64[],minimum=Float64[],feval=Float64[],niter=Int64[],convergence=Bool[])
    #list = ["parabola_-0.5_-9.0_0.0_650_143.csv"]
    for file in list
        prob = load_problem(file)
        display(prob)
        s = solve(prob,zeros(prob.dim),methodname)
        push!(df,[file,norm(prob.solution-s.solution,Inf),s.minimum,s.feval,s.niter,s.status])
    end
    CSV.write("$(methodname)Output.csv", df)
end
#end
