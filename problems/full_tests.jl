using DataFrames, DelimitedFiles, FitProblems
using LinearAlgebra, BenchmarkTools
using CSV
BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000
#function build_results(filename::String,method::String)
filename = "list.txt"
methods = ["LMlovo","NelderMead","SimulatedAnnealing","ParticleSwarm","RAFF"]
list = String.(readdlm(filename))
io = open("log_read_erros.txt", "w");
for methodname in methods
    # colocar em dataframe
    df = DataFrame(prob = String[], error= Float64[],minimum=Float64[],feval=Float64[],niter=Int64[],convergence=Bool[], minimum_PT = Float64[], maximum_PT = Float64[], mean_PT = Float64[], median_PT = Float64[])
    #list = ["parabola_-0.5_-9.0_0.0_650_143.csv"]
    for file in list
        local prob = load_problem(file)
        println(" 🏁 Testing problem $(file) using $(methodname)")
        try  s = solve(prob,zeros(prob.dim),methodname)
            println(" 🏆 Problem free from read errors!")
            println(" 🕥 Performing benchmark! Please wait!")
            local b = @benchmark solve($(prob),zeros($(prob.dim)),$(methodname))
            display(b)
            push!(df,[file,norm(prob.solution-s.solution,Inf),s.minimum,s.feval,s.niter,s.status,minimum(b.times),maximum(b.times),mean(b.times),median(b.times)])
        catch
            println(" 💣 Some read error in problem $(file) occurred")
            write(io," Error reading $(file) using $(methodname) occurred \n")
        end
    end
    CSV.write("$(methodname)Output.csv", df)
end
close(io)
#end
