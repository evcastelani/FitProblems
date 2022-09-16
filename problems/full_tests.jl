using Dataframes
run_probs() = begin
    # colocar em dataframe
    list = ["parabola_-0.5_-9.0_0.0_650_143.csv"]
    for file in list
        prob = load_problem(file)
        solve(prob,zeros(prob.dim),"LMlovo")
    end

end
