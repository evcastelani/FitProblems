using DelimitedFiles

function fit_res(filename)
    Data = Float64.(readdlm(filename)[2:end,:])
    model(x,θ) = θ[1]*x[1]^2+θ[2]*x[1]+θ[3]
    #model(x,θ) = θ[1]*exp(-θ[2]*x[1])
    s = raff(model,Data,3)
    plt = plot()
    scatter!(plt,Data[:,1],Data[:,2],markercolor=:blue)
    scatter!(plt,Data[s.outliers,1],Data[s.outliers,2],markercolor=:red)
    x = [2:0.01:4.5;]
    y = map(t->model(t,s.solution),x)
    #z = map(t->model(t,[2724.9,0.572]),x)
    plot!(plt,x,y,lw=2)
    #plot!(plt,x,z,lw=2)
end
