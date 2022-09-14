#=
Péndulo elástico usando integradores simplécticos por Jose Torrente Teruel
Constantes: σ = ω_z = ω_ϕ = 1 
=#
using Plots, Calculus

sf = true # save figures
sub="D:/Escritorio/U/Master/Practicum/Pendulum/figures/"

#q = (z, ϕ)
function V(q::Vector{<:Real})
        return 0.5*(q[1]+1)^2 - (q[1]+1)*cos(q[2])
end

function ∂zV(q::Vector{<:Real})
        return q[1]+1 -cos(q[2])
end

function ∂ϕV(q::Vector{<:Real})
        return (q[1]+1)*sin(q[2])
end

function E_d(q::Vector{<:Real},v::Vector{<:Real})
   T=0.5*(v[1]^2+v[2]^2 * (1+q[1])^2)
   U=V(q)
   return T-U
end

#Inicialización
epochs = 2^12;
#h=0.100462;
h=0.0035

q=Array{Float64,2}(undef,epochs,2)
μ=Array{Float64,2}(undef,epochs,2)
v=Array{Float64,2}(undef,epochs,2)
H=Array{Float64,1}(undef,epochs)

#Condiciones iniciales
z0=0.011
z1=0.011
ϕ0=-0.23
ϕ1=-0.225
q[1:2,:]=[z0 ϕ0
          z1 ϕ1];
v[1,:]=(q[2,:]-q[1,:])/h;  #1/0.009^2 ∼ 10000
μ[1,:]=[v[1,1]/h - v[1,2]^2 * (1+q[1,1]) + 1/2 * ∂zV(q[1,:]), 
        v[1,2]/h * (1+q[1,1])^2 + 0.5 * ∂ϕV(q[1,:])];
μ[2,:]=[v[1,1]/h - 1/2 * ∂zV(q[2,:]),
        v[1,2]/h * (1+q[1,1])^2 - 0.5* ∂ϕV(q[2,:])];

H[1]=E_d((q[1,:]+q[2,:])/2,v[1,:])

#Iteraciones
for k=2:epochs-1
    v[k,2]=h/(1+q[k,1])^2 * (μ[k,2] - 1/2*∂ϕV(q[k,:]));    #μ_k = -D_1 L_d^k -> v_k=Δq_k / h
    v[k,1]=h*(μ[k,1]+ v[k,2]^2 * (1+q[1,1]) - 1/2 * ∂zV(q[k,:]));

    #=if any(i -> isnan(i),v[k,:])
        print("v($k) is nan\n")
        break
    end=#

    q[k+1,:]=q[k,:]+h*v[k,:];

    μ[k+1,:]=[v[k,1]/h - 1/2 * ∂zV(q[k+1,:]),
              v[k,2]/h * (1+q[k,1])^2 - 0.5*∂ϕV(q[k+1,:])];

    H[k]=E_d((q[k,:]+q[k+1,:])/2,v[k,:])
end

##Plots and gifs
display(plot(q[1:end,1],label="",xlabel = "paso k",#=H[1:end-1,1],legend=:bottomleft,=#title="Pendulo elástico: elongación z")) 
sf && savefig(sub*"elastic_length_$h-$z0-$z1-$ϕ0-$ϕ1.svg")
display(plot(q[1:end,2],label="",xlabel = "paso k",#=H[1:end-1,1],legend=:bottomleft,=#title="Péndulo elástico: ángulo φ")) 
sf && savefig(sub*"elastic_angle_$h-$z0-$z1-$ϕ0-$ϕ1.svg")
display(plot(H[1:end-1],label="",xlabel = "paso k",#=H[1:end-1,1],legend=:bottomleft,=#title="Péndulo elástico: energía")) 
sf && savefig(sub*"elastic_energy_$h-$z0-$z1-$ϕ0-$ϕ1.svg")


xs=Vector{Float64}(undef,epochs)
ys=Vector{Float64}(undef,epochs)

xmin=Float64[]; xmax=Float64[]; ymin=Float64[]; ymax=Float64[]
xlim=Tuple[]; ylim=Tuple[]

xs=q[:,1] .* sin.(q[:,2])
ys=-q[:,1] .* cos.(q[:,2])

xmin = minimum(xs); xmax = maximum(xs)
ymin = minimum(ys); ymax = maximum(ys)
xlim = tuple(([xmin xmax] + (xmax-xmin)/20*[-1 1])...)
ylim = tuple(([ymin ymax] + (ymax-ymin)/20*[-1 1])...)

plt=plot(title="Péndulo elástico: trayectorias",legend=:topleft,xlim=xlim,ylim=ylim, size = (1024, 576))
plot([xs[1]], [ys[1]], #=xlim=xlim,ylim=ylim,=# label="")
anim= @animate for i in 1:10:epochs
    plot(title="Péndulo elástico: trayectorias",xlabel = "ϕ(rad)",legend=:topleft, size = (1024, 576))
    plot!([0],[0],xlim=xlim,ylim=ylim,label="")
    plot!(xs[1:i],ys[1:i],linewidth=2,color=1,label="") #trajectory 
    plot!([0,xs[i]],[0,ys[i]],color=1,label="") #rope 
    scatter!(xs[i,:],ys[i,:],color=3, label="") #weight of the pendulum

    i==epochs-1 && sf && savefig("trajectories.svg")
end
gif(anim,sub*"pendulum_elastic_$h-$z0-$z1-$ϕ0-$ϕ1.mp4",fps=20)