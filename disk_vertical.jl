#=
Disco vertical usando integradores simplécticos por Jose Torrente Teruel
Constantes: m=a=1
=#

using Plots, Calculus

sf = false # save figures
sub="D:/Escritorio/U/Master/Practicum/"

function L(v::Vector{<:Real})
    return 0.5*sum(v.^2)
end

function fc(ϕ0::Real, ϕ1::Real, Δθ::Real)
    ϕ12=(ϕ0+ϕ1)/2
    return [Δθ*cos(ϕ12),Δθ*sin(ϕ12)]
end

function fθ(ϕ1::Real,ϕ2::Real,Δx::Real,Δy::Real,Δθ::Real)
    return (Δθ + 2*(cos(ϕ1)*Δx + sin(ϕ1)*Δy))/(1+2*cos((ϕ2-ϕ1)/2))
end

#Inicializacion
epochs = 2^12
h=0.01

q=Array{Float64,2}(undef,epochs,4)   #q=[x, y, ϕ, θ]
Δq=Array{Float64,2}(undef,epochs,4)
H=Array{Float64}(undef,epochs-1)

#Condiciones iniciales
x0=-2.; y0=3.; ϕ0=2.0 
q[1:2,1:4]=[x0 y0 ϕ0 0.0
            -2. 3.1 1.9 0.1]
Δq[1,:]=q[2,:]-q[1,:] 
H[1]=L(q[1,:])

#Iteraciones
for k=1:epochs-2
    Δq[k+1,3]=Δq[k,3]   #Δϕ_k+1=Δϕ_k
    q[k+2,3]=Δq[k+1,3]+q[k+1,3]

    Δq[k+1,4]=fθ(q[k+1,3],q[k+2,3],Δq[k,1],Δq[k,2],Δq[k,4]) #Δθ
    q[k+2,4]=Δq[k+1,4]+q[k+1,4]

    Δq[k+1,1],Δq[k+1,2] = fc(q[k+1,3],q[k+2,3],Δq[k+1,4]) #Δx, Δy
    q[k+2,1]=Δq[k+1,1]+q[k+1,1]
    q[k+2,2]=Δq[k+1,2]+q[k+1,2]

    H[k+1]=L(Δq[k+1,:])
end

###Plots and gifs
display(plot(q[1:end-1,1],q[1:end-1,2],#=legend=:bottomleft,=#title="Disco vertical: trayectoria",label="")) 
sf && savefig(sub*"Euler disk/figures/vertical_trajectory_$h-$x0-$y0-$ϕ0.svg")
display(plot(q[1:end-1,3].%(2*pi),xlabel = "paso k",#=legend=:bottomleft,=#title="Disco vertical: ángulo φ",label=""))
sf && savefig(sub*"Euler disk/figures/vertical_angle_$h-$x0-$y0-$ϕ0.svg")
display(plot(q[1:end-1,4].%(2*pi),xlabel = "paso k",#=legend=:bottomleft,=#title="Disco vertical: ángulo interno θ",label="")) 
sf && savefig(sub*"Euler disk/figures/vertical_inner-angle_$h-$x0-$y0-$ϕ0.svg")
display(plot(H[2:end-1],xlabel = "paso k",#=legend=:bottomleft,=#title="Disco vertical: energía",label="")) 
sf && savefig(sub*"Euler disk/figures/vertical_energy_$h-$x0-$y0-$ϕ0.svg")
display(plot(Δq[2:end-1,:]./2,xlabel = "paso k",#=legend=:bottomleft,=#title="Disco vertical: velocidades",label=["v_x" "v_y" "v_φ" "v_θ"])) 
sf && savefig(sub*"Euler disk/figures/vertical_velocities_$h-$x0-$y0-$ϕ0.svg")


xs=Array{Float64,1}(undef,epochs-1)
ys=Array{Float64,1}(undef,epochs-1)

xmin=Float64[]; xmax=Float64[]; ymin=Float64[]; ymax=Float64[]
xlim=Tuple[]; ylim=Tuple[]

xs=q[:,1]
ys=q[:,2]

xmin = minimum(xs); xmax = maximum(xs)
ymin = minimum(ys); ymax = maximum(ys)
xlim = tuple(([xmin xmax] + (xmax-xmin)/20*[-1 1])...)
ylim = tuple(([ymin ymax] + (ymax-ymin)/20*[-1 1])...)


plot(title="Vertical Euler Disk: center",xlabel = "paso k",legend=:topleft, size = (1024, 576))
anim= @animate for i in 100:40:epochs-2
    plot(title="Vertical Euler Disk: center",xlabel = "paso k",legend=:topleft, size = (1024, 576))
    plot!(xs[i-70:i],ys[i-70:i],linewidth=2,color=1,xlim=xlim,ylim=ylim,label="") #trajectory center

    i==epochs-1 && sf && savefig(sub*"Euler Disk/figures/trajectories.svg")
end
gif(anim,sub*"Euler disk/figures/disk_vertical.mp4",fps=25)