#=
Doble péndulo usando integradores simplécticos por Jose Torrente Teruel
Constantes: g=1, m1=55/24, m2=3/8, L1=1, L2=8/3
=#
using Plots

sf = false # save figures
sub="D:/Escritorio/U/Master/Practicum/Pendulum/figures/"

# q = (θ1, θ2)

function E_d(q::Vector{<:Real},v::Vector{<:Real})
    T=1/3*(v[2]^2 + 4*v[1]^2 + 3*v[1]*v[2]*cos(q[1]-q[2]))
    V=-3*cos(q[1])-cos(q[2])
    return T+V
end

#Inicialización
epochs = 2^12;
#h=0.100462;
h=0.01

q=Array{Float64,2}(undef,epochs,2)
μ=Array{Float64,2}(undef,epochs,2)
v=Array{Float64,2}(undef,epochs,2)
H=Array{Float64,1}(undef,epochs)

#Condiciones iniciales
q0=0.1
q1=0.108
q[1:2,:]=[q0 q0
          q1 q1];
v[1,:]=(q[2,:]-q[1,:])/h
μ[1,:]=[8/3 * v[1,1]/h + v[1,2]/h * cos(q[1,1]-q[1,2]) + v[1,1]*v[1,2]*sin(q[1,1]-q[1,2]) + 3/2 * sin(q[1,1]),
        2/3 * v[1,2]/h + v[1,1]/h * cos(q[1,1]-q[1,2]) - v[1,1]*v[1,2]*sin(q[1,1]-q[1,2]) + 1/2 * sin(q[1,2])]
μ[2,:]=[8/3 * v[1,1]/h + v[1,2]/h * cos(q[1,1]-q[1,2]) - 3/2 * sin(q[2,1]),
        2/3 * v[1,2]/h + v[1,1]/h * cos(q[1,1]-q[1,2]) - 1/2 * sin(q[2,2])]
        
H[1]=E_d(0.5*(q[1,:]+q[2,:]),v[1,:])

#Iteraciones
for k=2:epochs-1
    local a1,a2,a,b1,b2,b,c,d1,d2,d,A,B,C,M,m
    c=sin(q[k,1]-q[k,2])
    a1=8/(3*h)
    a2=1/h *cos(q[k,1]-q[k,2])
    b1=a2
    b2=2/(3*h)
    d1=μ[k,1]-3/2 * sin(q[k,1])
    d2=μ[k,2]-1/2 * sin(q[k,2])

    if c==0.0
        v[k,1]=(b1*d2-b2*d1)/(a2*b1-b2*a1)
        v[k,2]=(d1 - a1*v[k,1])/(b1+c*v[k,1])
    else
        A=c*(a1+a2)
        B=a2*b1-b2*a1-c*(d1+d2)
        C=b2*d1-b1*d2
        v[k,1]=(-B-sqrt(B^2-4*A*C))/(2*A)
        v[k,2]=(d1 - a1*v[k,1])/(b1+c*v[k,1])
    end

    q[k+1,1]=q[k,1]+h*v[k,1]
    q[k+1,2]=q[k,2]+h*v[k,2]
    
    μ[k+1,:]=[8/3 * v[k,1]/h + v[k,2]/h * cos(q[k,1]-q[k,2]) - 3/2 * sin(q[k+1,1]),
              2/3 * v[k,2]/h + v[k,1]/h * cos(q[k,1]-q[k,2]) - 1/2 * sin(q[k+1,2])]

    H[k]=E_d(0.5*(q[k+1,:]+q[k,:]),v[k,:])
end

###Plots and gifs
display(plot(q[2:end,:].%(2*pi),label=["θ1" "θ2"], xlabel = "paso k",#=legend=:bottomleft,=#title="Péndulo doble: ángulos")) 
sf && savefig(sub*"double_angles_$h-$q0-$q1.svg")
display(plot(H[1:end-1],label="Energía", xlabel = "paso k",#=legend=:bottomleft,=#title="Péndulo doble: energía")) 
sf && savefig(sub*"double_energy_$h-$q0-$q1.svg")

xs=Array{Float64,2}(undef,epochs,2)
ys=Array{Float64,2}(undef,epochs,2)

xmin=Float64[]; xmax=Float64[]; ymin=Float64[]; ymax=Float64[]
xlim=Tuple[]; ylim=Tuple[]

xs[:,2]=sin.(q[:,2])
ys[:,2]=-cos.(q[:,2])
xs[:,1]=xs[:,2]+2.6*sin.(q[:,1])
ys[:,1]=ys[:,2]-2.6*cos.(q[:,1])

xmin = minimum(xs); xmax = maximum(xs)
ymin = minimum(ys); ymax = maximum(ys)
xlim = tuple(([xmin xmax] + (xmax-xmin)/20*[-1 1])...)
ylim = tuple(([ymin ymax] + (ymax-ymin)/20*[-1 1])...)

plot([0],[0],xlim=xlim,ylim=ylim)
anim= @animate for i in 1:20:epochs
    plot(title="Péndulo doble: trayectorias",legend=:topleft, xlabel = "paso k", size = (1024, 576))
    plot([0],[0],xlim=xlim,ylim=ylim,label="")
    plot!(xs[1:i,1],ys[1:i,1],linewidth=2,color=1,label="") #trajectory point 1
    plot!(xs[1:i,2],ys[1:i,2],linewidth=0.8,color=2,label="") #trajectory point 2
    plot!([0,xs[i,2]],[0,ys[i,2]],color=2,label="")    #rope 2
    plot!([xs[i,2],xs[i,1]],[ys[i,2],ys[i,1]],color=1,label="") #rope 1
    scatter!(xs[i,:],ys[i,:],color=3, label="")     #weight of the pendulums

    i==epochs-1 && sf && savefig("trajectories.svg")
end
gif(anim,sub*"pendulum_double_$h-$q0-$q1.mp4",fps=25)
