#=
Péndulo simple usando integradores simplécticos por Jose Torrente Teruel
Constantes: todas igual a la unidad
=#
using Plots, Calculus

sf = true # save figures
sub="D:/Escritorio/U/Master/Practicum/Pendulum/figures/" #choose a path to save the figures

function V(q::Real)
    return cos(q)
end

function dV(q::Real) #derivative of V
    return derivative(V)(q)
end

function E_d(q1::Real,q2::Real,f::Function)
    return 0.5*((q2-q1)/h)^2 - 0.5*(f(q1)+f(q2));
end

function g(q::Real,p::Real,f::Function) 
    return q+h^2*(p+0.5*f(q))
end

#Inicialización
epochs = 2^12
h=0.001

q=Array{Float64,2}(undef,epochs,2)
p=Vector{Float64}(undef,epochs)
v=Vector{Float64}(undef,epochs)
H=Array{Float64,2}(undef,epochs-1,2)

#Condiciones iniciales
q0=0.2
q[1,1:2]=[q0,q0]
v0=5e-4
v[1]=v0
p[1]=v0
#p[1]=(q[2,1]-q[1,1])/h^2 - 0.5*dV(q[1,1])

#Iteraciones
for k=1:epochs-1
    #symplectic method
    q[k+1,1]=g(q[k,1],p[k],dV);
    p[k+1]=p[k]+0.5*(dV(q[k,1])+dV(q[k+1,1]));
    H[k,1]=E_d(q[k,1],q[k+1,1],V);

    #Euler method in the Hamiltonian side
    q[k+1,2]=q[k,2]+h*v[k];
    v[k+1]=v[k]+0.5*h*(dV(q[k,2])+dV(q[k+1,2]));
    H[k,2]=E_d(q[k,2],q[k+1,2],V);
end


###Plots and gifs
#display(plot([q[1:end-1,1],H[1:end-1,1]],label=["Symplectic" "Energy"],#=legend=:bottomleft,=#title="Pendulum (h=$h)"))
#sf && savefig(sub*"simple_symplectic_energy_$h-$q0-$v0.svg")
#display(plot([q[1:end-1,2],H[1:end-1,2]],label=["Euler" "Energy"],#=legend=:bottomleft,=#title="Pendulum (h=$h)"))
#sf && savefig(sub*"simple_euler_energy_$h-$q0-$v0.svg")
display(plot(q[:,1:2],label=["Simpléctico" "Euler"],xlabel="paso k",legend=:topleft,title="Péndulo (h=$h): trayectorias")) 
sf && savefig(sub*"simple_trajectories_$h-$q0-$v0.svg")
display(plot(H[:,1:2],label=["Simpléctico" "Euler"],xlabel="paso k",legend=:topleft,title="Péndulo (h=$h): energía")) 
sf && savefig(sub*"simple_energy_$h-$q0-$v0.svg")



xs=Array{Float64,2}(undef,epochs,2)
ys=Array{Float64,2}(undef,epochs,2)
xmin=Float64[]; xmax=Float64[]; ymin=Float64[]; ymax=Float64[]
xlim=Tuple[]; ylim=Tuple[]

xs=sin.(q)
ys=-cos.(q)

xmin = minimum(xs); xmax = maximum(xs)
ymin = minimum(ys); ymax = maximum(ys)
xlim = tuple(([xmin xmax] + (xmax-xmin)/20*[-1 1])...)
ylim = tuple(([ymin ymax] + (ymax-ymin)/20*[-1 1])...)

plt=plot(title="Pendulum (h=$h): Trajectories",legend=:topleft,xlim=xlim,ylim=ylim, size = (1024, 576))
anim= @animate for i in 1:10:epochs
    plot(xs[1:i,1],ys[1:i,1],linewidth=2,color=1,label="") #trayectoria 1
    plot!(xs[1:i,2],ys[1:i,2],linewidth=0.8,color=2,label="") #trayectoria 2
    plot!([0,xs[i,1]],[0,ys[i,1]],color=1,label="Symplectic") #barra 1
    plot!([0,xs[i,2]],[0,ys[i,2]],color=2,label="Euler")    #barra 2
    scatter!(xs[i,:],ys[i,:],color=3, label="")     #peso de los péndulos
end
gif(anim,sub*"pendulum_simple_$h-$q0-$v0.mp4",fps=10)