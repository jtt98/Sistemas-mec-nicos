#=
Bola rodando en un plano que rota usando integradores simplécticos por Jose Torrente Teruel
Constantes: a elección
=#

using Plots, Calculus, NLsolve

include("SO.jl") #archivo auxiliar

sf = false # save figures
sub="D:/Escritorio/U/Master/Practicum/"

function proj(u::Vector{<:Real}) #projección de R^3 en R^2
    return [u[1], u[2]]
end

function L(v::Vector{<:Real},w::Vector{<:Real})
    return 0.5*m*sum(v.^2) + 0.5*J*sum(w.^2)
end

function fc(X::Vector{<:Real}, F::Matrix{<:Real})
    Y=X + a*(F-I)*γ
    return [Y[1],Y[2],0]
    #return X + a*(F-I)*γ
end

function fnl(x::Vector{<:Real}, α::Vector{<:Real})
    r2=sum(x.^2)
    aux=(x[3]^2-r2)/(1+r2)
    Ax=4+10*A-2*h*Ω*x[3]*(1+5*aux)-dot(x,α)
    Aγ=2*h*Ω*(1-5*aux)-10*A*x[3]
    return Ax*x + Aγ*γ-10*h*Ω*(1+aux)*cross(γ,x)-α-cross(α,x)
end

function fnl!(F::Vector{<:Real},x::Vector{<:Real},α::Vector{<:Real})
    F[:]=fnl(x,α)
end

#Inicializacion
epochs = 2^12;
global m=1; a=1; J=2/5; Ω=0.; γ=[0, 0, 1];
h=0.01
global A=-1 +0.25*h^2*Ω^2; B= 1+ 0.75*h^2*Ω^2

F=Array{Float64,3}(undef,epochs-1,3,3)   #q=[X=(x, y), R], F_k=R_k^t * R_{k+1}
w=Array{Float64,2}(undef,epochs,3)
X=Array{Float64,2}(undef,epochs,3)    #tratamos a X en R^3 y luego proyectamos
Y=Array{Float64,2}(undef,epochs,2) 
H=Array{Float64}(undef,epochs-1)

#Condiciones iniciales
θ1=0.001; n1=2
F[1,:,:]=rotations(n1,θ1)
w[1,:]=1/h*vee(F[1,:,:]-I)
#w[1,:]=[0.1,0.1,0]
#F[1,:,:]=I+h*hat(w[1,:])
X[1,:]=[0,1,0]
Y[1,:]=proj(X[1,:])
X[2,:]=fc(X[1,:],F[1,:,:])
Y[2,:]=proj(X[2,:])

v=1/h*(X[2,:]-X[1,:])
X12=0.5*(X[1,:]+X[2,:])
H[1]=L(v+Ω*cross(γ,X12),w[1,:]+Ω*γ) #energía es igual a energía cinética
x0=τ1(F[1,:,:])

#Iteraciones
for k=1:epochs-2
    local α=vee((I+h*Ω*hat(γ))*F[k,:,:]-F[k,:,:]'*(I-h*Ω*hat(γ))) - 5/a * (cross(γ,B*X[k+1,:]+A*X[k,:]) + h*Ω*(X[k+1,:]-X[k,:]))
    global x0=nlsolve((F,x)->fnl!(F,x,α), x0#=, ftol=1e-12=#).zero
    F[k+1,:,:]=cay(x0)

    X[k+2,:]=fc(X[k+1,:],F[k+1,:,:])
    #Y[k+2,:]=proj(X[k+2,:])
    w[k+1,:]=1/h*vee(F[k+1,:,:]-I)
    local v=1/h*(X[k+2,:]-X[k+1,:])
    local X12=0.5*(X[k+1,:]+X[k+2,:]) 
    H[k+1]=L(v+Ω*cross(γ,X12),w[k+1,:]+Ω*γ)
end

###Plots and gifs
display(plot(X[:,1],X[:,2],xlabel = "paso k",#=legend=:bottomleft,=#title="Bola rodando (h=$h,Ω=$Ω): punto de contacto",label="")) 
sf && savefig(sub*"Rolling Ball/ball_contact-point_$h-$Ω-$θ1.svg")

display(plot(H[2:end-1].-minimum(H[2:end-1]),xlabel = "paso k",#=legend=:bottomleft,=#title="Bola rodando (h=$h,Ω=$Ω): energía",label="")) 
sf && savefig(sub*"Rolling Ball/ball_energy_$h-$Ω-$θ1.svg")
