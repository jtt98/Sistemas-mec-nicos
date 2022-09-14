#=
Disco de Euler usando integradores simplécticos por Jose Torrente Teruel
Constantes: a elección
=#

using Plots, Calculus, NLsolve

sf = true # save figures
sub="D:/Escritorio/U/Master/Practicum/"

function Constraints(p0::Vector{<:Real},a0::Vector{<:Real},a1::Vector{<:Real})
    global r, h
    
    ϕm0=(a1[2]+a0[2])/2;
    ψd0=(a1[3]-a0[3])/h;
    
    return p0-h*r*ψd0*[cos(ϕm0), sin(ϕm0)];
end

function E(q::Vector{<:Real},v::Vector{<:Real})
    global M, g, r, I, J
    θ=q[1]
    ϕ = q[2]
    θd = v[1];
    ϕd = v[2];
    ψd = v[3];
    u = v[4];
    w = v[5];
    
    P = M*g*r*sin(θ);
    K = 0.5*( I*(θd^2+ϕd^2*(sin(θ))^2) + J*(ψd+ϕd*cos(θ))^2 + M*(u^2+w^2) + M*r^2*(θd^2+ϕd^2*(cos(θ))^2) ) + M*r*( θd*sin(θ)*(u*sin(ϕ)-w*cos(ϕ)) - ϕd*cos(θ)*(u*cos(ϕ)+w*sin(ϕ)) );
    
    return P + K
end

function Jacobian!(Jac::Matrix{<:Real},a1::Vector{<:Real},a2::Vector{<:Real})
    global M, g, r, I, J, h

    #local Jac = Array{Real,2}(undef,3,3)
    
    am1=(a2+a1)/2;
    ad1=(a2-a1)/h;
    
    ϕ1=a1[2]
    θm1=am1[1]
    ϕm1=am1[2]
    θd1=ad1[1]
    θd1=ad1[1]
    ϕd1=ad1[2]
    ψd1=ad1[3];
    
    Jac[1,1] = ( 4*(M*r^2+I)/h^2 - g*M*r*sin(θm1) + (M*r^2-I+J)*cos(2*θm1)*ϕd1^2 + (M*r^2+J)*cos(θm1)*ϕd1*ψd1 )/2;
    Jac[1,2] = ( (M*r^2-I+J)*sin(2*θm1)*ϕd1 + (M*r^2+J)*sin(θm1)*ψd1 )/h;
    Jac[1,3] = (M*r^2+J)*sin(θm1)*ϕd1/h;
    Jac[2,1] = (M*r^2-I+J)*sin(2*θm1)*ϕd1 + ( J*sin(θm1) + h*M*r^2/2*sin(θm1)*θd1 )*ψd1;
    Jac[2,2] = -( (M*r^2+I+J) + (M*r^2-I+J)*cos(2*θm1) )/h;
    Jac[2,3] = -2*(M*r^2+J)/h*cos(θm1) - M*r^2*sin(θm1)*θd1;
    Jac[3,1] = (J+M*r^2*cos(ϕm1-ϕ1)) * sin(θm1)*ϕd1/2 + M*r^2*sin(ϕm1-ϕ1) * (cos(θm1)*θd1/2+sin(θm1)/h);
    Jac[3,2] = M*r^2*( cos(ϕm1-ϕ1)*sin(θm1)*θd1 + sin(ϕm1-ϕ1)*(cos(θm1)*ϕd1+ψd1) )/2 - (J+M*r^2*cos(ϕm1-ϕ1))*cos(θm1)/h;
    Jac[3,3] = -(J+M*r^2*cos(ϕm1-ϕ1))/h;
end

function EulerLagrange!(EL::Vector{<:Real},a0::Vector{<:Real},a1::Vector{<:Real},a2::Vector{<:Real})
    global M, g, r, I, J, h
    
    am0=(a1+a0)/2;
    am1=(a2+a1)/2;
    ad0=(a1-a0)/h;
    ad1=(a2-a1)/h;
    
    ϕ1=a1[2];
    θm0=am0[1]
    ϕm0=am0[2]
    θm1=am1[1]
    ϕm1=am1[2]
    θd0=ad0[1]
    ϕd0=ad0[2]
    ψd0=ad0[3]
    θd1=ad1[1]
    ϕd1=ad1[2]
    ψd1=ad1[3]

    #Jac=zeros(3,3)
    
    # omitting factor (-1/(2h))
    EL[1] = 2*(M*r^2+I)*(θd1-θd0)/h + g*M*r*(cos(θm0) + cos(θm1)) + (M*r^2-I+J)/2*(sin(2*θm0)*ϕd0^2 + sin(2*θm1)*ϕd1^2) + (M*r^2+J)*(sin(θm0)*ϕd0*ψd0 + sin(θm1)*ϕd1*ψd1);
    # omitting factor (1/(2h))
    EL[2] = (M*r^2-I+J)*(cos(2*θm0)*ϕd0 - cos(2*θm1)*ϕd1) - (M*r^2+I+J)*(ϕd1-ϕd0) + 2*(M*r^2+J)*(cos(θm0)*ψd0 - cos(θm1)*ψd1) - h*M*r^2*(sin(θm0)*θd0*ψd0 + sin(θm1)*θd1*ψd1);
    # omitting factor (1/h)
    EL[3] = J*(cos(θm0)*ϕd0 - cos(θm1)*ϕd1 + ψd0-ψd1) + M*r^2*( sin(θm1)*sin(ϕm1-ϕ1)*θd1 - sin(θm0)*sin(ϕm0-ϕ1)*θd0 - cos(ϕm1-ϕ1)*(cos(θm1)*ϕd1+ψd1) + cos(ϕm0-ϕ1)*(cos(θm0)*ϕd0+ψd0) );
    #Jac = Jacobian(a1,a2);
end


# Inicializacion
global M, g, r, I, J, T, N, H

M = 1;   #Masa
g = 9.8; #Gravedad
r = 1;   #Radio
I = 1;   #Momento de inercia I_1=1_2
J = 2.5;   #Momento de inercia I_3!=I_1
#e=0

epochs = 2^14; #Numero de pasos
h = 0.001; #Tamaño del paso
q=Array{Float64,2}(undef,epochs,5)   
Δq=Array{Float64,2}(undef,epochs,5)
H=Array{Float64}(undef,epochs)

# q=[θ ϕ ψ x y]
# θ: inclinacion
# ϕ: giro
# ψ: rodamiento
# [x y]: posicion

#Condiciones iniciales
ejemplo = 3
if ejemplo==0
    #print("Iniciando simulacion en punto fijo y giro vertical\n")
    q[1,:]=[pi/2 0 0 0 0];
    q[2,:]=[pi/2 h 0 0 0];
elseif ejemplo==1
    #print("Iniciando simulacion en linea recta\n")
    q[1,:]=[pi/2 0 0    0 0];
    q[2,:]=[pi/2 0 h -h*r 0];
elseif ejemplo==2
    #print("Iniciando simulacion rara\n")
    q[1,:]=[pi/2   0 0 0 0];
    q[2,:]=[pi/2+h h h h h];
elseif ejemplo==3
    #print("Iniciando simulacion rara, pero vertical\n")
    q[1,:]=[pi/2 0 0 0 0];
    q[2,:]=[pi/2 h h h h];
end
a1=q[1,1:3];
a2=q[2,1:3];
p1=q[1,4:5];
p2=Constraints(p1,a1,a2);
if sum(q[2,4:5]≈p2)>0
    print("ATENCION: restriccion no satisfecha. Corrigiendo...")
    q[2,4:5]=p2;
end

# Iteraciones
for i=1:epochs-2
    local a0=q[i,1:3];
    local a1=q[i+1,1:3];
    local a2=nlsolve((EL,a2)->EulerLagrange!(EL,a0,a1,a2), (Jac,a2)->Jacobian!(Jac,a1,a2),a1#=,ftol=1e-8=#).zero;
    local p1=q[i+1,4:5];
    local p2=Constraints(p1,a1,a2);
    q[i+2,1:3]=a2
    q[i+2,4:5]=p2
end
qm = (q[2:epochs,:]+q[1:epochs-1,:])/2;
Δq = (q[2:epochs,:]-q[1:epochs-1,:])/h;
e=zeros(epochs-1);
for i=1:epochs-1
    e[i]=E(qm[i,:],Δq[i,:]);
end

###Plots y gifs
display(plot(q[1:end-1,4],q[1:end-1,5],#=legend=:bottomleft,=#title="Disco de Euler: trayectoria",label="")) 
sf && savefig(sub*"Euler disk/figures/euler_trajectory_$h.svg")
display(plot(e[2:end-1].-minimum(e[2:end-1]),xlabel="paso k",#=legend=:bottomleft,=#title="Disco de Euler: energía",label="")) 
sf && savefig(sub*"Euler disk/figures/euler_energy_$h.svg")
