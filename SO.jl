using LinearAlgebra

I=[1 0 0
    0 1 0
    0 0 1] #identity
#nskew=0 #0 to check all the matrices are skew-symmetric
nskew=1


function rotations(k::Integer,θ::Real)
    x=Array{Float64,2}(undef,(3,3))
    if k==1
        x=[1 0 0
            0 cos(θ) -sin(θ)
            0 sin(θ) cos(θ)]
    elseif k==2
        x=[cos(θ) 0 -sin(θ)
            0 1 0
            sin(θ) 0 cos(θ)]
    elseif k==3
        x=[cos(θ) -sin(θ) 0
            sin(θ) cos(θ) 0
            0 0 1]
    end
    return x
end

normrot(x::Matrix{<:Real}) = √max(0,min(4,3-tr(x)))
distrot(x::Matrix{<:Real},y::Matrix{<:Real}) = normrot(transpose(x)*y)
acosrot(x::Matrix{<:Real}) = acos(max(-1,min(1,(tr(x)-1)/2)))

function isskew(A::Matrix{<:Real})
    if A' != -A
        return false
    else 
        return true
    end
end

function skew(A::Matrix{<:Real})
    return 0.5*(A-A')
end

function hat(x::Vector{<:Real})
    aux=[0 -x[3] x[2]
        x[3] 0 -x[1]
        -x[2] x[1] 0]
    return aux
end
function vee(A::Matrix{<:Real})
    if nskew==0 && isskew(A)==false 
        return error("Not skew-symmetric")
    end
    return [A[3,2],A[1,3],A[2,1]]
end

function norm2(x::Vector{<:Real})
    return sum(x.^2)
end
function norm2(A::Matrix{<:Real})
    return norm2(vee(A))
end

#retractions
function cay(A::Matrix{<:Real})
    if nskew==0 && !isskew(A)
        return error("Not skew-symmetric")
    end 
    return I+ (2/(1+norm2(A))) * (A + A^2)
end
function cay(x::Vector{<:Real})
    return cay(hat(x))
end
function cay1(A::Matrix{<:Real})
    return 1/(1+tr(A))*(A-transpose(A))
end

function rexp(A::Matrix{<:Real})
    r=sqrt(norm2(A))
    return I+sin(r)/r * A +2*sin(0.5*r)^2/r^2 * A^2
end
function rlog(R::Matrix{<:Real})
    z = R-transpose(R)/2
    r = acosrot(R)
    return r == 0 ? zeros(3,3) : r/norm(vee(z))*z
end

function τ(x::Vector{<:Real})
    return cay(hat(x))
end
function τ1(A::Matrix{<:Real})
    return vee(cay1(A))
end

function ω(x::Vector{<:Real})
    return rexp(hat(x))
end
function ω1(A::Matrix{<:Real})
    return vee(rlog(A))
end