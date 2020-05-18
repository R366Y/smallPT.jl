module smallPT
import Base: +, -, *, % 

M_PI = 3.1415926535
M_1_PI = 1.0 / M_PI
@enum Refl_t DIFF SPEC REFR

struct Vec
    x::Float64
    y::Float64
    z::Float64

    Vec(x = 0., y = 0. , z = 0.) = new(x, y, z)
end

+(a::Vec, b::Vec) = return Vec(a.x + b.x, a.y + b.y, a.z + b.z)
-(a::Vec, b::Vec) = return Vec(a.x - b.x, a.y - b.y, a.z - b.z)
*(a::Vec, b) =  return Vec(a.x * b, a.y * b, a.z * b)
mult(a::Vec, b::Vec) = return Vec(a.x * b.x, a.y * b.y, a.z * b.z)
norm(a::Vec) = return a * (1 / sqrt(a.x^2 + a.y^2 + a.z^2))
dot(a::Vec, b::Vec) = return a.x * b.x + a.y * b.y + a.z * b.z 
%(a::Vec, b::Vec) = return Vec(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)

struct Ray
    o::Vec
    d::Vec
end

struct Sphere
    rad::Float64
    p::Vec  # position
    e::Vec  # emission
    c::Vec  # color
    refl    # reflection type: DIFFuse, SPECular , REFRactive
end

function intersect(s::Sphere, r::Ray)
    op = s.p - r.o
    t = epsi = 1e-04
    b = dot(op, r.d)
    det = b * b - dot(op, op) + s.rad * s.rad
    if det < 0
        return 0.
    else
        det = sqrt(det)
    end
    t1 = b - det
    t2 = b + det
    if t > epsi 
        return t1
    elseif t2 > epsi 
        return t2
    end
    return 0
end

spheres = Sphere[
    Sphere(1e5, Vec(1e5+1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF), # Left
    Sphere(1e5, Vec(-1e5+99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF) # Right
]

end 
