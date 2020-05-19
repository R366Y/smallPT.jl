module smallPT
import Base: +, -, *, % 

const M_PI = 3.1415926535
const M_1_PI = 1.0 / M_PI
const inf = 1e20
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
    epsi = 1e-04
    b = dot(op, r.d)
    det = b * b - dot(op, op) + s.rad * s.rad
    if det < 0.
        return inf
    else
        det = sqrt(det)
    end
    t1 = b - det
    t2 = b + det
    if t1 > epsi 
        return t1
    elseif t2 > epsi 
        return t2
    end
    return inf
end

spheres = Sphere[
    Sphere(1e5, Vec(1e5+1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),    # Left
    Sphere(1e5, Vec(-1e5+99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF),  # Right
    Sphere(1e5, Vec(50., 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),       # Back
    Sphere(1e5, Vec(50., 40.8, -1e5+170), Vec(), Vec(), DIFF),               # Front
    Sphere(1e5, Vec(50., 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),       # Bottom
    Sphere(1e5, Vec(50., -1e5+81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF), # Top
    Sphere(16.5, Vec(27., 16.5, 47.), Vec(), Vec(1., 1., 1.) *.999, SPEC),   # Mirror
    Sphere(16.5, Vec(73., 16.5, 78.), Vec(), Vec(1., 1., 1.)*.999, REFR),    # Glass
    Sphere(1.5, Vec(50., 81.6-16.5, 81.6), Vec(4.,4.,4.)*100, Vec(), DIFF)   # Lite
]

# Converts floats to integers to be saved in PPM file
@inline toInt(x::Float64) = return floor(Int, clip(x)^(1. / 2.2) * 255 + .5)

const num_spheres = length(spheres)

function intersect(r::Ray)
    t = 1e20
    id = nothing
    for sphere in spheres
        d = intersect(sphere, r)
        if d < t 
            t = d
            id = sphere
        end
    end
    return id, t
end

function main()
    w = 512
    h = 385
    samps = 64
    cam = Ray(Vec(50., 52., 295.6), norm(Vec(0., -0.042612, -1.))) # camera position, direction
    fov = .5135 # Field of view angle
    cx = Vec(w * fov/h, 0., 0.) # Horizontal (x) camera direction increment
    cy = norm(cx % cam.d) * fov # Vector up for the camera , y direction increment
    c = Array{Vec}(undef, w * h) # The image
    for i = 1:length(c)
        c[i] = Vec(0., 0., 0.)
    end

end

end 
