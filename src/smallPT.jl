import Base: +, -, *, % 
import Base.Threads: @threads

const inf = 1e20
const (DIFF, SPEC, REFR) = (1, 2, 3)

struct Vec
    x::Float64
    y::Float64
    z::Float64

    Vec(x = 0., y = 0. , z = 0.) = new(x, y, z)
end

const ZERO  = Vec(0., 0., 0.)
const XAXIS = Vec(1., 0., 0.)
const YAXIS = Vec(0., 1., 0.)

# vector operations 

+(a::Vec, b::Vec) = return Vec(a.x + b.x, a.y + b.y, a.z + b.z)
-(a::Vec, b::Vec) = return Vec(a.x - b.x, a.y - b.y, a.z - b.z)
*(a::Vec, b::Float64) =  return Vec(a.x * b, a.y * b, a.z * b)
*(a::Float64, b::Vec) = b * a
*(a::Vec, b::Vec) = return Vec(a.x * b.x, a.y * b.y, a.z * b.z)
norm(a::Vec) = return a * (1 / sqrt(a.x^2 + a.y^2 + a.z^2))
dot(a::Vec, b::Vec) = return a.x * b.x + a.y * b.y + a.z * b.z                                      # dot product 
%(a::Vec, b::Vec) = return Vec(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x) # cross product

struct Ray
    o::Vec
    d::Vec
end

struct Sphere
    rad::Float64
    p::Vec  # position
    e::Vec  # emission
    c::Vec  # color
    refl::Int   # reflection type: DIFFuse, SPECular , REFRactive
    maxRefl::Float64
    sqRad::Float64
    invC::Vec

    function Sphere(rad, p, e, c, refl)
        maxRefl = c.x > c.y && c.x > c.z ? c.x : c.y > c.z ? c.y : c.z    # max reflectivity
        sqRad = rad * rad
        new(rad, p, e, c, refl, maxRefl, sqRad, c*(1. / maxRefl) )
    end
end

function intersect(s::Sphere, r::Ray)
    op = s.p - r.o
    epsi = 1e-04
    b = dot(op, r.d)
    det = b * b - dot(op, op) + s.sqRad
    if det < 0.
        return inf
    else
        det = sqrt(det)
    end
    t1 = b - det
    t2 = b + det
    if  t1 > epsi 
        return t1
    elseif t2 > epsi 
        return t2
    end
    return inf
end


const left = Sphere(1e5, Vec(1e5+1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF)       # Left
const right = Sphere(1e5, Vec(-1e5+99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF)    # Right
const back =  Sphere(1e5, Vec(50., 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF)         # Back
const front = Sphere(1e5, Vec(50., 40.8, -1e5+170), Vec(), Vec(), DIFF)                 # Front
const bottom = Sphere(1e5, Vec(50., 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF)        # Bottom
const top = Sphere(1e5, Vec(50., -1e5+81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF)     # Top
const mirror = Sphere(16.5, Vec(27., 16.5, 47.), Vec(), Vec(1., 1., 1.) *.999, SPEC)    # Mirror
const glass = Sphere(16.5, Vec(73., 16.5, 78.), Vec(), Vec(1., 1., 1.)*.999, REFR)      # Glass
const lite = Sphere(600, Vec(50., 681.6 - .27, 81.6), Vec(12.,12.,12.), Vec(), DIFF)    # Lite

const spheres = [left, right, back, front, bottom, top, mirror, glass, lite]

# Converts floats to integers to be saved in PPM file
cl(x::Float64) = return max(min(x,one(x)),zero(x))
toInt(x::Float64) = return floor(Int, cl(x)^(1. / 2.2) * 255 + .5)

function intersectSpheres(r::Ray)
    t = inf
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

function radiance(r::Ray, depth::Int)
    obj,t = intersectSpheres(r)
    if isnothing(obj) 
        return ZERO
    end

    # Russian roulette: stop the recursion randomly based on the surface reflectivity

    depth+=1
    isMaxDepth = depth > 100
    useRR = depth > 5
    isRR = useRR && rand() < obj.maxRefl

    if isMaxDepth || (useRR && ~isRR)
        return obj.e
    end

    x = r.o + r.d * t                               # ray intersection position
    n = norm(x-obj.p)                               # sphere normal
    nl = dot(n, r.d) < 0. ? n : n * -1.             # properly oriented surface normal
    f = (useRR && isRR) ? obj.invC : obj.c    # object color (BRDF modulator)

    if obj.refl == DIFF
        r1 = 2 * pi * rand()
        r2 = rand()
        r2s = sqrt(r2)
        w = nl
        wo = w.x < -0.1 || w.x > 0.1 ? YAXIS : XAXIS
        u = norm(wo % w)
        v = w % u
        d = norm(u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1-r2))
        return obj.e + f * radiance(Ray(x, d), depth)
    elseif obj.refl == SPEC
        return obj.e + f * radiance(Ray(x, r.d - n * (2. * dot(n ,r.d))), depth)
    else # otherwise we have a dielectric (glass) surface
        reflRay = Ray(x, r.d - n * (2. * dot(n ,r.d)))  # ideal dielectric reflection
        into = dot(n, nl) > 0.                          # Ray from outside going in?
        nc = 1.
        nt = 1.5
        nnt = into ? nc/nt : nt / nc
        ddn = dot(r.d, nl)
        cos2t = 1 - nnt * nnt * (1 - ddn * ddn)
        
        # if total internal reflecttion, reflect
        if cos2t < 0
            return obj.e + f * radiance(reflRay, depth)
        else
            # otherwise, choose reflection or refraction
            tdir = norm(r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))))
            a = nt - nc
            b = nt + nc
            R0 = (a * a) / (b * b)
            c = 1. - (into ? -ddn : dot(tdir,n))
            Re = R0 + (1 - R0) * c * c * c * c * c
            Tr = 1 - Re
            P = .25 + .5 * Re
            RP = Re / P
            TP = Tr / (1. - P)

            result = ZERO
            if depth > 2
                if rand() < P
                    result = radiance(reflRay, depth) * RP
                else
                    result = radiance(Ray(x, tdir), depth) * TP
                end
            else
                result = radiance(reflRay, depth) * Re + radiance(Ray(x, tdir), depth) * Tr
            end

            return obj.e + f * result
        end
    end
end

function main()
    w = 1024
    h = 768
    samps = 16
    cam = Ray(Vec(50., 52., 295.6), norm(Vec(0., -0.042612, -1.)))  # camera position, direction
    fov = .5135                                                     # Field of view angle
    cx = Vec(w * fov/h, 0., 0.)                                     # Horizontal (x) camera direction increment
    cy = norm(cx % cam.d) * fov                                     # Vector up for the camera , y direction increment
    c = Array{Vec}(undef, w * h)                                    # The image
    for i = 1:length(c)
        c[i] = Vec(0., 0., 0.)
    end

    @threads for y in 1:h                                                    # Loop over image rows
        #print(stderr, "\rRendering ($(samps*4)) $(100. * y/(h-1))%") # Print progress
        for x in 1:w                                                # Loop columns 

            # For each pixel do 2x2 subsamples and "samps" samples per subsample
            for sy in 1:2
                i = (h - y) * w + x                                 # Calculate array index for pixel(x,y)
                for sx in 1:2
                    r = Vec()
                    for s in 1:samps
                        # Tent filter
                        r1 = 2 * rand()
                        r2 = 2 * rand()
                        dx = r1 < 1. ? sqrt(r1) - 1. : 1. - sqrt(2. - r1)
                        dy = r2 < 1. ? sqrt(r2) - 1. : 1. - sqrt(2. - r2)
                        # Compute ray direction
                        d = cx * (((sx-1 + .5 + dx) / 2 + x) / w - .5) +
                            cy * (((sy-1 + .5 + dy) / 2 + y) / h - .5) + cam.d
                        r = r + radiance(Ray(cam.o + d * 140., norm(d)), 0) * (1. / samps)
                    end
                    @inbounds c[i] = c[i] + Vec(cl(r.x), cl(r.y), cl(r.z)) * .25
                end
            end 

        end
    end

    open("image_julia.ppm", "w") do f
        print(f, "P3\n$w $h\n255\n")
        for i in 1:length(c)
            print(f, toInt(c[i].x))
            print(f, " ")
            print(f, toInt(c[i].y))
            print(f, " ")
            print(f, toInt(c[i].z))
            print(f, "\n")
        end
    end
end


#@time main()