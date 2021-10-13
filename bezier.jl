module Bezier

import Base: +, -, *, \, Base.zero

struct Point
    x::Float64
    y::Float64
    z::Float64
end

*(p::Point, d::Number) = Point(p.x * d, p.y * d, p.z * d)
*(d::Number, p::Point) = p * d

+(p::Point, q::Point) = Point(p.x + q.x, p.y + q.y, p.z + q.z)
-(p::Point, q::Point) = p + (-1 * q)

zero(::Type{Point}) = Point(0.0, 0.0, 0.0)
zero(::Point) = Point(0.0, 0.0, 0.0)

by_coordinates(obj) = [broadcast(p -> getfield(p, f), obj) for f in fieldnames(Point)]

function \(A::Matrix{Float64}, b::Vector{Point})
    bs = by_coordinates(b)
    sols = map(bs) do bi A \ bi end
    map(Point, sols...)
end

function read_bzr(filename)
    local cnet
    open(filename) do f
        n, m = map(split(readline(f))) do s parse(Int, s) + 1 end
        cnet = Matrix{Point}(undef, n, m)
        for i in 1:n, j in 1:m
            x, y, z = map(split(readline(f))) do s parse(Float64, s) end
            cnet[i,j] = Point(x, y, z)
        end
    end
    cnet
end

function write_bzr(filename, cnet)
    open(filename, "w") do f
        n, m = size(cnet)
        println(f, "$(n - 1) $(m - 1)")
        for i in 1:n, j in 1:m
            p = cnet[i,j]
            println(f, "$(p.x) $(p.y) $(p.z)")
        end
    end
end

# only uses the boundary
function optimize!(cnet, alpha = -0.25)
    beta = 0.25 - alpha
    n, m = size(cnet)
    nvars = (n - 2) * (m - 2)
    A = one(Matrix{Float64}(undef, nvars, nvars))
    b = zeros(Point, nvars)
    for i in 2:n-1, j in 2:m-1
        row = (i - 2) * (m - 2) + j - 1
        for (k, l) in [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]
            if k == 1 || l == 1 || k == n || l == m
                b[row] += cnet[k,l] * beta
            else
                col = (k - 2) * (m - 2) + l - 1
                A[row,col] = -beta
            end
        end
        for (k, l) in [(i-1,j-1), (i+1,j+1), (i+1,j-1), (i-1,j+1)]
            if k == 1 || l == 1 || k == n || l == m
                b[row] += cnet[k,l] * alpha
            else
                col = (k - 2) * (m - 2) + l - 1
                A[row,col] = -alpha
            end
        end
    end
    sol = A \ b
    for i in 2:n-1, j in 2:m-1
        row = (i - 2) * (m - 2) + j - 1
        cnet[i,j] = sol[row]
    end
    cnet
end

function c1_optimize!(cnet, alpha = 0)
    mask1 = 1/20 * [ 0  0  1  0  0;
                     0  2 -8  2  0;
                     1 -8 20 -8  1;
                     0  2 -8  2  0;
                     0  0  1  0  0]
    mask2 = 1/36 * [  1  -4   6  -4   1;
                     -4  16 -24  16  -4;
                      6 -24  36 -24   6;
                     -4  16 -24  16  -4;
                      1  -4   6  -4   1]
    mask = mask1 * (1 - alpha) + mask2 * alpha
    n, m = size(cnet)
    nvars = (n - 4) * (m - 4)
    A = zeros(Float64, nvars, nvars)
    b = zeros(Point, nvars)
    for i in 3:n-2, j in 3:m-2
        row = (i - 3) * (m - 4) + j - 2
        for k in 1:5, l in 1:5
            ik = i + k - 3
            jl = j + l - 3
            if ik <= 2 || jl <= 2 || ik >= n - 1 || jl >= m - 1
                b[row] -= cnet[ik,jl] * mask[k,l]
            else
                col = (ik - 3) * (m - 4) + jl - 2
                A[row,col] = mask[k,l]
            end
        end
    end
    sol = A \ b
    for i in 3:n-2, j in 3:m-2
        row = (i - 3) * (m - 4) + j - 2
        cnet[i,j] = sol[row]
    end
    cnet
end

function bernstein(n, k, u)
    res = zeros(n + 1)
    res[n-k+1] = 1
    for i in 1:n, j in n:-1:i
        res[j+1] = res[j] * u + res[j+1] * (1 - u)
    end
    res[n+1]
end

function bezier_eval(cps, u)
    n = length(cps) - 1
    p = Point(0, 0, 0)
    for i in 0:n
        p += cps[i+1] * bernstein(n, i, u)
    end
    p
end

# only uses the boundary
function coons_eval(cnet, u, v)
    corner = [cnet[1,1] cnet[1,end]; cnet[end,1] cnet[end,end]]
    bezier_eval(cnet[1,:], v) * (1 - u) +
        bezier_eval(cnet[end,:], v) * u +
        bezier_eval(cnet[:,1], u) * (1 - v) +
        bezier_eval(cnet[:,end], u) * v -
        ([1-u u] * corner * [1-v, v])[1,1]
end

function cubic_coons_eval(cnet, u, v)
    n, m = map(x -> x - 1, size(cnet))
    # S(0,v)
    p1 = bezier_eval(cnet[1,:], v)
    t1 = n * (bezier_eval(cnet[2,:], v) - p1)
    # S(1,v)
    p2 = bezier_eval(cnet[end,:], v)
    t2 = n * (p2 - bezier_eval(cnet[end-1,:], v))
    # S(u,0)
    p3 = bezier_eval(cnet[:,1], u)
    t3 = m * (bezier_eval(cnet[:,2], u) - p3)
    # S(u,1)
    p4 = bezier_eval(cnet[:,end], u)
    t4 = m * (p4 - bezier_eval(cnet[:,end-1], u))
    # S(0,0)
    cp1 = cnet[1,1]
    cu1 = n * (cnet[2,1] - cp1)
    cv1 = m * (cnet[1,2] - cp1)
    ct1 = m * n * (cnet[2,2] - cnet[1,2] - (cnet[2,1] - cp1))
    # S(1,0)
    cp2 = cnet[end,1]
    cu2 = n * (cp2 - cnet[end-1,1])
    cv2 = m * (cnet[end,2] - cp2)
    ct2 = m * n * (cnet[end,2] - cnet[end-1,2] - (cp2 - cnet[end-1,1]))
    # S(0,1)
    cp3 = cnet[1,end]
    cu3 = n * (cnet[2,end] - cp3)
    cv3 = m * (cp3 - cnet[1,end-1])
    ct3 = m * n * (cnet[2,end] - cp3 - (cnet[2,end-1] - cnet[1,end-1]))
    # S(1,1)
    cp4 = cnet[end,end]
    cu4 = n * (cp4 - cnet[end-1,end])
    cv4 = m * (cp4 - cnet[end,end-1])
    ct4 = m * n * (cp4 - cnet[end-1,end] - (cnet[end,end-1] - cnet[end-1,end-1]))
    corner = [cp1 cv1 cp3 cv3;
              cu1 ct1 cu3 ct3;
              cp2 cv2 cp4 cv4;
              cu2 ct2 cu4 ct4]
    f0 = x -> 2x^3 - 3x^2 + 1
    f1 = x -> -2x^3 + 3x^2
    g0 = x -> x^3 - 2x^2 + x
    g1 = x -> x^3 - x^2
    p1 * f0(u) + p2 * f1(u) + p3 * f0(v) + p4 * f1(v) +
        t1 * g0(u) + t2 * g1(u) + t3 * g0(v) + t4 * g1(v) -
        ([f0(u) g0(u) f1(u) g1(u)] * corner * [f0(v), g0(v), f1(v), g1(v)])[1,1]
end

function coons_mesh(filename, cnet, resolution, cubic = false)
    surf = cubic ? cubic_coons_eval : coons_eval
    open(filename, "w") do f
        for i in 0:resolution, j in 0:resolution
            u = i / resolution
            v = j / resolution
            p = surf(cnet, u, v)
            println(f, "v $(p.x) $(p.y) $(p.z)")
        end
        for i in 1:resolution, j in 1:resolution
            base = i * (resolution + 1) + j + 1
            println(f, "f $base $(base-1) $(base-resolution-1)")
            println(f, "f $(base-resolution-1) $(base-1) $(base-resolution-2)")
        end
    end
end

function test(root, alpha)
    cnet = read_bzr("$root.bzr")
    optimize!(cnet, -alpha)
    write_bzr("$root-$alpha.bzr", cnet)
end

function test_c1(root, alpha)
    cnet = read_bzr("$root.bzr")
    c1_optimize!(cnet, -alpha)
    write_bzr("$root-$alpha.bzr", cnet)
end

end
