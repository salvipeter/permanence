module DiscreteGB

import Base: getindex, setindex!
using LinearAlgebra

const Point = Vector{Float64}
const GBIndex = NTuple{3, Int}

struct BezierPatch
    n :: Int
    d :: Int
    center :: Point
    cpts :: Dict{GBIndex,Point}
end

getindex(s::BezierPatch, idx::GBIndex) = s.cpts[idx]
getindex(s::BezierPatch, i::Int, j::Int, k::Int) = s.cpts[i,j,k]
setindex!(s::BezierPatch, v::Point, idx::GBIndex) = s.cpts[idx] = v
setindex!(s::BezierPatch, v::Point, i::Int, j::Int, k::Int) = s.cpts[i,j,k] = v

function readGBP(filename)
    read_numbers(f, numtype) = map(s -> parse(numtype, s), split(readline(f)))
    local result
    open(filename) do f
        n, d = read_numbers(f, Int)
        l = Int(floor((d + 1) / 2))
        cp = 1 + Int(floor(d / 2))
        cp = n * cp * l + 1
        side, col, row = 1, 0, 0
        center = read_numbers(f, Float64)
        result = BezierPatch(n, d, center, Dict())
        for i in 1:cp-1
            if col >= d - row
                side += 1
                if side > n
                    side = 1
                    row += 1
                end
                col = row
            end
            p = read_numbers(f, Float64)
            result[side,col,row] = p
            if col < l
                result[mod1(side-1,n),d-row,col] = p
            elseif d - col < l
                result[mod1(side+1,n),row,d-col] = p
            end
            col += 1
        end
    end
    result
end

function writeOBJ(surface, filename)
    open(filename, "w") do f
        for (i, p) in surface.cpts
            println(f, "v $(p[1]) $(p[2]) $(p[3])")
        end
        n = length(surface.cpts)
        for i in 1:n
            println(f, "p $i")
        end
        for i in 1:surface.n, j in 0:surface.d, k in 1:(surface.d-1)÷2
            p = surface[i,j,k-1]
            println(f, "v $(p[1]) $(p[2]) $(p[3])")
            p = surface[i,j,k]
            println(f, "v $(p[1]) $(p[2]) $(p[3])")
            println(f, "l $(n+1) $(n+2)")
            n += 2
        end
        for i in 1:surface.n, k in 0:(surface.d-1)÷2, j in 1:surface.d
            p = surface[i,j-1,k]
            println(f, "v $(p[1]) $(p[2]) $(p[3])")
            p = surface[i,j,k]
            println(f, "v $(p[1]) $(p[2]) $(p[3])")
            println(f, "l $(n+1) $(n+2)")
            n += 2
        end
        if iseven(surface.d)
            p = surface.center
            println(f, "v $(p[1]) $(p[2]) $(p[3])")
            n += 1
            println(f, "p $n")
            l = surface.d ÷ 2
            for i in 1:surface.n
                p = surface[i,l,l-1]
                println(f, "v $(p[1]) $(p[2]) $(p[3])")
                println(f, "l $(n) $(n+i)")
            end
        end
    end
end

function smooth_controls(surface, α = -1/4)
    β = (1 - 4α) / 4
    boundary = surface.d * surface.n
    layers = (surface.d + 1) ÷ 2
    nvars = surface.n * (surface.d ÷ 2 + 1) * layers - boundary + (iseven(surface.d) ? 1 : 0)

    index = Dict{GBIndex,Int}()
    last = iseven(surface.d) ? 1 : 0
    function geti(i, j, k)
        iseven(surface.d) && j == k == layers && return 1 # central cp
        while j < k || j >= surface.d - k
            if j < k
                i = mod1(i - 1, surface.n)
                j, k = surface.d - k, j
            else
                i = mod1(i + 1, surface.n)
                j, k = k, surface.d - j
            end
        end
        haskey(index, (i,j,k)) && return index[i,j,k]
        last += 1
        index[i,j,k] = last
    end

    A = one(Matrix{Float64}(undef, nvars, nvars))
    b = zeros(nvars, 3)
    for i in 1:surface.n, j in 1:surface.d-1, k in 1:layers-1
        (j < k || j >= surface.d - k) && continue
        row = geti(i, j, k)
        for x in -1:1, y in -1:1
            x == 0 && y == 0 && continue
            weight = x * y == 0 ? β : α
            if j + x == 0 || k + y == 0 || j + x == surface.d
                b[row,:] += surface[i,j,k] * weight
            else
                col = geti(i, j + x, k + y)
                A[row,col] = -weight
            end
        end
    end
    if iseven(surface.d)
        for i in 1:surface.n
            col = geti(i, layers - 1, layers - 1)
            A[1,col] = -α * 4 / surface.n
            col = geti(i, layers, layers - 1)
            A[1,col] = -β * 4 / surface.n
        end
    end

    sol = A \ b
    result = BezierPatch(surface.n, surface.d,
                         iseven(surface.d) ? sol[1,:] : surface.center,
                         Dict{GBIndex,Point}())
    for i in 1:surface.n, j in 0:surface.d, k in 0:layers-1
        if j == 0 || k == 0 || j == surface.d
            result[i,j,k] = surface[i,j,k]
        else
            result[i,j,k] = sol[geti(i,j,k),:]
        end
    end
    result
end

function test(α = -0.25)
    s = readGBP("cagd86.gbp")
    writeOBJ(s, "/tmp/cagd86.obj")
    writeOBJ(smooth_controls(s, α), "/tmp/cagd86-smooth.obj")
    s = readGBP("setback.gbp")
    writeOBJ(s, "/tmp/setback.obj")
    writeOBJ(smooth_controls(s, α), "/tmp/setback-smooth.obj")
end

end # module
