"""
    icogrid(;N=nothing, sepAng=nothing, unit=:rad, type=:lla) -> Vector{Point{🌐,<:LatLon{WGS84Latest}}}

This function returns a vector `Nx2` of LAT, LON values for a `N` points grid
built with the Fibonacci Spiral method.

## Arguments:
- `N`: The number of points to generate.
- `sepAng`: The separation angle for the grid of points to be generated [rad].
- `unit`: `:rad` or `:deg`
- `type`: `:lla` or `:point`. Output type either `LLA` or `Point2`

## Output:
- `out`: an Vector{Point{🌐,<:LatLon{WGS84Latest}}} of points in the icosahedral grid.
"""
function icogrid(; N::Union{Int,Nothing}=nothing, sepAng::Union{ValidAngle,Nothing}=nothing)
    if isnothing(sepAng) && !isnothing(N)
        vec = _icogrid(N; coord=:sphe)
    elseif !isnothing(sepAng) && isnothing(N)
        # Inputs validation (sepAng must be positive Uniful value)    
        _sepAng = let
            x = sepAng isa Real ? sepAng * ° : sepAng |> u"°" # Convert to Uniful °
            abs(x) ≤ 360° || error(
#! format: off
"The sepAng provided as numbers must be expressed in radians and satisfy -360° ≤ x ≤ 360°. 
Consider using `°` (or `rad`) from `Unitful` if you want to pass numbers in degrees (or rad), by doing `x * °` (or `x * rad`)."
#! format: on   
            )
            if x < 0
                @warn "Input sepAng is negative, it will be converted to positive..."
                abs(x)
            else
                x
            end
        end
        N, _ = _points_required_for_separation_angle(_sepAng)
        vec = _icogrid(N; coord=:sphe)
    else
        error("Input one argument between N and sepAng...")
    end

    # Unit Conversion	
    out = map(x -> LatLon{WGS84Latest}(rad2deg.(x)...) |> Point, vec)

    return out
end

"""
    _icogrid(N::Int)
    
This function generates `N` uniformly distributed points on the surface of a
unitary sphere using the classic Fibonacci Spiral method as described in [1].
Contrary to the Ichosahedral grid generation process, with the Fibonacci Spiral
method it is possible to generate a grid of points uniformly distributed in the
area for a generic `N` value. As a drawback, the structure of the points do not
follow a "perfect simmetry" however, the density of points in the area is
preserved quasi-constant.

## Arguments:
- `N::Int`: The number of points to generate.
- `coord::Symbol`: The type of coordinates of generated points (`:sphe` | `:cart`).
- `radius`: the sphere radius in meters (unitary as default)

## Output:
- `pointsVec`: a `Vector{SVector}` of the generated points. Each element \
corresponds to a point on the surface of the sphere, the SVector contains \
either the x, y, and z (:cart) or lat, lon (:sphe) (LAT=x, LON=y) in rad \
coordinates of the point.

## References
1. http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
"""
function _icogrid(N::Int; coord::Symbol=:sphe, spheRadius=1.0)
    pointsVec = if coord == :sphe # :sphe | :cart
        map(0:N-1) do k
            θ, ϕ = _get_theta_phi(k, N)
            SVector(π / 2 - ϕ, rem(θ, 2π, RoundNearest)) # wrap (lat,lon)
        end
    else
        map(0:N-1) do k
            θ, ϕ = _get_theta_phi(k, N)
            sθ, cθ = sincos(θ)
            sϕ, cϕ = sincos(ϕ)
            spheRadius * SVector(sϕ * cθ, sϕ * sθ, cϕ) # Real, representing meters
        end
    end

    return pointsVec
end

"""
    _points_required_for_separation_angle(angle)

This function computes the minimum number of points required on the surface of a
sphere to achieve a desired separation angle between any two adjacent points.
The function uses the bisection method to find the minimum number of points
needed and returns the higher end of the precision. Instead of checking all the
possible pairs of the N points generated with the Fibonacci Spiral, only the
first `pointsToCheck` are checked in order to evaluate the minimum separation
angle.

## Arguments:
- `sepAng`: a float representing the desired separation angle between two \
adjacent points on the surface of the sphere.
- `spheRadius`: an optional float representing the radius of the sphere. If not \
provided, it defaults to 1.0.
- `pointsToCheck`: an optional integer representing the number of points to \
generate on the surface of the sphere. If not provided, it defaults to 50.
- `maxPrec`: an optional integer representing the maximum precision for the \
number of points generated on the surface of the sphere. If not provided, it \
defaults to 10^7.
- `tol`: an optional integer representing the tolerance for the bisection method \
used to find the minimum number of points needed to achieve the desired \
separation angle. If not provided, it defaults to 10.

## Output:
- `Ns[2]`: an integer representing the minimum number of points required on the \
surface of the sphere to achieve the desired separation angle.
- `thisSep`: an Uniful.Quantity in degrees representing the separation angle \
between two adjacent points on the surface of the sphere.
"""
function _points_required_for_separation_angle(sepAng::ValidAngle; spheRadius=1.0, pointsToCheck::Int=50, maxPrec=10^7, tol=10)
    # Input validation (make it deg)
    _sepAng = sepAng isa Real ? sepAng * ° : sepAng |> u"°" # Convert to Uniful °
    # Bisection
    Ns = [2, maxPrec]

    ## Define Inner Functions ---------------------------------------------------------------------
    # Find the separation angle for a given N
    f(N) = _find_min_separation_angle(_fibonaccisphere_classic_partial(N; spheRadius=spheRadius, pointsToCheck=pointsToCheck)) # radians
    # Distance between the tested points
    tolerance(v) = v[2] - v[1]

    # Compute the starting angles
    angs = map(f, Ns)
    _sepAng > angs[1] && return Ns[1], f(Ns[1]) # The lower is already sufficient
    _sepAng < angs[2] && error("$(maxPrec) points is not sufficient for the requested separation angle")
    thisSep = Inf
    thisTol = tolerance(Ns)

    while thisTol > tol
        N = Ns[1] + floor(thisTol / 2)
        thisSep = f(N)
        if thisSep > _sepAng
            # The separation angle is too big, the density of points has to be increased. We keep increasing the lower index
            Ns[1] = N
        else
            # # The separation angle is too small, the desity of points has to be decreased. We have to start decreasing the upper end
            Ns[2] = N
        end
        thisTol = tolerance(Ns)
    end

    return Ns[2], thisSep # We return the higher end of the precision and angle in degrees
end

"""
    _find_min_separation_angle(points)

This function takes an array of 3D Cartesian points as input and computes the
smallest angle between any two points in the array. It does this by iterating
over all unique pairs of points in the array and computing the angle between
them using the `angle`` function. The smallest angle encountered during the
iteration is stored in the variable `sep` and returned as the output of the
function.

## Arguments:
- `points`: an array of 3D points in the Cartesian plane represented as Tuples, \
Arrays, SVectors.

## Output:
- `sep`: the smallest angle between any two points in the input array, returned \
as Uniful.Quantity in degrees.
"""
function _find_min_separation_angle(points)
    sep = Inf
    L = length(points)

    if L < 2
        error("The input vector must contain at least 2 points...")
    end

    @inbounds for i in 1:L-1
        for j in i+1:L # Check triang sup (avoid repeating pairs of points, order is not important)
            # acos(dot(points[i],points[j]) / norm(points[i])*norm(points[j]))
            temp = angle(points[i], points[j]) # avoid errors with float approximations (<-1, >1) using AngleBetweenVectors.jl
            sep = temp < sep ? temp : sep
        end
    end

    return rad2deg(sep) * ° # return as ValidAngle (in deg)
end

"""
    _fibonaccisphere_classic_partial(N; spheRadius=1.0, pointsToCheck::Int=50)

## Arguments:
- `N`: an integer representing the number of points to generate on the surface \
of the sphere.
- `spheRadius`: (optional) a float representing the radius of the sphere.
- `pointsToCheck`: (optional) an integer representing the number of points to \
return starting from the first generated.

## Output:
- `points`: an array of 3D points on the surface of the sphere represented as \
SVector{3}.
"""
function _fibonaccisphere_classic_partial(N; spheRadius=1.0, pointsToCheck::Int=50)
    points = map(0:min(pointsToCheck, N)-1) do k
        θ, ϕ = _get_theta_phi(k, N)
        sθ, cθ = sincos(θ)
        sϕ, cϕ = sincos(ϕ)
        spheRadius * SVector(sϕ * cθ, sϕ * sθ, cϕ)
    end

    return points
end

# Helper Functions
"""
    _get_theta_phi(k::Number, N::Number) -> Tuple{Number, Number}

Calculate the spherical coordinates θ (theta) and ϕ (phi) for a given index `k`
and total number of points `N` using the Golden Ratio method.

## Arguments
- `k::Number`: The index of the point for which the spherical coordinates are to \
be calculated. This should typically be an integer between `0` and `N-1`.
- `N::Number`: The total number of points for which the spherical coordinates \
are to be calculated.

## Returns
- A tuple `(θ, ϕ)` where:
    - `θ::Number`: The longitude angle in rad, ranging from `[0, 360]`.
    - `ϕ::Number`: The latitude angle in rad, ranging from `[0, 180]` from \
the North Pole.
"""
function _get_theta_phi(k::Number, N::Number)
    goldenRatio = (1 + sqrt(5)) / 2
    θ = 2π * k / goldenRatio # [0,2π] [LON]
    ϕ = acos(1 - 2(k + 0.5) / N) # [0,π] from North Pole [LAT]

    return (θ, ϕ) # Return values as rad
end

## Alternative implementations ------------------------------------------------------------------------------
# These implementations are here only for keep track of the alternatives and for future developement
"""
    fibonaccisphere_optimization1(N)

Create a set of `N` uniformly distributed points on the surface of a sphere
using the Fibonacci spiral method optimized as described in [1].

This method is called Offset Fibonacci Lattice which is one method to optimize
the minimum nearest-neighbor distance. We need to move (offset) all the points
slightly farther away from the poles. This of course means, that almost all of
them become slightly closer together. Offsetting the points of the Fibonacci
lattice slightly away from the poles produces a packing that is up to 8.3%
tighter than the canonical Fibonacci lattice.

For `n>100`, an improvement can be made beyond this, by initially placing a
point at each pole, and then placing the remaining `n-2` points. This not only
(very sightly) improves minimal nearest packing, but it also prevents a large
gap at each pole.

## Arguments
- `N::Int`: The number of points to generate.

## Output
- `points::Matrix{Float64}`: A `N`x`3` matrix where each row corresponds to a \
point `(x,y,z)` on the surface of the unitary sphere.

## References
1. http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
"""
function fibonaccisphere_optimization1(N::Int)
    points = zeros(N, 3)

    if N >= 600000
        epsilon = 214
    elseif N >= 400000
        epsilon = 75
    elseif N >= 11000
        epsilon = 27
    elseif N >= 890
        epsilon = 10
    elseif N >= 177
        epsilon = 3.33
    elseif N >= 24
        epsilon = 1.33
    else
        epsilon = 0.33
    end

    goldenRatio = (1 + sqrt(5)) / 2

    points = map(0:N-1) do k
        θ = 2π * k / goldenRatio # Longitude
        ϕ = acos(1 - 2 * (k + epsilon) / (N - 1 + 2 * epsilon)) # Latitude

        SVector(sin(ϕ) * cos(θ), sin(ϕ) * sin(θ), cos(ϕ))
    end

    return points
end

"""
    fibonaccisphere_alternative1(N::Int)

This function generates points on the surface of a unit sphere using the
Fibonacci spiral method. The function takes an integer `N` as an input, which
specifies the number of points to be generated.

## Arguments
- `N::Int`: The number of points to generate. This is an integer value.
## Output
- `N x 3` matrix containing the generated points. Each row of the matrix \
corresponds to a point on the surface of the sphere, and the columns \
correspond to the x, y, and z coordinates of the point.
"""
function fibonaccisphere_alternative1(N::Int)
    # //FIX: to be further tested
    points = zeros(N, 3)
    goldenSphere = π * (3.0 - sqrt(5.0))
    off = 2.0 / N

    points = map(0:N-1) do k
        y = k * off - 1.0 + (off / 2.0)
        r = sqrt(1.0 - y * y)
        phi = k * goldenSphere

        SVector(cos(phi) * r, y, sin(phi) * r)
    end

    return points
end