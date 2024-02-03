"""
	fibonaccigrid(;N=nothing, sepAng=nothing)

This function returns a vector `Nx2` of LAT, LON values for a `N` points grid built with the Fibonacci Spiral method.

### Arguments:
- `N`: The number of points to generate.
- `sepAng`: The separation angle for the grid of points to be generated [rad].
- `unit`: `:rad` or `:deg`

### Output:
- `vec`: A vector of SVector{2}(lon,lat) objects of LAT-LON coordinates in rad (or deg) (LAT=y, LON=x).
"""
function fibonaccigrid(;N=nothing, sepAng=nothing, unit=:rad)	
	if N isa Nothing && sepAng isa Nothing
		error("Input one argument between N and sepAng...")
	elseif sepAng isa Nothing
		vec = fibonaccisphere_classic(N; coord=:sphe)
	elseif N isa Nothing
		N,sepAng = points_required_for_separation_angle(sepAng)
		vec = fibonaccisphere_classic(N; coord=:sphe)
	else
		error("Input one argument between N and sepAng...")
	end

	# Unit Conversion
	if unit == :deg
		vec = map(x -> rad2deg.(x), vec)
	end

	return vec
end

"""
	fibonaccisphere_classic(N::Int)
	
This function generates `N` uniformly distributed points on the surface of a unitary sphere using the classic Fibonacci Spiral method as described in [1].
Contrary to the Ichosahedral grid generation process, with the Fibonacci Spiral method it is possible to generate a grid of points uniformly distributed in the area for a generic `N` value. As a drawback, the structure of the points do not follow a "perfect simmetry" however, the density of points in the area is preserved quasi-constant.

### Arguments:
- `N::Int`: The number of points to generate.
- `coord::Symbol`: The type of coordinates of generated points (`:sphe` | `:cart`).
- `radius`: the sphere radius in meters (unitary as default)

### Output:
- `pointsVec``: `N x 1` array containing the SVector of the generated points. Each element corresponds to a point on the surface of the sphere, the SVector contains either the x, y, and z (:cart) or lon, lat (:sphe) (LAT=y, LON=x) coordinates of the point.

### References
1. http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
"""
function fibonaccisphere_classic(N::Int; coord::Symbol=:sphe, spheRadius=1.0)
	goldenRatio = (1 + sqrt(5))/2	
	if coord == :sphe # :sphe | :cart
		pointsVec = map(0:N-1) do k
			θ = 2π * k/ goldenRatio # [0,2π] [LON]
			ϕ = acos(1 - 2(k+0.5)/N) # [0,π] from North Pole [LAT]
			
			SVector(rem2pi(θ, RoundNearest), π/2 - ϕ) # wrap (lon,lat)
		end
	else
		pointsVec = map(0:N-1) do k
			θ = 2π * k/ goldenRatio # [0,2π] [LON]
			ϕ = acos(1 - 2(k+0.5)/N) # [0,π] from North Pole [LAT]
			
			SVector(spheRadius.*sin(ϕ)*cos(θ), spheRadius.*sin(ϕ)*sin(θ), spheRadius.*cos(ϕ))
		end
	end

	return pointsVec
end

"""
	points_required_for_separation_angle(angle)

This function computes the minimum number of points required on the surface of a sphere to achieve a desired separation angle between any two adjacent points. The function uses the bisection method to find the minimum number of points needed and returns the higher end of the precision.
Instead of checking all the possible pairs of the N points generated with the Fibonacci Spiral, only the first `pointsToCheck` are checked in order to evaluate the minimum separation angle.

### Arguments:
- `sepAng`: a float representing the desired separation angle between two adjacent points on the surface of the sphere.
- `spheRadius`: an optional float representing the radius of the sphere. If not provided, it defaults to 1.0.
- `pointsToCheck`: an optional integer representing the number of points to generate on the surface of the sphere. If not provided, it defaults to 50.
- `maxPrec`: an optional integer representing the maximum precision for the number of points generated on the surface of the sphere. If not provided, it defaults to 10^7.
- `tol`: an optional integer representing the tolerance for the bisection method used to find the minimum number of points needed to achieve the desired separation angle. If not provided, it defaults to 10.

### Output:
- `Ns[2]`: an integer representing the minimum number of points required on the surface of the sphere to achieve the desired separation angle.
"""
function points_required_for_separation_angle(sepAng; spheRadius=1.0, pointsToCheck::Int=50, maxPrec=10^7, tol=10)
	# Bisection
	Ns = [2, maxPrec]
	
	## Define Inner Functions ---------------------------------------------------------------------
	# Find the separation angle for a given N
	f(N) = find_separation_angle(fibonaccisphere_classic_partial(N; spheRadius=spheRadius, pointsToCheck=pointsToCheck))
	# Distance between the tested points
	tolerance(v) = v[2] - v[1]
	
	# Compute the startin angles
	angs = map(f, Ns)
	sepAng > angs[1] && return Ns[1] # The lower is already sufficient
	sepAng < angs[2] && error("$(maxPrec) points is not sufficient for the requested separation angle")
	thisSep = Inf
	thisTol = tolerance(Ns)
	
	while thisTol > tol
	  N = Ns[1] + floor(thisTol/2)
	  thisSep = f(N)
	  if thisSep > sepAng
		# The separation angle is too big, the density of points has to be increased. We keep increasing the lower index
		Ns[1] = N
	  else
		# # The separation angle is too small, the desity of points has to be decreased. We have to start decreasing the upper end
		Ns[2] = N
	  end
	  thisTol = tolerance(Ns)
	end
	
	return Ns[2],rad2deg(thisSep) # We return the higher end of the precision
  end

"""
	find_separation_angle(points)

This function takes an array of 3D Cartesian points as input and computes the smallest angle between any two points in the array. It does this by iterating over all unique pairs of points in the array and computing the angle between them using the `angle`` function. The smallest angle encountered during the iteration is stored in the variable `sep` and returned as the output of the function.

### Arguments:
- `points``: an array of 3D points in the Cartesian plane represented as Tuples, Arrays, SVectors.

### Output:
- `sep`: the smallest angle between any two points in the input array, measured in radians.
"""
function find_separation_angle(points)
	sep = Inf
	L = length(points)

	if L<2
		error("The input vector must contain at least 2 points...")
	end

	@inbounds for i in 1:L-1
		for j in i+1:L # Check triang sup (avoid repeating pairs of points, order is not important)
			# acos(dot(points[i],points[j]) / norm(points[i])*norm(points[j]))
			temp = angle(points[i],points[j]) # avoid errors with float approximations (<-1, >1)
			sep = temp < sep ? temp : sep
		end
	end
	
	return sep
end

"""
	fibonaccisphere_classic_partial(N; spheRadius=1.0, pointsToCheck::Int=50)

### Arguments:
- `N`: an integer representing the number of points to generate on the surface of the sphere.
- `spheRadius`: (optional) a float representing the radius of the sphere.
- `pointsToCheck`: (optional) an integer representing the number of points to return starting from the first generated.

### Output:
- `points`: an array of 3D points on the surface of the sphere represented as SVector{3}.
"""
function fibonaccisphere_classic_partial(N; spheRadius=1.0, pointsToCheck::Int=50)
	goldenRatio = (1 + sqrt(5))/2
	points = map(0:min(pointsToCheck,N)-1) do k
		θ = 2π * k/ goldenRatio # [0,2π] [LON]
		ϕ = acos(1 - 2(k+0.5)/N) # [0,π] from North Pole [LAT]

		SVector(spheRadius.*sin(ϕ)*cos(θ), spheRadius.*sin(ϕ)*sin(θ), spheRadius.*cos(ϕ))
	end

	return points
end

## Alternative implementations ------------------------------------------------------------------------------
# These implementations are here only for keep track of the alternatives and for future developement
"""
    fibonaccisphere_optimization1(N)

Create a set of `N` uniformly distributed points on the surface of a sphere using the Fibonacci spiral method optimized as described in [1].

This method is called Offset Fibonacci Lattice which is one method to optimize the minimum nearest-neighbor distance. 
We need to move (offset) all the points slightly farther away from the poles. This of course means, that almost all of them become slightly closer together. Offsetting the points of the Fibonacci lattice slightly away from the poles produces a packing that is up to 8.3% tighter than the canonical Fibonacci lattice.

For `n>100`, an improvement can be made beyond this, by initially placing a point at each pole, and then placing the remaining `n-2` points. This not only (very sightly) improves minimal nearest packing, but it also prevents a large gap at each pole.

### Arguments
- `N::Int`: The number of points to generate.

### Output
- `points::Matrix{Float64}`: A `N`x`3` matrix where each row corresponds to a point `(x,y,z)` on the surface of the unitary sphere.

### References
1. http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
"""
function fibonaccisphere_optimization1(N::Int)
	points = zeros(N,3)
	
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
	
	goldenRatio = (1 + sqrt(5))/2

	points = map(0:N-1) do k
		θ = 2π * k / goldenRatio # Longitude
		ϕ = acos(1 - 2*(k+epsilon)/(N-1+2*epsilon)) # Latitude

		SVector(sin(ϕ)*cos(θ), sin(ϕ)*sin(θ), cos(ϕ))
	end

	return points
end

"""
	fibonaccisphere_alternative1(N::Int)

This function generates points on the surface of a unit sphere using the Fibonacci spiral method. The function takes an integer `N` as an input, which specifies the number of points to be generated.

### Arguments
- `N::Int`: The number of points to generate. This is an integer value.
### Output
- `N x 3` matrix containing the generated points. Each row of the matrix corresponds to a point on the surface of the sphere, and the columns correspond to the x, y, and z coordinates of the point.
"""
function fibonaccisphere_alternative1(N::Int)
	# //FIX: to be further tested 
	points = zeros(N,3)
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

function points_required_for_separation_angle_var1(sepAng; pointToCheck::Int=20, altitude=1.0)
	# Brute Force
	N = 2
	goldenRatio = (1 + sqrt(5))/2
	while true	
		M = minimum([N,pointToCheck])
		testAng = 999 .* ones(M,M)
		# Fibonacci Sphere generation	
		points = map(0:M-1) do k
			θ = 2π * k/ goldenRatio # [0,2π] [LON]
			ϕ = acos(1 - 2(k+0.5)/N) # [0,π] from North Pole [LAT]
			
			SVector(altitude.*sin(ϕ)*cos(θ), altitude.*sin(ϕ)*sin(θ), altitude.*cos(ϕ))
		end

		# Angle between vectors computation
		for i in 1:M
			for j in 1:M
				if i!=j
					# acos(dot(points[i],points[j]) / norm(points[i])*norm(points[j]))
					testAng[i,j] = angle(points[i],points[j])
				end
			end
		end

		min_angle = minimum(testAng[:])
		if min_angle <= sepAng
			return N,rad2deg(min_angle)
		end
		N += 1
	end
end

function points_required_for_separation_angle_var2(sepAng; pointToCheck::Int=20, altitude=1.0)
	# Spherical Cap
	N = Int(round(2/(1-cos(sepAng)))) # Spherical cap
	goldenRatio = (1 + sqrt(5))/2
	while true	
		M = minimum([N,pointToCheck])
		testAng = 999 .* ones(M,M)
		# Fibonacci Sphere generation	
		points = map(0:M-1) do k
			θ = 2π * k/ goldenRatio # [0,2π] [LON]
			ϕ = acos(1 - 2(k+0.5)/N) # [0,π] from North Pole [LAT]
			
			SVector(altitude.*sin(ϕ)*cos(θ), altitude.*sin(ϕ)*sin(θ), altitude.*cos(ϕ))
		end

		# Angle between vectors computation
		for i in 1:M
			for j in 1:M
				if i!=j
					# acos(dot(points[i],points[j]) / norm(points[i])*norm(points[j]))
					testAng[i,j] = angle(points[i],points[j])
				end
			end
		end

        min_angle = minimum(testAng[:])
        if min_angle <= sepAng
            return N,rad2deg(min_angle)
        end
        N += 1
    end
end