"""
	color_greedy(points::AbstractVector{<:Union{LatLon,Point{🌐,<:LatLon{WGS84Latest}}}}, radius::Number, count::Int; cutoff_dist::Number=radius * 10.0, min_dist::Number=radius * 2)

Color a set of points with n colors with an extended greedy algorithm. The algorithm makes sure no points closer than min_dist share the same color.
Additionally, it tries to optimize the total spatial distribution of colors.

## Arguments
- `points`: A vector of points.
- `radius`: The radius of the individual cells.
- `count`: The number of colors to use.
- `cutoff_dist`: The cutoff distance for the greedy algorithm. A low value will result in a faster but less optimal coloring.
- `min_dist`: The minimum distance for the greedy algorithm. If the minimum distance between two colors is not met, the coloring will fail.

## Returns
- A vector of integers where each entry represents the "color" of the cell with the same index.
"""
function color_greedy(points::AbstractVector{<:Union{LatLon,Point{🌐,<:LatLon{WGS84Latest}}}}, radius::Number, count::Integer; cutoff_dist::Number=radius * 5.0, min_dist::Number=radius * 2.0)::Vector{Int}
	num_points = length(points)
	g = SimpleWeightedGraph(num_points)
	
	# Calculate the distance to adjacent points
	for i in 1:num_points
		for j in (i+1):num_points
			dist = great_circle_distance(points[i], points[j])
			if dist <= cutoff_dist
				Graphs.add_edge!(g, i, j, dist)
			end
		end
	end

	seq = convert(Vector, sortperm(Graphs.degree(g); rev=true))
	nvg = Graphs.nv(g)
	
	colors::Vector{Int} = zeros(Int, nvg)
	color_count = zeros(Int, count)
	
	for v in seq
		color_dist = fill(Inf, count)

		# Find the min distance for each color
		for nb in Graphs.neighbors(g, v)
			if colors[nb] == 0	
				continue
			end

			dist = SimpleWeightedGraphs.get_weight(g, v, nb)
            color_dist[colors[nb]] = min(color_dist[colors[nb]], dist)
		end

		inf_dist_cols = findall(c -> c == Inf, color_dist)
		
		# Select the color that is not used or has the largest distance
		color_idx = -1
		if length(inf_dist_cols) == 0
			color_idx = sortperm(color_dist; rev=true)[1]
			if color_dist[color_idx] < min_dist
				error("Could not color graph with $count colors")
			end
		elseif length(inf_dist_cols) == 1
			color_idx = inf_dist_cols[1]
		else
			color_idx = sort(inf_dist_cols; by=x -> color_count[x])[1]
		end

		colors[v] = color_idx
		color_count[color_idx] += 1
	end

	return colors
end
