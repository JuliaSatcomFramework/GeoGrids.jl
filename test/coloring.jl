@testitem "ColoringGreedy" begin
	using GeoGrids

	radius = 200e3
	min_dist = 2 * radius
	color_count = 5
	region = PolyRegion(domain=[
		LatLon(35°, -10°),
		LatLon(35°, 40°),
		LatLon(71°, 40°),
		LatLon(71°, -10°)
	])
	centers, tiles = generate_tesselation(region, radius, HEX(), EO())
	colors = color_greedy(centers, radius, color_count; min_dist)
	
	@test all(colors[i] != colors[j] || great_circle_distance(centers[i], centers[j]) > min_dist for i in eachindex(centers) for j in (i+1):length(centers))
end
