@testitem "NearestNeighborsExt" begin
    using GeoGrids
    using GeoGrids: authalic_radius
    using GeoGrids.Distances
    using GeoGrids.CoordRefSystems
    using GeoGrids.Meshes
    using NearestNeighbors

    @test GreatCircleMetric() isa Metric

    @test GreatCircleMetric(rand(Point; crs = LatLon)).radius === authalic_radius(rand(LatLon))

    authalic_wgs84 = authalic_radius(LatLon(0,0))

    @test great_circle_distance(LatLon(0, 0), LatLon(rad2deg(1), 0)) ≈ authalic_wgs84

    lat_vec = map(0:90) do lat
        LatLon(lat, 0)
    end
    btree = LatLonTree(lat_vec)

    @test btree.tree isa BallTree

    kdtree = LatLonTree{KDTree}(lat_vec)
    @test kdtree.tree isa KDTree

    # Basic tests for both trees
    for tree in (btree, kdtree)
        idxs, _ = knn(tree, LatLon(0,0), 10, true)
        @test idxs == 1:10

        idxs, _ = knn(tree, lat_vec, 1)
        @test first.(idxs) == eachindex(lat_vec)

        idxs = Int[-1]
        dists = Float64[Inf]
        idx = rand(eachindex(lat_vec))
        knn!(idxs, dists, tree, lat_vec[idx], 1)
        @test idxs == [idx]
        @test dists == [0.0]
    end

    # Tests for only BallTree
    gc = btree.tree.metric
    degm = gc.radius * deg2rad(1)
    _, dists = knn(btree, Point(LatLon(0,0)), 10, true)
    @test all(dists .≈ degm * collect(0:9))

    idxs = Int[]
    inrange!(idxs, btree, Point(LatLon(0,0)), degm * 10 + 1, true)
    @test idxs == collect(1:11)

    @test inrange(btree, LatLon(0,0), degm * .5) == Int[1]
    @test inrange(btree, lat_vec, degm * .5) == map(eachindex(lat_vec)) do i
        [i]
    end

    @test inrangecount(btree, LatLon(0,0), degm - 1) == 1
    @test inrangecount(btree, lat_vec, degm - 1) == ones(length(lat_vec))
end

