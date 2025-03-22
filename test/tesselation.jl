## HEX
@testitem "GeoRegion HEX Layout Tesselation and :hex Pattern (no EO)" tags = [:tesselation] begin
    # NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test with predefined points like in the cases with :circ.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    samples = [LatLon(43.2693, -8.83691), LatLon(43.6057, -8.11563), LatLon(36.808, -5.06421), LatLon(36.792, -2.3703)] # LatLon{WGS84Latest} coordinates
    corresponding_idxs = [1, 2, 121, 122]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers = generate_tesselation(reg, 40000, HEX())

    @test length(centers) == 122
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1
    end
end

@testitem "GeoRegion HEX Layout Tesselation and :circ Pattern (no EO)" tags = [:tesselation] begin
    samples = [LatLon(43.2693, -8.83691), LatLon(43.6057, -8.11563), LatLon(36.808, -5.06421), LatLon(36.792, -2.3703)] # LatLon{WGS84Latest} coordinates
    corresponding_idxs = [1, 2, 121, 122]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers = generate_tesselation(reg, 40000, HEX(; pattern=:circ))

    @test length(centers) == 122
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1e-4
    end
end

@testitem "GeoRegion HEX Layout Tesselation and :hex Pattern (with EO)" tags = [:tesselation] begin
    # //NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test like the previous cases using a circular tesselation.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers, ngon = generate_tesselation(reg, 40000, HEX(), EO())

    @test length(centers) == 122
    @test length(ngon) == 122

    samplePoints = [LatLon(43.2693, -8.83691), LatLon(43.6057, -8.11563), LatLon(36.808, -5.06421), LatLon(36.792, -2.3703)] # LatLon{WGS84Latest} coordinates
    corresponding_idxs_points = [1, 2, 121, 122]
    sampleNgons = [[LatLon(43.558,-9.1653), LatLon(43.2648,-9.27284), LatLon(42.937,-9.11636), LatLon(42.9797,-8.51291), LatLon(43.2764,-8.40112), LatLon(43.6025,-8.55324), LatLon(43.558,-9.1653)],
    [LatLon(36.4876,-2.64153), LatLon(36.4735,-2.11589), LatLon(36.7875,-1.96157), LatLon(37.0961,-2.09562), LatLon(37.1107,-2.62807), LatLon(36.7972,-2.77902), LatLon(36.4876,-2.64153)]] 
    corresponding_idxs_ngon = [1, 122]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
	centers, ngon = generate_tesselation(reg, 40000, HEX(), EO())

    @test length(centers) == 122
    for i in eachindex(samplePoints) 
        @test abs(get_lat(centers[corresponding_idxs_points[i]]) - samplePoints[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs_points[i]]) - samplePoints[i].lon) < 1    
    end

    @test length(ngon) == 122
    for i in eachindex(sampleNgons)
        for v in eachindex(sampleNgons[i])
            # @test get_lat(ngon[corresponding_idxs_ngon[i]][v]) == sampleNgons[i][v].lat
            @test abs(get_lat(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lat) < 1
            @test abs(get_lon(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lon) < 1
        end
    end
end

@testitem "GeoRegion HEX Layout Tesselation and :circ Pattern (with EO)" tags = [:tesselation] begin
    samplePoints = [LatLon(43.2693, -8.83691), LatLon(43.6057, -8.11563), LatLon(36.808, -5.06421), LatLon(36.792, -2.3703)] # LatLon{WGS84Latest} coordinates
    corresponding_idxs_points = [1, 2, 121, 122]
    sampleNgons = [[LatLon(42.9096, -8.83691), LatLon(42.9271, -8.68509), LatLon(42.9779, -8.5479), LatLon(43.0572, -8.43861), LatLon(43.1572, -8.36791), LatLon(43.2682, -8.34287), LatLon(43.3795, -8.36619), LatLon(43.48, -8.43583), LatLon(43.56, -8.54512), LatLon(43.6113, -8.68338), LatLon(43.629, -8.83691), LatLon(43.6113, -8.99044), LatLon(43.56, -9.12869), LatLon(43.48, -9.23798), LatLon(43.3795, -9.30762), LatLon(43.2682, -9.33094), LatLon(43.1572, -9.3059), LatLon(43.0572, -9.23521), LatLon(42.9779, -9.12591), LatLon(42.9271, -8.98872), LatLon(42.9096, -8.83691), LatLon(42.9096, -8.83691)],
        [LatLon(36.4322, -2.3703), LatLon(36.4498, -2.23211), LatLon(36.5006, -2.10727), LatLon(36.58, -2.00789), LatLon(36.68, -1.94371), LatLon(36.7911, -1.9211), LatLon(36.9024, -1.94247), LatLon(37.0029, -2.00588), LatLon(37.0827, -2.10526), LatLon(37.134, -2.23087), LatLon(37.1517, -2.3703), LatLon(37.134, -2.50974), LatLon(37.0827, -2.63535), LatLon(37.0029, -2.73472), LatLon(36.9024, -2.79814), LatLon(36.7911, -2.8195), LatLon(36.68, -2.7969), LatLon(36.58, -2.73272), LatLon(36.5006, -2.63334), LatLon(36.4498, -2.5085), LatLon(36.4322, -2.3703), LatLon(36.4322, -2.3703)]]
    corresponding_idxs_ngon = [1, 122]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers, ngon = generate_tesselation(reg, 40000, HEX(; pattern=:circ), EO())

    @test length(centers) == 122
    for i in eachindex(samplePoints)
        @test abs(get_lat(centers[corresponding_idxs_points[i]]) - samplePoints[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs_points[i]]) - samplePoints[i].lon) < 1e-4
    end

    @test length(ngon) == 122
    for i in eachindex(sampleNgons)
        for v in eachindex(sampleNgons[i])
            @test abs(get_lat(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lat) < 1e-4
            @test abs(get_lon(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lon) < 1e-4
        end
    end
end

@testitem "GeoRegion HEX Layout Tesselation and :hex Pattern :flat (no EO)" tags = [:tesselation] begin
    # NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test with predefined points like in the cases with :circ.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    samples = [LatLon(42.8428, -9.06506), LatLon(43.4251, -7.82139), LatLon(37.5057, -1.64231), LatLon(36.459, -4.77773),]
    corresponding_idxs = [1, 2, 117, 118]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers = generate_tesselation(reg, 40000, HEX(; direction=:flat))

    @test length(centers) == 118
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1
    end
end

@testitem "PolyRegion HEX Layout Tesselation and :hex Pattern (no EO)" tags = [:tesselation] begin
    # NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test with predefined points like in the cases with :circ.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.    
    samples = [LatLon(79.9691, 5.0), LatLon(79.6441, 7.98881), LatLon(60.0123, 7.14867), LatLon(60.3001, 8.25311)]
    corresponding_idxs = [1, 2, 318, 319]

    reg = PolyRegion(; domain=[LatLon(60, -5), LatLon(80, 0), LatLon(80, 10), LatLon(60, 15)])
    centers = generate_tesselation(reg, 40000, HEX())

    @test length(centers) == 319
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1
    end
end

## ICO
@testitem "GlobalRegion ICO Layout Tesselation and :circ Pattern (no EO)" tags = [:tesselation] begin
    samples = [LatLon(83.991, 0.0), LatLon(79.5826, -137.508), LatLon(-79.5826, 88.6025), LatLon(-83.991, -48.9053)]
    corresponding_idxs = [1, 2, 181, 182]

    reg = GlobalRegion()
    centers = generate_tesselation(reg, 1000000, ICO())

    @test length(centers) == 182
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1e-4
    end
end

@testitem "GlobalRegion ICO Layout Tesselation and :hex Pattern (with EO)" tags = [:tesselation] begin
     # NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test with predefined points like in the cases with :circ.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    samplePoints = [LatLon(83.991, 0.0), LatLon(79.5826, -137.508), LatLon(-79.5826, 88.6025), LatLon(-83.991, -48.9053)]
    corresponding_idxs_points = [1, 2, 181, 182]
    sampleNgons = [[LatLon(66.8411, -2.78182), LatLon(75.8026, 15.7204), LatLon(82.6627, 18.6819), LatLon(83.991, 0.0), LatLon(83.4003, -18.4249), LatLon(66.8411, -2.78182)],
        [LatLon(-82.3641, -67.4582), LatLon(-83.991, -48.9053), LatLon(-83.4003, -30.4803), LatLon(-66.8411, -46.1235), LatLon(-75.8026, -64.6257), LatLon(-82.3641, -67.4582)]]
    corresponding_idxs_ngon = [1, 182]

    reg = GlobalRegion()
    centers, ngon = generate_tesselation(reg, 1000000, ICO(; pattern=:hex), EO())

    @test length(centers) == 182
    for i in eachindex(samplePoints)
        @test abs(get_lat(centers[corresponding_idxs_points[i]]) - samplePoints[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs_points[i]]) - samplePoints[i].lon) < 1
    end

    @test length(ngon) == 182
    for i in eachindex(sampleNgons)
        for v in eachindex(sampleNgons[i])
            @test abs(get_lat(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lat) < 1
            @test abs(get_lon(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lon) < 1
        end
    end
end

@testitem "LatBeltRegion ICO Layout Tesselation and :circ Pattern (no EO)" tags = [:tesselation] begin
    samples = [LatLon(9.89685, 41.5062), LatLon(9.78996, -96.0016), LatLon(-9.78996, 25.0621), LatLon(-9.89685, -112.446)]
    corresponding_idxs = [1, 2, 187, 188]

    reg = LatBeltRegion(lim=(-10, 10))
    centers = generate_tesselation(reg, 400000, ICO())

    @test length(centers) == 188
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1e-4
    end
end

@testitem "LatBeltRegion ICO Layout Tesselation and :hex Pattern (with EO)" tags = [:tesselation] begin
    # //NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test like the previous cases using a circular tesselation.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    samplePoints = [LatLon(9.89685,41.5062), LatLon(9.78996,-96.0016), LatLon(-9.78996,25.0621), LatLon(-9.89685,-112.446)]
    corresponding_idxs_points = [1, 2, 187, 188]

    sampleNgons = [[LatLon(11.8183,45.104),LatLon(13.921,40.8336),LatLon(9.54199,37.4515),LatLon(8.05959,37.8844),LatLon(5.89612,42.1999),LatLon(10.2875,45.5542),LatLon(11.8183,45.104)],
    [LatLon(-11.8183,-116.044),LatLon(-13.921,-111.773),LatLon(-9.54199,-108.391),LatLon(-8.05959,-108.824),LatLon(-5.89611,-113.139),LatLon(-10.2875,-116.494),LatLon(-11.8183,-116.044)]]
    corresponding_idxs_ngon = [1, 188]

    reg = LatBeltRegion(lim=(-10,10))
	centers,ngon = generate_tesselation(reg, 400000, ICO(;pattern=:hex), EO())

    @test length(centers) == 188
    for i in eachindex(samplePoints) 
        @test abs(get_lat(centers[corresponding_idxs_points[i]]) - samplePoints[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs_points[i]]) - samplePoints[i].lon) < 1  
    end

    @test length(ngon) == 188
    for i in eachindex(sampleNgons)
        for v in eachindex(sampleNgons[i])
            @test abs(get_lat(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lat) < 1
            @test abs(get_lon(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lon) < 1
        end
    end
end

@testitem "PolyRegion ICO Layout Tesselation and :circ Pattern (no EO)" tags = [:tesselation] begin
    samples = [LatLon(79.9456, 4.99995), LatLon(79.6831, 9.73597), LatLon(60.0598, -3.47868), LatLon(60.0243, 8.92039)]
    corresponding_idxs = [1, 2, 218, 219]

    reg = PolyRegion(; domain=[LatLon(60, -5), LatLon(80, 0), LatLon(80, 10), LatLon(60, 15)])
    centers = generate_tesselation(reg, 40000, ICO(; correction=1.7))

    @test length(centers) == 219
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1e-4
    end
end

@testitem "PolyRegion ICO Layout Tesselation and :hex Pattern (with EO)" tags = [:tesselation] begin
  # //NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test like the previous cases using a circular tesselation.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    samplePoints = [LatLon(79.9456,4.99995), LatLon(79.6831,9.73597), LatLon(60.0598,-3.47868), LatLon(60.0243,8.92039)]
    corresponding_idxs_points = [1, 2, 218, 219]
    sampleNgons = [[LatLon(79.2839,5.7835),LatLon(80.0631,6.07365),LatLon(81.0731,4.99257),LatLon(80.8394,4.30809),LatLon(79.929,3.94491),LatLon(78.9239,4.90667),LatLon(79.2839,5.7835)],
        [LatLon(60.0209,9.55066),LatLon(60.4379,9.40303),LatLon(60.6574,8.86491),LatLon(60.0377,8.28825),LatLon(59.6336,8.42962),LatLon(59.4002,8.98413),LatLon(60.0209,9.55066)]]
    corresponding_idxs_ngon = [1, 219]
    
    reg = PolyRegion(;domain=[LatLon(60,-5), LatLon(80,0), LatLon(80,10), LatLon(60,15)])
	centers, ngon = generate_tesselation(reg, 40000, ICO(;correction=1.7, pattern=:hex), EO())

    @test length(centers) == 219
    for i in eachindex(samplePoints) 
        @test abs(get_lat(centers[corresponding_idxs_points[i]]) - samplePoints[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs_points[i]]) - samplePoints[i].lon) < 1  
    end

    @test length(ngon) == 219
    for i in eachindex(sampleNgons)
        for v in eachindex(sampleNgons[i])
            @test abs(get_lat(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lat) < 1
            @test abs(get_lon(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lon) < 1
        end
    end
end

@testitem "GeoRegion ICO Layout Tesselation and :circ Pattern (no EO)" tags = [:tesselation] begin
    samples = [LatLon(43.5772, -7.50902), LatLon(43.4468, -5.70002), LatLon(36.6657, -3.58858), LatLon(36.285, -5.82463)]
    corresponding_idxs = [1, 2, 103, 104]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers = generate_tesselation(reg, 40000, ICO())

    @test length(centers) == 104
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1e-4
    end
end

@testitem "GeoRegion ICO Layout Tesselation and :hex Pattern (with EO)" tags = [:tesselation] begin
    samplePoints = [LatLon(43.5772, -7.50902), LatLon(43.4468, -5.70002), LatLon(36.6657, -3.58858), LatLon(36.285, -5.82463)]
    corresponding_idxs_points = [1, 2, 103, 104]
    sampleNgons = [[LatLon(43.2175, -7.50902), LatLon(43.235, -7.35644), LatLon(43.2858, -7.21855), LatLon(43.3651, -7.1087), LatLon(43.4651, -7.03764), LatLon(43.5761, -7.01247), LatLon(43.6874, -7.03589), LatLon(43.7879, -7.10588), LatLon(43.8679, -7.21573), LatLon(43.9192, -7.35469), LatLon(43.9369, -7.50902), LatLon(43.9192, -7.66334), LatLon(43.8679, -7.8023), LatLon(43.7879, -7.91215), LatLon(43.6874, -7.98214), LatLon(43.5761, -8.00557), LatLon(43.4651, -7.9804), LatLon(43.3651, -7.90933), LatLon(43.2858, -7.79948), LatLon(43.235, -7.6616), LatLon(43.2175, -7.50902), LatLon(43.2175, -7.50902)],
        [LatLon(35.9252, -5.82463), LatLon(35.9428, -5.68733), LatLon(35.9937, -5.56329), LatLon(36.073, -5.46457), LatLon(36.1731, -5.40081), LatLon(36.2841, -5.37836), LatLon(36.3954, -5.3996), LatLon(36.4959, -5.46261), LatLon(36.5757, -5.56134), LatLon(36.627, -5.68612), LatLon(36.6447, -5.82463), LatLon(36.627, -5.96314), LatLon(36.5757, -6.08792), LatLon(36.4959, -6.18665), LatLon(36.3954, -6.24966), LatLon(36.2841, -6.27089), LatLon(36.1731, -6.24845), LatLon(36.073, -6.18469), LatLon(35.9937, -6.08597), LatLon(35.9428, -5.96193), LatLon(35.9252, -5.82463), LatLon(35.9252, -5.82463)]]
    corresponding_idxs_ngon = [1, 104]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers, ngon = generate_tesselation(reg, 40000, ICO(), EO())

    @test length(centers) == 104
    for i in eachindex(samplePoints)
        @test abs(get_lat(centers[corresponding_idxs_points[i]]) - samplePoints[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs_points[i]]) - samplePoints[i].lon) < 1e-4
    end

    @test length(ngon) == 104
    for i in eachindex(sampleNgons)
        for v in eachindex(sampleNgons[i])
            @test abs(get_lat(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lat) < 1e-4
            @test abs(get_lon(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lon) < 1e-4
        end
    end
end

@testitem "Misc coverage" begin
    using GeoGrids: gen_circle_pattern, get_lat, _wrap_latlon

    # We test wrapping around poles
    pts = gen_circle_pattern(LatLon(90, 0), 1000e3) |> first
    lat = get_lat(first(pts))
    @test all(p -> get_lat(p) ≈ lat, pts)

    pts = gen_circle_pattern(LatLon(-90, 0), 1000e3) |> first
    lat = get_lat(first(pts))
    @test all(p -> get_lat(p) ≈ lat, pts)

    @test _wrap_latlon(100, 0) == (80, 180)
    @test _wrap_latlon(-100, 0) == (-80, 180)
end

@testitem "Test missing Functions" tags = [:tesselation] begin
    reg = GlobalRegion()
    @test_throws "H3 tassellation is not yet implemented in this version..." generate_tesselation(reg, 10.0, H3())

    @test_throws ArgumentError generate_tesselation(GlobalRegion(), 10.0, HEX())
end