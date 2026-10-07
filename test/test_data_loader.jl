using Test
using TNM

@testset "Experimental Data Loader" begin
    tmp_dir = mktempdir()
    file1 = joinpath(tmp_dir, "exp_curve1.dat")
    file2 = joinpath(tmp_dir, "exp_curve2.csv")

    # Write synthetic DSC data with comments
    open(file1, "w") do io
        println(io, "# Temperature_K Cp_Normalized")
        println(io, "450.0  0.02")
        println(io, "500.0  0.50")
        println(io, "550.0  1.45")
        println(io, "600.0  1.00")
    end

    open(file2, "w") do io
        println(io, "# Time Temp Cp")
        println(io, "460.0  0.03")
        println(io, "510.0  0.60")
        println(io, "560.0  1.30")
        println(io, "610.0  1.00")
    end

    @testset "Load Single File" begin
        d1 = load_experimental_file(file1; temp_col=1, cp_col=2)
        @test length(d1.temperatures) == 4
        @test d1.temperatures[1] == 450.0
        @test d1.heat_capacities[3] == 1.45
        @test d1.is_normalized == true
    end

    @testset "Load Directory" begin
        datasets = load_experimental_directory(tmp_dir)
        @test length(datasets) == 2
    end
end
