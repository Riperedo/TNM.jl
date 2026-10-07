using Test
using TNM

@testset "Data Export and Gnuplot Generator" begin
    # Create simple simulation
    kin = ArrheniusKinetics(-35.0, 40000.0, 0.5)
    model = TNMModel(kin, KWWSpectrum(0.5))
    protocol = LinearRampProtocol(600.0, 400.0, -10.0)
    res = solve(model, protocol, ConvolutionSolver(n_steps=100))

    tmp_dir = mktempdir()
    dat_path = joinpath(tmp_dir, "test_output.dat")
    gnu_path = joinpath(tmp_dir, "test_plot.gnu")

    @testset "Export .dat file" begin
        export_dat(dat_path, res)
        @test isfile(dat_path)

        lines = readlines(dat_path)
        @test length(lines) == 101 # 1 header line + 100 data rows
        # Check that header begins with # prop1 prop2 ...
        @test startswith(lines[1], "# ")
        @test occursin("time_s", lines[1])
        @test occursin("temperature_K", lines[1])
        @test occursin("fictive_temp_K", lines[1])
        @test occursin("cp_normalized", lines[1])
    end

    @testset "Generate Gnuplot Script" begin
        generate_gnuplot_script(dat_path, gnu_path; title="Test Simulation", plot_type=:cp_normalized)
        @test isfile(gnu_path)

        gnu_content = read(gnu_path, String)
        @test occursin("set title 'Test Simulation'", gnu_content)
        @test occursin("plot 'test_output.dat'", gnu_content)
        @test occursin("using 2:4", gnu_content)
    end
end
