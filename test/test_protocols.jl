using Test
using TNM

@testset "Thermal Protocols" begin
    @testset "Linear Ramp" begin
        # Cooling ramp from 600 K to 300 K at 10 K/min = 10/60 K/s
        ramp = linear_ramp_k_min(600.0, 300.0, -10.0)
        @test duration(ramp) == 1800.0 # 30 min = 1800 s
        @test isapprox(temperature_at(ramp, 0.0), 600.0, atol=1e-10)
        @test isapprox(temperature_at(ramp, 1800.0), 300.0, atol=1e-10)
        @test isapprox(temperature_at(ramp, 900.0), 450.0, atol=1e-10)
        @test isapprox(heating_rate_at(ramp, 900.0), -10.0 / 60.0, atol=1e-10)
    end

    @testset "Multi Ramp" begin
        # 600K -> 500K at 20 K/s (dur 5s), then 500K -> 700K at 10 K/s (dur 20s)
        multi = MultiRampProtocol(600.0, [500.0, 700.0], [-20.0, 10.0])
        @test duration(multi) == 25.0
        @test isapprox(temperature_at(multi, 0.0), 600.0, atol=1e-10)
        @test isapprox(temperature_at(multi, 5.0), 500.0, atol=1e-10)
        @test isapprox(temperature_at(multi, 25.0), 700.0, atol=1e-10)
        @test isapprox(temperature_at(multi, 15.0), 600.0, atol=1e-10)
        @test isapprox(heating_rate_at(multi, 2.0), -20.0, atol=1e-10)
        @test isapprox(heating_rate_at(multi, 10.0), 10.0, atol=1e-10)
    end

    @testset "Isothermal Aging Protocol" begin
        # Cool 600K -> 500K at 10 K/s (10s), hold for 100s at 500K, heat 500K -> 650K at 5 K/s (30s)
        proto = IsothermalAgingProtocol(600.0, 500.0, 650.0, -10.0, 100.0, 5.0)
        @test duration(proto) == 140.0
        @test isapprox(temperature_at(proto, 0.0), 600.0, atol=1e-10)
        @test isapprox(temperature_at(proto, 10.0), 500.0, atol=1e-10)
        @test isapprox(temperature_at(proto, 60.0), 500.0, atol=1e-10)
        @test isapprox(temperature_at(proto, 110.0), 500.0, atol=1e-10)
        @test isapprox(temperature_at(proto, 140.0), 650.0, atol=1e-10)

        @test isapprox(heating_rate_at(proto, 5.0), -10.0, atol=1e-10)
        @test isapprox(heating_rate_at(proto, 50.0), 0.0, atol=1e-10)
        @test isapprox(heating_rate_at(proto, 120.0), 5.0, atol=1e-10)
    end

    @testset "Kovacs Crossover Protocol" begin
        kovacs = KovacsCrossoverProtocol(600.0, 480.0, 500.0, 100.0, 500.0; q_quench=100.0)
        @test duration(kovacs) > 600.0
        @test isapprox(temperature_at(kovacs, 0.0), 600.0, atol=1e-10)
        @test isapprox(temperature_at(kovacs, duration(kovacs)), 500.0, atol=1e-10)
    end

    @testset "TMDSC Protocol" begin
        # 500 K to 600 K at 2 K/min, amplitude 0.5 K, period 60 s
        tmdsc = tmdsc_protocol_k_min(500.0, 600.0, 2.0, 0.5, 60.0)
        dur = duration(tmdsc)
        @test isapprox(dur, 3000.0, atol=1e-10) # 100 K / (2/60 K/s) = 3000 s
        
        # Test sinusoidal fluctuation
        T0 = temperature_at(tmdsc, 0.0)
        @test isapprox(T0, 500.0, atol=1e-10)
        # At t = 15s (1/4 period), sin(omega*t) = 1 -> T = 500 + (2/60)*15 + 0.5*1 = 500.5 + 0.5 = 501.0
        T15 = temperature_at(tmdsc, 15.0)
        @test isapprox(T15, 501.0, atol=1e-10)
    end

    @testset "Tabular Protocol" begin
        t_pts = [0.0, 10.0, 30.0, 60.0]
        T_pts = [500.0, 520.0, 540.0, 600.0]
        tab = TabularProtocol(t_pts, T_pts)
        @test duration(tab) == 60.0
        @test isapprox(temperature_at(tab, 0.0), 500.0, atol=1e-10)
        @test isapprox(temperature_at(tab, 20.0), 530.0, atol=1e-10)
        @test isapprox(temperature_at(tab, 60.0), 600.0, atol=1e-10)
    end
end
