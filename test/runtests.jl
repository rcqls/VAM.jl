using Test
using VAM
@testset "Weibull mle" begin
	data = DataFrame(time=[3.36],type=[-1])
	m = mle(
        @vam(Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))), 
        data
    )
	θ = [0.3,1.8,0.6]
	lnL = -2.30449245951301
	dlnL = [-5.52619699756427,-1.45367181592636,0]
	d2lnL = vcat([[-11.1111111111111,-10.7372278181901,0],[-10.7372278181901,-4.21250787723964,0],[0,0,0]]...)
	c = -1.62415430907299
	dC = [0.555555555555556,0]
	d2C = reshape([-0.308641975308642,0,0,0],2,2)
	# @test logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	# expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	# expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	
    @test contrast(m,θ) == c atol = 0.0000000000001
	# expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	# expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
end