 """
```
pseudo_measurement(m::myModel1010{T}, TTT::Matrix{T}, RRR::Matrix{T},
                   CCC::Vector{T}) where {T<:AbstractFloat}
```

Assign pseudo-measurement equation (a linear combination of states):

```
x_t = ZZ_pseudo*s_t + DD_pseudo
```
"""
function pseudo_measurement(m::myModel1010{T},
                            TTT::Matrix{T},
                            RRR::Matrix{T},
                            CCC::Vector{T}) where {T<:AbstractFloat}
    endo      = m.endogenous_states
    endo_addl = m.endogenous_states_augmented
    pseudo    = m.pseudo_observables

    _n_states = n_states_augmented(m)
    _n_pseudo = n_pseudo_observables(m)

    # Compute TTT^10, used for ExpectedAvg10YearRateGap, ExpectedAvg10YearRate, and ExpectedAvg10YearNaturalRate
    TTT5  = (1/20)*((UniformScaling(1.) - TTT)\(UniformScaling(1.) - TTT^20))
    TTT10 = (1/40)*((UniformScaling(1.) - TTT)\(UniformScaling(1.) - TTT^40))
    TTT20 = (1/80)*((UniformScaling(1.) - TTT)\(UniformScaling(1.) - TTT^80))

	# Initialize pseudo ZZ and DD matrices
    ZZ_pseudo = zeros(_n_pseudo, _n_states)
    DD_pseudo = zeros(_n_pseudo)

	##########################################################
	## PSEUDO-OBSERVABLE EQUATIONS
	##########################################################

	## Output
	ZZ_pseudo[pseudo[:y_t],endo[:y_t]] = 1.

	## Flexible Output
	ZZ_pseudo[pseudo[:y_f_t],endo[:y_f_t]] = 1.

	## Output Gap
	ZZ_pseudo[pseudo[:OutputGap],endo[:y_t]] = 1
	ZZ_pseudo[pseudo[:OutputGap],endo[:y_f_t]] = -1

    if subspec(m) in ["ss2", "ss4", "ss22"]
	    ## π_t
	    ZZ_pseudo[pseudo[:π_t],endo[:π_t]] = 1.
	    DD_pseudo[pseudo[:π_t]]            = 100*(m[:π_star]-1)

	    ## Long Run Inflation
	    ZZ_pseudo[pseudo[:LongRunInflation],endo[:π_star_t]] = 1.
	    DD_pseudo[pseudo[:LongRunInflation]]                 = 100.0*(m[:π_star]-1.)

        ## Marginal Cost
        ZZ_pseudo[pseudo[:MarginalCost],endo[:mc_t]] = 1.

	    ## Wages
	    ZZ_pseudo[pseudo[:Wages],endo[:w_t]] = 1.

	    ## Flexible Wages
	    ZZ_pseudo[pseudo[:FlexibleWages],endo[:w_f_t]] = 1.

	    ## z_t
	    ZZ_pseudo[pseudo[:z_t], endo[:z_t]] = 1.
    end

	## Hours
	ZZ_pseudo[pseudo[:Hours],endo[:L_t]] = 1.

    if subspec(m) in ["ss2", "ss4", "ss22"]
	    ## Flexible Hours
	    ZZ_pseudo[pseudo[:FlexibleHours],endo[:L_f_t]] = 1.
    end

	## Natural Rate
	ZZ_pseudo[pseudo[:RealNaturalRate],endo[:r_f_t]] = 1.
	DD_pseudo[pseudo[:RealNaturalRate]]              = 100. * (m[:rstar]-1.)

	## Ex Ante Real Rate
	ZZ_pseudo[pseudo[:ExAnteRealRate],endo[:R_t]]  = 1
	ZZ_pseudo[pseudo[:ExAnteRealRate],endo[:Eπ_t]] = -1
	DD_pseudo[pseudo[:ExAnteRealRate]]             = m[:Rstarn] - 100. * (m[:π_star]-1.)

	## Nominal risk-free Rate
	ZZ_pseudo[pseudo[:NominalRate],endo[:R_t]]  = 1
	DD_pseudo[pseudo[:NominalRate]]             = m[:Rstarn]

	## Ex Ante Real Deposit Rate
	ZZ_pseudo[pseudo[:RealDepositRate],endo[:Rd_t]] = 1
	ZZ_pseudo[pseudo[:RealDepositRate],endo[:Eπ_t]] = -1
	DD_pseudo[pseudo[:RealDepositRate]]             = m[:Rstarn] - 100. * (m[:π_star]-1.)

	## Nominal Deposit Rate
	ZZ_pseudo[pseudo[:NominalDepRate],endo[:Rd_t]]  = 1
	DD_pseudo[pseudo[:NominalDepRate]]             = m[:Rstarn]

	## Nominal FFR
	ZZ_pseudo[pseudo[:NominalFFR], endo[:R_t]] = 1.
	DD_pseudo[pseudo[:NominalFFR]] = m[:Rstarn]

	## Output: GDP
	ZZ_pseudo[pseudo[:gdp], endo[:y_t]]          = 1.0
    ZZ_pseudo[pseudo[:gdp], endo_addl[:y_t1]]    = -1.0
    ZZ_pseudo[pseudo[:gdp], endo[:z_t]]          = 1.0
    DD_pseudo[pseudo[:gdp]]                      = 100*(exp(m[:z_star])-1)

	## Output: GDI
	ZZ_pseudo[pseudo[:gdi], endo[:y_t]]          = m[:γ_gdi]
    ZZ_pseudo[pseudo[:gdi], endo_addl[:y_t1]]    = -m[:γ_gdi]
    ZZ_pseudo[pseudo[:gdi], endo[:z_t]]          = m[:γ_gdi]
    DD_pseudo[pseudo[:gdi]]                      = 100*(exp(m[:z_star])-1) + m[:δ_gdi]

	## PCE Inflation
	ZZ_pseudo[pseudo[:InflationPCE], endo[:π_t]] = 1.0
    DD_pseudo[pseudo[:InflationPCE]]             = 100*(m[:π_star]-1)

	## Inflation (GDP Deflator)
	ZZ_pseudo[pseudo[:GDPdeflator], endo[:π_t]] = m[:Γ_gdpdef]
	DD_pseudo[pseudo[:GDPdeflator]]             = 100*(m[:π_star]-1) + m[:δ_gdpdef]

	## 10-year nominal bond yield
	ZZ_pseudo[pseudo[:LongRate], :] = ZZ_pseudo[pseudo[:NominalFFR], :]' * TTT10
	DD_pseudo[pseudo[:LongRate]]    = m[:Rstarn]

	## 10y inflaton expectations
	ZZ_pseudo[pseudo[:LongInflation], :] = TTT10[endo[:π_t], :]
    DD_pseudo[pseudo[:LongInflation]]    = 100*(m[:π_star]-1)

	## 10-year real yield
	ZZ_pseudo[pseudo[:LongRealRate], :] = ZZ_pseudo[pseudo[:LongRate], :] - ZZ_pseudo[pseudo[:LongInflation], :]
	DD_pseudo[pseudo[:LongRealRate]] = m[:Rstarn] - 100*(m[:π_star]-1)

	## 20 yrs forward transition matrix
    TTT20 = if subspec(m) == "ss16"
        eye(size(TTT,1))
    else
        (1/80)*((UniformScaling(1.) - TTT) \ (UniformScaling(1.) - TTT^80))
    end

    ## BAA or BBB Spread
    ZZ_pseudo[pseudo[:BBBspread], endo[:ERtil_k_t]]    = 1.0
    ZZ_pseudo[pseudo[:BBBspread], endo[:R_t]]          = -1.0
    ZZ_pseudo[pseudo[:BBBspread], :]                   = ZZ_pseudo[pseudo[:BBBspread], :]' * TTT20
    DD_pseudo[pseudo[:BBBspread]]                      = 100*log(m[:spr]*m[:lnb_safe]*m[:lnb_liq])

    ## AAA Spread
    ZZ_pseudo[pseudo[:AAAspread], endo[:b_liq_t]]       = -(1.0-m[:λ_AAA])*(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))/(1 - m[:h]*exp(-m[:z_star]))
    ZZ_pseudo[pseudo[:AAAspread], endo[:ERtil_k_t]]     = m[:λ_AAA]
    ZZ_pseudo[pseudo[:AAAspread], endo[:R_t]]           = -m[:λ_AAA]
    ZZ_pseudo[pseudo[:AAAspread], :]                    = ZZ_pseudo[pseudo[:AAAspread], :]' * TTT20
    DD_pseudo[pseudo[:AAAspread]]                       = 100*log(m[:lnb_liq]) + 100*m[:λ_AAA]*log(m[:spr]*m[:lnb_safe])


    if subspec(m) in ["ss2", "ss4", "ss22"]
        ## Expected Average Nominal Natural Rate
        ZZ_pseudo[pseudo[:ExpectedAvgNominalNaturalRate], endo[:r_f_t]] = 1.
        ZZ_pseudo[pseudo[:ExpectedAvgNominalNaturalRate], endo[:Eπ_t]]  = 1.
        DD_pseudo[pseudo[:ExpectedAvgNominalNaturalRate]]               = m[:Rstarn]

        ## Nominal Rate Gap
        ZZ_pseudo[pseudo[:NominalRateGap], endo[:R_t]]   = 1.
        ZZ_pseudo[pseudo[:NominalRateGap], endo[:r_f_t]] = -1.
        ZZ_pseudo[pseudo[:NominalRateGap], endo[:Eπ_t]]  = -1.
    end

    ## Expected Average 10-Year Real Interest Rate
    ZZ_pseudo[pseudo[:ExpectedAvg10YearRealRate], :] = TTT10[endo[:R_t], :] - TTT10[endo[:Eπ_t], :]
    DD_pseudo[pseudo[:ExpectedAvg10YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

    ## Expected Average 10-Year Real Natural Rate
    ZZ_pseudo[pseudo[:ExpectedAvg10YearRealNaturalRate], :] = TTT10[endo[:r_f_t], :]
    DD_pseudo[pseudo[:ExpectedAvg10YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

    if subspec(m) in ["ss2", "ss4", "ss22"]
        ## Expected Average 10-Year Nominal Interest Rate
        ZZ_pseudo[pseudo[:ExpectedAvg10YearNominalRate], :] = TTT10[endo[:R_t], :]
        DD_pseudo[pseudo[:ExpectedAvg10YearNominalRate]]    = m[:Rstarn]

        ## Expected Average 10-Year Nominal Natural Rate
        ZZ_pseudo[pseudo[:ExpectedAvg10YearNominalNaturalRate], :] = TTT10[endo[:r_f_t], :] + TTT10[endo[:Eπ_t], :]
        DD_pseudo[pseudo[:ExpectedAvg10YearNominalNaturalRate]]    = m[:Rstarn]

        ## Expected Average 10-Year Rate Gap
        ZZ_pseudo[pseudo[:ExpectedAvg10YearRateGap], :] = TTT10[endo[:R_t], :] - TTT10[endo[:r_f_t], :] - TTT10[endo[:Eπ_t], :]
    end

    ## Expected Average 5-Year Real Interest Rate
    ZZ_pseudo[pseudo[:ExpectedAvg5YearRealRate], :] = TTT5[endo[:R_t], :] - TTT5[endo[:Eπ_t], :]
    DD_pseudo[pseudo[:ExpectedAvg5YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

    ## Expected Average 5-Year Real Natural Rate
    ZZ_pseudo[pseudo[:ExpectedAvg5YearRealNaturalRate], :] = TTT5[endo[:r_f_t], :]
    DD_pseudo[pseudo[:ExpectedAvg5YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

    if subspec(m) in ["ss2", "ss4", "ss22"]
        ## Expected Average 5-Year Nominal Interest Rate
        ZZ_pseudo[pseudo[:ExpectedAvg5YearNominalRate], :] = TTT5[endo[:R_t], :]
        DD_pseudo[pseudo[:ExpectedAvg5YearNominalRate]]    = m[:Rstarn]

        ## Expected Average 5-Year Nominal Natural Rate
        ZZ_pseudo[pseudo[:ExpectedAvg5YearNominalNaturalRate], :] = TTT5[endo[:r_f_t], :] + TTT5[endo[:Eπ_t], :]
        DD_pseudo[pseudo[:ExpectedAvg5YearNominalNaturalRate]]    = m[:Rstarn]

        ## Expected Average 5-Year Rate Gap
        ZZ_pseudo[pseudo[:ExpectedAvg5YearRateGap], :] = TTT5[endo[:R_t], :] -
        	TTT5[endo[:r_f_t], :] - TTT5[endo[:Eπ_t], :]
    end

    ## Expected Average 20-Year Real Interest Rate
    ZZ_pseudo[pseudo[:ExpectedAvg20YearRealRate], :] = TTT20[endo[:R_t], :] - TTT20[endo[:Eπ_t], :]
    DD_pseudo[pseudo[:ExpectedAvg20YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

    ## Expected Average 20-Year Real Natural Rate
    ZZ_pseudo[pseudo[:ExpectedAvg20YearRealNaturalRate], :] = TTT20[endo[:r_f_t], :]
    DD_pseudo[pseudo[:ExpectedAvg20YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

    if subspec(m) in ["ss2", "ss4", "ss22"]
        ## Expected Average 20-Year Nominal Interest Rate
        ZZ_pseudo[pseudo[:ExpectedAvg20YearNominalRate], :] = TTT20[endo[:R_t], :]
        DD_pseudo[pseudo[:ExpectedAvg20YearNominalRate]]    = m[:Rstarn]

        ## Expected Average 20-Year Nominal Natural Rate
        ZZ_pseudo[pseudo[:ExpectedAvg20YearNominalNaturalRate], :] = TTT20[endo[:r_f_t], :] + TTT20[endo[:Eπ_t], :]
        DD_pseudo[pseudo[:ExpectedAvg20YearNominalNaturalRate]]    = m[:Rstarn]

        ## Expected Average 20-Year Rate Gap
        ZZ_pseudo[pseudo[:ExpectedAvg20YearRateGap], :] = TTT20[endo[:R_t], :] -
        	TTT20[endo[:r_f_t], :] - TTT20[endo[:Eπ_t], :]

        ## 5-year forward rate gap
        TTT5_fwd = TTT^20
        ZZ_pseudo[pseudo[:Forward5YearRateGap], :] = TTT5_fwd[endo[:R_t], :] -
        	TTT5_fwd[endo[:r_f_t], :] - TTT5_fwd[endo[:Eπ_t], :]

        ## 10-year forward rate gap
        TTT10_fwd = TTT^40
        ZZ_pseudo[pseudo[:Forward10YearRateGap], :] = TTT10_fwd[endo[:R_t], :] - TTT10_fwd[endo[:r_f_t], :] - TTT10_fwd[endo[:Eπ_t], :]
    end

    ## 20-Year Forward Real Rate
    TTT20_fwd = TTT^80
    ZZ_pseudo[pseudo[:Forward20YearRealRate], :] = ZZ_pseudo[pseudo[:ExAnteRealRate], :]' * TTT20_fwd
    DD_pseudo[pseudo[:Forward20YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

    ## 20-Year Real Natural Forward Rate
    ZZ_pseudo[pseudo[:Forward20YearRealNaturalRate], :] = ZZ_pseudo[pseudo[:RealNaturalRate], :]' * TTT20_fwd
    DD_pseudo[pseudo[:Forward20YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

    ## 30-Year Real Forward Rate
    TTT30_fwd = TTT^120
    ZZ_pseudo[pseudo[:Forward30YearRealRate], :] = ZZ_pseudo[pseudo[:ExAnteRealRate], :]' * TTT30_fwd
    DD_pseudo[pseudo[:Forward30YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

    ## 30-Year Real Natural Forward Rate
    ZZ_pseudo[pseudo[:Forward30YearRealNaturalRate], :] = ZZ_pseudo[pseudo[:RealNaturalRate], :]' * TTT30_fwd
    DD_pseudo[pseudo[:Forward30YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

    if !(subspec(m) in ["ss2", "ss4", "ss22"])
        ## Nominal Natural Rate
        ZZ_pseudo[pseudo[:NominalNaturalRate], endo[:r_f_t]] = 1.
        ZZ_pseudo[pseudo[:NominalNaturalRate], endo[:Eπ_t]]  = 1.
        DD_pseudo[pseudo[:NominalNaturalRate]]               = m[:Rstarn]

        ## 5-Year Forward Real Rate
        TTT5_fwd = TTT^20
        ZZ_pseudo[pseudo[:Forward5YearRealRate], :] = ZZ_pseudo[pseudo[:ExAnteRealRate], :]' * TTT5_fwd
        DD_pseudo[pseudo[:Forward5YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

        ## 5-Year Real Natural Forward Rate
        ZZ_pseudo[pseudo[:Forward5YearRealNaturalRate], :] = ZZ_pseudo[pseudo[:RealNaturalRate], :]' * TTT5_fwd
        DD_pseudo[pseudo[:Forward5YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

        ## 10-Year Forward Real Rate
        TTT10_fwd = TTT^40
        ZZ_pseudo[pseudo[:Forward10YearRealRate], :] = ZZ_pseudo[pseudo[:ExAnteRealRate], :]' * TTT10_fwd
        DD_pseudo[pseudo[:Forward10YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)

        ## 10-Year Real Natural Forward Rate
        ZZ_pseudo[pseudo[:Forward10YearRealNaturalRate], :] = ZZ_pseudo[pseudo[:RealNaturalRate], :]' * TTT10_fwd
        DD_pseudo[pseudo[:Forward10YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1)
    end

    if subspec(m) in ["ss13", "ss18"]
        ZZ_pseudo[pseudo[:LiquidityConvenienceYield], endo[:b_liq_t]] = 1.
        DD_pseudo[pseudo[:LiquidityConvenienceYield]] = 100*log(m[:lnb_liq])

        ZZ_pseudo[pseudo[:SafetyConvenienceYield], endo[:b_safe_t]] = 1.
        DD_pseudo[pseudo[:SafetyConvenienceYield]] = 100*log(m[:lnb_safe])

        ZZ_pseudo[pseudo[:SDF], endo[:r_f_t]]    = 1.
        ZZ_pseudo[pseudo[:SDF], endo[:b_liq_t]]  = 1.
        ZZ_pseudo[pseudo[:SDF], endo[:b_safe_t]] = 1.
        DD_pseudo[pseudo[:SDF]] =  m[:Rstarn] - 100*(m[:π_star]-1) + 100*log(m[:lnb_liq]) + 100*log(m[:lnb_safe])

        ZZ_pseudo[pseudo[:Forward20YearLiquidityConvenienceYield], :] = ZZ_pseudo[pseudo[:LiquidityConvenienceYield], :]' * TTT20_fwd
        DD_pseudo[pseudo[:Forward20YearLiquidityConvenienceYield]] = 100*log(m[:lnb_liq])

        ZZ_pseudo[pseudo[:Forward20YearSafetyConvenienceYield], :] = ZZ_pseudo[pseudo[:SafetyConvenienceYield], :]' * TTT20_fwd
        DD_pseudo[pseudo[:Forward20YearSafetyConvenienceYield]] = 100*log(m[:lnb_safe])

        ZZ_pseudo[pseudo[:Forward20YearSDF], :] = ZZ_pseudo[pseudo[:SDF], :]' * TTT20_fwd
        DD_pseudo[pseudo[:Forward20YearSDF]]    = m[:Rstarn] - 100*(m[:π_star]-1) + 100*log(m[:lnb_liq]) + 100*log(m[:lnb_safe])
    end

    ## Exogenous processes
    if subspec(m) in ["ss22"]
        to_add = [:g_t, :b_liq_t, :b_safe_t, :μ_t, :z_t, :λ_f_t, :λ_w_t, :rm_t, :σ_ω_t, :μ_e_t,
                  :γ_t, :π_star_t]
        to_add_addl = [:lr_t, :tfp_t, :e_gdpdef_t, :e_corepce_t, :e_gdp_t, :e_gdi_t, :e_BBB_t, :e_AAA_t]
        for i in to_add
            ZZ_pseudo[pseudo[i], endo[i]] = 1.
        end
        for i in to_add_addl
            ZZ_pseudo[pseudo[i], endo_addl[i]] = 1.
        end
    end

    # Collect indices and transforms
    return PseudoMeasurement(ZZ_pseudo, DD_pseudo)
end
