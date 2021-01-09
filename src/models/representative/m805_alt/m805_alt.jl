"""
Model805_alt{T} <: AbstractDSGEModel{T}

The Model805_alt type defines the structure of the Smets & Wouters (2007) model augmented
with long run inflation expetations. We can then concisely pass around a Model object
to the remaining steps of the model (solve, estimate, and forecast).

### Fields

#### Parameters and Steady-States
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model parameters.

* `steady_state::Vector`: Model steady-state values, computed as a
function of elements of `parameters`.

* `keys::OrderedDict{Symbol,Int}`: Maps human-readable names for all model
parameters and steady-states to their indices in `parameters` and
`steady_state`.

#### Inputs to Measurement and Equilibrium Condition Equations

The following fields are dictionaries that map human-readible names to
row and column indices in the matrix representations of of the
measurement equation and equilibrium conditions.

* `endogenous_states::OrderedDict{Symbol,Int}`: Maps each state to a column
in the measurement and equilibrium condition matrices.

* `exogenous_shocks::OrderedDict{Symbol,Int}`: Maps each shock to a column in
the measurement and equilibrium condition matrices.

* `expected_shocks::OrderedDict{Symbol,Int}`: Maps each expected shock to a
column in the measurement and equilibrium condition matrices.

* `equilibrium_conditions::OrderedDict{Symbol,Int}`: Maps each equlibrium
condition to a row in the model's equilibrium condition matrices.

* `endogenous_states_augmented::OrderedDict{Symbol,Int}`: Maps lagged states
to their columns in the measurement and equilibrium condition
equations. These are added after Gensys solves the model.

* `observables::OrderedDict{Symbol,Int}`: Maps each observable to a row in
the model's measurement equation matrices.

* `pseudo_observables::OrderedDict{Symbol,Int}`: Maps each pseudo-observable to
  a row in the model's pseudo-measurement equation matrices.

#### Model Specifications and Settings

* `spec::String`: The model specification identifier, \"smets_wouters\",
cached here for filepath computation.

* `subspec::String`: The model subspecification number,
indicating that some parameters from the original model spec (\"ss0\")
are initialized differently. Cached here for filepath computation.

* `settings::Dict{Symbol,Setting}`: Settings/flags that affect
computation without changing the economic or mathematical setup of
the model.

* `test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

#### Other Fields

* `rng::MersenneTwister`: Random number generator. By default, it is
seeded to ensure replicability in algorithms that involve randomness
(such as Metropolis-Hastings).

* `testing::Bool`: Indicates whether the model is in testing mode. If
`true`, settings from `m.test_settings` are used in place of those in
`m.settings`.

* `observable_mappings::OrderedDict{Symbol,Observable}`: A dictionary that
  stores data sources, series mnemonics, and transformations to/from model units.
  DSGE.jl will fetch data from the Federal Reserve Bank of
  St. Louis's FRED database; all other data must be downloaded by the
  user. See `load_data` and `Observable` for further details.

* `pseudo_observable_mappings::OrderedDict{Symbol,PseudoObservable}`: A
  dictionary that stores names and transformations to/from model units. See
  `PseudoObservable` for further details.
"""
mutable struct Model805_alt{T} <: AbstractDSGEModel{T}
    parameters::ParameterVector{T}                  # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                # model steady-state values
    keys::OrderedDict{Symbol,Int}                   # human-readable names for all the model
                                                    # parameters and steady-n_states

    endogenous_states::OrderedDict{Symbol,Int}            # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}             # measurement and equilibrium condition equations.
    expected_shocks::OrderedDict{Symbol,Int}              #
    equilibrium_conditions::OrderedDict{Symbol,Int}       #
    endogenous_states_augmented::OrderedDict{Symbol,Int}  #
    observables::OrderedDict{Symbol,Int}                  #
    pseudo_observables::OrderedDict{Symbol,Int}           #

    spec::String                                    # Model specification number (eg "m805")
    subspec::String                                 # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                  # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}             # Settings/flags for testing mode
    rng::MersenneTwister                            # Random number generator
    testing::Bool                                   # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
    pseudo_observable_mappings::OrderedDict{Symbol, PseudoObservable}
end

description(m::Model805_alt) = "Smets-Wouters model, with long run inflation expectations"

"""
`init_model_indices!(m::Model805_alt)`

Arguments:
`m:: Model805_alt`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::Model805_alt)
    # Endogenous states
    endogenous_states = [[
        :y_t, :c_t, :i_t, :qk_t, :k_t, :kbar_t, :u_t, :rk_t, :mc_t,
        :π_t, :μ_ω_t, :w_t, :L_t, :R_t, :g_t, :b_t, :μ_t, :z_t,
        :λ_f_t, :λ_f_t1, :λ_w_t, :λ_w_t1, :rm_t, :π_star_t, :Ec_t, :Eqk_t, :Ei_t,
        :Eπ_t, :EL_t, :Erk_t, :Ew_t, :y_f_t, :c_f_t,
        :i_f_t, :qk_f_t, :k_f_t, :kbar_f_t, :u_f_t, :rk_f_t, :w_f_t,
        :L_f_t, :r_f_t, :Ec_f_t, :Eqk_f_t, :Ei_f_t, :EL_f_t, :Erk_f_t, :ztil_t,
        :b_til_t, :b_p_t, :zp_t, :Rd_t, :Ez_t];
        [Symbol("rm_tl$i") for i = 1:n_anticipated_shocks(m)]]

    # Exogenous shocks
    exogenous_shocks = [[
        :g_sh, :b_til_sh, :b_p_sh, :μ_sh, :z_sh, :λ_f_sh, :λ_w_sh, :rm_sh,
        :π_star_sh, :zp_sh, :tfp_sh, :gdpdef_sh, :corepce_sh, :gdp_sh, :gdi_sh, :lr_sh];
        [Symbol("rm_shl$i") for i = 1:n_anticipated_shocks(m)]]

    # Expectations shocks
    expected_shocks = [
        :Ec_sh, :Eqk_sh, :Ei_sh, :Eπ_sh, :EL_sh, :Erk_sh, :Ew_sh, :Ec_f_sh,
        :Eqk_f_sh, :Ei_f_sh, :EL_f_sh, :Erk_f_sh]

    # Equilibrium conditions
    equilibrium_conditions = [[
        :eq_euler, :eq_inv, :eq_capval, :eq_output, :eq_caputl, :eq_capsrv, :eq_capev,
        :eq_mkupp, :eq_phlps, :eq_caprnt, :eq_msub, :eq_wage, :eq_mp, :eq_res, :eq_g, :eq_b, :eq_μ, :eq_z,
        :eq_λ_f, :eq_λ_w, :eq_rm, :eq_λ_f1, :eq_λ_w1, :eq_Ec,
        :eq_Eqk, :eq_Ei, :eq_Eπ, :eq_EL, :eq_Erk, :eq_Ew, :eq_euler_f, :eq_inv_f,
        :eq_capval_f, :eq_output_f, :eq_caputl_f, :eq_capsrv_f, :eq_capev_f, :eq_mkupp_f, :eq_caprnt_f, :eq_msub_f,
        :eq_res_f, :eq_Ec_f, :eq_Eqk_f, :eq_Ei_f, :eq_EL_f, :eq_Erk_f, :eq_ztil, :eq_π_star,
        :eq_b_til, :eq_b_p, :eq_zp, :eq_Ez, :eq_dep];
        [Symbol("eq_rml$i") for i=1:n_anticipated_shocks(m)]]

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = [
        :y_t1, :c_t1, :i_t1, :w_t1, :π_t1, :L_t1, :Et_π_t, :lr_t, :tfp_t, :e_gdpdef_t,
        :e_corepce_t, :e_gdp_t, :e_gdi_t, :e_gdp_t1, :e_gdi_t1, :u_t1]

    # Observables
    observables = keys(m.observable_mappings)

    # Pseudo-observables
    pseudo_observables = keys(m.pseudo_observable_mappings)

    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(exogenous_shocks);            m.exogenous_shocks[k]            = i end
    for (i,k) in enumerate(expected_shocks);             m.expected_shocks[k]             = i end
    for (i,k) in enumerate(equilibrium_conditions);      m.equilibrium_conditions[k]      = i end
    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(observables);                 m.observables[k]                 = i end
    for (i,k) in enumerate(pseudo_observables);          m.pseudo_observables[k]          = i end
end

function Model805_alt(subspec::String="ss20";
                      custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                      testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)        # Random Number Generator

    # initialize empty model
    m = Model805_alt{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(), OrderedDict{Symbol,Int}(),

            # model indices
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}(),
            OrderedDict{Symbol,PseudoObservable}())

    # Set settings
    model_settings!(m)
    default_test_settings!(m)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Set observable and pseudo-observable transformations
    init_observable_mappings!(m)

    # Initialize parameters
    init_parameters!(m)
    init_model_indices!(m)
    init_subspec!(m)
    steadystate!(m)
    return m
end

"""
```
init_parameters!(m::Model805_alt)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::Model805_alt)
    m <= parameter(:α,      0.1687, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     Normal(0.30, 0.05),         fixed=false,
                   description="α: Capital elasticity in the intermediate goods sector's Cobb-Douglas production function.",
                   tex_label="\\alpha")

    m <= parameter(:ζ_p,   0.7467, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.1),          fixed=false,
                   description="ζ_p: The Calvo parameter. In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")

    m <= parameter(:ι_p,   0.2684, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.15),         fixed=false,
                   description="ι_p: The weight attributed to last period's inflation in price indexation. (1-ι_p) is the weight attributed to steady-state inflation.",

                   tex_label="\\iota_p")

    m <= parameter(:δ,      0.025,  fixed=true,
                   description="δ: The capital depreciation rate.", tex_label="\\delta" )


    m <= parameter(:Upsilon,  1.000,  (0., 10.),     (1e-5, 0.),      ModelConstructors.Exponential(),    GammaAlt(1., 0.5),          fixed=true,
                   description="Υ: The trend evolution of the price of investment goods relative to consumption goods. Set equal to 1.",
                   tex_label="\\mathcal{\\Upsilon}")

    m <= parameter(:Φ,   1.4933, (1., 10.),     (1.00, 10.00),   ModelConstructors.Exponential(),    Normal(1.25, 0.12),         fixed=false,
                   description="Φ: Fixed costs.",
                   tex_label="\\Phi")

    m <= parameter(:S′′,       3.3543, (-15., 15.),   (-15., 15.),     Untransformed(),  Normal(4., 1.5),            fixed=false,
                   description="S'': The second derivative of households' cost of adjusting investment.", tex_label="S\\prime\\prime")

    m <= parameter(:h,        0.4656, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.7, 0.1),          fixed=false,
                   description="h: Consumption habit persistence.", tex_label="h")

    m <= parameter(:ppsi,     0.7614, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.15),         fixed=false,
                   description="ppsi: Utilization costs.", tex_label="ppsi")

    m <= parameter(:ν_l,     1.0647, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    Normal(2, 0.75),            fixed=false,
                   description="ν_l: The coefficient of relative risk aversion on the labor term of households' utility function.", tex_label="\\nu_l")

    m <= parameter(:ζ_w,   0.7922, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.1),          fixed=false,
                   description="ζ_w: (1-ζ_w) is the probability with which households can freely choose wages in each period. With probability ζ_w, wages increase at a geometrically weighted average of the steady state rate of wage increases and last period's productivity times last period's inflation.",
                   tex_label="\\zeta_w")

    m <= parameter(:ι_w,   0.5729, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.15),         fixed=false,
                   description="ι_w: No description available.",
                   tex_label="\\iota_w")

    m <= parameter(:λ_w,      1.5000,                                                                               fixed=true,
                   description="λ_w: The wage markup, which affects the elasticity of substitution between differentiated labor services.",
                   tex_label="\\lambda_w")

    m <= parameter(:β, 0.7420, (1e-7, 10.),   (1e-7, 10.),     ModelConstructors.Exponential(),    GammaAlt(0.25, 0.1),        fixed=false,  scaling = x -> 1/(1 + x/100),
                   description="β: Discount rate.",
                   tex_label="\\beta ")

    m <= parameter(:ψ1,  1.8678, (1e-5, 10.),   (1e-5, 10.00),   ModelConstructors.Exponential(),    Normal(1.5, 0.25),          fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2,  0.0715, (-0.5, 0.5),   (-0.5, 0.5),     Untransformed(),  Normal(0.12, 0.05),         fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3, 0.2131, (-0.5, 0.5),   (-0.5, 0.5),     Untransformed(),  Normal(0.12, 0.05),         fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:π_star,   0.6231, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    GammaAlt(0.75, 0.4),        fixed=false,  scaling = x -> 1 + x/100,
                   description="π_star: The steady-state rate of inflation.",
                   tex_label="\\pi_*")

    m <= parameter(:σ_c, 1.5073, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    Normal(1.5, 0.37),          fixed=false,
                   description="σ_c: No description available.",
                   tex_label="\\sigma_{c}")

    m <= parameter(:ρ,      .8519, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.75, 0.10),        fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho")

    m <= parameter(:ϵ_p,     10.000,                                                                               fixed=true,
                   description="ϵ_p: No description available.",
                   tex_label="\\varepsilon_{p}")

    m <= parameter(:ϵ_w,     10.000,                                                                               fixed=true,
                   description="ϵ_w: No description available.",
                   tex_label="\\varepsilon_{w}")


    # exogenous processes - level
    m <= parameter(:γ,      0.3085, (-5.0, 5.0),     (-5., 5.),     Untransformed(), Normal(0.4, 0.1),            fixed=false, scaling = x -> x/100,
                   description="γ: The log of the steady-state growth rate of technology.",
                   tex_label="\\gamma")

    m <= parameter(:Lmean,  -45., (-1000., 1000.), (-1e3, 1e3),   Untransformed(), Normal(-45, 5),   fixed=false,
                   description="Lmean: No description available.",
                   tex_label="Lmean")

    m <= parameter(:g_star,    0.1800,                                                                               fixed=true,
                   description="g_star: No description available.",
                   tex_label="g_*")

    # exogenous processes - autocorrelation
    m <= parameter(:ρ_g,      0.9930, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label="\\rho_g")

    m <= parameter(:ρ_b_til, 0.9410, (0.0, 1.0), (0.0, 1.0), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                  description="ρ_b_til: AR(1) coefficient in the non-permanent component of the intertemporal preference shifter process for safe/liquid assets.",
                  tex_label="\\rho_{\\tilde{b}}")

    m <= parameter(:ρ_b_p, 0.99, (0.0, 1.0), (0.0, 1.0), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=true,
                  description="ρ_b_p: AR(1) coefficient in the permanent component of the intertemporal preference shifter process for safe/liquid assets.",
                  tex_label="\\rho_{b^p}")

    m <= parameter(:ρ_μ,     0.7790, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   tex_label="\\rho_{\\mu}")

    m <= parameter(:ρ_z,      0.9676, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_z")

    m <= parameter(:ρ_λ_f,    0.8372, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_λ_f: AR(1) coefficient in the price mark-up shock process.",
                   tex_label="\\rho_{\\lambda_f}")

    m <= parameter(:ρ_λ_w,    0.9853, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.", # CHECK THIS
                   tex_label="\\rho_{\\lambda_w}")

    m <= parameter(:ρ_rm,     0.9000, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_rm: AR(1) coefficient in the monetary policy shock process.", # CHECK THIS
                   tex_label="\\rho_{rm}")

    m <= parameter(:ρ_π_star,     0.9900, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),
                  fixed=true,
                  description="ρ_π_star: AR(1) coefficient in the π* process.",
                  tex_label="\\rho_{\\pi^*}")

                  m <= parameter(:ρ_lr, 0.6936, (0.0, 1.0), (0.0, 1.0), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                                 tex_label="\\rho_{10y}")

    m <= parameter(:ρ_z_p, 0.8910, (0.0, 1.0), (0.0, 1.0), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                  description="ρ_z_p: AR(1) coefficient in the process describing the permanent component of productivity.",
                  tex_label="\\rho_{z^p}")

    m <= parameter(:ρ_tfp, 0.1953, (0.0, 1.0), (0.0, 1.0), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                                 tex_label="\\rho_{tfp}")

    m <= parameter(:ρ_gdpdef, 0.5379, (0.0, 1.0), (0.0, 1.0), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                                 tex_label="\\rho_{gdpdef}")

    m <= parameter(:ρ_corepce, 0.2320, (0.0, 1.0), (0.0, 1.0), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                                 tex_label="\\rho_{pce}")

    m <= parameter(:ρ_gdp, 0., (-1.0, 1.0), (-0.999, 0.999), ModelConstructors.SquareRoot(), Normal(0.0, 0.2), fixed=false,
                                 tex_label="\\rho_{gdp}")

    m <= parameter(:ρ_gdi, 0., (-0.999, 0.999), (-0.999, 0.999), ModelConstructors.SquareRoot(), Normal(0.0, 0.2), fixed=false,
                                 tex_label="\\rho_{gdi}")

    m <= parameter(:ρ_gdpvar, 0., (-0.999, 0.999), (-0.999, 0.999), ModelConstructors.SquareRoot(), Normal(0.0, 0.4), fixed=false,
                                 tex_label="\\varrho_{gdp}")

    m <= parameter(:me_level, 1., (-0.999, 0.999), (-0.999, 0.999), ModelConstructors.Untransformed(), Normal(0.0, 0.4), fixed=true,
                                 description="me_level: Indicator of cointegration of GDP and GDI.",
                                 tex_label="\\mathcal{C}_{me}")


    # exogenous processes - standard deviation
    m <= parameter(:σ_g,      0.6090, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_g: The standard deviation of the government spending process.",
                   tex_label="\\sigma_{g}")

    m <= parameter(:σ_b_til, 0.0292, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                                  description="σ_b_til: Standard deviation of stationary component of asset preference shifter process.",
                                  tex_label="\\sigma_{\\tilde{b}}")

    m <= parameter(:σ_b_p, 0.0269, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(6, 0.03), fixed=false,
                                  description="σ_b_p: Standard deviation of non-stationary component of asset preference shifter process.",
                                  tex_label="\\sigma_{b^p}")

    m <= parameter(:σ_μ,     0.4601, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label="\\sigma_{\\mu}")

    m <= parameter(:σ_z,      0.4618, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_z: No description available.",
                   tex_label="\\sigma_{z}")

    m <= parameter(:σ_λ_f,    0.1455, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_λ_f: The mean of the process that generates the price elasticity of the composite good.  Specifically, the elasticity is (1+λ_{f,t})/(λ_{f_t}).",
                   tex_label="\\sigma_{\\lambda_f}")

    m <= parameter(:σ_λ_w,    0.2089, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_λ_w: No description available.",
                   tex_label="\\sigma_{\\lambda_w}")

    m <= parameter(:σ_rm,     0.2397, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_rm: No description available.",
                   tex_label="\\sigma_{rm}")

    m <= parameter(:σ_π_star,      0.0300, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_π_star: The standard deviation of the π* process.",
                   tex_label="\\sigma_{\\pi^*}")

                   m <= parameter(:σ_lr, 0.1766, (1e-8,10.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.75), fixed=false,
                                  tex_label="\\sigma_{10y}")

    m <= parameter(:σ_z_p, 0.1662, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                                  description="σ_z_p: The standard deviation of the shock to the permanent component of productivity.",
                                  tex_label="\\sigma_{z^p}")

    m <= parameter(:σ_tfp, 0.9391, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                                  tex_label="\\sigma_{tfp}")

    m <= parameter(:σ_gdpdef, 0.1575, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                                  tex_label="\\sigma_{gdpdef}")

    m <= parameter(:σ_corepce, 0.0999, (1e-8, 5.),(1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                                  tex_label="\\sigma_{pce}")

    m <= parameter(:σ_gdp, 0.1, (1e-8, 5.),(1e-8, 5.),ModelConstructors.Exponential(),RootInverseGamma(2, 0.10), fixed=false,
                                  tex_label="\\sigma_{gdp}")

    m <= parameter(:σ_gdi, 0.1, (1e-8, 5.),(1e-8, 5.),ModelConstructors.Exponential(),RootInverseGamma(2, 0.10), fixed=false,
                                  tex_label="\\sigma_{gdi}")

                                  # standard deviations of the anticipated policy shocks
                                  for i = 1:n_anticipated_shocks_padding(m)
                                      if i < 13
                                          m <= parameter(Symbol("σ_r_m$i"), .2, (1e-7, 100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, .2), fixed=false,
                                                         description="σ_r_m$i: Standard deviation of the $i-period-ahead anticipated policy shock.",
                                                         tex_label=@sprintf("\\sigma_{%d,r}",i))
                                      else
                                          m <= parameter(Symbol("σ_r_m$i"), .0, (1e-7, 100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, .2), fixed=true,
                                                         description="σ_r_m$i: Standard deviation of the $i-period-ahead anticipated policy shock.",
                                                         tex_label=@sprintf("\\sigma_{%d,r}",i))
                                      end
                                  end

    m <= parameter(:η_gz,       0.5632, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.50, 0.20), fixed=false,
                   description="η_gz: Correlate g and z shocks.",
                   tex_label="\\eta_{gz}")

    m <= parameter(:η_λ_f,      0.7652, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.50, 0.20),         fixed=false,
                   description="η_λ_f: No description available.",
                   tex_label="\\eta_{\\lambda_f}")

    m <= parameter(:η_λ_w,      0.8936, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.50, 0.20),         fixed=false,
                   description="η_λ_w: AR(2) coefficient on wage markup shock process.",
                   tex_label="\\eta_{\\lambda_w}")

                   m <= parameter(:Γ_gdpdef, 1.0354, (-10., 10.), (-10., -10.), ModelConstructors.Untransformed(), Normal(1.00, 2.), fixed=false,
                                  tex_label="\\gamma_{gdpdef}")

                   m <= parameter(:δ_gdpdef, 0.0181, (-10., 10.), (-10., -10.), ModelConstructors.Untransformed(), Normal(0.00, 2.), fixed=false,
                                  tex_label="\\delta_{gdpdef}")

                   m <= parameter(:γ_gdi, 1., (-10., 10.), (-10., -10.), ModelConstructors.Untransformed(), Normal(1., 2.), fixed=true,
                                  tex_label="\\gamma_{gdi}")

                   m <= parameter(:δ_gdi, 0., (-10., 10.), (-10., -10.), ModelConstructors.Untransformed(), Normal(0.00, 2.), fixed=true,
                                  tex_label="\\delta_{gdi}")

    # steady states
    m <= SteadyStateParameter(:zstar,  NaN, description="steady-state growth rate of productivity", tex_label="\\z_*")
    m <= SteadyStateParameter(:rstar,   NaN, description="steady-state something something", tex_label="\\r_*")
    m <= SteadyStateParameter(:Rstarn,  NaN, description="steady-state something something", tex_label="\\R_*_n")
    m <= SteadyStateParameter(:rkstar,  NaN, description="steady-state something something", tex_label="\\BLAH")
    m <= SteadyStateParameter(:wstar,   NaN, description="steady-state something something", tex_label="\\w_*")
    m <= SteadyStateParameter(:Lstar,   NaN, description="steady-state something something", tex_label="\\L_*")
    m <= SteadyStateParameter(:kstar,   NaN, description="Effective capital that households rent to firms in the steady state.", tex_label="\\k_*")
    m <= SteadyStateParameter(:kbarstar, NaN, description="Total capital owned by households in the steady state.", tex_label="\\bar{k}_*")
    m <= SteadyStateParameter(:istar,  NaN, description="Detrended steady-state investment", tex_label="\\i_*")
    m <= SteadyStateParameter(:ystar,  NaN, description="steady-state something something", tex_label="\\y_*")
    m <= SteadyStateParameter(:cstar,  NaN, description="steady-state something something", tex_label="\\c_*")
    m <= SteadyStateParameter(:wl_c,   NaN, description="steady-state something something", tex_label="\\wl_c")
end


"""
```
steadystate!(m::Model805_alt)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever the parameters of `m` are updated.
"""
function steadystate!(m::Model805_alt)
    m[:zstar]    = log(1+m[:γ]) + m[:α]/(1-m[:α])*log(m[:Upsilon])
    m[:rstar]    = exp(m[:σ_c]*m[:zstar]) / m[:β]
    m[:Rstarn]   = 100*(m[:rstar]*m[:π_star] - 1)
    m[:rkstar]   = m[:rstar]*m[:Upsilon] - (1-m[:δ])
    m[:wstar]    = (m[:α]^m[:α] * (1-m[:α])^(1-m[:α]) * m[:rkstar]^(-m[:α]) / m[:Φ])^(1/(1-m[:α]))
    m[:Lstar]    = 1.
    m[:kstar]    = (m[:α]/(1-m[:α])) * m[:wstar] * m[:Lstar] / m[:rkstar]
    m[:kbarstar] = m[:kstar] * (1+m[:γ]) * m[:Upsilon]^(1 / (1-m[:α]))
    m[:istar]    = m[:kbarstar] * (1-((1-m[:δ])/((1+m[:γ]) * m[:Upsilon]^(1/(1-m[:α])))))
    m[:ystar]    = (m[:kstar]^m[:α]) * (m[:Lstar]^(1-m[:α])) / m[:Φ]
    m[:cstar]    = (1-m[:g_star])*m[:ystar] - m[:istar]
    m[:wl_c]     = (m[:wstar]*m[:Lstar])/(m[:cstar]*m[:λ_w])

    return m
end


function model_settings!(m::Model805_alt)

    default_settings!(m)

    # Anticipated shocks
    m <= Setting(:n_anticipated_shocks, 6)
    m <= Setting(:n_anticipated_shocks_padding, 20)

    # Data
    m <= Setting(:data_id, 4, "Dataset identifier")

    # Conditional data variables
    m <= Setting(:cond_semi_names, [:obs_nominalrate])
    m <= Setting(:cond_full_names, [:obs_gdp, :obs_nominalrate])

    # Forecast
    m <= Setting(:forecast_pseudoobservables, false)
end



"""
```
parameter_groupings(m::Model805_alt)
```

Returns an `OrderedDict{String, Vector{Parameter}}` mapping descriptions of
parameter groupings (e.g. \"Policy Parameters\") to vectors of
`Parameter`s. This dictionary is passed in as a keyword argument to
`prior_table`.
"""
function parameter_groupings(m::Model805_alt)
    steadystate = [:γ, :α, :β, :σ_c, :h, :ν_l, :δ, :Φ, :S′′, :ppsi,
                   :δ_gdpdef, :Lmean, :λ_w, :π_star, :g_star]

    sticky      = [:ζ_p, :ζ_w, :ι_p, :ι_w, :ϵ_p, :ϵ_w]

    policy      = [:ψ1, :ψ2, :ψ3, :ρ, :ρ_rm]

    processes   = [[:ρ_g, :ρ_μ, :ρ_z_p, :ρ_z, :ρ_b_p, :ρ_b_til, :ρ_σ_w,
                  :ρ_π_star,  :ρ_λ_f, :ρ_λ_w, :η_λ_f, :η_λ_w, :η_gz, :σ_g, :σ_μ,
                  :σ_z_p, :σ_z, :σ_b_til, :σ_b_p, :σ_π_star,
                  :σ_λ_f, :σ_λ_w, :σ_r_m];
                  [Symbol("σ_r_m$i") for i in 1:n_anticipated_shocks(m)]]

    error       = [:δ_gdpdef, :Γ_gdpdef, :ρ_gdp, :ρ_gdi, :ρ_gdpvar, :ρ_gdpdef, :ρ_corepce,
                   :ρ_lr, :ρ_tfp, :σ_gdp, :σ_gdi, :σ_gdpdef, :σ_corepce,
                   :σ_lr, :σ_tfp]

    all_keys     = Vector[steadystate, sticky, policy, processes, error]
    all_params   = map(keys -> [m[θ]::Parameter for θ in keys], all_keys)
    descriptions = ["Steady State",  "Nominal Rigidities", "Policy",
                    "Exogenous Processes",
                    "Measurement"]

    groupings = OrderedDict{String, Vector{Parameter}}(zip(descriptions, all_params))

    # Ensure no parameters missing
    incl_params = vcat(collect(values(groupings))...)
    excl_params = [m[θ] for θ in vcat([:Upsilon, :ρ_μ_e, :ρ_γ, :σ_μ_e, :σ_γ, :γ_gdi, :δ_gdi, :me_level],
                                      [Symbol("σ_r_m$i") for i=2:20])]
    @assert isempty(setdiff(m.parameters, vcat(incl_params, excl_params)))

    return groupings
end
