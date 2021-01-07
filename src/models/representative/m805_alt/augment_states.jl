"""
```
augment_states{T<:AbstractFloat}(m::AbstractDSGEModel, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T})
```

### Arguments

-`m`: the m object
-`TTT`, `RRR`, and `CCC`: matrices of the state transition equation

### Return values

- `TTT_aug`, `RRR_aug`, and `CCC_aug`: extend the corresponding input matrices to include
  jobservables which are growth rates.

### Description

Some observables in the model are growth rates, which are calculated as a linear combination
of a present and lagged state (which is not yet accounted for in the `TTT`, `RRR`,and `CCC`
matrices). To improve the performance of `gensys`, these additional states are added after
the model is solved. `augment_states` assigns an index to each lagged state, and extends the
input `TTT`, `RRR`, and `CCC` matrices to accommodate the additional states and capture the
lagged state value in the current state vector. `RRR` and `CCC` are mostly augmented with
zeros.

The diagram below shows how `TTT` is extended to `TTT_aug`.

                TTT_aug
     (m.endogenous_states_additional
                  x
     m.endogenous_states_additional)
     _________________________________
    |                     |           |
    |          TTT        | endog_    |
    | (endogenous_states  | states_   |
    |          x          | augmented |
    |  endogenous_states) |           |
    |_____________________|           |
    |                                 |
    |    endogenous_states_augmented  |
    |_________________________________|

"""
function augment_states(m::Model805_alt{T}, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T}) where {T<:AbstractFloat}
    endo = m.endogenous_states
    endo_addl = m.endogenous_states_augmented
    exo = m.exogenous_shocks

    n_endo = n_states(m)
    n_exo = n_shocks_exogenous(m)
    @assert (n_endo, n_endo) == size(TTT)
    @assert (n_endo, n_exo) == size(RRR)
    @assert n_endo == length(CCC)

    # Initialize augmented matrices
    n_addl_states = length(m.endogenous_states_augmented)
    n_addl_eqs = n_addl_states
    TTT_aug = zeros(n_endo + n_addl_eqs, n_endo + n_addl_states)
    TTT_aug[1:n_endo, 1:n_endo] = TTT
    RRR_aug = [RRR; zeros(n_addl_eqs, n_exo)]
    CCC_aug = [CCC; zeros(n_addl_eqs)]

    ### TTT modifications

    # Track Lags
    TTT_aug[endo_addl[:y_t1], endo[:y_t]] = 1.0
    TTT_aug[endo_addl[:c_t1], endo[:c_t]] = 1.0
    TTT_aug[endo_addl[:i_t1], endo[:i_t]] = 1.0
    TTT_aug[endo_addl[:w_t1], endo[:w_t]] = 1.0
    TTT_aug[endo_addl[:π_t1], endo[:π_t]] = 1.0
    TTT_aug[endo_addl[:L_t1], endo[:L_t]]  = 1.0
    TTT_aug[endo_addl[:e_gdp_t1], endo_addl[:e_gdp_t]] = 1.0
    TTT_aug[endo_addl[:e_gdi_t1], endo_addl[:e_gdi_t]] = 1.0

    ## We construct state for expected inflation using the fact
    ##
    ##   E_t[s_{t+1}] = CCC + TTT*s_{t} = CCC + TTT * (CCC + TTT*s_{t-1} + RRR)
    ##                = (CCC+TTT*CCC) + (TTT^2)s_{t-1} + TTT*RRR
    ##
    ## So to construct the state for E_t[p_{t+1}], we need to make all these
    ## objects, and just index out the row relevant to pi_t

    # Expected inflation
    TTT_aug[endo_addl[:Et_π_t], 1:n_endo] = (TTT^2)[endo[:π_t], :]

    TTT_aug[endo_addl[:lr_t], endo_addl[:lr_t]]               = m[:ρ_lr]
    TTT_aug[endo_addl[:tfp_t], endo_addl[:tfp_t]]             = m[:ρ_tfp]
    TTT_aug[endo_addl[:e_gdpdef_t], endo_addl[:e_gdpdef_t]]   = m[:ρ_gdpdef]
    TTT_aug[endo_addl[:e_corepce_t], endo_addl[:e_corepce_t]] = m[:ρ_corepce]
    TTT_aug[endo_addl[:e_gdp_t], endo_addl[:e_gdp_t]]         = m[:ρ_gdp]
    TTT_aug[endo_addl[:e_gdi_t], endo_addl[:e_gdi_t]]         = m[:ρ_gdi]


    ### RRR modfications

    # Expected inflation
    RRR_aug[endo_addl[:Et_π_t], :] = (TTT*RRR)[endo[:π_t], :]

    # Measurement Error on long rate
    RRR_aug[endo_addl[:lr_t], exo[:lr_sh]] = 1.0

    # Measurement Error on TFP
    RRR_aug[endo_addl[:tfp_t], exo[:tfp_sh]] = 1.0

    # Measurement Error on GDP Deflator
    RRR_aug[endo_addl[:e_gdpdef_t], exo[:gdpdef_sh]] = 1.0

    # Measurement Error on Core PCE
    RRR_aug[endo_addl[:e_corepce_t], exo[:corepce_sh]] = 1.0

    # Measurement Error on GDP Deflator
    RRR_aug[endo_addl[:e_gdp_t], exo[:gdp_sh]] = 1.0
    RRR_aug[endo_addl[:e_gdp_t], exo[:gdi_sh]] = m[:ρ_gdpvar] * m[:σ_gdp] ^ 2

    RRR_aug[endo_addl[:e_gdi_t], exo[:gdi_sh]] = 1.0

    ### CCC Modifications

    # Expected inflation
    CCC_aug[endo_addl[:Et_π_t]] = (CCC + TTT*CCC)[endo[:π_t]]

    return TTT_aug, RRR_aug, CCC_aug
end
