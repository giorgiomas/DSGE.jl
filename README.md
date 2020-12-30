# New York Fed DSGE Model (Version 1002)
[![Build Status](https://travis-ci.org/FRBNY-DSGE/DSGE.jl.svg)](https://travis-ci.org/FRBNY-DSGE/DSGE.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://frbny-dsge.github.io/DSGE.jl/dev)
[![Coverage Status](https://coveralls.io/repos/github/FRBNY-DSGE/DSGE.jl/badge.svg?branch=master)](https://coveralls.io/github/FRBNY-DSGE/DSGE.jl?branch=master)

The `DSGE.jl` package implements the New York Fed dynamic stochastic general equilibrium (DSGE) model and provides general code to estimate many user-specified DSGE models. The package is introduced in the Liberty Street Economics blog post
[The FRBNY DSGE Model Meets Julia](http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html).
(We previously referred to our model as the "FRBNY DSGE Model.")

This Julia-language implementation mirrors the MATLAB code included in the
Liberty Street Economics blog post
[The FRBNY DSGE Model Forecast](http://libertystreeteconomics.newyorkfed.org/2015/05/the-frbny-dsge-model-forecast-april-2015.html).

Documentation for the *code* can be accessed by clicking on the `docs|dev` button above. For documentation about the most recent *model version*, read this [pdf](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/docs/DSGE_Model_Documentation_1002.pdf).

The New York Fed DSGE team is currently extending the code to solve and estimate heterogeneous agent models. Filtering and smoothing algorithms are available in the registered package [StateSpaceRoutines.jl](https://github.com/FRBNY-DSGE/StateSpaceRoutines.jl).
An implementation of Sequential Monte Carlo (SMC) sampling, used for the estimation of DSGE models, can be found in the registered package [SMC.jl](https://github.com/FRBNY-DSGE/SMC.jl). The foundational `AbstractModel` type, from which the `AbstractDSGEModel` type derives, is defined in the registered package [ModelConstructors.jl](https://github.com/FRBNY-DSGE/ModelConstructors.jl).

Further extensions of the DSGE model code may be released at the discretion of the New York Fed.

## Installation

`DSGE.jl` is a registered Julia package in the [`General`](https://github.com/JuliaRegistries/General) registry. To install it, open your Julia REPL, type `]` to enter the package manager, and run

```julia
pkg> add DSGE
```

If you use any code that loads data (e.g. the example script `run_default.jl` and `make_packet.jl`), then you need make sure you have a FRED API key by following these [instructions for the FredData.jl package](https://github.com/micahjsmith/FredData.jl).

If you are using Windows OS and you encounter the error `AssertionError: length(dirs) == 1`, please see this [issue](https://github.com/JuliaLang/Pkg.jl/issues/1943). Additionally, please do not run the `plot.jl` test if you are using Windows OS
because the generated output will violate the default filename length restriction on Windows. If you want to run this test, then
you need to enable [long paths](https://docs.microsoft.com/en-us/windows/win32/fileio/naming-a-file#enable-long-paths-in-windows-10-version-1607-and-later).

*Note we do not test our code in Windows OS, so we cannot guarantee the code works properly in Windows.*

## Versioning

`DSGE.jl` is currently compatible with Julia `v1.x` (as of `v1.1.6`).

To use `DSGE.jl` with Julia `v0.7`, please check out tag `0.8.1`. To do this, click on the drop-down menu that reads `branch:master` on the left-hand side of the page. Select `tags`, then `v0.8.1`.  If you've already cloned the repo, you can simply run `git checkout v0.8.1`.

To use `DSGE.jl` with Julia `v0.6`, please check out tag `0.4.1`.

## Precompilation

The `DSGE.jl` package is not precompiled by default because when running code in parallel, we want to re-compile
the copy of `DSGE.jl` on each processor to guarantee the right version of the code is being used. If users do not
anticipate using parallelism, then users ought to change the first line of `src/DSGE.jl` from

```
isdefined(Base, :__precompile__) && __precompile__(false)
```

to

```
isdefined(Base, :__precompile__) && __precompile__(true)
```

Disclaimer
------
Copyright Federal Reserve Bank of New York. You may reproduce, use, modify, make derivative works of, and distribute and this code in whole or in part so long as you keep this notice in the documentation associated with any distributed works. Neither the name of the Federal Reserve Bank of New York (FRBNY) nor the names of any of the authors may be used to endorse or promote works derived from this code without prior written permission. Portions of the code attributed to third parties are subject to applicable third party licenses and rights. By your use of this code you accept this license and any applicable third party license.

THIS CODE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT ANY WARRANTIES OR CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, EXCEPT TO THE EXTENT THAT THESE DISCLAIMERS ARE HELD TO BE LEGALLY INVALID. FRBNY IS NOT, UNDER ANY CIRCUMSTANCES, LIABLE TO YOU FOR DAMAGES OF ANY KIND ARISING OUT OF OR IN CONNECTION WITH USE OF OR INABILITY TO USE THE CODE, INCLUDING, BUT NOT LIMITED TO DIRECT, INDIRECT, INCIDENTAL, CONSEQUENTIAL, PUNITIVE, SPECIAL OR EXEMPLARY DAMAGES, WHETHER BASED ON BREACH OF CONTRACT, BREACH OF WARRANTY, TORT OR OTHER LEGAL OR EQUITABLE THEORY, EVEN IF FRBNY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES OR LOSS AND REGARDLESS OF WHETHER SUCH DAMAGES OR LOSS IS FORESEEABLE.


# My additions

## Models
I created a model called my1010, adding some elements to m1010:
- equilibrium conditions (to keep track of some not specified variables):
  - spread between the return on capital and the risk-free rate
  - entrepreneurs' leverage
  - convenience yield
  - ex-ante real rate
  - credit
  - nominal deposit rate
- pseudo-observables and related measurement equations:
  - to compare model estimate vs actual path of variables observed with measurement error: GDP, GDI, PCE inflation, GDP deflator, Aaa spread, Baa spread, 10-year yield
  - nominal deposit rate, real deposit rate, 10-year real yield

Further, I added two models, called m1010depo and my1010depo, which are identical to m1010 and my1010 respectively, apart from the addition of the interest rate on deposits as observable (also the real deposit rate is included in the pseudo-observables).
Only in my1010 and my1010depo, the observable and measurement equations for the nominal interest rate on deposits are subject to an if-condition that checks the model setting `:add_deposits`: this is set to `true` in `defaults.jl`, and can be modified in `my1010.jl` and `my1010depo.jl`.

`m1010depo_alt` is another version of `m1010`, where
- the deposit rate is added as observable,
- the Baa-Treasury spread and the Aaa-Treasury spread are dropped from the observables, in order to make the model comparable, in terms of data, to Smets and Wouters (2007).
The distinction between liquidity and safety is therefore omitted, while the distinction between transient and permanent convenience yield shocks is maintained; accordingly, the measurement errors on the two spreads are also deleted, while the steady state convenience yield is set equal to the sum of the steady state liquidity and safety components.

## Input data
So far, the data considered for interest rate on deposits are:
- the secondary market interest rate on 3-month certificates of deposits, available since 1964:III (source is the OECD but data are available on FRED as `IR3TCD01USQ156N`)
- the arithmetic average of secondary market interest rates on 1-, 3-, and 6-month certificates of deposits, together available from 1966:I to 2013:II (`CD1M`, `CD3M`, and `CD6M` on FRED); these are the data used by Hollander and Liu (2016), and Pesaran and Xu (2016)
Actually, `IR3TCD01USQ156N` and `CD3M` coincide.


## Estimation
I fixed the bug in the Metropolis-Hastings step by building a new function in `estimate/smc/helpers.jl`, named `my_generate_param_blocks`; this is called in `metropolis-hastings.jl`. In order for the function to be recognized, uncomment `include("estimate/smc/helpers.jl")` in `DSGE.jl`.
Whenever, we run estimation without reoptimizing, that is we compute the posterior distribution starting from the pre-computed mode and Hessian, we should check that the correct files are in `input_data/user/model_name`, and we should check that the correct path is specified in the main file.
