this repository collects the codes I typically use in the Advanced Macroeconomics II or Macroeconomic Models and Methods courses (2nd semester MA course at WU (Vienna University of Economics and Business))

NKmodel_loglin.mod: Dynare code that solves the already log-linearized equations of the Galí (2008), chapter 3 baseline New Keynesian model; gives solution in terms of variables in percentage deviations from steady state (hat variables)
(NKmodel_loglinFLEX.mod: same but in the special case of flexible prices)

NKmodel_nonlin.mod: Dynare code that codes up the nonlinear system of (my teaching version) of the baseline New Keynesian  model. This is almost the same as Galí chapter 3, just that I use Rotemberg quadratic price adjustment cost instead of Calvo pricing and a linear-in-labor production function. Uses variable transformation (x_t = exp(log(x_t))) to obtain a solution to the stochastic growth model in terms of percentage deviations from steady state (log-variables)
(NKmodel_nonlinFLEX.mod: same but in the special case of flexible prices)
