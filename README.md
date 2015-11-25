# liquid-helium-physical-properties-interpolator
This small programm allows you to get physical properties of liquid helium 
(especially superfluid densities and viscosities) at different pressure and temperature 
It performs interpolation on several experimental datasets, the source of the data are in the files.

It doesn't require any external library that are not already installed with python when downloading python XY.

The pressures have to be input in psi (my gauges are in psi, so I use psi and I have no incentive to convert is to mbar 
or a better unit). The temperature you can input varies depending on the property you are interested in (which is utlimately 
depending on the availiable reliable datasets). 

You can go from 0.1 K to 5 K in temperature and from standard vapor pressure SVP to 5 atmosphere in pressure for the density.

You can go from 0.8 K to 4.2 K in temperature and from SVP to 4.4 atmosphere in pressure for the viscosity.

You can also get the lambda transition temperature (not extremely accurate) as a function of pressure.

I made it to be as flexible as possible in terms of inputs, ie if you input numpy arrays for pressure (n elements) 
and temperature (m elements) it will return a m by n matrix.
