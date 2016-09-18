efficiency_simulation
======
Calculates the detection efficiency for the MiDM analysis of XENON100.
A parameter scan is done for a fixed mass and varying values of the mass splitting and magnetic moment.
## Usage
```
./efficiency_simulation <efficiency_simulation_input.txt>
```
## Input
* efficiency_simulation_input.txt: contains input parameters
    * number of simulated particles per point in parameter space
    * fixed WIMP mass (eV)
    * lower value of mass splitting (eV)
    * upper value of mass splitting (eV)
    * steps between lower and upper value of mass splitting
    * lower value of magentic moment (mu_N)
    * upper value of magnetic moment (mu_N)
    * steps between lower and upper value of magentic moment

## Output
* efficiency_m_deltalow-deltahigh_mulow-muhigh.root 
    * m: WIMP mass
    * deltalow: lower limit of mass splitting
    * deltahigh: upper limit of mass splitting
    * mulow: lower limit of magnetic moment
    * muhigh: upper limit of magentic moment
    * the file contains a TH2D "hist": xbin center = mu, ybin center = delta, bin content = efficiency
