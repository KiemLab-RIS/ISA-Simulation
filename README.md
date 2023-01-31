# ISA-Simulation
Clonal Diversity Modelling and Simulations

The simulation StemCellProb2 is meant to expand on our simple HSC stem cell model
and expand the HSC model to cover many different blood types (gran/mono/nk etc).
This model expands these blood types over time to create a very simplistic approximation
of what blood would look like. We then simulate taking a 500,000 genome sample of this blood.
We use this blood sample to simulate viral Insertion Site Analysis (ISA). We repeat these samples
over time to simulate (very approximately) what ISA samples will look like over time with our model.

The parameters for blood cell expansion and time delay are not based on any actual animal, they are instead
just values chosen to make a simulation easy and possible.

The file main.cpp contains an explanation of how to run this program.
