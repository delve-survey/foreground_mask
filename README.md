# foreground_mask
Foreground mask for the delve footprint


# Y6A2 Foreground mask for the LSS sample (from David and Nacho)

Scripts to generate a foreground mask for large and bright objects affecting the LSS
sample.

Please, be aware that there are several hardwire paths spread in the code (mainly pointing
to necessary catalogues of stars) which need to be modified accordingly to your setup.

Raw description of the generation of the mask:

1. Compute the avoidance radii per each object in the different stars catalogues using
   main.py 

2. Have a look at the computed radii with radius_code.py (the definition is based in the
   relative density of good and bad objects near the bright object)

3. Modify y6a2_mask_bits.py according to the previously estimated radii

4. Define the final mask considering all the previous information with
   y6a2_make_foreground_mask_28042023.py 

If in doubt, do not hesitate to contact david.sanchez@ciemat.es
