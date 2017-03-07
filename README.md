# Rob's Adventure into the Effect of Various Maxfilter Parameters on MEG Source Localisation

I recorded 4 alien datasets with increasing amounts of measurement noise and movement:

1. Very Still (0.1cm)     (noisy1) 
2. Some Movement (0.5cm)      (noisy2)
3. “Child-like” movement (0.5-1cm)      (noisy3)
4. Excessive Movement (>1cm)        (noisy4)

.. and 5 different Maxfilter settings:

1. noSSS (removing coil noise)
2. SSS
3. tSSS
4. MC
5. MC+tSSS

.. and then applied my Fieldtrip source analysis & statistics pipeline.

## Analysis Steps

- Preprocessing was performed using the preprocessing_maxfilter_adventure.m script

- Visual Source Localisation was was performed using the lcmv_analysis.m script

- Visual Source Localisation was was performed lcmv_analysis_auditory.m

