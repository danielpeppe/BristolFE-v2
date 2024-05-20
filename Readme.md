
# BristolFE v2 - Daniel Peppe (Fork)

This is a fork of BristolFE v2 used in my dissertation: **Generating Synthetic Ultrasonic Testing Datasets to Advance Data-Driven Porosity Detection in Fibre-Reinforced Composites**.

## Examples

The `examples-DP` directory contains two main example scripts:

1. **run_single**: 
    - Run a single iteration of the simulation.
    - Useful for debugging and animating.

2. **run_many**:
    - Main data generation script.

## Default Options

Default options are defined in the `fn_set_options` function.

## Helper Functions

- **plot_batch_from_save**:
    - Plots all responses in a batch.
  
- **fn_process_mat_data**:
    - Converts batches of .mat files to a CSV format.

# BristolFE v2 - Paul Wilcox
=======================

This repository contains a number of Matlab functions and example scripts for performing basic Finite Element (FE) simulations, in particular explicit dynamic ones. Currently only 2D models are supported and the only element types are 3-noded triangular ones for elastic solids (CPE3), fluids (AC2D3), and  fluid-solid interface elements (ASI2D2).

To use, clone (or download and unzip) the repository and add the BristolFE-v2/code folder to the Matlab path.

The entry-point for FE analysis is the function fn_BristolFE_v2 in the BristolFE-v2/code folder. The BristolFE-v2/code folder also contains numerous helper functions for creating meshes and displaying results.

The core FE code is in various functions in the BristolFE-v2/code/private folder, which are not called directly by the user.

Some example scripts are provided in the BristolFE-v2/examples. Most likely you will start with one of these and modify it according to your requirements.

OVERVIEW
========

The entry point function is res = fn_BristolFE_v2(mod, matls, steps, fe_options). 

When this function is called, a complete mesh must be specified (in mod), the materials used must be defined (in matls), and one or more loading steps and the required outputs defined (in the cell array steps).

The requested results for the corresponding loading step are returned in the cell array res. Typical outputs are one or both of: 
    1. History outputs -complete time histories of the displacement (or pressure in fluids) at one or mode nodes, typically plotted as time-domain signals.
    2. Field output - snapshots of the complete wavefield (its local kinetic energy) at intervals in time, typically displayed as a movie and used as a visualisation tool.

Most of the code in the example scripts is concerned with preparing mod, matls, and steps before fn_BristolFE_v2 is called and then displaying the outputs.

EXAMPLES
========

In BristolFE-v2/examples you will find the following scripts which provide simple examples how to set up different features in models:
    1. fluid_example.m - simulate pressure waves in a fluid domain
    2. solid_example.m - simulate longitudinal and shear waves in a solid domain
    3. coupled_solid_fluid_example.m - simulate waves in a fluid domain coupled to a solid one, showing mode conversions at the interface
    4. absorbing_layer_example.m - same as 3 but this time with an absorbing layer on 3 sides of the domain to prevent reflections
