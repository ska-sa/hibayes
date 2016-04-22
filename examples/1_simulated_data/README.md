## Simulated data

Firstly, run nested sampling routine from the project root directory:

      mpiexec -n 8 ./hi_multinest.py examples/1_simulated_data/config.ini

Next, we can create a triangle plot:

      ./hi_plot.py examples/1_simulated_data/config.ini

And we can run the reconstruction with:

      ./hi_recon.py examples/1_simulated_data/config.ini