
# Generate EoR lightcone
# uses modules: py21cmfast, tools21cm
python generate_lightcone.py

# Generate Exgal source catalogue
trecs -c -w -p input_files/trecs_SDC3a.ini

#Convert the Exgal source catalogue to an image cube
python dependencies/SKAO-SDC/run_SDC3_pipeline.py input_files/ska-sdc_SDC3a.ini

#call ska-low_simulator
python dependencies/ska-low-simulator/sim_low.py input_files/ska-low-sim_SDC3a.ini








