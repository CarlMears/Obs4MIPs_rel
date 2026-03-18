from Obs4MIPS_forward_operator.check_model_data import check_model_data
from calc_tb_sounder_simple.sounder_tbs_simple import SounderForwardOperatorSimple
from era5 import era5_monthly_files,read_era5_data_monthly
import numpy as np
import xarray as xr
from pathlib import Path
# for graphical debugging
import matplotlib.pyplot as plt

# This not strictlyneeded, but a lot problems can arise if the wrong Python environment is being used,
# so this is just a sanity check to print the Python executable being used.
import sys
print()
print('------------------------------------------')
print("Python executable being used: ")
print(sys.executable)
print()
print('------------------------------------------')
print()

if __name__ == "__main__":
    month = 1
    year = 2024

    path_to_era5 = Path('/mnt/m/Obs4MIPs_rel/examples/input_data/ERA5/monthly')  # Change this to the path where your ERA5 monthly data is stored.

    # find a list of the ERA5 files needed
    era5_files = era5_monthly_files(year_to_do=year, month_to_do=month, path_to_era5=path_to_era5)

    # read the ERA5 data for the specified month and year.  Model data is a dictionary a 2D and 3D numpy arrays
    model_data = read_era5_data_monthly(era5_files) 

    # some simple checks on the data
    result = check_model_data(model_data, verbose=False)

    # intialize the AMSU forward operator for choosen channel
    amsu_op = SounderForwardOperatorSimple(sat='AMSU', channel='TMT')

    #compute the brightness temperatures for the specified month and year
    tb = amsu_op.compute_tbs(model_data)


    

    plt.figure(figsize=(12,6))
    plt.imshow(np.flipud(tb), vmin=200.0, vmax=300.0, cmap='viridis')
    plt.colorbar(label='Brightness Temperature (K)')
    plt.title(f'{amsu_op.sat}, {amsu_op.channel} TB')
    plt.xlabel('Longitude Index')
    plt.ylabel('Latitude Index')

    plt.show()

    print()