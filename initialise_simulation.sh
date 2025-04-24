# ----------------------------------------------------------------|
#                 SETUP NEW SIMULATIONS                           |
# ----------------------------------------------------------------|
#                                                                 |
# This is a bash code to setup a new diablo simulation. It        |
# executes the follwing steps                                     |
#                                                                 |
# 1. Copy the grid files ygrid01.txt, ..., ygrid16.txt            |
#                                                                 |
#                   y01.txt                                       |
#                      .                                          |
#                      .                                          |
#                      .                                          |
#                   y16.txt                                       |
#                                                                 |
#                                                                 |
# 2. Removes and creates again (empty) the files:                 |
#        * diss.txt                                               |
#        * kin.txt                                                |
#        * output.dat                                             |
#        * pot.txt                                                |
#        * tot.txt                                                |
#                                                                 |
#                                                                 |
# 3. Creates the path.dat file with the current location of       |
#    the workin folder                                            |      
#                                                                 |
#       * path.dat                                                |
#                                                                 |   
# 4. Initialise the file "save_flow.ind.txt" with a 0 in then     |
#    first line                                                   |
#                                                                 |
#       * save_flow.ind                                           |   
#                                                                 |
# ----------------------------------------------------------------|
# Written by      : Jorge Sandoval (jsandoval001.dundee.ac.uk)    |
# Role            : PDRA, University of Dundee                    |
# Project         : Saving energy via drag reduction: a           |
#                   mathematical description of oscillatory flows |
# PI              : Tom Eaves (teaves001@dunde.ac.uk)             |
#                                                                 |
# Previous update : 03/01/2024                                    |
# Last update     : 20/03/2024                                    |
#   DETAILS       : The necessary instructions for the Bisection  |
#                   Edge Tracking mode were added                 |
# ----------------------------------------------------------------|

# INPUT FOLDERS
# change this whenever is needed
grid_folder="/home/jorge/wall_stokes_ET_2"
source_folder="/home/jorge/wall_stokes_ET_2"

# ----------------------------------------------------------------|
# 1. COPY GRID FILES AND EXECUTABLE FILE

# Check if the source folder exists
if [ ! -d "$grid_folder" ]; then
    echo "grid folder not found. Exiting."
    exit 1
fi

# Loop over the grid files
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16; do
    file="ygrid$i.txt"
    cp "$grid_folder/$file" .
done

file="diablo"
cp "$source_folder/$file" .

# ----------------------------------------------------------------|

# 2. REMOVE OUTPUT FILES

rm -f *_out*

if [ -d last_ET_it_flow_field ]; then
    # Delete all files inside the folder
    rm -f last_ET_it_flow_field/*
else
    mkdir last_ET_it_flow_field
fi

if [ -d latest_turbulent_solutions ]; then
    # Delete all files inside the folder
    rm -f latest_turbulent_solutions/*
else
    mkdir latest_turbulent_solutions
fi

if [ -d latest_laminar_solutions ]; then
    # Delete all files inside the folder
    rm -f latest_laminar_solutions/*
else
    mkdir latest_laminar_solutions
fi

if [ -f energy_shift_history.dat ]; then
  rm -f energy_shift_history.dat
fi

if [ -f lambda_shift_history.dat ]; then
  rm -f lambda_shift_history.dat
fi

if [ -f sp_shifts_history.dat ]; then
  rm -f sp_shifts_history.dat
fi

if [ -f stop.now ]; then
  rm -f stop.now
fi

if [ -f diss.txt ]; then
  rm -f diss.txt
fi

if [ -f kin.txt ]; then
  rm kin*
fi

if [ -f output.dat  ]; then
  rm -f output.dat 
fi

if [ -f pot.txt  ]; then
  rm -f pot.txt 
fi

if [ -f tot.txt  ]; then
  rm -f tot.txt 
fi

touch diss.txt
touch kin.txt
touch pot.txt
touch tot.txt
touch energy_shift_history.dat
touch lambda_shift_history.dat
touch sp_shifts_history.dat

# ----------------------------------------------------------------|
# 2. CREATE path.dat FILE

if [ -f path.dat  ]; then
  rm -f path.dat 
fi

touch path.dat

# Get the current directory.
CURRENT_DIR=$(pwd)/

# Write the current directory to the file "path.dat" in two
# consecutive lines (input and output directories)

echo $CURRENT_DIR >> path.dat
echo $CURRENT_DIR >> path.dat

# ----------------------------------------------------------------|
# 3. INITIALISE THE FILE "save_flow.ind.txt"

# Check if the file exists.
if [ -f save_flow.ind ]; then
  # Remove the file.
  rm -f save_flow.ind
fi

# Create the file again.
touch save_flow.ind

# Write a 0 to the first line and the fourth column.
echo "    0" | cat >> save_flow.ind

# ----------------------------------------------------------------|
# 4. Run the program with 16 cores and write the output onto the
# file output.dat

# mpirun -n 16 diablo > output.dat &

# ----------------------------------------------------------------|



