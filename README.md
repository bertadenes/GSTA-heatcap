# Installation of requirements

## Python 3
Download one version of python 3 environment, available [here](https://www.python.org/downloads/).
This program was tested on python 3.5.2

*Help for installation on unix systems*

```commandline
tar xf Python-3.5.2.tar.xz #extract tarball
cd Python-3.5.2
./configure --prefix=$HOME ##configure python
make #compile
make test
make install
echo "export PYTHONPATH="${PYTHONPATH}:$HOME/local/lib/python3.5/site-packages"" >> $HOME/.bashrc #define pythonpath NOTE: local can be .local instead
. $HOME/.bashrc
```

## Help for installation of python packages
We recommend easy_istall and pip, if you do not have root access.
```commandline
easy_install-3.5 --prefix=$HOME/local cclib-1.4 #should be 1.4 *
easy_install-3.5 --prefix=$HOME/local pip
pip3.5 install --user numpy #make sure to install numpy first
pip3.5 install --user scipy
```

Unfortunately cclib-1.5 does not work properly with our code, so version 1.4 is recommended.
If you'd still like to use 1.5 you'll need to modify the source code of cclib.
The simplest (and least elegant) way to do this is to find and modify the following file:
```commandline
/home/user/lib/python3.5/site-packages/cclib/parser/gaussianparser.py
```
where lines 677 and 678 have to be commented out with a '#' symbol at the beginning of each line.
The code should look something like this:
```python
                    if len(parts) > 4:
                        energy = parts[2].split('=')[1]
                        if energy == "":
                            energy = self.float(parts[3])
                        else:
                            energy = self.float(energy)
                   # Next 2 lines commented by User
                   # if len(self.scftargets[0]) == 3:  # Only add the energy if it's a target criteria
                   #     newlist.append(energy)
                    scfvalues.append(newlist)
```


# Installation

Download this repository and note where you decide to store it. For proper execution we strongly suggest to create a
symbolic link to the main file: 
```commandline
ln -s /wherever/this/program/is/main.py $HOME/bin/GSTA-hc
```
Make sure to use the absolute path, otherwise the link might be broken.

After reloading the environment, the program should be ready to use.
```commandline
. $HOME/.bashrc && . $HOME/.bash_profile
```

# Usage

The utility will run interactively simply by executing it:
```commandline
GSTA-hc
```
For more advance usage, flags are available:
```commandline
GSTA-hc -f [freq.out file] -c [config file] -t [type]
```

**Config file**

The parameters given in a config file are optional. The program may request interaction for necessary information.

*example:*
parameter = value

parameter|var type|value
---|---|---
command| (str) | *command to run Gaussian 09*
options| (str) | *all options you wish to use with command*
temp| (float) | *temperature where data should be calculated*
NTraj| (int) | *number of independet trajectories*
NCalc| (int) | *number of Gaussian calculations allowed to be run at one time*
rotation| (bool) | *include or exclude rotation*
end| (int) | *value*

Note, that even though the number of independent trajectories can be one,
for the linear regression and determination of heat capacity it must be at least two.

*value* can be:
1. create velocity generating inputs.
2. run velocity generating inputs.
3. create BOMD inputs.
4. run BOMD inputs.
5. process BOMD data. *Not available automatically.*

## Process the data

After the simulations all terminated, the data is ready to be processed to obtain heat capacities. Unfortunately the
automatization of the process is not yet stable, so we recommend to you the postprocess function:
 ```commandline
GSTA-hc -t pp -f .Calcfile*
```
where .Calcfile* is the binary file used as a checkpoint during the run. This results in an output *cv.out.

## Further remarks

The directory structure of one run is important. Do not move or remove files from the associated directories.
Once Gaussian calculations have been submitted, a process monitors them in the background. This can be followed by
the watchlog files in the directories. This daemon recalls the main program when a calculation is finished. This can
be done manually if needed:
```commandline
GSTA-hc -f [file to process] -t [type] -c [config]
```
Configuration file is only needed if the initial run was set to an end already achieved (*see above*).
Type can be *MOVE*(only if rotation was excluded) *VELGEN* and *MD*.
The simulations' length is set to 1 ps, the timestep is 0.1 fs at the moment. We wish to optimize these parameters to
actual cases and implement the algorithm.

Due to differences in output syntax simulation with Gaussian16 is not stable. At the moment we suggest to use Gaussian09
instead.