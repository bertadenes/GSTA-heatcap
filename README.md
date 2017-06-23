#Installation of requirements

##Python 3
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

##Help for installation of python packages
We recommend easy_istall and pip, if you do not have root access.
```commandline
easy_install-3.5 --prefix=$HOME/local cclib-1.5 #should be 1.4 or over
easy_install-3.5 --prefix=$HOME/local pip
pip3.5 install --user numpy #make sure to install numpy first
pip3.5 install --user scipy
```

#Installation

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

#Usage

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

|parameter|var type|value|
command| (str) | *command to run Gaussian 09*
options| (str) | *all options you wish to use with command*
temp| (float) | *temperature where data should be calculated*
NTraj| (int) | *number of independet trajectories*
NCalc| (int) | *number of Gaussian calculations allowed to be run at one time*
rotation| (bool) | *include or exclude rotation*
end| (int) | *value*

*value* can be:
1. create velocity generating inputs.
2. run velocity generating inputs.
3. create BOMD inputs.
4. run BOMD inputs.
5. process BOMD data. *Not available automatically.*

##Further remarks

The directory structure of one run is important. Do not move or remove files from the associated directories.
Once Gaussian calculations have been submitted, a process monitors them in the background. This can be followed by
the watchlog files in the directories. This daemon recalls the main program when a calculation is finished. This can
be done manually if needed:
```commandline
GSTA-hc -f [file to process] -t [type] -c [config]
```
Configuration file is only needed if the initial run was set to an end already achieved (*see above*).
Type can be *MOVE*(only if rotation was excluded) *VELGEN* and *MD*.