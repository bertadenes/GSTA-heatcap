Installation of requirements

download the latest python environment, currently 3.5.2
also available at /home/berta/Python-3.5.2.tar.xz
extract tarball:
tar xf Python-3.5.2.tar.xz
configure python:
cd Python-3.5.2.tar.xz
./configure --prefix=$HOME
compile:
make
make test
make install
define pythonpath:
NOTE: local can be .local instead
echo "export PYTHONPATH="${PYTHONPATH}:/home/berta/local/lib/python3.5/site-packages"" >> $HOME/.bashrc
. $HOME/.bashrc
now you should have easy_install-3.5
install cclib 1.4(!) and pip:
easy_install-3.5 --prefix=$HOME/local cclib-1.4
easy_install-3.5 --prefix=$HOME/local pip
install numpy and scipy, order is important:
pip3.5 install --user numpy
pip3.5 install --user scipy

NOT NEEDED after v2.0.0
pip3.5 install --user portalocker
pip3.5 install --user --upgrade portalocker

Installation

currently the names of the program and the daemon are fixed to "entropy" and "g09d" therefore 
commands need to be added:
THIS SHOULD BE DONE USING ABSOLUTE PATH:
ln -s /wherever/this/program/is/main.py $HOME/bin/entropy
ln -s /wherever/this/program/is/daemon/g09daemon.py $HOME/bin/g09d

refresh
. $HOME/.bashrc && . $HOME/.bash_profile

Usage

basics v3.0.0:

arguments are optional:
entropy -f [freq.out file] -c [config file] -t [type for continuing session]

will check if it is a freq output

config file (all optional):
type=value
1 for vibration wise approach, 2 for NVE sampling
command=value
value is the command to start g09, will ask if not given
options=value
value is all options you wish to use with commands
maxvib=value
number of modes to calculate for, accepts integers and all keyword, will ask if not given
temp=value
temperature where data should be calculated, will ask if not given
end=value
value can be:
1 -- create velocity generating inputs.
2 -- run velocity generating inputs.
3 -- create BOMD inputs.
4 -- run BOMD inputs.
5 -- process BOMD data.
target=value
for setting up the a variable for filtering, default is 1 for energy-type variables, not yet implemented

continuing calculation:
entropy velgen.out -f velgen/MD.out -t velgen/MD -c must be given to refresh end variable (see above)
only one CalcFile_ can be in the working directory

To continue a session, -t flag can be used with "velgen" or "MD" values. The checkpoint CalcFile must be in the parent directory.
