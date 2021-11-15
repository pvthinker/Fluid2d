#!/bin/bash

pydir=$HOME/.fluid2d
default=core/defaults.json
srcdir=`pwd`
myexpdir=$srcdir/myexp

echo "--------------------------------------------------------------------------------"
echo ""
echo " Installing Fluid2d"
echo ""
echo "It is recommended to create a virtual environnement. "
echo "If you are using anaconda then you may create the 'pyrsw' environment with"
echo ""
echo "> conda env create -f environment.yml"
echo ""
echo "then whenever you want to use Fluid2d, switch to this environement with"
echo ""
echo "> conda activate fluid2d"
echo ""
echo "Are you ok with you environment? (y/n)"
if [ $ok = "n" ]; then
   exit 42
fi

if [ ! -d "$pydir" ]; then
    echo "  Create $pydir"
    mkdir $pydir
fi
if [ ! -f "$pydir/$default" ]; then
    echo "  Copy $default in $pydir"
    cp $default $pydir/
fi
if [ ! -d "$myexpdir" ]; then
    echo "  Create $myexpdir"
    mkdir $myexpdir
fi
echo "  Copy reference experiments in $myexpdir"
cp -pR $srcdir/experiments/* $myexpdir

# for bash users
cat > $pydir/activate.sh << EOF
export PYTHONPATH=`pwd`/core
echo Python now knows that Fluid2d is in `pwd`
EOF

# for csh, tcsh users
cat > $pydir/activate.csh << EOF
setenv PYTHONPATH `pwd`/core
echo Python now knows that Fluid2d is in `pwd`
EOF

# for fish users
cat > $pydir/activate.fish << EOF
set -gx PYTHONPATH `(pwd)`/core
echo Python now knows that Fluid2d is in `(pwd)`
EOF

# compile the modules with module
echo "--------------------------------------------------------------------------------"
echo ""
echo " Compile modules with f2py"
echo ""
{
    # try
    make
    
} || {
    #catch
    echo "Unable to compile"
    echo "Are you sure f2py is installed?"
    exit
}

# copy the experiment into

echo ""
echo "  Before starting, please read this note carefully"
echo ""
echo "  As it configured, Fluid2d will store the results in"
echo ""
echo "      *** $HOME/data/fluid2d ***"
echo ""
echo "  If you don't run the code from your laptop then it is likely that"
echo "  you are not allowed to store large binary files on $HOME/data"
echo "  because this is your home."
echo ""
echo "  In this case, edit $HOME/.fluid2d/defaults.json"
echo "  and set 'datadir' default value (in the output section) to a"
echo "  place where you are authorized to store large binary files."
echo "  In the jargon, this place is usually a 'work' directory."
echo ""
echo "  If you are unsure where to store your results check that"
echo "  with your system administrator."
echo ""
echo "  Once you know where to store the results, then you're good to go"
echo ""
echo "  Each time you open a new terminal you need to"
echo "     source ~/.fluid2d/activate.sh  if you're under bash"
echo "     source ~/.fluid2d/activate.csh  if you're under csh/tcsh"
echo "     source ~/.fluid2d/activate.fish  if you're under fish"
echo ""
echo "  To run your first experiment"
echo "      cd $myexpdir/Vortex"
echo "      python3 vortex.py"
echo ""
echo "  Write your new experiments in $myexpdir"
echo "  or wherever you want, but not in $srcdir/experiments !"
echo ""

