# Fluid2d
Fluid2d is a versatile Python-Fortran CFD code that solves a large
class of 2D flows. It is designed to be used easily by Students
learning Fluid Mechanics or Geophysical Fluid Dynamics and by their
professors willing to illustrate their course and to organize
numerical practicals. The idea is to visualize flows on the fly, as
they are computed. The effect of parameter changes can be seen
immediately. The key quantity of fluid2D is the vorticity. If you feel
weak on vorticity dynamics, this code is for you. You should rapidly
become as expert as the experts.

You can learn how basic processes work because of the power of
animations. It is quite easy to go beyond textbooks and to reach
research questions.


# Install
Fluid2d installation does not rely on the standard `pip install Fluid2d`. The procedure is a bit more convoluted

  1) git clone Fluid2d. To copy-paste from your computer to the virtual desktop, use the copy-paste interactive window on the left

> git clone https://github.com/pvthinker/Fluid2d.git

  2) create the conda environment (this step may fail in that case, execute the next step). It takes a few minutes to download everything. It’s a bit a shame to have to download all these packages in your $HOME but that’s how it works.

> cd Fluid2d

> conda create --name fluid2d --file requirements.txt

if conda create fails

> conda init bash

then close the terminal, reopen a new one, module load anaconda3 and repeat the conda

> create --name fluid2d --file requirements.txt

That should work.

  3) activate this environment

> conda activate fluid2d

  4) you can now build Fluid2d

> ./install.sh

  5) last step, make Python aware of where Fluid2d is

> source ~/.fluid2d/activate.sh

you’re good to run Fluid2d

  6) run your first experiment

> cd myexp/Vortex

> python vortex.py


# Run an experiment
  Once the code is installed, you don't need to repeat stage 4).

  What you need to do though, every time you run Fluid2d in a new
  terminal is to repeat stages 3) and 5), then do something like 6)
  but with another experiment.

# Most frequent error

  When trying to run an experiment you may get the following error message

  ModuleNotFoundError: No module named 'fluid2d'

  it is most likely because you forgot to activate Fluid2d. Fix it with

> source ~/.fluid2d/activate.sh

  and it should work