VibroSim Simulator
==================

VibroSim Simulator is a set of tools to organize and facilitate
simulating Vibrothermography testing using VibroSim.

PLEASE NOTE THAT VIBROSIM MAY NOT HAVE BEEN ADEQUATELY 
VALIDATED, THE NUMBERS BUILT INTO IT ALMOST CERTAINLY 
DO NOT APPLY TO YOUR VIBROTHERMOGRAPHY PROCESS. ITS OUTPUT 
CANNOT BE TRUSTED AND IS NOT SUITABLE FOR ENGINEERING 
REQUIREMENTS WITHOUT APPLICATION- AND PROCESS-SPECIFIC 
VALIDATION. 

VibroSim Simulator relies directly on the following packages:
  * VibroSim_COMSOL
  * angled_friction_model
  * crackclosuresim2
  * VibroSim_WelderModel (if ultrasonic welder-base excitation
    is to be simulated)
  * Limatix

As such, it requires MATLAB, COMSOL, and Python. The COMSOL
Structural Mechanics Module and COMSOL LiveLink for MATLAB
are also necessary. The Python version must be at least 2.7 
and should include the full IPython, Matplotlib, Numpy, Scipy 
stack as well as Pandas v0.17.1 or later. To build 
crackclosuresim2 you will also need the platform compiler 
for your Python version (see the crackclosuresim2 documentation 
for more information).

The Git version control system and the GitPython bindings
are strongly recommended.

While the current implementation uses COMSOL for vibration
calculation and for heat flow evaluation, because of the
modular nature of VibroSim it would be reasonably
straightforward to re-implement those steps using other
tools.


Installation
------------

Like most the other VibroSim components, VibroSim simulator is 
a Python package. Use the usual

    ``python setup.py build``
    
    ``python setup.py install``

commands to install it. Then look in the ``examples/`` folder.

Windows Installation
--------------------

Additional steps must be followed for this package to work on 
Windows.
If COMSOL is to be used, the command line executables for COMSOL 
need to be added to the path. For COMSOL 5.4 these executables 
were installed to the following directory:

	 ``C:\Program Files\COMSOL\COMSOL54\Multiphysics\bin\win64``

This directory needs to be added to the end of the ``path`` 
environment variable. Searching "environment variables" in the
start menu is a good way to find where to make this edit.

VibroSim Simulator Workflow
---------------------------

The VibroSim simulator workflow splits the process of performing a
VibroSim simulation into a sequence of steps. Each step can be run
manually, but as the manual steps can be rather complicated we
recommend the use of the ProcessTrak tool from the Limatix
package to automate the execution of the steps.

The conceptual steps involved in a VibroSim simulation are:
  1. Creation of a geometric model
  2. Vibrational analysis
  3. Modeling of the vibrothermography excitation 
  4. Prediction of heating power
  5. Modeling of the heat from the crack conducting through the
     material


ProcessTrak
-----------

ProcessTrak is a commandline tool from the Limatix package that is
used to keep track of what has been performed in a multistep
process. It is executed by typing ``processtrak`` at the command line.
ProcessTrak is configured by an XML listing of input file steps in a
``.prx`` file.  Usually ProcessTrak is run referencing that ``.prx`` file
followed by additional instructions for what is desired. For example

     ``processtrak cantilever.prx --status``

will list out the status of each step in the process for each input
file.

ProcessTrak is designed to process input data into output results. The
input data is specified in the form of an XML "experiment log" (``.xlg``
file). The experiment log specifies or references the various inputs.
The first ProcessTrak step is always an implicit ``copyinput`` step
which copies the input ``.xlg`` to an output "processed experiment log"
(``.xlp``) file.  The processed experiment log is annotated with
Provenance information, log output from the various processing steps,
and the result data from each processing step. For example,

     ``processtrak cantilever.prx -s copyinput``

will run the implicit ``copyinput`` step on the input files listed in
``cantilever.prx``, generating an output ``.xlp`` file (the input ``.xlg``
is never touched).


Git and Limatix-Git
--------------------

Having confidence in simulation output requires confidence that you
executed the code you intended and confidence that you have a
repeatable process. We recommend the use of Git and Limatix-Git
to perform version management both on the scripts and parameters
of the simulation and on the generated output from the simulation.
Specifically, entering

     ``git init``

in your simulation directory will create a new Git repository there. 

We recommend managing your simulation process with two branches:
"master" which contains the scripts and instructions but no output,
and "master_processed" which also includes processed output.
(These two branches can of course themselves be branched as desired).

The ``limatix-git`` program exists to help automate the process of
committing changed scripts and simulation output to the proper
branches. It is based on the assumption that the name of any
branch intended to contain processed output ends with ``_processed``.
It operates on the principle that scripts, input data, etc. should
be committed to the master branch, and processed output should be
cross-referenced in the ``.xlp`` files.

To add files to the unprocessed branch, check out that branch,
run ``limatix-git add -a`` to stage files for commit, ``git status``
to verify only input files have been staged, and ``git commit``
to perform the commit. 

To add files to the processed branch, check out that branch, run
``limatix-git add-processed -a`` to stage files for commit, ``git status``
to verify only processed output has been staged, and ``git commit`` to
perform the commit.


COMSOL-based VibroSim Workflow
------------------------------

The COMSOL-based VibroSim workflow follows roughly the conceptual
steps listed above, but the model creation is nominally all done
up-front (in reality the first few steps will be iterated to get
the model where it needs to be). 

The steps involved in a COMSOL-based VibroSim simulation are:
  1. Scripting COMSOL to create a geometric and physics model,
     including mounting, excitation position/couplant,
     vibration monitoring, and a healed internal boundary
     representing the crack, 
  2. Vibration analysis of sample including:

    a. Modal analysis
    b. Spectrum verification
    c. Frequency response calculation
    d. Generation of time-domain response. 

  3. Modeling of the vibrothermography excitation to evaluate
     response at the crack
  4. Prediction of heating power from response at the crack.
  5. Modeling of the heat from the crack conducting through the
     material to the surface. 

Building the Docs
------------------------------

This documentation was built using `Sphinx
<https://www.sphinx-doc.org/en/master/>`_. Documentation source code can be
found in the ``docs`` folder. If you are using Fedora, Sphinx can be installed
using the following command:

``dnf install python-sphinx``

Or similarly for Ubuntu:

``apt-get install python3-sphinx``

Once Sphinx is installed an html version of the documentation can be built
using the makefile in the ``docs`` folder.

``make html``

On Windows Sphinx can be installed using ``pip``.

``pip install sphinx``

Sphinx can also be used to create ``.tex`` source files, which can be converted
to pdf using Latex.
