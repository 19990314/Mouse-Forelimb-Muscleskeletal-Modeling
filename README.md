# Mouse-Forelimb-Muscleskeletal-Modeling Running guidline


Step 1. Obtain an OpenSim installation (GUI optional)

You need an OpenSim installation folder that contains sdk/, sdk/Java/, and the native libraries (bin/ or lib/).
The OpenSim GUI is optional;
A standalone opensim-core installation is sufficient.

Follow the link to download the latest Opensim 4.5:     https://simtk.org/frs/?group_id=91


Step 2. Run configureOpenSim.m once (provided along the folder)

In MATLAB, run:

configureOpenSim


This script registers the JAR and native libraries with MATLAB’s Java class path / Java library path (and also detects and comments out older OpenSim entries to avoid conflicts). After it finishes, restart MATLAB for the changes to take effect.

Step 3. Verify
import org.opensim.modeling.*;
disp(char(org.opensim.modeling.Model.getVersion()));


Note for the instructor: The GUI is not required. The code only needs the OpenSim core libraries and Java bindings; configureOpenSim.m registers org-opensim-modeling.jar and the native libraries with MATLAB.


Folders introduction:
Workspace: Randomly generates workspace using real time clock seed, and saved as mouse_workspace_samples.m
The other three folders: Muscle Fiber Velocity, Polar Graph, and Sensitivity test all used a provided mouse_workspace_samples.m which was generated during our simulations.


Formula deduction was contained in workspace folders for interpretating Minimum Jerk Trajectory.
