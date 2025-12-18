# During 3D center-out reaching movement of the mouse forelimb, fascicular velocities are greatly affected by target location and can at times exhibit both eccentric and concentric contractions within a movement

# Introduction 
The fusimotor system is critical for limb movement [3]. Murine models help understand mismodulation of stretch reflexes relevant to spasticity in neurological conditions such as stroke [1]. Therefore, there is a need to quantify fascicular velocities for the mouse forelimb to inform experimental protocols for muscle mechanics [2].

# Methods 
Our computational model of the murine forelimb quantifies fascicular velocities for 3D center-out movements.  We extended a 3D, 6-DOF, 21-muscle forelimb model [2] by adding wrist flexion-extension and radio-ulnar deviation in OpenSim 4.5, and reading out fascicular velocities via the Java API in MATLAB. 100 Monte Carlo simulations of 500-ms center-out movements spanned the 3-D workspace via quintic minimum-jerk S-curve trajectories from reference to target. Fascicular velocities were normalized to each muscles’ optimal fiber length L0. Additionally, we randomly perturbed each S-curve and L0 by ±10%.

# Results 
Most simulated peak fascicle velocities across all muscles and targets were within ±1 L0/s, and only 5 muscles surpassed 1 L0/s. The fastest fascicle velocities consistently occurred within Pectoralis Major, Latissimus Dorsi, Medial Deltoid, and Flexor Carpi Radialis.Our Monte Carlo analysis showed that <10% of perturbed trajectories exceeded 1 L0/s. Varying L0 by ±10% yielded modest variability in peak velocities (coefficients of variation ~6-7% for most sensitive muscles) suggesting our predictions are robust to muscle architecture. Most importantly, we find that even center-out movements induce some muscles to have fascicular velocities that switch between eccentric and concentric velocities, as has also been predicted in humans and non-human primates [3].

# Discussion
Murine experimental paradigms to understand fusimotor regulation of reflexes in the forelimb must identify both movement targets to use, and muscles to record from. Our study identifies both, and also shows that fascicular velocities are relatively slow for 500 ms movements. We also confirm predictions that 3D forelimb movements induce target-dependent complex fascicular velocities that are useful to study the fusimotor system.

![alt text](https://github.com/19990314/Mouse-Forelimb-Muscleskeletal-Modeling.git/PolarGraph/fig_export/polar_graph.png)

Figure 1. Maximum normalized velocity (Lo/s) for top 8 muscles for 100 Monte Carlo simulated trajectories. 


## Running guidline


#Step 1. Obtain an OpenSim installation (GUI optional)

You need an OpenSim installation folder that contains sdk/, sdk/Java/, and the native libraries (bin/ or lib/).
The OpenSim GUI is optional;
A standalone opensim-core installation is sufficient.

Follow the link to download the latest Opensim 4.5:     https://simtk.org/frs/?group_id=91


#Step 2. Run configureOpenSim.m once (provided along the folder)

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
