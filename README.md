# Pedestrian-Walk
Simulation of Pedestrian Random Walk using Brownian Motion coded in Matlab.

## Introduction

Pedestrian movement is a complex phenomenon influenced by various factors, such as the environment, social context, and personal preferences. To better understand how people move in crowded areas and avoid collisions, I studied the pedestrian flow on a CSULB campus in front of the bookstore.

This project, I used [Euler–Maruyama method](https://en.wikipedia.org/wiki/Euler%E2%80%93Maruyama_method) for solving stochastic differential equations numerically and to analyze how people adjust their intended path to avoid obstacles and maintain a comfortable distance from others. By examining factors like walking speed, direction, and decision-making, I gained insights into pedestrian behavior during different times of the day.

While studying pedestrian movement can be challenging due to the numerous factors involved, it is essential to improving safety and efficiency in crowded areas. These findings can inform the design and management of public spaces, helping to create more comfortable and safer environments for everyone.

<p align="center">
  <img width="600" height="300" src="https://github.com/LiliyaSemenenko/Pedestrian-Walk/blob/main/space_animation.gif">
</p>
<p align="center">
Figure of the space including the bookstore (black rectangle on the top), trashcans (small squares on the top colored black and blue), and tents (big gray squares on the bottom) and pedestrian flow within the space. Pedestrians are represented as red circles.
</p>

# Features

[Click here to watch a simulation video](https://youtu.be/EuF8drNEohI)

- Pedestrians can walk both directions (right and left). This was achieved by a direction inversion, which is the change of the longitudinal velocity $u → −u$ of the pedestrian’s walking direction.
- Pedestrians can turn to the bookstore and enter it. Once the pedestrian is inside a bookstore, they are not leaving it.
- Pedestrians can turn down and walk in between the tents, eventually exiting the space.
- Pedestrians can go under the tents and stay there for some amount of time, then possibly leave and continue walking along their mean path.
- If a pedestrian walks to one of the trashcans, it could bump into it and immediatelly come back to the center of the space, continuing walking on their mean path.

# Simulation description
## Parameters used in the simulation

This simulation involves the use of several parameters to control the behavior of pedestrians in a virtual environment. These parameters are as follows:

- $t_{0}$: Initial time, which specifies the starting point of the simulation.
- $x_{0}$ and $y_{0}$: Initial positions of the pedestrians.
- $u_{0}$: Longitudinal velocity, representing the initial speed in the x-direction.
- $v_{0}$: Transversal velocity, representing the initial speed in the y-direction.
- $u_{1}$ and $u_{2}$: Speed constraints, defining the minimum and maximum allowable speeds for the pedestrians.
- $dt$: Time step, which plays a crucial role in the simulation by defining the interval at which the pedestrians' positions and velocities are updated.
- $\alpha$, $\beta$, and $\gamma$: Coefficients for controlling the pedestrians' movement, influencing their acceleration and deceleration patterns in the x and y directions.
- $\sigma_{x}$ and $\sigma_{y}$: Standard deviations of the random noise added to the pedestrians' positions, introducing randomness and variability into the simulation.
- $\kappa$: Influences the size of the pedestrians in the simulation, allowing for a diverse crowd with varying physical characteristics.
- $\theta_{1}$ and $\theta_{2}$: Define the range of transparency values for the pedestrians, allowing for the representation of partially opaque or translucent individuals.
- $offset$: Enables fine-tuning of the pedestrians' positions in the simulation, introducing small displacements to avoid collisions or align the pedestrians with specific areas or paths within the virtual environment.

It is important to note that adjusting these parameters can greatly impact the behavior and movement patterns of the pedestrians in the simulation. Careful consideration and testing should be undertaken to ensure a realistic and consistent virtual environment.

## Mean paths

We can define an average path by considering all pedestrian trajectories that connect L to R and vice versa.



First part of the code is used to create a mean path for pedestrians in a simulation. The mean path is created using splines and is defined by the functions y_star, dy_star, and ddy_star. These functions are used to generate a trajectory that a group of pedestrians can follow.

The code initializes the values of several parameters and generates a random initial position for the first pedestrian. It then creates an array 'ped' that will store the positions and velocities of all the pedestrians at each time step. The first row of 'ped' is initialized with the initial position and velocity of the first pedestrian.

## Method 

Simulating pedestrian movement is a complex task that requires careful consideration of various factors, including the environment, pedestrian behavior, and crowd dynamics. To achieve this, I used the Euler-Maruyama method, a numerical scheme for approximating stochastic differential equations (SDEs).

I started by modeling pedestrian movement as a random process, where the direction and speed of each pedestrian are influenced by various factors, such as the presence of obstacles, proximity to other pedestrians, and the location of the destination. Using the Euler-Maruyama method, I created a simulation that approximates the position and velocity of pedestrians as they navigate through the virtual environment.

To achieve this, I discretized time and used a time step $dt$ to simulate the movement of pedestrians at each time interval. I used a set of SDEs from Euler–Maruyama method to model the stochastic behavior of pedestrian movement, including the velocity and direction of each pedestrian. A noise term was also incorporated to account for the random movements of pedestrians in real-life scenarios using random variables $\Delta W_{n}$.

<p align="center">
<b>Euler–Maruyama method</b>
</p>

Consider the stochastic differential equation

$$
d X_{t} = a\left(X_{t}, d t\right)d t + b\left(X_{t}, t\right) dW_{t}
$$

with initial condition $X_{0}=x_{0}$, where $W_{t}$ stands for the Wiener process, and suppose that we wish to solve this SDE on some interval of time. Then the Euler-Maruyama approximation to the true solution $X$ is the Markov chain $Y$ defined as follows:


- partition the interval $0, \tau$ into $N$ equal subintervals of width $\Delta t>0$ :

$$
0=\tau_{0}<\tau_{1}<\cdots<\tau_{N}=T \text { and } \Delta t=T / N
$$

- set $Y_{0}=x_{0}$

- recursively define $Y_{n}$ for $0 \leq n \leq N-1$ by

$$
Y_{n+1}=Y_{n}+a\left(Y_{n}, \tau_{n}\right) \Delta t+b\left(Y_{n}, \tau_{n}\right) \Delta W_{n}
$$

where

$$
\Delta W_{n}=W_{\tau_{n+1}}-W_{\tau_{n}} .
$$

The random variables $\Delta W_{n}$ are independent and identically distributed normal random variables with expected value zero and variance $\Delta t$.

<br />
<br />

<p align="center">
<b>Normal vector, radius of curvature, and centripetal force</b>
</p>

Before making any SDEs, we need to calulate the normal vector ($\vec{n}$), rho ($\rho$), and centripetal force ($F_{\text {cent}}$) terms to model the curvature of the path and the effects of centripetal forces on pedestrian motion.

- $\vec{n}$: represents the normal vector pointing towards the center of curvature of the path that the pedestrian is following. It is calculated based on the slope of the path at the current position of the pedestrian.

$$
\vec{n} = \frac{1}{\sqrt{1+\left(\frac{dy_{\text{star}}}{dX}\right)^2}}\begin{bmatrix} \frac{dy_{\text{star}}}{dX} & -1 \end{bmatrix} \
$$

- $\rho$: represents the radius of curvature of the path at the current position of the pedestrian. It is calculated based on the second derivative of the path equation.

$$
\rho = \frac{\left|1 + \left(\frac{dy_{\text{star}}}{dX}\right)^2\right|^{1.5}}{\left|\frac{d^2y_{\text{star}}}{dX^2}\right|} \
$$

- $Centripetal Force$: represents the force acting on the pedestrian towards the center of curvature of the path. It is calculated based on the magnitude of the velocity of the pedestrian and the radius of curvature of the path at the current position of the pedestrian.

$$
F_{\text {cent}} = \vec{n} \cdot \frac{U^2 + V^2}{\rho}
$$



<br />
<br />


<p align="center">
<b>Stochastic differential equation we used and applied within the code</b>
</p>

$$
d X=U \cdot d t
$$

$$
d U=-\alpha U\left(U-u_{\text {intial }}\right) \cdot d t+F_{\text {xcent }} \cdot d t+F_{\text {xrepul }} \cdot d t+\sigma_{x} \cdot d W
$$

$$
d Y=V \cdot d t
$$

$$
d V=-\beta \left(Y-y_{\text {star }}(x)\right) \cdot d t-V_{\gamma} \cdot d t+F_{ycent} \cdot d t+F_{yrepul} \cdot d t+\sigma_{y} \cdot d W
$$


<br />
<br />


Next, we created a virtual environment using a rectangular corridor with a bookstore, tents, and trash cans to simulate pedestrian flow. We created these objects using specifications and metrics to ensure they were accurately represented in the virtual environment. To control pedestrian movement, we used a set of algorithms to ensure that pedestrians move in the right direction and avoid collisions with each other and objects.

We then ran the simulation, adjusting the parameters and refining the code to ensure smooth and realistic pedestrian movement. By using the Euler-Maruyama method, we were able to create a robust simulation of pedestrian movement that captures the stochastic nature of pedestrian behavior in crowded environments. This knowledge can help urban planners and designers to design safer and more efficient public spaces.



## Limitations

Avoiding Pedestrian Merging: Difficulties arose in preventing pedestrians from merging together as they traveled from the left and right sides of the space. The code is not using offset and repulsive force to keep pedestrians separated and maintain appropriate distances between them.
