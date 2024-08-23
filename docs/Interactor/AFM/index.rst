Atomic Force Microscopy
=======================

To simulate Atomic Force Microscopy (AFM), we follow the scheme used by previous studies. In this model, the AFM tip is represented as a sphere connected to a virtual point, referred to as the chip, through an anisotropic spring:

.. math::

   U_{\text{tip-chip}} = \frac{1}{2}K(z_{\text{tip}}-z_{\text{chip}})^2 + \frac{1}{2}K_{xy}\left[ (x_{\text{tip}}-x_{\text{chip}})^2 + (y_{\text{tip}}-y_{\text{chip}})^2\right]

Adjusting the position of the chip allows us to manipulate the location of the tip. Typically, K_xy is significantly larger than K to restrict lateral movement. To simulate indentation, we decrease the position of the chip, dragging the AFM tip primarily along the z axis due to the large value of K_xy. The deflection force of the tip can be measured by calculating the magnitude K(z_tip-z_chip). Figure 1 illustrates this setup. AFM potential does not include the surface interaction, which must be added separately using an external potential.

.. figure:: /img/afm_tmv.png
   :alt: Virtual AFM

   The figure shows a model of the tobacco mosaic virus. Also included are the AFM tip and the chip, represented by a blue dot. The particle representing the tip is connected to the chip through the potential described by the equation above. As the position of the chip is lowered, the AFM tip indents the sample. A surface is also added, which is not defined by the AFM interaction and must be included separately.

Once we have added the tip particle and the sample we can add the AFM potential using the interaction of type ["AFM","SphericalTip"]:

.. code-block:: yaml

   topology:
     forceField:
       AFM:
         type: ["AFM", "SphericalTip"]
         labels: ["idSet_i", "idSet_j",
                         "epsilon", "sigma",
                         "K", "Kxy",
                         "tipVelocity", "startChipPosition",
                         "indentationStartStep"]
         data:
           - [[0], [1, 2, 3, "..."],
              -1.0, 1.0,
              1.0, 10.0,
              -0.5, [0.0, 0.0, 50.0], 0]

Initially, we must identify which particle represents the AFM tip and which group constitutes the sample. The AFM tip is a normal particle, requiring definition in Sections like "state" and "structure" within the Data Entry using specific labels/data fields. "idSet_i" indicates the particle for the tip (particle 0 in our example), while "idSet_j" designates the particle group of the sample. Parameters such as "epsilon," "sigma," "K," and "Kxy" define the characteristics of the tip. The interaction between the tip and the sample is modeled using a modified Lennard-Jones potential for larger, non-point-like particles:

.. math::

   U_{\text{tip-sample}} = \sum_{i \in \text{sample}}{\varepsilon_{\text{tip}} \left [\left( \frac{\sigma_{\text{tip}}}{|\vec{r}_i - \vec{r}_{\text{tip}}| - R_{\text{tip}}} \right)^{12} - 2 \left( \frac{\sigma_{\text{tip}}}{|\vec{r}_i - \vec{r}_{\text{tip}}| - R_{\text{tip}}} \right)^6 \right ]}

Here, ε_tip and σ_tip are the predefined values, and R_tip is the radius of the particle representing the tip. The position of this particle is denoted by vec{r}_tip.

Finally, we define the velocity of the chip ("tipVelocity") along the z-axis, its initial position ("startChipPosition"), and the moment indentation begins ("indentationStartStep"). The position of the chip at simulation step n > nStart (where nStart is the step when indentation starts) is given by:

.. math::

   z_{\text{chip}} = {z_{\text{chip}}}_{0} + v_{\text{chip}}\cdot (n - nStart) \cdot \text{dt}.

The position of the chip in the xy plane remains fixed. Consequently, the chip will progressively lower (for v_chip < 0), dragging the tip due to the potential U_tip-chip. If the tip encounters an obstacle, such as the particles of the sample, it will exert increasing force upon them.

This configuration allows us to simulate the AFM indentation process on a sample. Additional AFMs can be included as rows in the Data Entry, operating in a similar manner.

----

There are different AFM tip models:

.. toctree::
   :maxdepth: 1

   SphericalTip
   SphericallyBluntedConicTip
