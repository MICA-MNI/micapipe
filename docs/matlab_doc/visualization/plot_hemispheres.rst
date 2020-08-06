.. _plot_hemispheres_matlab:

.. role:: matlab(code)
   :language: matlab

==================
plot_hemispheres
==================

------------------
Synopsis
------------------

Plots data on the cortical surface (`source code
<https://github.com/MICA-MNI/BrainStat/blob/master/matlab/visualization/%40plot_hemispheres/plot_hemispheres.m>`_).


------------------
Usage
------------------

::

   obj = plot_hemispheres(data,surface,varargin);

- *data*: an n-by-m vector of data to plot where n is the number of vertices or parcels, and m the number of markers to plot.
- *surface*: a surface readable by :ref:`convert_surface_matlab` or a two-element cell array containing left and right hemispheric surfaces in that order. 
- *varargin*: Name-Value Pairs (see below).
- *obj*: an object that can be used to further modify the figure.

------------------
Description
------------------

Plots any data vector onto cortical surfaces. Data is always provided as a
single vector; if two surfaces are provided then the *n* vertices of the first
surface will be assigned datapoints 1:*n* and the second surface is assigned the
remainder. If a parcellation scheme is provided, data should have as many
datapoints as there are parcels.  

BrainStat only provides basic figure building functionality. For more
information on how to use MATLAB to create publication-ready figures we
recommend delving into `graphics object properties
<https://www.mathworks.com/help/matlab/graphics-object-properties.html>`_ (e.g.
`figure
<https://www.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html>`_,
`axes
<https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.axes-properties.html>`_,
`surface
<https://www.mathworks.com/help/matlab/ref/matlab.graphics.primitive.surface-properties.html>`_).
Also see the `source code
<https://github.com/MICA-MNI/BrainSpace/blob/master/matlab/plot_data/plot_hemispheres.m>`_
for basic graphic object property modifications.

Note: This class is a copy of a `BrainSpace <https://brainspace.readthedocs.io/>`_ class by the same name. 

Name-Value Pairs
^^^^^^^^^^^^^^^^^
- *parcellation*: an k-by-1 vector containing the parcellation scheme, where k is the number of vertices. 
- *labeltext*: A cell array of m elements containing labels for each column of data. These will be printed next to the hemispheres. 
- *views*: a character vector containing the angles from which to view the brain. Valid characters are: l(ateral), m(edial), i(nferior), s(uperior), a(nterior), and p(osterior). When supplying a single surface, the lateral/medial designations must be inverted for the right hemisphere. Default is 'lm'. 

-------------------
Methods
-------------------

The plot_hemispheres class has three public methods that may be used to further
modify the figure. 

- *colorlimits*: Used as :matlab:`obj.colorlimits(M)` where M is either a 1-by-2 vector containing the lower and upper limit of the colorbar, or a a n-by-2 matrix containing the lower and upper limits of the colorbar for each plotted row. 
- *colormaps*: Used as :matlab:`obj.colormaps(M)` where M is either a c-by-3 colormap, where c is the number of different colors, or a n-element cell vector where the ith cell contains the colormap for the ith plotted row. 
- *labels*: Used as :matlab:`obj.labels(varargin)` where varargin are name-value pairs containing the new properties of the labeltext e.g. :matlab:`obj.labels('FontSize',10,'FontName',DroidSans)`. Any properties left unset will default to the BrainStat defaults i.e. Rotation: 90, Units: Normalized, HorizontalAlignment: 'Center', FontName: DroidSans, FontSize 18, otherwise default text() parameters. 