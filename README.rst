=========
Soil Heat
=========


.. image:: https://img.shields.io/pypi/v/soil_heat.svg
        :target: https://pypi.python.org/pypi/soil_heat

.. image:: https://readthedocs.org/projects/soil-heat/badge/?version=latest
        :target: https://soil-heat.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status


A Python library for calculating soil heat flux using various scientific models.

This package provides a collection of Python functions for calculating soil heat flux and related soil thermal properties. The implementations are based on various models and parameterizations published in the peer-reviewed scientific literature, making it a useful tool for researchers and practitioners in hydrology, meteorology, and agricultural sciences.

The library is designed for ease of use, with a modular structure and vectorized functions that leverage NumPy for high performance.

Features
--------

*   **Multiple Models**: Implements a wide range of soil heat flux estimation methods from the scientific literature.
*   **Modular Design**: Each scientific paper or set of related equations is organized into its own module for clarity and easy cross-referencing.
*   **Vectorized Calculations**: Leverages NumPy for efficient, element-wise calculations on time-series data.
*   **Comprehensive Documentation**: All public functions include high-quality, Numpy-style docstrings with mathematical formulas, parameter descriptions, and examples.

Installation
------------

You can install `soil-heat` using pip:

.. code-block:: bash

    pip install soil-heat

Usage
-----

Here is a basic example of how to use one of the functions from the library. This example uses the `gao2010_gz` function, which estimates soil heat flux based on a sinusoidal ground-temperature forcing.

.. code-block:: python

    import numpy as np
    from soil_heat import gao2010_gz

    # --- Parameters for a daily cycle ---
    # Amplitude of the sinusoidal ground-surface temperature (K)
    AT = 8.0
    # Soil thermal conductivity (W m-1 K-1)
    lambda_s = 1.2
    # Soil thermal diffusivity (m2 s-1)
    kappa = 1.0e-6
    # Time array for one day (in seconds), sampled every 15 minutes
    t_day = np.linspace(0, 86400, 97)

    # --- Calculate heat flux ---
    Gz = gao2010_gz(AT, lambda_s, kappa, t_day)

    print("Calculated heat flux at the first 5 time steps:")
    print(Gz[:5])


Implemented Models
------------------

The library includes implementations from the following key publications:

*   **Gao, Z., et al. (2017)**: Soil thermal conductivity parameterization.
*   **Liebethal, C., & Foken, T. (2006)**: Evaluation of six parameterization approaches for the ground heat flux.
*   **Wang, Z.-H., & Bou-Zeid, E. (2012)**: A novel approach for the estimation of soil ground heat flux.
*   **Yang, K., & Wang, J. (2008)**: A temperature prediction-correction method for estimating surface soil heat flux.

Please refer to the documentation for each module for a detailed list of the implemented equations.

Running Tests
-------------

To run the test suite, first clone the repository and install the development dependencies:

.. code-block:: bash

    git clone https://github.com/inkenbrandt/soil-heat.git
    cd soil-heat
    pip install -r requirements_dev.txt

Then, run pytest from the root directory:

.. code-block:: bash

    pytest

Credits
-------

This package was created by Paul Inkenbrandt. It was originally created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
