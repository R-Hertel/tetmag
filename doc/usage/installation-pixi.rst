Installation using pixi
======================

This section describes how to build and run tetmag using `pixi`_.
This is an optional method that provides a self-contained and reproducible
environment without relying on system-wide installations of compilers and
libraries.

What is pixi?
-------------

`pixi`_ is a cross-platform package manager and workflow tool based on the
conda-forge ecosystem. It allows users to define project-specific environments
that include compilers, build tools, and libraries. These environments are
isolated from the system and help avoid dependency conflicts.

Install pixi
------------

pixi can be installed using the official installer:

.. code-block:: bash

   curl -fsSL https://pixi.sh/install.sh | bash

After installation, ensure that pixi is available:

.. code-block:: bash

   pixi --version

If necessary, restart your shell or ensure that ``~/.pixi/bin`` is in your
``PATH``.

Set up the environment
----------------------

Clone the tetmag repository and install the required dependencies:

.. code-block:: bash

   git clone https://github.com/R-Hertel/tetmag.git
   cd tetmag
   pixi install

This will create a local environment in the ``.pixi/`` directory containing all
required tools and libraries.

Compile tetmag
--------------

Configure and build the project using pixi:

.. code-block:: bash

   pixi run configure  # runs "cmake -S . -B build"
   pixi run build  # runs "cmake --build build -j"

The executable will be created in the ``build/`` directory.

Run tetmag
----------

You can run tetmag either through pixi:

.. code-block:: bash

   pixi run build/tetmag <input-file>

or directly:

.. code-block:: bash

   ./build/tetmag <input-file>

Notes
-----

- Pixi provides its own compiler toolchain and libraries, so no system-wide
  installation of dependencies is required.
- This makes the installation more portable across different Linux distributions.
- Changing compiler versions or versions of dependencies may be easier
  than for a given Linux distribution.
- It is recommended to use ``pixi run`` for configuration and build steps to
  ensure the pixi environment is active.
- If you encounter build issues, a clean rebuild can help:

  .. code-block:: bash

     rm -rf build
     pixi install

.. _pixi: https://pixi.sh
