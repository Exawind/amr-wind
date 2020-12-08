.. _dev-documenting:

Documentation - user manual and source code docs
================================================

AMR-Wind comes with two different types of documentation:

- User & developer manuals, such as the document you are reading now, that are
  written using `Sphinx <https://www.sphinx-doc.org/en/master/index.html>`_, and

- Inline documentation within C++ source code that are written in a format that can be
  processed automatically by `Doxygen <http://www.doxygen.nl/manual/index.html>`_

User documentation
------------------

AMR-Wind user documentation is written using a special format called
ReStructured Text (ReST) and is converted into HTML and PDF formats using a
python package Sphinx. Since the manuals are written in simple text files, they
can be version controlled alongside the source code. Documentation is
automatically generated with new updates to the GitHub repository and deployed
at `AMR-Wind documentation site <https://exawind.github.io/amr-wind>`_.

Building documentation
``````````````````````

To update and build documentation locally on your system you will need Sphinx
installed on your system. Please consult the `Sphinx installation page
<https://www.sphinx-doc.org/en/master/usage/installation.html>`_ for
instructions to install Sphinx. Users of `Anaconda python distribution
<https://www.anaconda.com/>`_ can install Sphinx by executing the following
command within their desired environment.

.. code-block:: console

   conda install sphinx

Once Sphinx is installed properly, you should have access to the
:program:`sphinx-build` executable. To build documentation follow the instructions below:

.. code-block:: console

   # Switch to AMR-Wind build directory
   cd ~/exawind/source/amr-wind/build/

   # Build docs in HTML format
   sphinx-build -M html ../docs/sphinx .

The above command will build docs in HTML format and the output is placed in
:file:`html` directory within the :file:`build` directory. Documentation can
also be generated in other formats, consult `Sphinx docs
<https://www.sphinx-doc.org/en/master/usage/builders/index.html>`_ for available
formats and their usage.

Writing user documentation
``````````````````````````

As mentioned previously, documentation is written using a special text format
called reStructuredText. Sphinx user manual provides a `reST Primer
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_ that
provides an overview of this format and how to write documentation using this format.


Documenting source code
-------------------------

Source code (C++ files) are commented using a special format that allows Doxygen
to extract the annotated comments and create source code documentation as well
as inheritance diagrams. API documentation for the latest snapshot of the
codebase can be browsed online `here
<https://exawind.github.io/amr-wind/api_docs>`_. To build the documentation
locally, first install ``doxygen`` and ``graphviz`` executables on your system.
Once they are successfully installed, execute the following command from the
root directory of ``amr-wind``

.. code-block:: console

   doxygen ./docs/doxygen/Doxyfile

The default format is HTML, and upon successful completion of the above command,
the documentation files are available in :file:`build/html` directory. Open
:file:`build/html/index.html` on your browser to browse the locally generated
documentation.

`Doxygen manual <http://www.doxygen.nl/manual/index.html>`_ provides an overview
of the syntax that must be used. Please follow the Doxygen style of commenting
code when commenting AMR-Wind sources.

When commenting code, try to use self-documenting code, i.e., descriptive names
for variables and functions that eliminate the need to describe what is going on
in comments. In general, comments should address *why* something is being coded
in a particular way, rather than how the code does things. Try to write the code
in a clear manner so that it is obvious from reading the code instead of having
to rely on comments to follow the code structure.
