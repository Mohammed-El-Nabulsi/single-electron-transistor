Getting started
===============

``HaPPPy`` is available as a Git repository at https://trac.physnet.uni-hamburg.de/git/Happpy

Step 1. Get the code
--------------------

If you are seeing this text, you probably have found the Git repository. Clone the repository by typing

.. code-block:: bash

    git clone https://trac.physnet.uni-hamburg.de/git/Happpy

in a terminal and authenticate with your Physnet credentials. If access is denied, please send an email to lfrahm@physnet.uni-hamburg.de and ask for permission.

If cloning was successful, you have now a directory named *Happpy* on your computer, which holds the source code of the program. Type :code:`ls` to check the location of the folder. The terminal should prompt something like

.. code-block:: console

    user@machine:~/$ ls
    SomeFolderA  Happpy  SomeFolderB  SomeFolderC  SomeFolderD

Step 2. Run the first tests
---------------------------
Let's move with the terminal into that directory and run the very first test. The comand :code:`cd Happpy` changes the directory of the terminal into *Happpy*. Typing :code:`ls` again shows all files in the directory. The output should look like

.. code-block:: console

    user@machine:~/Happpy$ ls
    docs  HaPPPy  LICENSE  MANIFEST.in  README.md  setup.py  tests

To run the first tests, run :code:`python3 -m tests.AllTests`. This runs all tests preset tests. The output looks like

.. code-block:: console

    user@machine:~/Happpy$ python3 -m tests.test_main
    tests_HaPPPy_main_exists (__main__.HappyBasicTestSuite) ... ok
    tests_HaPPPy_main_runs (__main__.HappyBasicTestSuite) ... ok

    ----------------------------------------------------------------------
    Ran 2 tests in 0.000s

    OK

which means all tests ran successfully. Get in contact with Lars or write in the Mattermost channel if the tests fail on your machine.

.. note :: If Python 3 is not installed on your machine, get it from https://www.python.org/
Step 3. Run HaPPPy for the first time
-------------------------------------

Now let's run HaPPPy. To run the program itself (not the tests), type 

.. code-block:: shell

    python3 -m HaPPPy

You will see a welcome message and information about your setup.

You have now successfully setup HaPPPy and are able to dive into the code.

Step 4. Build the Documentation 
-------------------------------
Next you should build the documentation. You need *sphinx* and *spinx_rtd_theme* in order to do that. 
If you installed a recent Python version you can just run 

.. code-block:: shell

    pip3 install sphinx
    pip3 install sphinx_rtd_theme

If installation was successful, you are ready to build the docs. 

You find a directory called *docs* int the project, which holds the documentation.
To build the documentation move into the directory and type run `make html`

.. code-block:: shell

    cd docs
    make html

.. note:: You can also run other targets like `make latexpdf` or `make text`, but in most cases we will use the html output.

This creates a directory called *build* in your *docs* directory. In that *build* directory you find a directory called *html*
which holds a file called *index.html*. So the path to the file is something like 

.. code-block:: shell

    /Happpy/docs/build/html/index.html

Double click it. A browser should open showing you the HaPPPy documentation. Have a look.
