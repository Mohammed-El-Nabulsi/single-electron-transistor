The HaPPPy program
------------------

The **Ha**\ mburg **P**\ rogrammier **P**\ rojekt **Py**\ thon is the program written by the students of the course *66-527 Proseminar: Programmierprojekt* at the University of Hamburg. The program serves to simulate various properties of quantum dots.

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

To run the first tests, run :code:`python3 -m test.test_main`. This runs all tests preset tests. The output looks like

.. code-block:: console

    user@machine:~/Happpy$ python3 -m test.test_main
    test_HaPPPy_main_exists (__main__.HappyBasicTestSuite) ... ok
    test_HaPPPy_main_runs (__main__.HappyBasicTestSuite) ... ok

    ----------------------------------------------------------------------
    Ran 2 tests in 0.000s

    OK

which means all tests ran successfully. Get in contact with Lars or write in the Mattermost channel if the tests fail on your machine.

Step 3. Run HaPPPy for the first time
-------------------------------------

Now let's run HaPPPy. To run the program itself (not the tests), type 

.. code-block:: shell

    python3 -m HaPPPy

You will see a welcome message and information about your setup.

You have now successfully setup HaPPPy and are able to dive into the code.