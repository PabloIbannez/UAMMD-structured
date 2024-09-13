Extending Interactor
====================

To add a new Interactor from zero, it is recommendable to read first the original UAMMD's `Interactor <https://uammd.readthedocs.io/en/latest/Interactor/index.html>`_
documentation. There, you will find what the Interactor class is, how to construct a new one, and the most basic example of an Interactor. 
This is useful if you want to implement something significantly different from anything that already exists in UAMMD-structured. 
However, UAMMD-structured offers a variety of Interactors that you can use to easily implement your own Potentials. For example, 
if you just need a new Bond Potential, it is not necessary to create your own Interactor from scratch; you can use a template. 
You can find the Interactors that already exist in the Interactor section of this documentation—take a look at them. 

All interactors have a ``Family``, ``Type``, and ``subType``. These three keys are
necessary to organize the folder structure of UAMMD-structured.
To add a new Interactor, you must first create the file ``MyInteractor.cu``
in the directory ``/structured/src/Interactor/Family/Type/MyInteractor.cu``
and then add it to the file ``/structured/Components.json``, this will tell
the compiler to take into account your new file.

.. code-block:: json

    "Interactor":
        "Family":[
                    ["..."],
                    ["Type","MyInteractor","MyInteractor.cu"],
                    ["..."]
                    ],


In this case the ``subType`` is our new interactor ``MyInteractor``.
Family and Type are simply the organization of folders within
UAMMD-structured. The user is free to create new directories for Family
and Type as long as they remain consistent. If you are adding a new
Interactor to an existing type, it is important to respect the hierarchy.

REGISTER functions at the end of your ``MyInteractor.cu`` will add the
Interactor to the list of available Interactor. Once you recompile
UAMMDlauncher you can use your new Interactor. Just add to the topology/forceField
section of your input a dataEntry with ``type : ["Type", "MyInteractor"]``

-----

.. toctree::
   :maxdepth: 3

   Generic/index
   Potentials/index

